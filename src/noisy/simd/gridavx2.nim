## Only scalar arguments cross this module's boundary. All vector operations
## and generic kernel instantiations stay in AVX2-enabled translation units.
import avx2base, valueavx2, valuecubicavx2, opensimplex2avx2,
  opensimplex2savx2, perlinavx2, cellularavx2
include avx2flags

type
  AvxKind* = enum
    avxPerlin, avxOpenSimplex2, avxOpenSimplex2S, avxValue, avxValueCubic,
    avxCellular
  AvxSettings* = object
    seed*: int32
    octaves*: int
    amplitude*, frequency*, lacunarity*, gain*, jitter*: float32
    perm*: ptr array[512, int32]
    distanceFunction*, returnType*: int
  AvxBatch* = object
    kind*: AvxKind
    settings*: AvxSettings

func finishCellular(settings: AvxSettings, features: NearestFeatures): Float8 {.avxInline.} =
  ## Enum ordinals match CellularDistance and CellularReturn in cellulars.nim.
  let one = float8(1)
  var
    distance0 = features.first
    distance1 = features.second
  if settings.distanceFunction == 0 and settings.returnType >= 1:
    distance0 = sqrt8(distance0)
    if settings.returnType >= 2:
      distance1 = sqrt8(distance1)
  case settings.returnType
  of 0: toFloats(features.hash) * float8(1.0'f / 2147483648.0'f)
  of 1: distance0 - one
  of 2: distance1 - one
  of 3: (distance1 + distance0) * float8(0.5'f) - one
  of 4: distance1 - distance0 - one
  of 5: distance1 * distance0 * float8(0.5'f) - one
  else: distance0 / distance1 - one

func sampleCellular[N: static[int]](
  settings: AvxSettings, point: array[N, Float8]
): Float8 {.avxInline.} =
  template metric(trackSecond, trackHash: static[bool]): Float8 =
    template search(squared, manhattan: static[bool]): Float8 =
      finishCellular(settings, searchFeatures8(settings.seed, settings.jitter,
        point, squared, manhattan, trackSecond, trackHash))
    case settings.distanceFunction
    of 0, 1: search(true, false)
    of 2: search(false, true)
    else: search(true, true)
  case settings.returnType
  of 0: metric(false, true)
  of 1: metric(false, false)
  else: metric(true, false)

func sample[N: static[int]](
  settings: AvxSettings, point: array[N, Float8], kind: static[AvxKind]
): Float8 {.avxInline.} =
  var
    amplitude = settings.amplitude
    frequency = settings.frequency
    total = float8(0)
  for _ in 0 ..< settings.octaves:
    var scaled: array[N, Float8]
    for axis in 0 ..< N:
      scaled[axis] = point[axis] * float8(frequency)
    let value =
      when kind == avxPerlin:
        const Scales = [0.188'f, 0.507'f, 0.936'f, 0.87'f]
        float8(Scales[N - 1]) * samplePerlinGrid8(settings.perm[], scaled)
      elif kind == avxCellular:
        sampleCellular(settings, scaled)
      else:
        template evaluate(kernel: untyped): Float8 =
          when N == 2: kernel(settings.seed, scaled[0], scaled[1])
          else: kernel(settings.seed, scaled[0], scaled[1], scaled[2])
        when kind == avxOpenSimplex2: evaluate(sample8)
        elif kind == avxOpenSimplex2S: evaluate(sampleSmooth8)
        elif kind == avxValue: evaluate(sampleValue8)
        else: evaluate(sampleCubic8)
    total += value * float8(amplitude)
    amplitude *= settings.gain
    frequency *= settings.lacunarity
  total / float8(settings.octaves.float32)

func fillRow[N: static[int]](
  settings: AvxSettings, start: array[4, float32], length: int,
  output: ptr UncheckedArray[float32], kind: static[AvxKind]
) =
  const axis = min(N - 1, 2)
  var point: array[N, Float8]
  for i in 0 ..< N:
    point[i] = float8(start[i])
  for offset in countup(0, length - 8, 8):
    point[axis] = coordinates8(start[axis], offset)
    mm256_storeu_ps(output[offset].addr, sample(settings, point, kind))

func fillRowAvx2*(
  batch: AvxBatch, start: array[4, float32], dimensions, length: int,
  output: ptr float32
) =
  ## Fills complete groups of eight; the baseline caller handles scalar tails.
  template fill(kind: static[AvxKind]) =
    template run(n: static[int]) =
      fillRow[n](batch.settings, start, length,
        cast[ptr UncheckedArray[float32]](output), kind)
    when kind == avxPerlin:
      case dimensions
      of 1: run(1)
      of 2: run(2)
      of 3: run(3)
      else: run(4)
    else:
      if dimensions == 2: run(2)
      else: run(3)
  case batch.kind
  of avxPerlin: fill(avxPerlin)
  of avxOpenSimplex2: fill(avxOpenSimplex2)
  of avxOpenSimplex2S: fill(avxOpenSimplex2S)
  of avxValue: fill(avxValue)
  of avxValueCubic: fill(avxValueCubic)
  of avxCellular: fill(avxCellular)

func sampleTileAvx2*(
  batch: AvxBatch, points: array[4, array[8, float32]], dimensions: int,
  output: ptr float32
) =
  ## Internal: samples eight circle/torus positions through the scalar ABI.
  ## The caller validates settings, selects 2D/4D, and supplies eight outputs.
  template run(n: static[int], kind: static[AvxKind]) =
    var point: array[n, Float8]
    for axis in 0 ..< n:
      point[axis] = mm256_loadu_ps(points[axis][0].unsafeAddr)
    mm256_storeu_ps(output, sample(batch.settings, point, kind))
  case batch.kind
  of avxPerlin:
    if dimensions == 2: run(2, avxPerlin)
    else: run(4, avxPerlin)
  of avxOpenSimplex2: run(2, avxOpenSimplex2)
  of avxOpenSimplex2S: run(2, avxOpenSimplex2S)
  of avxValue: run(2, avxValue)
  of avxValueCubic: run(2, avxValueCubic)
  of avxCellular: run(2, avxCellular)
