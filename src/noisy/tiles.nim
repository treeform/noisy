## Seamless loops and tiles via circle/torus mapping of existing noise samplers.
import std/math, common, grids, simplexcommon, perlins, opensimplexes, values,
  cellulars

type
  TileNoise = Simplex | Perlin | OpenSimplex2 | OpenSimplex2S | ValueNoise |
    ValueCubic | Cellular
  TileNoise2D = Simplex | Perlin

when useGridAvx2:
  import simd/[avx2, gridavx2, runtimecheck]

  func tileBatch[T: TileNoise](noise: T): AvxBatch =
    when T isnot Simplex:
      result.settings = AvxSettings(octaves: noise.octaves,
        amplitude: noise.amplitude, frequency: noise.frequency,
        lacunarity: noise.lacunarity, gain: noise.gain)
      when T is Perlin:
        result.kind = avxPerlin
      else:
        result.settings.seed = cast[int32](noise.seed)
        when T is OpenSimplex2: result.kind = avxOpenSimplex2
        elif T is OpenSimplex2S: result.kind = avxOpenSimplex2S
        elif T is ValueNoise: result.kind = avxValue
        elif T is ValueCubic: result.kind = avxValueCubic
        elif T is Cellular:
          result.kind = avxCellular
          result.settings.jitter = noise.jitter
          result.settings.distanceFunction = ord(noise.distanceFunction)
          result.settings.returnType = ord(noise.returnType)

  template sampleBatch(noise, batch, points, dimensions, output: untyped) =
    when noise is Simplex:
      sampleTileAvx2(noise, points, dimensions, output)
    else:
      sampleTileAvx2(batch, points, dimensions, output)

func validateTile[T: TileNoise](
  noise: T, start, period: openArray[float64]
) =
  ## Validate before writes, including the integer domain of every octave.
  if noise.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  for value in start:
    if not (abs(value) < Inf):
      raise newException(NoisyError, "Tile origins must be finite")
  var radiusSum = 0.0
  for value in period:
    if not (value > 0 and value < Inf):
      raise newException(NoisyError, "Tile periods must be finite and > 0")
    radiusSum += value / (2 * PI)
  # Circle coordinates are stored as float32 before octave scaling.
  if not (radiusSum < 3.4028234e38):
    raise newException(NoisyError, "Tile periods are too large")
  for value in [noise.amplitude, noise.frequency, noise.lacunarity, noise.gain]:
    if not (abs(value) < Inf):
      raise newException(NoisyError, "Tile noise settings must be finite")
  when T is Cellular:
    if not (noise.jitter >= 0 and noise.jitter <= 1):
      raise newException(NoisyError, "Cellular jitter must be in [0, 1]")
  var frequency = noise.frequency
  for octave in 0 ..< noise.octaves:
    # Allow margin for the simplex skew and all neighboring lattice cells.
    if not (abs(frequency.float64) * radiusSum <= 268435456.0):
      raise newException(NoisyError, "Tiled octave coordinates are too large")
    if octave + 1 < noise.octaves:
      frequency = frequency * noise.lacunarity

func circle(coordinate, period: float64): array[2, float32] {.inline, raises: [].} =
  ## Reduce phase first so representable period shifts use the same trig input.
  var phase = coordinate mod period
  if phase < 0:
    phase += period
  if phase == 0 or phase >= period:
    phase = 0
  let
    angle = (phase / period) * (2 * PI)
    radius = period / (2 * PI)
  [(cos(angle) * radius).float32, (sin(angle) * radius).float32]

func fillTile[T: TileNoise](
  noise: T, start, period: float64, width: int, buffer: var openArray[float32]
) =
  mixin value
  var scalarStart = 0
  when nimvm:
    discard
  else:
    when useGridAvx2:
      if hasAvx2():
        var batch = tileBatch(noise)
        when T is Perlin:
          batch.settings.perm = noise.perm.unsafeAddr
        scalarStart = width - width mod 8
        var points: array[4, array[8, float32]]
        for x in countup(0, scalarStart - 8, 8):
          for lane in 0 ..< 8:
            let point = circle(start + float64(x + lane), period)
            points[0][lane] = point[0]
            points[1][lane] = point[1]
          sampleBatch(noise, batch, points, 2, buffer[x].addr)
  for x in scalarStart ..< width:
    let point = circle(start + x.float64, period)
    buffer[x] = noise.value(point[0], point[1])

func fillTile[T: TileNoise2D](
  noise: T, start, period: (float64, float64), dimens: (int, int),
  buffer: var openArray[float32]
) =
  mixin value
  let (width, height) = dimens
  var scalarStart = 0
  when useGridAvx2:
    var batch: AvxBatch
  when nimvm:
    discard
  else:
    when useGridAvx2:
      if hasAvx2():
        scalarStart = height - height mod 8
        batch = tileBatch(noise)
        when T is Perlin:
          batch.settings.perm = noise.perm.unsafeAddr
  for x in 0 ..< width:
    let xPoint = circle(start[0] + x.float64, period[0])
    when nimvm:
      discard
    else:
      when useGridAvx2:
        var points: array[4, array[8, float32]]
        for lane in 0 ..< 8:
          points[0][lane] = xPoint[0]
          points[1][lane] = xPoint[1]
        for y in countup(0, scalarStart - 8, 8):
          for lane in 0 ..< 8:
            let yPoint = circle(start[1] + float64(y + lane), period[1])
            points[2][lane] = yPoint[0]
            points[3][lane] = yPoint[1]
          sampleBatch(noise, batch, points, 4, buffer[x * height + y].addr)
    for y in scalarStart ..< height:
      let yPoint = circle(start[1] + y.float64, period[1])
      buffer[x * height + y] = noise.value(
        xPoint[0], xPoint[1], yPoint[0], yPoint[1])

func tileGrid*[T: TileNoise; S, P: SomeNumber](
  noise: T, start: S, width: int, period: P, buffer: var openArray[float32]
) =
  ## Fills a seamless 1D loop using a circle through the 2D noise sampler.
  ## Period is in world units; samples advance by one and exclude the endpoint.
  ## Extra buffer elements are untouched. Invalid input raises before writes.
  let length = checkGridBuffer(buffer.len, [width])
  validateTile(noise, [start.float64], [period.float64])
  if length > 0:
    fillTile(noise, start.float64, period.float64, width, buffer)

func tileGrid*[T: TileNoise; S, P: SomeNumber](
  noise: T, start: S, width: int, period: P
): Grid =
  ## Allocates a 1D tiled grid and calls the buffer-taking overload.
  validateTile(noise, [start.float64], [period.float64])
  result = allocateGrid(width)
  noise.tileGrid(start, width, period, result.values)

func tileGrid*[T: TileNoise; S: SomeNumber](
  noise: T, start: S, width: int, buffer: var openArray[float32]
) =
  ## Fills one complete 1D repeat; the period defaults to width (one if empty).
  noise.tileGrid(start, width, max(width, 1), buffer)

func tileGrid*[T: TileNoise; S: SomeNumber](noise: T, start: S, width: int): Grid =
  ## Allocates one complete 1D repeat; the period defaults to width.
  noise.tileGrid(start, width, max(width, 1))

func tileGrid*[T: TileNoise2D; S, P: SomeNumber](
  noise: T, start: (S, S), dimens: (int, int), period: (P, P),
  buffer: var openArray[float32]
) =
  ## Fills a seamless 2D tile by mapping both axes onto a torus in 4D noise.
  ## Simplex and Perlin support this mapping. Periods are in world units.
  ## Extra buffer elements are untouched. Invalid input raises before writes.
  let length = checkGridBuffer(buffer.len, [dimens[0], dimens[1]])
  validateTile(noise, [start[0].float64, start[1].float64],
    [period[0].float64, period[1].float64])
  if length > 0:
    fillTile(noise, (start[0].float64, start[1].float64),
      (period[0].float64, period[1].float64), dimens, buffer)

func tileGrid*[T: TileNoise2D; S, P: SomeNumber](
  noise: T, start: (S, S), dimens: (int, int), period: (P, P)
): Grid =
  ## Allocates a 2D tiled grid and calls the buffer-taking overload.
  validateTile(noise, [start[0].float64, start[1].float64],
    [period[0].float64, period[1].float64])
  result = allocateGrid(dimens)
  noise.tileGrid(start, dimens, period, result.values)

func tileGrid*[T: TileNoise2D; S: SomeNumber](
  noise: T, start: (S, S), dimens: (int, int), buffer: var openArray[float32]
) =
  ## Fills one complete 2D repeat; periods default to the grid dimensions.
  noise.tileGrid(start, dimens, (max(dimens[0], 1), max(dimens[1], 1)), buffer)

func tileGrid*[T: TileNoise2D; S: SomeNumber](
  noise: T, start: (S, S), dimens: (int, int)
): Grid =
  ## Allocates one complete 2D repeat; periods default to the dimensions.
  noise.tileGrid(start, dimens, (max(dimens[0], 1), max(dimens[1], 1)))
