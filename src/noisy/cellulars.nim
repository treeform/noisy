import
  std/math,
  common, grids, lattices

when useGridNeon:
  import simd, simd/cellularneon

export common, Grid, `[]`

type
  CellularDistance* = enum
    Euclidean, EuclideanSquared, Manhattan, Hybrid

  CellularReturn* = enum
    CellValue, Distance, Distance2, Distance2Add, Distance2Sub,
    Distance2Mul, Distance2Div

  Cellular* = object
    ## 2D/3D cellular noise using FastNoiseLite's conventions.
    seed*, octaves*: int
    amplitude*, frequency*, lacunarity*, gain*: float32
    distanceFunction*: CellularDistance
    returnType*: CellularReturn
    jitter*: float32

func initCellular*(seed: int): Cellular =
  ## Initializes cellular noise with squared distance and unit jitter.
  result.seed = seed
  result.octaves = 1
  result.amplitude = 1
  result.frequency = 1
  result.lacunarity = 2
  result.gain = 0.5
  result.distanceFunction = EuclideanSquared
  result.returnType = Distance
  result.jitter = 1

func metric[N: static[int]](
  distanceFunction: CellularDistance,
  offset: array[N, float32]
): float32 {.inline.} =
  ## Measures a feature point using the selected distance convention.
  var squared, manhattan: float32
  for axis in 0 ..< N:
    squared += offset[axis] * offset[axis]
    manhattan += abs(offset[axis])
  case distanceFunction
  of Euclidean, EuclideanSquared: squared
  of Manhattan: manhattan
  of Hybrid: manhattan + squared

func finish(
  cellular: Cellular,
  closestHash: uint32,
  nearest, second: float32
): float32 =
  ## Converts the two nearest feature distances to the selected return mode.
  var
    distance0 = nearest
    distance1 = second
  if cellular.distanceFunction == Euclidean and
    cellular.returnType >= Distance:
      distance0 = sqrt(distance0)
      if cellular.returnType >= Distance2:
        distance1 = sqrt(distance1)
  case cellular.returnType
  of CellValue:
    cast[int32](closestHash).float32 * (1.0'f / 2147483648.0'f)
  of Distance:
    distance0 - 1
  of Distance2:
    distance1 - 1
  of Distance2Add:
    (distance1 + distance0) * 0.5'f - 1
  of Distance2Sub:
    distance1 - distance0 - 1
  of Distance2Mul:
    distance1 * distance0 * 0.5'f - 1
  of Distance2Div:
    distance0 / distance1 - 1

func noise[N: static[int]](
  cellular: Cellular,
  point: array[N, float32]
): float32 =
  ## Searches neighboring feature points in FastNoiseLite's traversal order.
  const
    Primes = [PrimeX, PrimeY, PrimeZ]
    Cells = when N == 2: 9 else: 27
    JitterScale = when N == 2: 0.43701595'f else: 0.39614353'f
  var
    center: array[N, int32]
    distance0 = 1.0e10'f
    distance1 = 1.0e10'f
    closestHash: uint32
  let
    seed = cast[int32](cellular.seed)
    jitter = cellular.jitter * JitterScale
  for axis in 0 ..< N:
    center[axis] = fastRound(point[axis])
  for candidate in 0 ..< Cells:
    var
      ordinal = candidate
      coordinates: array[N, int32]
      primed: array[N, uint32]
      offset: array[N, float32]
    for axis in countdown(N - 1, 0):
      coordinates[axis] = center[axis] + (ordinal mod 3).int32 - 1
      primed[axis] = prime(coordinates[axis], Primes[axis])
      ordinal = ordinal div 3
    let hashValue =
      when N == 2: hash(seed, primed[0], primed[1])
      else: hash(seed, primed[0], primed[1], primed[2])
    when N == 2:
      let index = (hashValue and (255'u32 shl 1)).int
      for axis in 0 ..< N:
        offset[axis] = (coordinates[axis].float32 - point[axis]) +
          RandVecs2D[index or axis] * jitter
    else:
      let index = (hashValue and (255'u32 shl 2)).int
      for axis in 0 ..< N:
        offset[axis] = (coordinates[axis].float32 - point[axis]) +
          RandVecs3D[index or axis] * jitter
    let distance = metric(cellular.distanceFunction, offset)
    distance1 = max(min(distance1, distance), distance0)
    if distance < distance0:
      distance0 = distance
      closestHash = hashValue
  cellular.finish(closestHash, distance0, distance1)

func validate(cellular: Cellular) {.inline, raises: [NoisyError].} =
  ## Checks the shared requirements for point and grid sampling.
  if cellular.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  if not (cellular.jitter >= 0 and cellular.jitter <= 1):
    raise newException(NoisyError, "Cellular jitter must be in [0, 1]")

func sample[N: static[int]](
  cellular: Cellular,
  point: array[N, float32]
): float32 {.raises: [NoisyError].} =
  ## Layers scalar samples using the same octave convention as Simplex.
  cellular.validate()
  var
    amplitude = cellular.amplitude
    frequency = cellular.frequency
  for octave in 0 ..< cellular.octaves:
    var scaled: array[N, float32]
    for axis in 0 ..< N:
      scaled[axis] = point[axis] * frequency
    result += cellular.noise(scaled) * amplitude
    amplitude *= cellular.gain
    frequency *= cellular.lacunarity
  result /= cellular.octaves.float32

func value*(cellular: Cellular, x, y: float32): float32 =
  ## Samples native 2D cellular noise with the configured octaves.
  cellular.sample([x, y])

func value*(cellular: Cellular, x, y, z: float32): float32 =
  ## Samples native 3D cellular noise with the configured octaves.
  cellular.sample([x, y, z])

func value*(cellular: Cellular, x, y: int): float32 {.inline.} =
  ## Samples 2D cellular noise at integer coordinates.
  cellular.value(x.float32, y.float32)

func value*(cellular: Cellular, x, y, z: int): float32 {.inline.} =
  ## Samples 3D cellular noise at integer coordinates.
  cellular.value(x.float32, y.float32, z.float32)

when useGridNeon:
  func finish4(cellular: Cellular, features: NearestFeatures): Float4 {.inline.} =
    ## Converts the nearest feature distances using the scalar return rules.
    let one = float4(1)
    var
      distance0 = features.first
      distance1 = features.second
    if cellular.distanceFunction == Euclidean and cellular.returnType >= Distance:
      distance0 = sqrt4(distance0)
      if cellular.returnType >= Distance2:
        distance1 = sqrt4(distance1)
    case cellular.returnType
    of CellValue:
      toFloats(cast[Int4](features.hash)) * float4(1.0'f / 2147483648.0'f)
    of Distance:
      distance0 - one
    of Distance2:
      distance1 - one
    of Distance2Add:
      (distance1 + distance0) * float4(0.5'f) - one
    of Distance2Sub:
      distance1 - distance0 - one
    of Distance2Mul:
      distance1 * distance0 * float4(0.5'f) - one
    of Distance2Div:
      distance0 / distance1 - one

  func sample4Impl[N: static[int]](
    cellular: Cellular, point: array[N, Float4],
    squared, manhattan, trackSecond, trackHash: static[bool]
  ): Float4 =
    ## Layers four independent samples with the same seed in each octave.
    let seed = cast[int32](cellular.seed)
    var
      amplitude = cellular.amplitude
      frequency = cellular.frequency
      total = float4(0)
    for _ in 0 ..< cellular.octaves:
      let scale = float4(frequency)
      var scaled: array[N, Float4]
      for axis in 0 ..< N:
        scaled[axis] = point[axis] * scale
      let features = searchFeatures4(seed, cellular.jitter, scaled,
        squared, manhattan, trackSecond, trackHash)
      total += cellular.finish4(features) * float4(amplitude)
      amplitude *= cellular.gain
      frequency *= cellular.lacunarity
    total / float4(cellular.octaves.float32)

  func sample4[N: static[int]](
    cellular: Cellular, point: array[N, Float4]
  ): Float4 =
    ## Selects the required search calculations before visiting any features.
    template sampleMetric(trackSecond, trackHash: static[bool]): Float4 =
      case cellular.distanceFunction
      of Euclidean, EuclideanSquared:
        sample4Impl(cellular, point, true, false, trackSecond, trackHash)
      of Manhattan:
        sample4Impl(cellular, point, false, true, trackSecond, trackHash)
      of Hybrid:
        sample4Impl(cellular, point, true, true, trackSecond, trackHash)
    case cellular.returnType
    of CellValue: sampleMetric(false, true)
    of Distance: sampleMetric(false, false)
    else: sampleMetric(true, false)

  func sample2x4(cellular: Cellular, x, y: Float4): Float4 =
    cellular.sample4([x, y])

  func sample3x4(cellular: Cellular, x, y, z: Float4): Float4 =
    cellular.sample4([x, y, z])

func grid*(
  cellular: Cellular, start: (float32, float32), dimens: (int, int)
): Grid =
  ## Generates a 2D Cellular grid at unit intervals.
  cellular.validate()
  when useGridNeon:
    cellular.fillGrid(start, dimens, sample2x4)
  else:
    cellular.scalarGrid(start, dimens)

func grid*(
  cellular: Cellular,
  start: (float32, float32, float32),
  dimens: (int, int, int)
): Grid =
  ## Generates a 3D Cellular grid at unit intervals.
  cellular.validate()
  when useGridNeon:
    cellular.fillGrid(start, dimens, sample3x4)
  else:
    cellular.scalarGrid(start, dimens)

func grid*(
  cellular: Cellular, start: (int, int), dimens: (int, int)
): Grid {.inline.} =
  ## Generates a 2D Cellular grid from an integer origin.
  cellular.grid((start[0].float32, start[1].float32), dimens)

func grid*(
  cellular: Cellular, start: (int, int, int), dimens: (int, int, int)
): Grid {.inline.} =
  ## Generates a 3D Cellular grid from an integer origin.
  cellular.grid((start[0].float32, start[1].float32, start[2].float32), dimens)
