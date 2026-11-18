import
  std/math,
  common, lattices

export common

type
  CellularDistance* = enum
    Euclidean, EuclideanSquared, Manhattan, Hybrid

  CellularReturn* = enum
    CellValue, Distance, Distance2, Distance2Add, Distance2Sub,
    Distance2Mul, Distance2Div

  Cellular* = object
    ## Scalar 2D/3D cellular noise using FastNoiseLite's conventions.
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

func sample[N: static[int]](
  cellular: Cellular,
  point: array[N, float32]
): float32 {.raises: [NoisyError].} =
  ## Layers scalar samples using the same octave convention as Simplex.
  if cellular.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  if not (cellular.jitter >= 0 and cellular.jitter <= 1):
    raise newException(NoisyError, "Cellular jitter must be in [0, 1]")
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
