import common, lattices

export common

type
  ValueNoise* = object
    ## Scalar 2D and 3D Value noise using FastNoiseLite's lattice conventions.
    seed*, octaves*: int
    amplitude*, frequency*, lacunarity*, gain*: float32

func initValue*(seed: int): ValueNoise =
  ## Initializes Value noise with the same octave defaults as Simplex.
  result.seed = seed
  result.octaves = 1
  result.amplitude = 1
  result.frequency = 1
  result.lacunarity = 2
  result.gain = 0.5

func sample(seed: int32, x, y: float32): float32 =
  ## Interpolates hashed corner values using FastNoiseLite's cubic fade.
  let
    ix = fastFloor(x)
    iy = fastFloor(y)
    xs = hermite(x - ix.float32)
    ys = hermite(y - iy.float32)
    x0 = prime(ix, PrimeX)
    y0 = prime(iy, PrimeY)
    x1 = x0 + PrimeX
    y1 = y0 + PrimeY
    row0 = lerp(latticeValue(seed, x0, y0), latticeValue(seed, x1, y0), xs)
    row1 = lerp(latticeValue(seed, x0, y1), latticeValue(seed, x1, y1), xs)
  lerp(row0, row1, ys)

func sample(seed: int32, x, y, z: float32): float32 =
  ## Interpolates hashed corner values in native three-dimensional space.
  let
    ix = fastFloor(x)
    iy = fastFloor(y)
    iz = fastFloor(z)
    xs = hermite(x - ix.float32)
    ys = hermite(y - iy.float32)
    zs = hermite(z - iz.float32)
    x0 = prime(ix, PrimeX)
    y0 = prime(iy, PrimeY)
    z0 = prime(iz, PrimeZ)
    x1 = x0 + PrimeX
    y1 = y0 + PrimeY
    z1 = z0 + PrimeZ
    row00 = lerp(
      latticeValue(seed, x0, y0, z0), latticeValue(seed, x1, y0, z0), xs
    )
    row10 = lerp(
      latticeValue(seed, x0, y1, z0), latticeValue(seed, x1, y1, z0), xs
    )
    row01 = lerp(
      latticeValue(seed, x0, y0, z1), latticeValue(seed, x1, y0, z1), xs
    )
    row11 = lerp(
      latticeValue(seed, x0, y1, z1), latticeValue(seed, x1, y1, z1), xs
    )
  lerp(lerp(row00, row10, ys), lerp(row01, row11, ys), zs)

func value*(noise: ValueNoise, x, y: float32): float32 =
  ## Samples native 2D Value noise with the configured octaves.
  if noise.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  var
    amplitude = noise.amplitude
    frequency = noise.frequency
  for octave in 0 ..< noise.octaves:
    result += sample(cast[int32](noise.seed), x * frequency, y * frequency) *
      amplitude
    amplitude *= noise.gain
    frequency *= noise.lacunarity
  result /= noise.octaves.float32

func value*(noise: ValueNoise, x, y, z: float32): float32 =
  ## Samples native 3D Value noise with the configured octaves.
  if noise.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  var
    amplitude = noise.amplitude
    frequency = noise.frequency
  for octave in 0 ..< noise.octaves:
    result += sample(
      cast[int32](noise.seed), x * frequency, y * frequency, z * frequency
    ) * amplitude
    amplitude *= noise.gain
    frequency *= noise.lacunarity
  result /= noise.octaves.float32

func value*(noise: ValueNoise, x, y: int): float32 {.inline.} =
  ## Samples 2D Value noise at integer coordinates.
  noise.value(x.float32, y.float32)

func value*(noise: ValueNoise, x, y, z: int): float32 {.inline.} =
  ## Samples 3D Value noise at integer coordinates.
  noise.value(x.float32, y.float32, z.float32)
