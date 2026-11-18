import common, lattices

export common

type
  ValueNoise* = object
    ## Scalar 2D and 3D Value noise using FastNoiseLite's lattice conventions.
    seed*, octaves*: int
    amplitude*, frequency*, lacunarity*, gain*: float32

  ValueCubic* = object
    ## Scalar 2D and 3D cubic Value noise matching FastNoiseLite.
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

func initValueCubic*(seed: int): ValueCubic =
  ## Initializes cubic Value noise with the standard seeded octave settings.
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

func cubic(a, b, c, d, t: float32): float32 {.inline.} =
  ## Interpolates four samples in FastNoiseLite's arithmetic order.
  let p = (d - c) - (a - b)
  t * t * t * p + t * t * ((a - b) - p) + t * (c - a) + b

func sampleCubic(seed: int32, x, y: float32): float32 =
  ## Interpolates sixteen lattice values, first along x and then y.
  let
    ix = fastFloor(x)
    iy = fastFloor(y)
    xs = x - ix.float32
    ys = y - iy.float32
    x1 = prime(ix, PrimeX)
    y1 = prime(iy, PrimeY)
    columns = [x1 - PrimeX, x1, x1 + PrimeX, x1 + 2'u32 * PrimeX]
    rows = [y1 - PrimeY, y1, y1 + PrimeY, y1 + 2'u32 * PrimeY]
  var interpolated: array[4, float32]
  for row in 0 ..< 4:
    interpolated[row] = cubic(
      latticeValue(seed, columns[0], rows[row]),
      latticeValue(seed, columns[1], rows[row]),
      latticeValue(seed, columns[2], rows[row]),
      latticeValue(seed, columns[3], rows[row]),
      xs
    )
  cubic(
    interpolated[0],
    interpolated[1],
    interpolated[2],
    interpolated[3],
    ys
  ) * (1.0'f / (1.5'f * 1.5'f))

func sampleCubic(seed: int32, x, y, z: float32): float32 =
  ## Interpolates sixty-four lattice values in x, y, then z order.
  let
    ix = fastFloor(x)
    iy = fastFloor(y)
    iz = fastFloor(z)
    xs = x - ix.float32
    ys = y - iy.float32
    zs = z - iz.float32
    x1 = prime(ix, PrimeX)
    y1 = prime(iy, PrimeY)
    z1 = prime(iz, PrimeZ)
    columns = [x1 - PrimeX, x1, x1 + PrimeX, x1 + 2'u32 * PrimeX]
    rows = [y1 - PrimeY, y1, y1 + PrimeY, y1 + 2'u32 * PrimeY]
    layers = [z1 - PrimeZ, z1, z1 + PrimeZ, z1 + 2'u32 * PrimeZ]
  var interpolated: array[4, float32]
  for layer in 0 ..< 4:
    var values: array[4, float32]
    for row in 0 ..< 4:
      values[row] = cubic(
        latticeValue(seed, columns[0], rows[row], layers[layer]),
        latticeValue(seed, columns[1], rows[row], layers[layer]),
        latticeValue(seed, columns[2], rows[row], layers[layer]),
        latticeValue(seed, columns[3], rows[row], layers[layer]),
        xs
      )
    interpolated[layer] = cubic(values[0], values[1], values[2], values[3], ys)
  cubic(
    interpolated[0],
    interpolated[1],
    interpolated[2],
    interpolated[3],
    zs
  ) * (1.0'f / (1.5'f * 1.5'f * 1.5'f))

func value*(noise: ValueCubic, x, y: float32): float32 =
  ## Samples native 2D cubic Value noise with the configured octaves.
  if noise.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  var
    amplitude = noise.amplitude
    frequency = noise.frequency
  for octave in 0 ..< noise.octaves:
    result += sampleCubic(
      cast[int32](noise.seed),
      x * frequency,
      y * frequency
    ) * amplitude
    amplitude *= noise.gain
    frequency *= noise.lacunarity
  result /= noise.octaves.float32

func value*(noise: ValueCubic, x, y, z: float32): float32 =
  ## Samples native 3D cubic Value noise with the configured octaves.
  if noise.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  var
    amplitude = noise.amplitude
    frequency = noise.frequency
  for octave in 0 ..< noise.octaves:
    result += sampleCubic(
      cast[int32](noise.seed),
      x * frequency,
      y * frequency,
      z * frequency
    ) * amplitude
    amplitude *= noise.gain
    frequency *= noise.lacunarity
  result /= noise.octaves.float32

func value*(noise: ValueCubic, x, y: int): float32 {.inline.} =
  ## Samples 2D cubic Value noise at integer coordinates.
  noise.value(x.float32, y.float32)

func value*(noise: ValueCubic, x, y, z: int): float32 {.inline.} =
  ## Samples 3D cubic Value noise at integer coordinates.
  noise.value(x.float32, y.float32, z.float32)
