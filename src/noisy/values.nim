import common, grids, lattices

when useGridNeon:
  import simd, simd/[valuecubicneon, valueneon]

export common, Grid, `[]`

type
  ValueNoise* = object
    ## 2D and 3D Value noise using FastNoiseLite's lattice conventions.
    seed*, octaves*: int
    amplitude*, frequency*, lacunarity*, gain*: float32

  ValueCubic* = object
    ## 2D and 3D cubic Value noise matching FastNoiseLite.
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

include kernelchecks

func sample(seed: int32, x, y: float32): float32 {.noiseKernel.} =
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

func sample(seed: int32, x, y, z: float32): float32 {.noiseKernel.} =
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

func value*(noise: ValueNoise, x, y: float32): float32 {.raises: [NoisyError].} =
  ## Samples native 2D Value noise with the configured octaves.
  if noise.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  var
    amplitude = noise.amplitude
    frequency = noise.frequency
  for octave in 0 ..< noise.octaves:
    result = result + sample(cast[int32](noise.seed), x * frequency, y * frequency) *
      amplitude
    amplitude = amplitude * noise.gain
    frequency = frequency * noise.lacunarity
  result = result / noise.octaves.float32

func value*(noise: ValueNoise, x, y, z: float32): float32 {.raises: [NoisyError].} =
  ## Samples native 3D Value noise with the configured octaves.
  if noise.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  var
    amplitude = noise.amplitude
    frequency = noise.frequency
  for octave in 0 ..< noise.octaves:
    result = result + sample(
      cast[int32](noise.seed), x * frequency, y * frequency, z * frequency
    ) * amplitude
    amplitude = amplitude * noise.gain
    frequency = frequency * noise.lacunarity
  result = result / noise.octaves.float32

func value*(noise: ValueNoise, x, y: int): float32 {.inline.} =
  ## Samples 2D Value noise at integer coordinates.
  noise.value(x.float32, y.float32)

func value*(noise: ValueNoise, x, y, z: int): float32 {.inline.} =
  ## Samples 3D Value noise at integer coordinates.
  noise.value(x.float32, y.float32, z.float32)

func cubic(a, b, c, d, t: float32): float32 {.noiseKernel, inline.} =
  ## Interpolates four samples in FastNoiseLite's arithmetic order.
  let p = (d - c) - (a - b)
  t * t * t * p + t * t * ((a - b) - p) + t * (c - a) + b

func sampleCubic(seed: int32, x, y: float32): float32 {.noiseKernel.} =
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

func sampleCubic(seed: int32, x, y, z: float32): float32 {.noiseKernel.} =
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

func value*(noise: ValueCubic, x, y: float32): float32 {.raises: [NoisyError].} =
  ## Samples native 2D cubic Value noise with the configured octaves.
  if noise.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  var
    amplitude = noise.amplitude
    frequency = noise.frequency
  for octave in 0 ..< noise.octaves:
    result = result + sampleCubic(
      cast[int32](noise.seed),
      x * frequency,
      y * frequency
    ) * amplitude
    amplitude = amplitude * noise.gain
    frequency = frequency * noise.lacunarity
  result = result / noise.octaves.float32

func value*(noise: ValueCubic, x, y, z: float32): float32 {.raises: [NoisyError].} =
  ## Samples native 3D cubic Value noise with the configured octaves.
  if noise.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  var
    amplitude = noise.amplitude
    frequency = noise.frequency
  for octave in 0 ..< noise.octaves:
    result = result + sampleCubic(
      cast[int32](noise.seed),
      x * frequency,
      y * frequency,
      z * frequency
    ) * amplitude
    amplitude = amplitude * noise.gain
    frequency = frequency * noise.lacunarity
  result = result / noise.octaves.float32

func value*(noise: ValueCubic, x, y: int): float32 {.inline.} =
  ## Samples 2D cubic Value noise at integer coordinates.
  noise.value(x.float32, y.float32)

func value*(noise: ValueCubic, x, y, z: int): float32 {.inline.} =
  ## Samples 3D cubic Value noise at integer coordinates.
  noise.value(x.float32, y.float32, z.float32)

when useGridNeon:
  func sample2x4[T: ValueNoise | ValueCubic](noise: T, x, y: Float4): Float4 =
    ## Samples four 2D points with the configured octaves.
    let seed = cast[int32](noise.seed)
    var
      amplitude = noise.amplitude
      frequency = noise.frequency
      total = float4(0)
    for _ in 0 ..< noise.octaves:
      let scale = float4(frequency)
      let sample =
        when T is ValueNoise: sampleValue4(seed, x * scale, y * scale)
        else: sampleCubic4(seed, x * scale, y * scale)
      total = total + sample * float4(amplitude)
      amplitude = amplitude * noise.gain
      frequency = frequency * noise.lacunarity
    total / float4(noise.octaves.float32)

  func sample3x4[T: ValueNoise | ValueCubic](noise: T, x, y, z: Float4): Float4 =
    ## Samples four 3D points with the configured octaves.
    let seed = cast[int32](noise.seed)
    var
      amplitude = noise.amplitude
      frequency = noise.frequency
      total = float4(0)
    for _ in 0 ..< noise.octaves:
      let scale = float4(frequency)
      let sample =
        when T is ValueNoise: sampleValue4(seed, x * scale, y * scale, z * scale)
        else: sampleCubic4(seed, x * scale, y * scale, z * scale)
      total = total + sample * float4(amplitude)
      amplitude = amplitude * noise.gain
      frequency = frequency * noise.lacunarity
    total / float4(noise.octaves.float32)

{.pop.}

func grid*(
  noise: ValueNoise, start: (float32, float32), dimens: (int, int),
  buffer: var openArray[float32]
) =
  ## Fills 2D samples into buffer without allocating or resizing it.
  ## Requires room for the full grid; extra elements are left untouched.
  when useGridAvx2:
    noise.avxGrid(start, dimens, buffer, avxValue)
  elif useGridNeon:
    noise.fillGrid(start, dimens, buffer, sample2x4[ValueNoise])
  else:
    noise.scalarGrid(start, dimens, buffer)

func grid*(
  noise: ValueNoise, start: (float32, float32), dimens: (int, int)
): Grid =
  ## Allocates and fills a 2D grid at unit intervals from start.
  if noise.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  result = allocateGrid(dimens)
  noise.grid(start, dimens, result.values)

func grid*(
  noise: ValueNoise, start: (int, int), dimens: (int, int),
  buffer: var openArray[float32]
) {.inline.} =
  ## Fills a 2D grid into caller-owned storage from an integer start.
  noise.grid((start[0].float32, start[1].float32), dimens, buffer)

func grid*(
  noise: ValueNoise, start: (int, int), dimens: (int, int)
): Grid {.inline.} =
  ## Allocates a 2D grid from an integer start.
  noise.grid((start[0].float32, start[1].float32), dimens)

func grid*(
  noise: ValueNoise, start: (float32, float32, float32), dimens: (int, int, int),
  buffer: var openArray[float32]
) =
  ## Fills 3D samples into buffer without allocating or resizing it.
  ## Requires room for the full grid; extra elements are left untouched.
  when useGridAvx2:
    noise.avxGrid(start, dimens, buffer, avxValue)
  elif useGridNeon:
    noise.fillGrid(start, dimens, buffer, sample3x4[ValueNoise])
  else:
    noise.scalarGrid(start, dimens, buffer)

func grid*(
  noise: ValueNoise, start: (float32, float32, float32), dimens: (int, int, int)
): Grid =
  ## Allocates and fills a 3D grid at unit intervals from start.
  if noise.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  result = allocateGrid(dimens)
  noise.grid(start, dimens, result.values)

func grid*(
  noise: ValueNoise, start: (int, int, int), dimens: (int, int, int),
  buffer: var openArray[float32]
) {.inline.} =
  ## Fills a 3D grid into caller-owned storage from an integer start.
  noise.grid((start[0].float32, start[1].float32, start[2].float32), dimens, buffer)

func grid*(
  noise: ValueNoise, start: (int, int, int), dimens: (int, int, int)
): Grid {.inline.} =
  ## Allocates a 3D grid from an integer start.
  noise.grid((start[0].float32, start[1].float32, start[2].float32), dimens)

func grid*(
  noise: ValueCubic, start: (float32, float32), dimens: (int, int),
  buffer: var openArray[float32]
) =
  ## Fills 2D samples into buffer without allocating or resizing it.
  ## Requires room for the full grid; extra elements are left untouched.
  when useGridAvx2:
    noise.avxGrid(start, dimens, buffer, avxValueCubic)
  elif useGridNeon:
    noise.fillGrid(start, dimens, buffer, sample2x4[ValueCubic])
  else:
    noise.scalarGrid(start, dimens, buffer)

func grid*(
  noise: ValueCubic, start: (float32, float32), dimens: (int, int)
): Grid =
  ## Allocates and fills a 2D grid at unit intervals from start.
  if noise.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  result = allocateGrid(dimens)
  noise.grid(start, dimens, result.values)

func grid*(
  noise: ValueCubic, start: (int, int), dimens: (int, int),
  buffer: var openArray[float32]
) {.inline.} =
  ## Fills a 2D grid into caller-owned storage from an integer start.
  noise.grid((start[0].float32, start[1].float32), dimens, buffer)

func grid*(
  noise: ValueCubic, start: (int, int), dimens: (int, int)
): Grid {.inline.} =
  ## Allocates a 2D grid from an integer start.
  noise.grid((start[0].float32, start[1].float32), dimens)

func grid*(
  noise: ValueCubic, start: (float32, float32, float32), dimens: (int, int, int),
  buffer: var openArray[float32]
) =
  ## Fills 3D samples into buffer without allocating or resizing it.
  ## Requires room for the full grid; extra elements are left untouched.
  when useGridAvx2:
    noise.avxGrid(start, dimens, buffer, avxValueCubic)
  elif useGridNeon:
    noise.fillGrid(start, dimens, buffer, sample3x4[ValueCubic])
  else:
    noise.scalarGrid(start, dimens, buffer)

func grid*(
  noise: ValueCubic, start: (float32, float32, float32), dimens: (int, int, int)
): Grid =
  ## Allocates and fills a 3D grid at unit intervals from start.
  if noise.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  result = allocateGrid(dimens)
  noise.grid(start, dimens, result.values)

func grid*(
  noise: ValueCubic, start: (int, int, int), dimens: (int, int, int),
  buffer: var openArray[float32]
) {.inline.} =
  ## Fills a 3D grid into caller-owned storage from an integer start.
  noise.grid((start[0].float32, start[1].float32, start[2].float32), dimens, buffer)

func grid*(
  noise: ValueCubic, start: (int, int, int), dimens: (int, int, int)
): Grid {.inline.} =
  ## Allocates a 3D grid from an integer start.
  noise.grid((start[0].float32, start[1].float32, start[2].float32), dimens)
