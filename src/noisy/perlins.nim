import
  std/[math, random],
  common, grids

when useGridNeon:
  import simd, simd/perlinneon

export common, Grid, `[]`

const
  PerlinScales = [0.188'f, 0.507'f, 0.936'f, 0.87'f]

type
  Perlin* = object
    ## Improved Perlin noise with native 1D through 4D sampling.
    octaves*: int
    amplitude*, frequency*, lacunarity*, gain*: float32
    perm*: array[512, int32]

func initPerlin*(seed: int): Perlin =
  ## Initializes Improved Perlin noise with a seeded permutation table.
  result.octaves = 1
  result.amplitude = 1
  result.frequency = 1
  result.lacunarity = 2
  result.gain = 0.5
  var
    permutation: array[256, int32]
    rng = initRand(seed)
  for i in 0 ..< permutation.len:
    permutation[i] = i.int32
  shuffle(rng, permutation)
  for i in 0 ..< result.perm.len:
    result.perm[i] = permutation[i and 255]

func fade(t: float32): float32 {.inline.} =
  ## Uses Perlin's quintic interpolant with continuous second derivatives.
  t * t * t * (t * (t * 6 - 15) + 10)

func interpolate(t, a, b: float32): float32 {.inline.} =
  ## Interpolates the two lattice contributions.
  a + t * (b - a)

func gradient[N: static[int]](
  hash: int,
  offset: array[N, float32]
): float32 {.inline.} =
  ## Evaluates Gustavson's dimension-specific Improved Perlin gradients.
  when N == 1:
    let
      h = hash and 15
      magnitude = (1 + (h and 7)).float32
      signed = if (h and 8) == 0: magnitude else: -magnitude
    signed * offset[0]
  elif N == 2:
    let
      h = hash and 7
      u = if h < 4: offset[0] else: offset[1]
      v = if h < 4: offset[1] else: offset[0]
    (if (h and 1) == 0: u else: -u) +
      (if (h and 2) == 0: 2 * v else: -2 * v)
  elif N == 3:
    let
      h = hash and 15
      u = if h < 8: offset[0] else: offset[1]
      v =
        if h < 4:
          offset[1]
        elif h == 12 or h == 14:
          offset[0]
        else:
          offset[2]
    (if (h and 1) == 0: u else: -u) +
      (if (h and 2) == 0: v else: -v)
  else:
    let
      h = hash and 31
      u = if h < 24: offset[0] else: offset[1]
      v = if h < 16: offset[1] else: offset[2]
      w = if h < 8: offset[2] else: offset[3]
    (if (h and 1) == 0: u else: -u) +
      (if (h and 2) == 0: v else: -v) +
      (if (h and 4) == 0: w else: -w)

func noise[N: static[int]](
  perlin: Perlin,
  point: array[N, float32]
): float32 =
  ## Samples a lattice using Gustavson's public-domain noise1234 conventions.
  var
    lattice: array[N, int]
    fraction, weights: array[N, float32]
    corners: array[1 shl N, float32]
  for axis in 0 ..< N:
    let cell = floor(point[axis])
    lattice[axis] = cell.int and 255
    fraction[axis] = point[axis] - cell
    weights[axis] = fade(fraction[axis])
  for corner in 0 ..< corners.len:
    var
      hash = 0
      offset: array[N, float32]
    for axis in countdown(N - 1, 0):
      let
        upper = (corner shr (N - 1 - axis)) and 1
        index = (lattice[axis] + upper) and 255
      hash = perlin.perm[index + hash].int
      offset[axis] = fraction[axis] - upper.float32
    corners[corner] = gradient(hash, offset)
  var count = corners.len
  for axis in countdown(N - 1, 0):
    count = count shr 1
    for corner in 0 ..< count:
      corners[corner] = interpolate(
        weights[axis], corners[corner * 2], corners[corner * 2 + 1]
      )
  PerlinScales[N - 1] * corners[0]

func sample[N: static[int]](
  perlin: Perlin,
  point: array[N, float32]
): float32 {.raises: [NoisyError].} =
  ## Layers octaves with the same amplitude convention as Simplex.
  if perlin.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  var
    amplitude = perlin.amplitude
    frequency = perlin.frequency
  for octave in 0 ..< perlin.octaves:
    var scaled: array[N, float32]
    for axis in 0 ..< N:
      scaled[axis] = point[axis] * frequency
    result += perlin.noise(scaled) * amplitude
    amplitude *= perlin.gain
    frequency *= perlin.lacunarity
  result /= perlin.octaves.float32

func value*(perlin: Perlin, x: float32): float32 =
  ## Samples native 1D Improved Perlin noise with the configured octaves.
  perlin.sample([x])

func value*(perlin: Perlin, x, y: float32): float32 =
  ## Samples native 2D Improved Perlin noise with the configured octaves.
  perlin.sample([x, y])

func value*(perlin: Perlin, x, y, z: float32): float32 =
  ## Samples native 3D Improved Perlin noise with the configured octaves.
  perlin.sample([x, y, z])

func value*(perlin: Perlin, x, y, z, w: float32): float32 =
  ## Samples native 4D Improved Perlin noise with the configured octaves.
  perlin.sample([x, y, z, w])

func value*(perlin: Perlin, x: int): float32 {.inline.} =
  ## Samples 1D Improved Perlin noise at an integer coordinate.
  perlin.value(x.float32)

func value*(perlin: Perlin, x, y: int): float32 {.inline.} =
  ## Samples 2D Improved Perlin noise at integer coordinates.
  perlin.value(x.float32, y.float32)

func value*(perlin: Perlin, x, y, z: int): float32 {.inline.} =
  ## Samples 3D Improved Perlin noise at integer coordinates.
  perlin.value(x.float32, y.float32, z.float32)

func value*(perlin: Perlin, x, y, z, w: int): float32 {.inline.} =
  ## Samples 4D Improved Perlin noise at integer coordinates.
  perlin.value(x.float32, y.float32, z.float32, w.float32)

when useGridNeon:
  func sampleBatch[N: static[int]](
    perlin: Perlin, point: array[N, Float4]
  ): Float4 =
    ## Layers four independently positioned samples with the scalar settings.
    var
      amplitude = perlin.amplitude
      frequency = perlin.frequency
      total = float4(0)
    for _ in 0 ..< perlin.octaves:
      let scale = float4(frequency)
      var scaled: array[N, Float4]
      for axis in 0 ..< N:
        scaled[axis] = point[axis] * scale
      let sample = float4(PerlinScales[N - 1]) *
        samplePerlinGrid4(perlin.perm, scaled)
      total += sample * float4(amplitude)
      amplitude *= perlin.gain
      frequency *= perlin.lacunarity
    total / float4(perlin.octaves.float32)

  func sample1x4(perlin: Perlin, x: Float4): Float4 =
    perlin.sampleBatch([x])

  func sample2x4(perlin: Perlin, x, y: Float4): Float4 =
    perlin.sampleBatch([x, y])

  func sample3x4(perlin: Perlin, x, y, z: Float4): Float4 =
    perlin.sampleBatch([x, y, z])

  func sample4x4(perlin: Perlin, x, y, z, w: Float4): Float4 =
    perlin.sampleBatch([x, y, z, w])

func grid*(perlin: Perlin, start: float32, width: int): Grid =
  ## Samples a 1D Improved Perlin grid at unit intervals from start.
  when useGridNeon:
    perlin.fillGrid(start, width, sample1x4)
  else:
    perlin.scalarGrid(start, width)

func grid*(perlin: Perlin, start, width: int): Grid {.inline.} =
  ## Samples a 1D Improved Perlin grid from an integer start.
  perlin.grid(start.float32, width)

func grid*(
  perlin: Perlin, start: (float32, float32), dimens: (int, int)
): Grid =
  ## Samples a 2D Improved Perlin grid at unit intervals from start.
  when useGridNeon:
    perlin.fillGrid(start, dimens, sample2x4)
  else:
    perlin.scalarGrid(start, dimens)

func grid*(
  perlin: Perlin, start: (int, int), dimens: (int, int)
): Grid {.inline.} =
  ## Samples a 2D Improved Perlin grid from an integer start.
  perlin.grid((start[0].float32, start[1].float32), dimens)

func grid*(
  perlin: Perlin,
  start: (float32, float32, float32),
  dimens: (int, int, int)
): Grid =
  ## Samples a 3D Improved Perlin grid at unit intervals from start.
  when useGridNeon:
    perlin.fillGrid(start, dimens, sample3x4)
  else:
    perlin.scalarGrid(start, dimens)

func grid*(
  perlin: Perlin, start: (int, int, int), dimens: (int, int, int)
): Grid {.inline.} =
  ## Samples a 3D Improved Perlin grid from an integer start.
  perlin.grid((start[0].float32, start[1].float32, start[2].float32), dimens)

func grid*(
  perlin: Perlin,
  start: (float32, float32, float32, float32),
  dimens: (int, int, int, int)
): Grid =
  ## Samples a 4D Improved Perlin grid at unit intervals from start.
  when useGridNeon:
    perlin.fillGrid(start, dimens, sample4x4)
  else:
    perlin.scalarGrid(start, dimens)

func grid*(
  perlin: Perlin, start: (int, int, int, int), dimens: (int, int, int, int)
): Grid {.inline.} =
  ## Samples a 4D Improved Perlin grid from an integer start.
  perlin.grid((
    start[0].float32, start[1].float32, start[2].float32, start[3].float32
  ), dimens)
