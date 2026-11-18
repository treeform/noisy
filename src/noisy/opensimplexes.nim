import common, grids, lattices

export common, Grid, `[]`

when useGridNeon:
  import simd, simd/[opensimplex2neon, opensimplex2sneon]

type
  OpenSimplex2* = object
    ## OpenSimplex2 using FastNoiseLite's 2D and 3D conventions.
    seed*, octaves*: int
    amplitude*, frequency*, lacunarity*, gain*: float32

  OpenSimplex2S* = object
    ## OpenSimplex2S using FastNoiseLite's 2D and 3D conventions.
    seed*, octaves*: int
    amplitude*, frequency*, lacunarity*, gain*: float32

  LatticeAxis = tuple[cell: uint32, delta: float32]

func initOpenSimplex2*(seed: int): OpenSimplex2 =
  ## Initializes OpenSimplex2 with a seed and one octave.
  result.seed = seed
  result.octaves = 1
  result.amplitude = 1
  result.frequency = 1
  result.lacunarity = 2
  result.gain = 0.5

func initOpenSimplex2S*(seed: int): OpenSimplex2S =
  ## Initializes the smoother OpenSimplex2 variant with one octave.
  result.seed = seed
  result.octaves = 1
  result.amplitude = 1
  result.frequency = 1
  result.lacunarity = 2
  result.gain = 0.5

func sample(seed: int32, x, y: float32): float32 =
  ## Applies FastNoiseLite's 2D skew and triangular lattice contributions.
  # Retain the reference's float32 rounding before deriving skew constants.
  const Sqrt3 = 1.73205077648162841796875'f
  let
    f2 = 0.5'f * (Sqrt3 - 1)
    g2 = (3 - Sqrt3) / 6
    skew = (x + y) * f2
    xs = x + skew
    ys = y + skew
    ix = fastFloor(xs)
    iy = fastFloor(ys)
    xi = xs - ix.float32
    yi = ys - iy.float32
    unskew = (xi + yi) * g2
    x0 = xi - unskew
    y0 = yi - unskew
    ip = prime(ix, PrimeX)
    jp = prime(iy, PrimeY)
    a = 0.5'f - x0 * x0 - y0 * y0
    c = (2 * (1 - 2 * g2) * (1 / g2 - 2)) * unskew +
      (-2 * (1 - 2 * g2) * (1 - 2 * g2) + a)
  var n0, n1, n2: float32
  if a > 0:
    n0 = (a * a) * (a * a) * gradient(seed, ip, jp, x0, y0)
  if c > 0:
    let
      x2 = x0 + (2 * g2 - 1)
      y2 = y0 + (2 * g2 - 1)
    n2 = (c * c) * (c * c) *
      gradient(seed, ip + PrimeX, jp + PrimeY, x2, y2)
  if y0 > x0:
    let
      x1 = x0 + g2
      y1 = y0 + (g2 - 1)
      b = 0.5'f - x1 * x1 - y1 * y1
    if b > 0:
      n1 = (b * b) * (b * b) *
        gradient(seed, ip, jp + PrimeY, x1, y1)
  else:
    let
      x1 = x0 + (g2 - 1)
      y1 = y0 + g2
      b = 0.5'f - x1 * x1 - y1 * y1
    if b > 0:
      n1 = (b * b) * (b * b) *
        gradient(seed, ip + PrimeX, jp, x1, y1)
  (n0 + n1 + n2) * 99.83685446303647'f

func sample(seed: int32, x, y, z: float32): float32 =
  ## Samples two offset cube lattices after the default OpenSimplex2 rotation.
  const R3 = (2.0 / 3.0).float32
  let
    rotation = (x + y + z) * R3
    xr = rotation - x
    yr = rotation - y
    zr = rotation - z
    ix = fastRound(xr)
    iy = fastRound(yr)
    iz = fastRound(zr)
  var
    x0 = xr - ix.float32
    y0 = yr - iy.float32
    z0 = zr - iz.float32
    xSign = (-1.0'f - x0).int32 or 1
    ySign = (-1.0'f - y0).int32 or 1
    zSign = (-1.0'f - z0).int32 or 1
    ax0 = xSign.float32 * -x0
    ay0 = ySign.float32 * -y0
    az0 = zSign.float32 * -z0
    ip = prime(ix, PrimeX)
    jp = prime(iy, PrimeY)
    kp = prime(iz, PrimeZ)
    a = (0.6'f - x0 * x0) - (y0 * y0 + z0 * z0)
    latticeSeed = seed
  for lattice in 0 .. 1:
    if a > 0:
      result += (a * a) * (a * a) *
        gradient(latticeSeed, ip, jp, kp, x0, y0, z0)
    var
      b = a + 1
      i1 = ip
      j1 = jp
      k1 = kp
      x1 = x0
      y1 = y0
      z1 = z0
    if ax0 >= ay0 and ax0 >= az0:
      x1 += xSign.float32
      b -= (xSign * 2).float32 * x1
      i1 -= cast[uint32](xSign) * PrimeX
    elif ay0 > ax0 and ay0 >= az0:
      y1 += ySign.float32
      b -= (ySign * 2).float32 * y1
      j1 -= cast[uint32](ySign) * PrimeY
    else:
      z1 += zSign.float32
      b -= (zSign * 2).float32 * z1
      k1 -= cast[uint32](zSign) * PrimeZ
    if b > 0:
      result += (b * b) * (b * b) *
        gradient(latticeSeed, i1, j1, k1, x1, y1, z1)
    if lattice == 1:
      break
    ax0 = 0.5'f - ax0
    ay0 = 0.5'f - ay0
    az0 = 0.5'f - az0
    x0 = xSign.float32 * ax0
    y0 = ySign.float32 * ay0
    z0 = zSign.float32 * az0
    a += (0.75'f - ax0) - (ay0 + az0)
    if xSign < 0:
      ip += PrimeX
    if ySign < 0:
      jp += PrimeY
    if zSign < 0:
      kp += PrimeZ
    xSign = -xSign
    ySign = -ySign
    zSign = -zSign
    latticeSeed = not latticeSeed
  result *= 32.69428253173828125'f

func value*(noise: OpenSimplex2, x, y: float32): float32 =
  ## Samples native 2D OpenSimplex2 with the configured octave settings.
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

func value*(noise: OpenSimplex2, x, y, z: float32): float32 =
  ## Samples native 3D OpenSimplex2 with the configured octave settings.
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

func value*(noise: OpenSimplex2, x, y: int): float32 {.inline.} =
  ## Samples 2D OpenSimplex2 at integer coordinates.
  noise.value(x.float32, y.float32)

func value*(noise: OpenSimplex2, x, y, z: int): float32 {.inline.} =
  ## Samples 3D OpenSimplex2 at integer coordinates.
  noise.value(x.float32, y.float32, z.float32)

func contribution(
  seed: int32,
  x, y: uint32,
  dx, dy, attenuation: float32
): float32 {.inline.} =
  ## Weights a 2D OpenSimplex2S gradient by its compact support.
  (attenuation * attenuation) * (attenuation * attenuation) *
    gradient(seed, x, y, dx, dy)

func contribution(
  seed: int32, x, y, z: LatticeAxis, attenuation: float32
): float32 {.inline.} =
  ## Weights a 3D OpenSimplex2S gradient by its compact support.
  (attenuation * attenuation) * (attenuation * attenuation) *
    gradient(seed, x.cell, y.cell, z.cell, x.delta, y.delta, z.delta)

func sampleSmooth(seed: int32, x, y: float32): float32 =
  ## Samples FastNoiseLite's modified triangular OpenSimplex2S lattice.
  const Sqrt3 = 1.73205077648162841796875'f
  let
    f2 = 0.5'f * (Sqrt3 - 1)
    g2 = (3 - Sqrt3) / 6
    skew = (x + y) * f2
    xs = x + skew
    ys = y + skew
    ix = fastFloor(xs)
    iy = fastFloor(ys)
    xi = xs - ix.float32
    yi = ys - iy.float32
    ip = prime(ix, PrimeX)
    jp = prime(iy, PrimeY)
    i1 = ip + PrimeX
    j1 = jp + PrimeY
    unskew = (xi + yi) * g2
    x0 = xi - unskew
    y0 = yi - unskew
    a0 = (2.0'f / 3.0'f) - x0 * x0 - y0 * y0
    a1 = (2 * (1 - 2 * g2) * (1 / g2 - 2)) * unskew +
      (-2 * (1 - 2 * g2) * (1 - 2 * g2) + a0)
    x1 = x0 - (1 - 2 * g2)
    y1 = y0 - (1 - 2 * g2)
    difference = xi - yi
  result = contribution(seed, ip, jp, x0, y0, a0)
  result += contribution(seed, i1, j1, x1, y1, a1)
  if unskew > g2:
    if xi + difference > 1:
      let
        x2 = x0 + (3 * g2 - 2)
        y2 = y0 + (3 * g2 - 1)
        a2 = (2.0'f / 3.0'f) - x2 * x2 - y2 * y2
      if a2 > 0:
        result += contribution(
          seed, ip + (PrimeX shl 1), jp + PrimeY, x2, y2, a2
        )
    else:
      let
        x2 = x0 + g2
        y2 = y0 + (g2 - 1)
        a2 = (2.0'f / 3.0'f) - x2 * x2 - y2 * y2
      if a2 > 0:
        result += contribution(seed, ip, jp + PrimeY, x2, y2, a2)
    if yi - difference > 1:
      let
        x3 = x0 + (3 * g2 - 1)
        y3 = y0 + (3 * g2 - 2)
        a3 = (2.0'f / 3.0'f) - x3 * x3 - y3 * y3
      if a3 > 0:
        result += contribution(
          seed, ip + PrimeX, jp + (PrimeY shl 1), x3, y3, a3
        )
    else:
      let
        x3 = x0 + (g2 - 1)
        y3 = y0 + g2
        a3 = (2.0'f / 3.0'f) - x3 * x3 - y3 * y3
      if a3 > 0:
        result += contribution(seed, ip + PrimeX, jp, x3, y3, a3)
  else:
    if xi + difference < 0:
      let
        x2 = x0 + (1 - g2)
        y2 = y0 - g2
        a2 = (2.0'f / 3.0'f) - x2 * x2 - y2 * y2
      if a2 > 0:
        result += contribution(seed, ip - PrimeX, jp, x2, y2, a2)
    else:
      let
        x2 = x0 + (g2 - 1)
        y2 = y0 + g2
        a2 = (2.0'f / 3.0'f) - x2 * x2 - y2 * y2
      if a2 > 0:
        result += contribution(seed, ip + PrimeX, jp, x2, y2, a2)
    if yi < difference:
      let
        x2 = x0 - g2
        y2 = y0 - (g2 - 1)
        a2 = (2.0'f / 3.0'f) - x2 * x2 - y2 * y2
      if a2 > 0:
        result += contribution(seed, ip, jp - PrimeY, x2, y2, a2)
    else:
      let
        x2 = x0 + g2
        y2 = y0 + (g2 - 1)
        a2 = (2.0'f / 3.0'f) - x2 * x2 - y2 * y2
      if a2 > 0:
        result += contribution(seed, ip, jp + PrimeY, x2, y2, a2)
  result *= 18.24196194486065'f

func sampleSmooth(seed: int32, x, y, z: float32): float32 =
  ## Samples two offset cube lattices after the default OpenSimplex2 rotation.
  const
    Rotation = (2.0 / 3.0).float32
    RadiusSquared = 0.75'f
    AxisPrimes = [PrimeX, PrimeY, PrimeZ]
  let
    rotation = (x + y + z) * Rotation
    rotated = [rotation - x, rotation - y, rotation - z]
    secondSeed = cast[int32](cast[uint32](seed) + 1293373'u32)
  var
    nearest, opposite, centered, outer: array[3, LatticeAxis]
    firstFlip, secondFlip: array[3, float32]
  for axis, coordinate in rotated:
    let
      cell = fastFloor(coordinate)
      fraction = coordinate - cell.float32
      axisPrime = AxisPrimes[axis]
      hashed = prime(cell, axisPrime)
      mask = (-0.5'f - fraction).int32
      bits = cast[uint32](mask)
      direction = (mask or 1).float32
      firstDelta = fraction + mask.float32
      secondDelta = fraction - 0.5'f
    nearest[axis] = (hashed + (bits and axisPrime), firstDelta)
    opposite[axis] = (hashed + ((not bits) and axisPrime), firstDelta - direction)
    centered[axis] = (hashed + axisPrime, secondDelta)
    outer[axis] = (hashed + (bits and (axisPrime shl 1)), direction + secondDelta)
    firstFlip[axis] = (2 * direction) * secondDelta
    secondFlip[axis] = (-2 - 4 * mask).float32 * secondDelta - 1.0'f
  let
    x0 = nearest[0].delta
    y0 = nearest[1].delta
    z0 = nearest[2].delta
    x1 = centered[0].delta
    y1 = centered[1].delta
    z1 = centered[2].delta
    a0 = RadiusSquared - x0 * x0 - y0 * y0 - z0 * z0
    a1 = RadiusSquared - x1 * x1 - y1 * y1 - z1 * z1
  result = contribution(seed, nearest[0], nearest[1], nearest[2], a0)
  result += contribution(secondSeed, centered[0], centered[1], centered[2], a1)

  # For each axis, consider the adjacent first-lattice corner or its
  # complementary pair and the corresponding second-lattice corner.
  var usedOuterX, usedOuterY, usedOuterZ: bool
  let alongX = firstFlip[0] + a0
  if alongX > 0:
    result += contribution(seed, opposite[0], nearest[1], nearest[2], alongX)
  else:
    let acrossYZ = firstFlip[1] + firstFlip[2] + a0
    if acrossYZ > 0:
      result += contribution(seed, nearest[0], opposite[1], opposite[2], acrossYZ)
    let outerX = secondFlip[0] + a1
    if outerX > 0:
      result += contribution(secondSeed, outer[0], centered[1], centered[2], outerX)
      usedOuterX = true

  let alongY = firstFlip[1] + a0
  if alongY > 0:
    result += contribution(seed, nearest[0], opposite[1], nearest[2], alongY)
  else:
    let acrossXZ = firstFlip[0] + firstFlip[2] + a0
    if acrossXZ > 0:
      result += contribution(seed, opposite[0], nearest[1], opposite[2], acrossXZ)
    let outerY = secondFlip[1] + a1
    if outerY > 0:
      result += contribution(secondSeed, centered[0], outer[1], centered[2], outerY)
      usedOuterY = true

  let alongZ = firstFlip[2] + a0
  if alongZ > 0:
    result += contribution(seed, nearest[0], nearest[1], opposite[2], alongZ)
  else:
    let acrossXY = firstFlip[0] + firstFlip[1] + a0
    if acrossXY > 0:
      result += contribution(seed, opposite[0], opposite[1], nearest[2], acrossXY)
    let outerZ = secondFlip[2] + a1
    if outerZ > 0:
      result += contribution(secondSeed, centered[0], centered[1], outer[2], outerZ)
      usedOuterZ = true

  # A second-lattice pair is outside the support if its opposite axis
  # already contributed a second-lattice corner.
  if not usedOuterX:
    let acrossYZ = secondFlip[1] + secondFlip[2] + a1
    if acrossYZ > 0:
      result += contribution(secondSeed, centered[0], outer[1], outer[2], acrossYZ)
  if not usedOuterY:
    let acrossXZ = secondFlip[0] + secondFlip[2] + a1
    if acrossXZ > 0:
      result += contribution(secondSeed, outer[0], centered[1], outer[2], acrossXZ)
  if not usedOuterZ:
    let acrossXY = secondFlip[0] + secondFlip[1] + a1
    if acrossXY > 0:
      result += contribution(secondSeed, outer[0], outer[1], centered[2], acrossXY)
  result *= 9.046026385208288'f

func value*(noise: OpenSimplex2S, x, y: float32): float32 =
  ## Samples native 2D OpenSimplex2S with the configured octave settings.
  if noise.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  var
    amplitude = noise.amplitude
    frequency = noise.frequency
  for _ in 0 ..< noise.octaves:
    result += sampleSmooth(cast[int32](noise.seed), x * frequency, y * frequency) *
      amplitude
    amplitude *= noise.gain
    frequency *= noise.lacunarity
  result /= noise.octaves.float32

func value*(noise: OpenSimplex2S, x, y, z: float32): float32 =
  ## Samples native 3D OpenSimplex2S with the configured octave settings.
  if noise.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  var
    amplitude = noise.amplitude
    frequency = noise.frequency
  for _ in 0 ..< noise.octaves:
    result += sampleSmooth(
      cast[int32](noise.seed), x * frequency, y * frequency, z * frequency
    ) * amplitude
    amplitude *= noise.gain
    frequency *= noise.lacunarity
  result /= noise.octaves.float32

func value*(noise: OpenSimplex2S, x, y: int): float32 {.inline.} =
  ## Samples 2D OpenSimplex2S at integer coordinates.
  noise.value(x.float32, y.float32)

func value*(noise: OpenSimplex2S, x, y, z: int): float32 {.inline.} =
  ## Samples 3D OpenSimplex2S at integer coordinates.
  noise.value(x.float32, y.float32, z.float32)

when useGridNeon:
  func sample2x4[T: OpenSimplex2 | OpenSimplex2S](noise: T, x, y: Float4): Float4 =
    ## Samples four 2D positions through the same octave settings as value.
    let seed = cast[int32](noise.seed)
    var
      amplitude = noise.amplitude
      frequency = noise.frequency
    result = float4(0)
    for _ in 0 ..< noise.octaves:
      let scale = float4(frequency)
      let sample =
        when T is OpenSimplex2: sample4(seed, x * scale, y * scale)
        else: sampleSmooth4(seed, x * scale, y * scale)
      result += sample * float4(amplitude)
      amplitude *= noise.gain
      frequency *= noise.lacunarity
    result /= float4(noise.octaves.float32)

  func sample3x4[T: OpenSimplex2 | OpenSimplex2S](
    noise: T, x, y, z: Float4
  ): Float4 =
    ## Samples four 3D positions through the same octave settings as value.
    let seed = cast[int32](noise.seed)
    var
      amplitude = noise.amplitude
      frequency = noise.frequency
    result = float4(0)
    for _ in 0 ..< noise.octaves:
      let scale = float4(frequency)
      let sample =
        when T is OpenSimplex2: sample4(seed, x * scale, y * scale, z * scale)
        else: sampleSmooth4(seed, x * scale, y * scale, z * scale)
      result += sample * float4(amplitude)
      amplitude *= noise.gain
      frequency *= noise.lacunarity
    result /= float4(noise.octaves.float32)

func grid*(
  noise: OpenSimplex2, start: (float32, float32), dimens: (int, int)
): Grid =
  ## Samples a 2D OpenSimplex2 grid at unit intervals from start.
  when useGridNeon:
    noise.fillGrid(start, dimens, sample2x4[OpenSimplex2])
  else:
    noise.scalarGrid(start, dimens)

func grid*(
  noise: OpenSimplex2, start: (int, int), dimens: (int, int)
): Grid {.inline.} =
  ## Samples a 2D OpenSimplex2 grid from an integer start.
  noise.grid((start[0].float32, start[1].float32), dimens)

func grid*(
  noise: OpenSimplex2,
  start: (float32, float32, float32),
  dimens: (int, int, int)
): Grid =
  ## Samples a 3D OpenSimplex2 grid at unit intervals from start.
  when useGridNeon:
    noise.fillGrid(start, dimens, sample3x4[OpenSimplex2])
  else:
    noise.scalarGrid(start, dimens)

func grid*(
  noise: OpenSimplex2, start: (int, int, int), dimens: (int, int, int)
): Grid {.inline.} =
  ## Samples a 3D OpenSimplex2 grid from an integer start.
  noise.grid((start[0].float32, start[1].float32, start[2].float32), dimens)

func grid*(
  noise: OpenSimplex2S, start: (float32, float32), dimens: (int, int)
): Grid =
  ## Samples a 2D OpenSimplex2S grid at unit intervals from start.
  when useGridNeon:
    noise.fillGrid(start, dimens, sample2x4[OpenSimplex2S])
  else:
    noise.scalarGrid(start, dimens)

func grid*(
  noise: OpenSimplex2S, start: (int, int), dimens: (int, int)
): Grid {.inline.} =
  ## Samples a 2D OpenSimplex2S grid from an integer start.
  noise.grid((start[0].float32, start[1].float32), dimens)

func grid*(
  noise: OpenSimplex2S,
  start: (float32, float32, float32),
  dimens: (int, int, int)
): Grid =
  ## Samples a 3D OpenSimplex2S grid at unit intervals from start.
  when useGridNeon:
    noise.fillGrid(start, dimens, sample3x4[OpenSimplex2S])
  else:
    noise.scalarGrid(start, dimens)

func grid*(
  noise: OpenSimplex2S, start: (int, int, int), dimens: (int, int, int)
): Grid {.inline.} =
  ## Samples a 3D OpenSimplex2S grid from an integer start.
  noise.grid((start[0].float32, start[1].float32, start[2].float32), dimens)
