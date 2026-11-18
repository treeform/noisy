import common, lattices

export common

type
  OpenSimplex2* = object
    ## Scalar OpenSimplex2 using FastNoiseLite's 2D and 3D conventions.
    seed*, octaves*: int
    amplitude*, frequency*, lacunarity*, gain*: float32

func initOpenSimplex2*(seed: int): OpenSimplex2 =
  ## Initializes OpenSimplex2 with a seed and one octave.
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
