import nimsimd/neon, ../[lattices, simd], latticesneon

func positiveWeight(attenuation: Float4): Float4 {.inline.} =
  ## Keeps the fourth power of attenuation within the corner's support.
  let
    weight = vmaxq_f32(attenuation, float4(0))
    squared = weight * weight
  squared * squared

func sample4*(seed: int32, x, y: Float4): Float4 {.inline.} =
  ## Evaluates four lanes on the triangular OpenSimplex2 lattice.
  const
    Sqrt3 = 1.73205077648162841796875'f
    F2 = 0.5'f * (Sqrt3 - 1)
    G2 = (3 - Sqrt3) / 6
  let
    skew = (x + y) * float4(F2)
    xs = x + skew
    ys = y + skew
    ix = fastFloor4(xs)
    iy = fastFloor4(ys)
    xi = xs - toFloats(ix)
    yi = ys - toFloats(iy)
    unskew = (xi + yi) * float4(G2)
    x0 = xi - unskew
    y0 = yi - unskew
    xp = vmovq_n_u32(PrimeX)
    yp = vmovq_n_u32(PrimeY)
    ip = multiply4(cast[uint32x4](ix), xp)
    jp = multiply4(cast[uint32x4](iy), yp)
    i2 = vaddq_u32(ip, xp)
    j2 = vaddq_u32(jp, yp)
    alongY = vcltq_f32(x0, y0)
    x1 = x0 + vbslq_f32(alongY, float4(G2), float4(G2 - 1))
    y1 = y0 + vbslq_f32(alongY, float4(G2 - 1), float4(G2))
    i1 = vbslq_u32(alongY, ip, i2)
    j1 = vbslq_u32(alongY, j2, jp)
    x2 = x0 + float4(2 * G2 - 1)
    y2 = y0 + float4(2 * G2 - 1)
    a = float4(0.5) - x0 * x0 - y0 * y0
    b = float4(0.5) - x1 * x1 - y1 * y1
    c = float4(2 * (1 - 2 * G2) * (1 / G2 - 2)) * unskew +
      (float4(-2 * (1 - 2 * G2) * (1 - 2 * G2)) + a)
    n0 = positiveWeight(a) * gradient4(seed, ip, jp, x0, y0)
    n1 = positiveWeight(b) * gradient4(seed, i1, j1, x1, y1)
    n2 = positiveWeight(c) * gradient4(seed, i2, j2, x2, y2)
  (n0 + n1 + n2) * float4(99.83685446303647'f)

type LatticeAxis = object
  delta, direction, magnitude: Float4
  cell, negative: uint32x4

func nearestAxis(coordinate: Float4, axisPrime: static[uint32]): LatticeAxis {.inline.} =
  ## Computes the nearest first-lattice cell and its signed displacement.
  let
    bias = vbslq_f32(vcltq_f32(coordinate, float4(0)), float4(-0.5), float4(0.5))
    cell = truncate(coordinate + bias)
    delta = coordinate - toFloats(cell)
    direction = toFloats(cast[Int4](vorrq_u32(
      cast[uint32x4](truncate(float4(-1) - delta)), vmovq_n_u32(1)
    )))
  LatticeAxis(
    delta: delta,
    direction: direction,
    magnitude: (float4(0) - direction) * delta,
    cell: multiply4(cast[uint32x4](cell), vmovq_n_u32(axisPrime)),
    negative: vcltq_f32(direction, float4(0))
  )

func secondAxis(axis: LatticeAxis, axisPrime: static[uint32]): LatticeAxis {.inline.} =
  ## Shifts to the corresponding corner of the offset cubic lattice.
  let magnitude = float4(0.5) - axis.magnitude
  LatticeAxis(
    delta: axis.direction * magnitude,
    direction: float4(0) - axis.direction,
    magnitude: magnitude,
    cell: vaddq_u32(axis.cell, vandq_u32(axis.negative, vmovq_n_u32(axisPrime))),
    negative: vmvnq_u32(axis.negative)
  )

func neighborCell(
  axis: LatticeAxis, selected: uint32x4, axisPrime: static[uint32]
): uint32x4 {.inline.} =
  ## Moves the selected lanes one cell along the indicated axis.
  let step = vbslq_u32(axis.negative, vmovq_n_u32(0'u32 - axisPrime),
    vmovq_n_u32(axisPrime))
  vsubq_u32(axis.cell, vandq_u32(selected, step))

func addLattice(
  total: var Float4, seed: int32, x, y, z: LatticeAxis, attenuation: Float4
) {.inline.} =
  ## Accumulates the nearest corner and one axis neighbor in each lane.
  total += positiveWeight(attenuation) *
    gradient4(seed, x.cell, y.cell, z.cell, x.delta, y.delta, z.delta)
  let
    xAtLeastY = vcgeq_f32(x.magnitude, y.magnitude)
    alongX = vandq_u32(xAtLeastY, vcgeq_f32(x.magnitude, z.magnitude))
    alongY = vbicq_u32(vcgeq_f32(y.magnitude, z.magnitude), xAtLeastY)
    alongZ = vmvnq_u32(vorrq_u32(alongX, alongY))
    dx = vbslq_f32(alongX, x.delta + x.direction, x.delta)
    dy = vbslq_f32(alongY, y.delta + y.direction, y.delta)
    dz = vbslq_f32(alongZ, z.delta + z.direction, z.delta)
    maxMagnitude = vmaxq_f32(x.magnitude, vmaxq_f32(y.magnitude, z.magnitude))
    # Along the selected axis, sign * (delta + sign) is 1 - magnitude.
    neighborWeight = (attenuation + float4(1)) -
      float4(2) * (float4(1) - maxMagnitude)
  total += positiveWeight(neighborWeight) * gradient4(seed,
    neighborCell(x, alongX, PrimeX), neighborCell(y, alongY, PrimeY),
    neighborCell(z, alongZ, PrimeZ), dx, dy, dz)

func sample4*(seed: int32, x, y, z: Float4): Float4 {.inline.} =
  ## Evaluates four lanes on the two offset cubic OpenSimplex2 lattices.
  let
    rotation = (x + y + z) * float4((2.0 / 3.0).float32)
    ax = nearestAxis(rotation - x, PrimeX)
    ay = nearestAxis(rotation - y, PrimeY)
    az = nearestAxis(rotation - z, PrimeZ)
    bx = secondAxis(ax, PrimeX)
    by = secondAxis(ay, PrimeY)
    bz = secondAxis(az, PrimeZ)
    a = (float4(0.6) - ax.delta * ax.delta) -
      (ay.delta * ay.delta + az.delta * az.delta)
    b = a + ((float4(0.75) - bx.magnitude) - (by.magnitude + bz.magnitude))
  result = float4(0)
  addLattice(result, seed, ax, ay, az, a)
  addLattice(result, not seed, bx, by, bz, b)
  result = result * float4(32.69428253173828125'f)
