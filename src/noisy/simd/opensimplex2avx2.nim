import ../lattices, avx2base, latticesavx2
include avx2flags

func positiveWeight(attenuation: Float8): Float8 {.avxInline.} =
  ## Keeps the fourth power of attenuation within the corner's support.
  let
    weight = mm256_max_ps(attenuation, float8(0))
    squared = weight * weight
  squared * squared

func sample8*(seed: int32, x, y: Float8): Float8 {.avxInline.} =
  ## Evaluates eight lanes on the triangular OpenSimplex2 lattice.
  const
    Sqrt3 = 1.73205077648162841796875'f
    F2 = 0.5'f * (Sqrt3 - 1)
    G2 = (3 - Sqrt3) / 6
  let
    skew = (x + y) * float8(F2)
    xs = x + skew
    ys = y + skew
    ix = fastFloor8(xs)
    iy = fastFloor8(ys)
    xi = xs - toFloats(ix)
    yi = ys - toFloats(iy)
    unskew = (xi + yi) * float8(G2)
    x0 = xi - unskew
    y0 = yi - unskew
    xp = uints8(PrimeX)
    yp = uints8(PrimeY)
    ip = mm256_mullo_epi32(ix, xp)
    jp = mm256_mullo_epi32(iy, yp)
    i2 = mm256_add_epi32(ip, xp)
    j2 = mm256_add_epi32(jp, yp)
    alongY = less(x0, y0)
    x1 = x0 + select(alongY, float8(G2), float8(G2 - 1))
    y1 = y0 + select(alongY, float8(G2 - 1), float8(G2))
    i1 = select(alongY, ip, i2)
    j1 = select(alongY, j2, jp)
    x2 = x0 + float8(2 * G2 - 1)
    y2 = y0 + float8(2 * G2 - 1)
    a = float8(0.5) - x0 * x0 - y0 * y0
    b = float8(0.5) - x1 * x1 - y1 * y1
    c = float8(2 * (1 - 2 * G2) * (1 / G2 - 2)) * unskew +
      (float8(-2 * (1 - 2 * G2) * (1 - 2 * G2)) + a)
    n0 = positiveWeight(a) * gradient8(seed, ip, jp, x0, y0)
    n1 = positiveWeight(b) * gradient8(seed, i1, j1, x1, y1)
    n2 = positiveWeight(c) * gradient8(seed, i2, j2, x2, y2)
  (n0 + n1 + n2) * float8(99.83685446303647'f)

type LatticeAxis = object
  delta, direction, magnitude: Float8
  cell, negative: Int8

func nearestAxis(coordinate: Float8, axisPrime: static[uint32]): LatticeAxis {.avxInline.} =
  ## Computes the nearest first-lattice cell and its signed displacement.
  let
    bias = select(less(coordinate, float8(0)), float8(-0.5), float8(0.5))
    cell = truncate(coordinate + bias)
    delta = coordinate - toFloats(cell)
    direction = toFloats(mm256_or_si256(
      truncate(float8(-1) - delta), int8(1)))
  LatticeAxis(
    delta: delta,
    direction: direction,
    magnitude: (float8(0) - direction) * delta,
    cell: mm256_mullo_epi32(cell, uints8(axisPrime)),
    negative: less(direction, float8(0))
  )

func secondAxis(axis: LatticeAxis, axisPrime: static[uint32]): LatticeAxis {.avxInline.} =
  ## Shifts to the corresponding corner of the offset cubic lattice.
  let magnitude = float8(0.5) - axis.magnitude
  LatticeAxis(
    delta: axis.direction * magnitude,
    direction: float8(0) - axis.direction,
    magnitude: magnitude,
    cell: mm256_add_epi32(axis.cell, mm256_and_si256(axis.negative, uints8(axisPrime))),
    negative: bitNot(axis.negative)
  )

func neighborCell(
  axis: LatticeAxis, selected: Int8, axisPrime: static[uint32]
): Int8 {.avxInline.} =
  ## Moves the selected lanes one cell along the indicated axis.
  let step = select(axis.negative, uints8(0'u32 - axisPrime),
    uints8(axisPrime))
  mm256_sub_epi32(axis.cell, mm256_and_si256(selected, step))

func addLattice(
  total: var Float8, seed: int32, x, y, z: LatticeAxis, attenuation: Float8
) {.avxInline.} =
  ## Accumulates the nearest corner and one axis neighbor in each lane.
  total += positiveWeight(attenuation) *
    gradient8(seed, x.cell, y.cell, z.cell, x.delta, y.delta, z.delta)
  let
    xAtLeastY = greaterEqual(x.magnitude, y.magnitude)
    alongX = mm256_and_si256(xAtLeastY, greaterEqual(x.magnitude, z.magnitude))
    alongY = bitClear(greaterEqual(y.magnitude, z.magnitude), xAtLeastY)
    alongZ = bitNot(mm256_or_si256(alongX, alongY))
    dx = select(alongX, x.delta + x.direction, x.delta)
    dy = select(alongY, y.delta + y.direction, y.delta)
    dz = select(alongZ, z.delta + z.direction, z.delta)
    maxMagnitude = mm256_max_ps(x.magnitude, mm256_max_ps(y.magnitude, z.magnitude))
    # Along the selected axis, sign * (delta + sign) is 1 - magnitude.
    neighborWeight = (attenuation + float8(1)) -
      float8(2) * (float8(1) - maxMagnitude)
  total += positiveWeight(neighborWeight) * gradient8(seed,
    neighborCell(x, alongX, PrimeX), neighborCell(y, alongY, PrimeY),
    neighborCell(z, alongZ, PrimeZ), dx, dy, dz)

func sample8*(seed: int32, x, y, z: Float8): Float8 {.avxInline.} =
  ## Evaluates eight lanes on the two offset cubic OpenSimplex2 lattices.
  let
    rotation = (x + y + z) * float8((2.0 / 3.0).float32)
    ax = nearestAxis(rotation - x, PrimeX)
    ay = nearestAxis(rotation - y, PrimeY)
    az = nearestAxis(rotation - z, PrimeZ)
    bx = secondAxis(ax, PrimeX)
    by = secondAxis(ay, PrimeY)
    bz = secondAxis(az, PrimeZ)
    a = (float8(0.6) - ax.delta * ax.delta) -
      (ay.delta * ay.delta + az.delta * az.delta)
    b = a + ((float8(0.75) - bx.magnitude) - (by.magnitude + bz.magnitude))
  result = float8(0)
  addLattice(result, seed, ax, ay, az, a)
  addLattice(result, not seed, bx, by, bz, b)
  result = result * float8(32.69428253173828125'f)
