import ../lattices, avx2base, latticesavx2
include avx2flags

type
  Corner = tuple[cell: Int8, delta: Float8]
  Axis = object
    nearest, opposite, centered, outer: Corner
    firstFlip, secondFlip: Float8

func choose(mask: Int8, yes, no: Corner): Corner {.avxInline.} =
  ## Chooses a lattice coordinate and displacement independently in each lane.
  (select(mask, yes.cell, no.cell), select(mask, yes.delta, no.delta))

func fourthPower(a: Float8): Float8 {.avxInline.} =
  ## Computes the compact-support attenuation weight.
  let squared = a * a
  squared * squared

func contribution(
  seed: int32, x, y: Corner, attenuation: Float8
): Float8 {.avxInline.} =
  ## Weights eight 2D gradients at their selected corners.
  fourthPower(attenuation) * gradient8(seed, x.cell, y.cell, x.delta, y.delta)

func contribution(
  seed: int32, x, y, z: Corner, attenuation: Float8
): Float8 {.avxInline.} =
  ## Weights eight 3D gradients at their selected corners.
  fourthPower(attenuation) *
    gradient8(seed, x.cell, y.cell, z.cell, x.delta, y.delta, z.delta)

func optionalContribution(
  seed: int32, x, y: Corner
): Float8 {.avxInline.} =
  ## Evaluates optional 2D corners with zero weight outside their support.
  let attenuation = float8(2.0'f / 3.0'f) - x.delta * x.delta - y.delta * y.delta
  contribution(seed, x, y, mm256_max_ps(attenuation, float8(0)))

func sampleSmooth8*(seed: int32, x, y: Float8): Float8 {.avxInline.} =
  ## Evaluates eight lanes on the modified triangular lattice.
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
    xp = uints8(PrimeX)
    yp = uints8(PrimeY)
    ip = mm256_mullo_epi32(ix, xp)
    jp = mm256_mullo_epi32(iy, yp)
    i1 = mm256_add_epi32(ip, xp)
    j1 = mm256_add_epi32(jp, yp)
    unskew = (xi + yi) * float8(G2)
    x0 = xi - unskew
    y0 = yi - unskew
    a0 = float8(2.0'f / 3.0'f) - x0 * x0 - y0 * y0
    a1 = float8(2 * (1 - 2 * G2) * (1 / G2 - 2)) * unskew +
      (float8(-2 * (1 - 2 * G2) * (1 - 2 * G2)) + a0)
    x1 = x0 - float8(1 - 2 * G2)
    y1 = y0 - float8(1 - 2 * G2)
    difference = xi - yi
    upper = less(float8(G2), unskew)
    xFar = mm256_and_si256(upper, less(float8(1), xi + difference))
    yFar = mm256_and_si256(upper, less(float8(1), yi - difference))
    xBelow = bitClear(less(xi + difference, float8(0)), upper)
    yBelow = bitClear(less(yi, difference), upper)
    xNeighborAlongX = (i1, x0 + float8(G2 - 1))
    xNeighborAlongY = (jp, y0 + float8(G2))
    yNeighborAlongX = (ip, x0 + float8(G2))
    yNeighborAlongY = (j1, y0 + float8(G2 - 1))

  result = contribution(seed, (ip, x0), (jp, y0), a0)
  result += contribution(seed, (i1, x1), (j1, y1), a1)

  # Each lane chooses two optional corners, in the scalar addition order.
  let
    firstX = choose(xFar, (mm256_add_epi32(i1, xp), x0 + float8(3 * G2 - 2)),
      choose(xBelow, (mm256_sub_epi32(ip, xp), x0 + float8(1 - G2)),
        choose(upper, yNeighborAlongX, xNeighborAlongX)))
    firstY = choose(xFar, (j1, y0 + float8(3 * G2 - 1)),
      choose(xBelow, (jp, y0 - float8(G2)),
        choose(upper, yNeighborAlongY, xNeighborAlongY)))
    secondX = choose(yFar, (i1, x0 + float8(3 * G2 - 1)),
      choose(yBelow, (ip, x0 - float8(G2)),
        choose(upper, xNeighborAlongX, yNeighborAlongX)))
    secondY = choose(yFar, (mm256_add_epi32(j1, yp), y0 + float8(3 * G2 - 2)),
      choose(yBelow, (mm256_sub_epi32(jp, yp), y0 - float8(G2 - 1)),
        choose(upper, xNeighborAlongY, yNeighborAlongY)))
  result += optionalContribution(seed, firstX, firstY)
  result += optionalContribution(seed, secondX, secondY)
  result = result * float8(18.24196194486065'f)

func latticeAxis(coordinate: Float8, axisPrime: uint32): Axis {.avxInline.} =
  ## Computes corner choices and attenuation changes on one 3D axis.
  let
    cell = fastFloor8(coordinate)
    fraction = coordinate - toFloats(cell)
    primeLanes = uints8(axisPrime)
    hashed = mm256_mullo_epi32(cell, primeLanes)
    mask = truncate(float8(-0.5) - fraction)
    bits = mask
    direction = toFloats(mm256_or_si256(bits, int8(1)))
    firstDelta = fraction + toFloats(mask)
    secondDelta = fraction - float8(0.5)
  result.nearest = (mm256_add_epi32(hashed, mm256_and_si256(bits, primeLanes)), firstDelta)
  result.opposite = (mm256_add_epi32(hashed, bitClear(primeLanes, bits)),
    firstDelta - direction)
  result.centered = (mm256_add_epi32(hashed, primeLanes), secondDelta)
  result.outer = (mm256_add_epi32(hashed, mm256_and_si256(bits, uints8(axisPrime shl 1))),
    direction + secondDelta)
  result.firstFlip = (float8(2) * direction) * secondDelta
  result.secondFlip = (float8(-2) - float8(4) * toFloats(mask)) * secondDelta -
    float8(1)

func sampleSmooth8*(seed: int32, x, y, z: Float8): Float8 {.avxInline.} =
  ## Evaluates eight lanes on two offset 3D lattices.
  let
    rotation = (x + y + z) * float8((2.0 / 3.0).float32)
    ax = latticeAxis(rotation - x, PrimeX)
    ay = latticeAxis(rotation - y, PrimeY)
    az = latticeAxis(rotation - z, PrimeZ)
    secondSeed = cast[int32](cast[uint32](seed) + 1293373'u32)
    zero = float8(0)
    radius = float8(0.75)
    x0 = ax.nearest.delta
    y0 = ay.nearest.delta
    z0 = az.nearest.delta
    x1 = ax.centered.delta
    y1 = ay.centered.delta
    z1 = az.centered.delta
    a0 = radius - x0 * x0 - y0 * y0 - z0 * z0
    a1 = radius - x1 * x1 - y1 * y1 - z1 * z1
  result = contribution(seed, ax.nearest, ay.nearest, az.nearest, a0)
  result += contribution(secondSeed, ax.centered, ay.centered, az.centered, a1)

  template addAxis(
    axis, firstOther, secondOther: Axis,
    alongX, alongY, alongZ, acrossX, acrossY, acrossZ: Corner,
    outerX, outerY, outerZ, pairX, pairY, pairZ: Corner,
    deferred: var Float8
  ) =
    # Select one optional corner per lattice independently in each lane.
    let
      along = axis.firstFlip + a0
      across = firstOther.firstFlip + secondOther.firstFlip + a0
      useAlong = less(zero, along)
      firstWeight = mm256_max_ps(select(useAlong, along, across), zero)
      outer = axis.secondFlip + a1
      pair = firstOther.secondFlip + secondOther.secondFlip + a1
      useOuter = bitClear(less(zero, outer), useAlong)
      secondWeight = mm256_max_ps(select(useOuter, outer, pair), zero)
    result += contribution(seed,
      choose(useAlong, alongX, acrossX),
      choose(useAlong, alongY, acrossY),
      choose(useAlong, alongZ, acrossZ),
      firstWeight)
    let second = contribution(secondSeed,
      choose(useOuter, outerX, pairX),
      choose(useOuter, outerY, pairY),
      choose(useOuter, outerZ, pairZ),
      secondWeight)
    result += select(useOuter, second, zero)
    # Accumulate second-lattice pairs after all axes, as the scalar kernel does.
    deferred = select(useOuter, zero, second)

  var pairX, pairY, pairZ: Float8
  addAxis(ax, ay, az,
    ax.opposite, ay.nearest, az.nearest, ax.nearest, ay.opposite, az.opposite,
    ax.outer, ay.centered, az.centered, ax.centered, ay.outer, az.outer,
    pairX)
  addAxis(ay, ax, az,
    ax.nearest, ay.opposite, az.nearest, ax.opposite, ay.nearest, az.opposite,
    ax.centered, ay.outer, az.centered, ax.outer, ay.centered, az.outer,
    pairY)
  addAxis(az, ax, ay,
    ax.nearest, ay.nearest, az.opposite, ax.opposite, ay.opposite, az.nearest,
    ax.centered, ay.centered, az.outer, ax.outer, ay.outer, az.centered,
    pairZ)
  result += pairX
  result += pairY
  result += pairZ
  result = result * float8(9.046026385208288'f)
