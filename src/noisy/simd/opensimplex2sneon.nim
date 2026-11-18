import nimsimd/neon, ../lattices, ../simd, latticesneon

type
  Corner = tuple[cell: uint32x4, delta: Float4]
  Axis = object
    nearest, opposite, centered, outer: Corner
    firstFlip, secondFlip: Float4

func choose(mask: uint32x4, yes, no: Corner): Corner {.inline.} =
  ## Chooses a lattice coordinate and displacement independently in each lane.
  (vbslq_u32(mask, yes.cell, no.cell), vbslq_f32(mask, yes.delta, no.delta))

func fourthPower(a: Float4): Float4 {.inline.} =
  ## Computes the compact-support attenuation weight.
  let squared = a * a
  squared * squared

func contribution(
  seed: int32, x, y: Corner, attenuation: Float4
): Float4 {.inline.} =
  ## Weights four 2D gradients at their selected corners.
  fourthPower(attenuation) * gradient4(seed, x.cell, y.cell, x.delta, y.delta)

func contribution(
  seed: int32, x, y, z: Corner, attenuation: Float4
): Float4 {.inline.} =
  ## Weights four 3D gradients at their selected corners.
  fourthPower(attenuation) *
    gradient4(seed, x.cell, y.cell, z.cell, x.delta, y.delta, z.delta)

func optionalContribution(
  seed: int32, x, y: Corner
): Float4 {.inline.} =
  ## Evaluates optional 2D corners with zero weight outside their support.
  let attenuation = float4(2.0'f / 3.0'f) - x.delta * x.delta - y.delta * y.delta
  contribution(seed, x, y, vmaxq_f32(attenuation, float4(0)))

func sampleSmooth4*(seed: int32, x, y: Float4): Float4 {.inline.} =
  ## Evaluates four lanes on the modified triangular lattice.
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
    xp = vmovq_n_u32(PrimeX)
    yp = vmovq_n_u32(PrimeY)
    ip = multiply4(cast[uint32x4](ix), xp)
    jp = multiply4(cast[uint32x4](iy), yp)
    i1 = vaddq_u32(ip, xp)
    j1 = vaddq_u32(jp, yp)
    unskew = (xi + yi) * float4(G2)
    x0 = xi - unskew
    y0 = yi - unskew
    a0 = float4(2.0'f / 3.0'f) - x0 * x0 - y0 * y0
    a1 = float4(2 * (1 - 2 * G2) * (1 / G2 - 2)) * unskew +
      (float4(-2 * (1 - 2 * G2) * (1 - 2 * G2)) + a0)
    x1 = x0 - float4(1 - 2 * G2)
    y1 = y0 - float4(1 - 2 * G2)
    difference = xi - yi
    upper = vcltq_f32(float4(G2), unskew)
    xFar = vandq_u32(upper, vcltq_f32(float4(1), xi + difference))
    yFar = vandq_u32(upper, vcltq_f32(float4(1), yi - difference))
    xBelow = vbicq_u32(vcltq_f32(xi + difference, float4(0)), upper)
    yBelow = vbicq_u32(vcltq_f32(yi, difference), upper)
    xNeighborAlongX = (i1, x0 + float4(G2 - 1))
    xNeighborAlongY = (jp, y0 + float4(G2))
    yNeighborAlongX = (ip, x0 + float4(G2))
    yNeighborAlongY = (j1, y0 + float4(G2 - 1))

  result = contribution(seed, (ip, x0), (jp, y0), a0)
  result += contribution(seed, (i1, x1), (j1, y1), a1)

  # Each lane chooses two optional corners, in the scalar addition order.
  let
    firstX = choose(xFar, (vaddq_u32(i1, xp), x0 + float4(3 * G2 - 2)),
      choose(xBelow, (vsubq_u32(ip, xp), x0 + float4(1 - G2)),
        choose(upper, yNeighborAlongX, xNeighborAlongX)))
    firstY = choose(xFar, (j1, y0 + float4(3 * G2 - 1)),
      choose(xBelow, (jp, y0 - float4(G2)),
        choose(upper, yNeighborAlongY, xNeighborAlongY)))
    secondX = choose(yFar, (i1, x0 + float4(3 * G2 - 1)),
      choose(yBelow, (ip, x0 - float4(G2)),
        choose(upper, xNeighborAlongX, yNeighborAlongX)))
    secondY = choose(yFar, (vaddq_u32(j1, yp), y0 + float4(3 * G2 - 2)),
      choose(yBelow, (vsubq_u32(jp, yp), y0 - float4(G2 - 1)),
        choose(upper, xNeighborAlongY, yNeighborAlongY)))
  result += optionalContribution(seed, firstX, firstY)
  result += optionalContribution(seed, secondX, secondY)
  result = result * float4(18.24196194486065'f)

func latticeAxis(coordinate: Float4, axisPrime: uint32): Axis {.inline.} =
  ## Computes corner choices and attenuation changes on one 3D axis.
  let
    cell = fastFloor4(coordinate)
    fraction = coordinate - toFloats(cell)
    primeLanes = vmovq_n_u32(axisPrime)
    hashed = multiply4(cast[uint32x4](cell), primeLanes)
    mask = truncate(float4(-0.5) - fraction)
    bits = cast[uint32x4](mask)
    direction = toFloats(cast[Int4](vorrq_u32(bits, vmovq_n_u32(1))))
    firstDelta = fraction + toFloats(mask)
    secondDelta = fraction - float4(0.5)
  result.nearest = (vaddq_u32(hashed, vandq_u32(bits, primeLanes)), firstDelta)
  result.opposite = (vaddq_u32(hashed, vbicq_u32(primeLanes, bits)),
    firstDelta - direction)
  result.centered = (vaddq_u32(hashed, primeLanes), secondDelta)
  result.outer = (vaddq_u32(hashed, vandq_u32(bits, vmovq_n_u32(axisPrime shl 1))),
    direction + secondDelta)
  result.firstFlip = (float4(2) * direction) * secondDelta
  result.secondFlip = (float4(-2) - float4(4) * toFloats(mask)) * secondDelta -
    float4(1)

func sampleSmooth4*(seed: int32, x, y, z: Float4): Float4 {.inline.} =
  ## Evaluates four lanes on two offset 3D lattices.
  let
    rotation = (x + y + z) * float4((2.0 / 3.0).float32)
    ax = latticeAxis(rotation - x, PrimeX)
    ay = latticeAxis(rotation - y, PrimeY)
    az = latticeAxis(rotation - z, PrimeZ)
    secondSeed = cast[int32](cast[uint32](seed) + 1293373'u32)
    zero = float4(0)
    radius = float4(0.75)
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
    deferred: var Float4
  ) =
    # Select one optional corner per lattice independently in each lane.
    let
      along = axis.firstFlip + a0
      across = firstOther.firstFlip + secondOther.firstFlip + a0
      useAlong = vcltq_f32(zero, along)
      firstWeight = vmaxq_f32(vbslq_f32(useAlong, along, across), zero)
      outer = axis.secondFlip + a1
      pair = firstOther.secondFlip + secondOther.secondFlip + a1
      useOuter = vbicq_u32(vcltq_f32(zero, outer), useAlong)
      secondWeight = vmaxq_f32(vbslq_f32(useOuter, outer, pair), zero)
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
    result += vbslq_f32(useOuter, second, zero)
    # Accumulate second-lattice pairs after all axes, as the scalar kernel does.
    deferred = vbslq_f32(useOuter, zero, second)

  var pairX, pairY, pairZ: Float4
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
  result = result * float4(9.046026385208288'f)
