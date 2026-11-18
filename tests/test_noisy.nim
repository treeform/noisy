import
  std/[math, strformat],
  noisy

proc checkClose(actual, expected: float32, context: string) =
  ## Allows small rounding differences between scalar and SIMD samples.
  const
    AbsoluteTolerance = 0.0001'f
    RelativeTolerance = 0.00001'f
  doAssert abs(actual) < Inf and abs(expected) < Inf,
    context & &" non-finite sample: SIMD={actual}, scalar={expected}"
  let
    error = abs(actual - expected)
    tolerance = AbsoluteTolerance +
      RelativeTolerance * max(abs(actual), abs(expected))
  doAssert error <= tolerance,
    context & &" SIMD={actual}, scalar={expected}, " &
      &"error={error}, tolerance={tolerance}"

let simplex = initSimplex(1988) # 1988 is the random seed, generate this however

block continuity_2d:
  const
    F2 = ((sqrt(3.0) - 1) / 2).float32
    Epsilon = 0.000001'f
    Tolerance = 0.0001'f
  for seed in [1, 42, 1988]:
    let noise = initSimplex(seed)
    for y in [-0.25'f, 0.25'f]:
      # Cross i = 0 in the skewed lattice, where masking wraps -1 to 255.
      let
        x = -y * F2 / (1 + F2)
        left = noise.value(x - Epsilon, y)
        right = noise.value(x + Epsilon, y)
        leftGrid = noise.grid((x - Epsilon, y), (4, 4))
        rightGrid = noise.grid((x + Epsilon, y), (4, 4))
      doAssert abs(left - right) < Tolerance
      doAssert abs(leftGrid[0, 0] - rightGrid[0, 0]) < Tolerance
      doAssert abs(leftGrid[0, 0] - left) < Tolerance
      doAssert abs(rightGrid[0, 0] - right) < Tolerance

block continuity_3d:
  const
    Epsilon = 0.000001'f
    Tolerance = 0.0001'f
  for seed in [1, 42, 1988]:
    let
      noise = initSimplex(seed)
      y = 0.4'f
      z = -0.2'f
      # The 3D skew places i = 0 at x = -(y + z) / 4.
      x = -(y + z) / 4
      left = noise.value(x - Epsilon, y, z)
      right = noise.value(x + Epsilon, y, z)
      leftGrid = noise.grid((x - Epsilon, y, z), (4, 4, 4))
      rightGrid = noise.grid((x + Epsilon, y, z), (4, 4, 4))
    doAssert abs(left - right) < Tolerance
    doAssert abs(leftGrid[0, 0, 0] - rightGrid[0, 0, 0]) < Tolerance
    doAssert abs(leftGrid[0, 0, 0] - left) < Tolerance
    doAssert abs(rightGrid[0, 0, 0] - right) < Tolerance

block basic_2d:
  let grid = simplex.grid((0, 0), (19, 17))
  for x in 0 ..< 19:
    for y in 0 ..< 17:
      let value = simplex.value(x, y)
      doAssert value >= -1.0 and value <= 1.0
      checkClose(grid[x, y], value, &"basic 2D index=({x}, {y})")

block basic_3d:
  let grid = simplex.grid((0, 0, 0), (19, 19, 17))
  for x in 0 ..< 19:
    for y in 0 ..< 19:
      for z in 0 ..< 17:
        let value = simplex.value(x, y, z)
        doAssert value >= -1.0 and value <= 1.0
        checkClose(grid[x, y, z], value, &"basic 3D index=({x}, {y}, {z})")

block range_2d:
  var total: float32
  for offset in countup(-1200, 1200-1, 4):
    let simd = simplex.grid((offset, offset), (4, 4))
    for x in 0 ..< 4:
      for y in 0 ..< 4:
        let n = simplex.value(x + offset, y + offset)
        total += n
        checkClose(simd[x, y], n, &"range 2D offset={offset}, index=({x}, {y})")
  echo total

block range_3d:
  var total: float32
  for offset in countup(-1200, 1200-1, 4):
    let simd = simplex.grid((offset, offset, offset), (4, 4, 4))
    for x in 0 ..< 4:
      for y in 0 ..< 4:
        for z in 0 ..< 4:
          let n = simplex.value(x + offset, y + offset, z + offset)
          total += n
          checkClose(
            simd[x, y, z],
            n,
            &"range 3D offset={offset}, index=({x}, {y}, {z})"
          )
  echo total
