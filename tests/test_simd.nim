import
  std/[math, random, strformat],
  noisy

const
  # Allow float32 rounding without accepting visible differences in the noise.
  AbsoluteTolerance = 0.0001'f
  RelativeTolerance = 0.00001'f

type Comparison = object
  samples, simdSamples: int
  maxError: float32

proc checkSample(
  actual, expected: float32,
  context: string,
  comparison: var Comparison
) =
  ## Compares finite samples with an absolute and relative error allowance.
  doAssert abs(actual) < Inf and abs(expected) < Inf,
    context & &" non-finite sample: SIMD={actual}, scalar={expected}"
  let
    error = abs(actual - expected)
    tolerance = AbsoluteTolerance +
      RelativeTolerance * max(abs(actual), abs(expected))
  doAssert error <= tolerance,
    context & &" SIMD={actual}, scalar={expected}, " &
    &"error={error}, tolerance={tolerance}"
  inc comparison.samples
  comparison.maxError = max(comparison.maxError, error)

proc compareGrid(
  noise: Simplex,
  start: (float32, float32),
  dimens: (int, int),
  context: string,
  comparison: var Comparison
) =
  ## Compares all SIMD lanes and scalar tails against point sampling in 2D.
  let grid = noise.grid(start, dimens)
  doAssert grid.width == dimens[0] and grid.height == dimens[1]
  doAssert grid.depth == 1 and grid.values.len == dimens[0] * dimens[1]
  comparison.simdSamples += (dimens[0] div 4 * 4) * (dimens[1] div 4 * 4)
  for x in 0 ..< dimens[0]:
    for y in 0 ..< dimens[1]:
      let expected = noise.value(start[0] + x.float32, start[1] + y.float32)
      checkSample(
        grid[x, y],
        expected,
        context & &" 2D start={start}, dimens={dimens}, index=({x}, {y})",
        comparison
      )

proc compareGrid(
  noise: Simplex,
  start: (float32, float32, float32),
  dimens: (int, int, int),
  context: string,
  comparison: var Comparison
) =
  ## Compares all SIMD lanes and scalar tails against point sampling in 3D.
  let grid = noise.grid(start, dimens)
  doAssert grid.width == dimens[0] and grid.height == dimens[1]
  doAssert grid.depth == dimens[2]
  doAssert grid.values.len == dimens[0] * dimens[1] * dimens[2]
  comparison.simdSamples += (dimens[0] div 4 * 4) *
    (dimens[1] div 4 * 4) * (dimens[2] div 4 * 4)
  for x in 0 ..< dimens[0]:
    for y in 0 ..< dimens[1]:
      for z in 0 ..< dimens[2]:
        let expected = noise.value(
          start[0] + x.float32,
          start[1] + y.float32,
          start[2] + z.float32
        )
        checkSample(
          grid[x, y, z],
          expected,
          context & &" 3D start={start}, dimens={dimens}, " &
            &"index=({x}, {y}, {z})",
          comparison
        )

proc describe(noise: Simplex, seed: int): string =
  ## Identifies the full noise configuration in a failed comparison.
  &"seed={seed}, octaves={noise.octaves}, frequency={noise.frequency}, " &
    &"amplitude={noise.amplitude}, gain={noise.gain}, " &
    &"lacunarity={noise.lacunarity}"

var comparison2d, comparison3d: Comparison

echo "Testing scalar/SIMD parity across fixed configurations"
block:
  for seed in [-12345, 0, 1, 42, 1988]:
    for settings in [
      (1, 1.0'f, 1.0'f, 2.0'f, 0.5'f),
      (1, 0.1'f, 0.7'f, 2.0'f, 0.5'f),
      (3, 0.07'f, 1.8'f, 1.9'f, 0.65'f),
      (5, 0.013'f, 0.8'f, 2.3'f, 0.45'f),
      (3, -0.25'f, -1.5'f, 1.5'f, -0.4'f),
      (3, 0.3'f, 0.2'f, 1.5'f, 1.2'f),
      (2, 0.0'f, 1.0'f, 2.0'f, 0.5'f),
      (2, 0.1'f, 0.0'f, 2.0'f, 0.5'f)
    ]:
      var noise = initSimplex(seed)
      noise.octaves = settings[0]
      noise.frequency = settings[1]
      noise.amplitude = settings[2]
      noise.lacunarity = settings[3]
      noise.gain = settings[4]
      let context = describe(noise, seed)
      for start in [
        (0.0'f, 0.0'f, 0.0'f),
        (-0.125'f, -0.75'f, -0.25'f),
        (-5.37'f, 9.13'f, -2.61'f),
        (-16.0'f, 0.0'f, 16.0'f),
        (-257.25'f, 255.125'f, -128.75'f),
        (0.000001'f, -0.000001'f, 0.999999'f)
      ]:
        for dimens in [
          (1, 1), (3, 7), (7, 3), (4, 4),
          (5, 7), (8, 8), (9, 11), (11, 9)
        ]:
          compareGrid(
            noise,
            (start[0], start[1]),
            dimens,
            context,
            comparison2d
          )
        for dimens in [
          (1, 1, 1), (3, 5, 7), (5, 3, 7), (5, 7, 3),
          (4, 4, 4), (5, 7, 9), (8, 8, 8), (9, 10, 11)
        ]:
          compareGrid(noise, start, dimens, context, comparison3d)

echo "Testing lattice boundaries in every SIMD lane and scalar tails"
block:
  const F2 = ((sqrt(3.0) - 1) / 2).float32
  for seed in [1, 42, 1988]:
    let noise = initSimplex(seed)
    for boundary in [-256.0'f, -1.0'f, 0.0'f, 255.0'f]:
      for delta in [-0.0001'f, 0.0'f, 0.0001'f]:
        for lane in 0 .. 4:
          let
            x = -0.125'f
            y = 0.25'f
            # Place a skewed lattice boundary in each lane, then in the tail.
            boundary2d = (boundary - x * F2) / (1 + F2)
            boundary3d = (3 * boundary - x - y) / 4
            context = &"seed={seed}, boundary={boundary}, delta={delta}, " &
              &"lane={lane}"
          compareGrid(
            noise,
            (x, boundary2d + delta - lane.float32),
            (4, 5),
            context,
            comparison2d
          )
          compareGrid(
            noise,
            (x, y, boundary3d + delta - lane.float32),
            (4, 4, 5),
            context,
            comparison3d
          )

echo "Testing scalar/SIMD parity across reproducible random configurations"
block:
  var rng = initRand(20261117)
  for i in 0 ..< 128:
    let seed = rng.rand(-1_000_000 .. 1_000_000)
    var noise = initSimplex(seed)
    noise.octaves = rng.rand(1 .. 5)
    noise.frequency = (0.01 + rng.rand(0.49)).float32
    noise.amplitude = (0.1 + rng.rand(1.9)).float32
    noise.lacunarity = (1.1 + rng.rand(1.3)).float32
    noise.gain = (0.2 + rng.rand(0.6)).float32
    let
      start = (
        (rng.rand(128.0) - 64).float32,
        (rng.rand(128.0) - 64).float32,
        (rng.rand(128.0) - 64).float32
      )
      dimens = (rng.rand(4 .. 11), rng.rand(4 .. 11), rng.rand(4 .. 11))
      context = &"random case={i}, " & describe(noise, seed)
    compareGrid(
      noise,
      (start[0], start[1]),
      (dimens[0], dimens[1]),
      context,
      comparison2d
    )
    compareGrid(noise, start, dimens, context, comparison3d)

doAssert comparison2d.simdSamples > 0 and comparison3d.simdSamples > 0
echo &"2D: {comparison2d.samples} comparisons, " &
  &"{comparison2d.simdSamples} SIMD samples, max error={comparison2d.maxError}"
echo &"3D: {comparison3d.samples} comparisons, " &
  &"{comparison3d.simdSamples} SIMD samples, max error={comparison3d.maxError}"
