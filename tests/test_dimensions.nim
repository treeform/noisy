import
  std/[math, random, strformat],
  noisy

const
  AbsoluteTolerance = 0.0001'f
  RelativeTolerance = 0.00001'f

type Comparison = object
  samples, simdSamples: int
  maxError: float32

proc checkClose(actual, expected: float32, context: string) =
  ## Compares finite results while allowing normal float32 rounding.
  doAssert abs(actual) < Inf and abs(expected) < Inf, context
  let
    error = abs(actual - expected)
    tolerance = AbsoluteTolerance +
      RelativeTolerance * max(abs(actual), abs(expected))
  doAssert error <= tolerance,
    context & &" actual={actual}, expected={expected}, error={error}"

template expectError(error: typedesc, body: untyped) =
  ## Requires the operation to raise the specified exception.
  block:
    var raised = false
    try:
      body
    except error:
      raised = true
    doAssert raised, "Expected " & astToStr(error)

proc compareGrid(
  noise: Simplex,
  start: float32,
  width: int,
  context: string,
  comparison: var Comparison
) =
  ## Checks 1D SIMD lanes, scalar tails, and public storage dimensions.
  let grid = noise.grid(start, width)
  doAssert grid.width == width and grid.height == 1
  doAssert grid.depth == 1 and grid.frames == 1
  doAssert grid.values.len == width
  comparison.simdSamples += width div 4 * 4
  for x in 0 ..< width:
    let
      expected = noise.value(start + x.float32)
      actual = grid[x]
    checkClose(actual, expected, context & &" 1D start={start}, x={x}")
    doAssert grid.values[x] == actual
    inc comparison.samples
    comparison.maxError = max(comparison.maxError, abs(actual - expected))

proc compareGrid(
  noise: Simplex,
  start: (float32, float32, float32, float32),
  dimens: (int, int, int, int),
  context: string,
  comparison: var Comparison
) =
  ## Checks 4D sampling and frame-major storage with z as the fastest axis.
  let grid = noise.grid(start, dimens)
  doAssert grid.width == dimens[0] and grid.height == dimens[1]
  doAssert grid.depth == dimens[2] and grid.frames == dimens[3]
  doAssert grid.values.len == dimens[0] * dimens[1] * dimens[2] * dimens[3]
  comparison.simdSamples += dimens[0] * dimens[1] *
    (dimens[2] div 4 * 4) * dimens[3]
  var index = 0
  for w in 0 ..< dimens[3]:
    for x in 0 ..< dimens[0]:
      for y in 0 ..< dimens[1]:
        for z in 0 ..< dimens[2]:
          let
            expected = noise.value(
              start[0] + x.float32,
              start[1] + y.float32,
              start[2] + z.float32,
              start[3] + w.float32
            )
            actual = grid[x, y, z, w]
          checkClose(
            actual,
            expected,
            context & &" 4D start={start}, dimens={dimens}, " &
              &"index=({x}, {y}, {z}, {w})"
          )
          doAssert grid.values[index] == actual
          inc index
          inc comparison.samples
          comparison.maxError = max(comparison.maxError, abs(actual - expected))

proc describe(noise: Simplex, seed: int): string =
  ## Identifies a configuration when a parity comparison fails.
  &"seed={seed}, octaves={noise.octaves}, frequency={noise.frequency}, " &
    &"amplitude={noise.amplitude}, gain={noise.gain}, " &
    &"lacunarity={noise.lacunarity}"

echo "Testing samples against an independent C reference"
block:
  # Samples come from Stefan Gustavson's public-domain simplexnoise1234.c.
  # Reference: stegu/perlin-noise.
  # Commit: a624f5a18395c0ebc04738191efddf1dbd2cad9d.
  # The C permutation was replaced by 0..255 repeated, and its 1D scale was
  # changed from 0.25 to 0.395. The 4D reference uses support 0.5 and scale 62.
  var noise = initSimplex(1)
  for i in 0 ..< noise.perm.len:
    noise.perm[i] = (i and 255).int32
    noise.permMod12[i] = noise.perm[i] mod 12
  for sample in [
    (0.125'f, 0.0442748033'f), (0.25'f, 0.0545749888'f),
    (0.5'f, -0.062490236'f), (0.75'f, -0.141710535'f),
    (1.25'f, 0.120003492'f), (-0.125'f, -0.0547039844'f),
    (-0.5'f, -0.562412143'f), (-1.25'f, 0.534281552'f),
    (255.75'f, -0.163110182'f), (256.25'f, 0.0545749888'f),
    (-256.25'f, -0.163110182'f)
  ]:
    checkClose(noise.value(sample[0]), sample[1], "1D C reference")
    checkClose(noise.grid(sample[0], 4)[0], sample[1], "1D SIMD C reference")
  for sample in [
    ((0.1'f, 0.2'f, 0.3'f, 0.4'f), -0.148842052'f),
    ((0.125'f, 0.25'f, 0.5'f, 0.75'f), -0.170217291'f),
    ((-0.1'f, -0.2'f, -0.3'f, -0.4'f), 0.102250859'f),
    ((-0.75'f, 0.25'f, -0.125'f, 0.5'f), -0.0745871663'f),
    ((0.25'f, 0.25'f, 0.25'f, 0.25'f), 0.00956536643'f),
    ((-0.25'f, -0.25'f, -0.25'f, -0.25'f), -0.00956536643'f),
    ((1.0'f, 2.0'f, 3.0'f, 4.0'f), -0.148286268'f),
    ((-1.25'f, 2.5'f, -3.75'f, 4.125'f), 0.171521381'f),
    ((255.25'f, -256.5'f, 0.125'f, 0.75'f), -0.00690730661'f),
    ((0.75'f, 0.5'f, 0.25'f, 0.125'f), 0.171521291'f),
    ((0.125'f, 0.75'f, 0.25'f, 0.5'f), 0.0483709686'f),
    ((0.25'f, 0.125'f, 0.75'f, 0.5'f), -0.238764793'f)
  ]:
    let point = sample[0]
    checkClose(
      noise.value(point[0], point[1], point[2], point[3]),
      sample[1],
      "4D C reference"
    )
    checkClose(
      noise.grid(point, (1, 1, 4, 1))[0, 0, 0, 0],
      sample[1],
      "4D SIMD C reference"
    )

echo "Testing native dimensions, determinism, and output range"
block:
  var
    rng = initRand(7419)
    varied1d, varied4d, native1d, fourthAxis: bool
  for seed in [-12345, 0, 1, 42, 1988]:
    let
      noise = initSimplex(seed)
      repeated = initSimplex(seed)
      other = initSimplex(seed + 1)
    doAssert noise.value(0) == 0
    doAssert noise.value(0, 0, 0, 0) == 0
    for i in 0 ..< 256:
      let
        x = (rng.rand(64.0) - 32).float32
        y = (rng.rand(64.0) - 32).float32
        z = (rng.rand(64.0) - 32).float32
        w = (rng.rand(64.0) - 32).float32
        value1d = noise.value(x)
        value4d = noise.value(x, y, z, w)
      doAssert value1d >= -1 and value1d <= 1
      doAssert value4d >= -1 and value4d <= 1
      doAssert value1d == repeated.value(x)
      doAssert value4d == repeated.value(x, y, z, w)
      varied1d = varied1d or value1d != other.value(x)
      varied4d = varied4d or value4d != other.value(x, y, z, w)
      native1d = native1d or value1d != noise.value(x, 0.0'f)
      fourthAxis = fourthAxis or value4d != noise.value(x, y, z, w + 0.25'f)
    for x in -8 .. 8:
      doAssert noise.value(x) == noise.value(x.float32)
      doAssert noise.value(x, 2, -3, 5) ==
        noise.value(x.float32, 2.0'f, -3.0'f, 5.0'f)
    doAssert noise.grid(-3, 9).values == noise.grid(-3.0'f, 9).values
    doAssert noise.grid((-3, 1, -2, 4), (2, 3, 5, 2)).values ==
      noise.grid((-3.0'f, 1.0'f, -2.0'f, 4.0'f), (2, 3, 5, 2)).values
  doAssert varied1d and varied4d and native1d and fourthAxis

echo "Testing continuity at lattice and tied-rank boundaries"
block:
  const
    Epsilon = 0.00001'f
    F4 = ((sqrt(5.0) - 1) / 4).float32
  for seed in [1, 42, 1988]:
    let noise = initSimplex(seed)
    for x in [-257, -256, -1, 0, 1, 255, 256, 257]:
      checkClose(
        noise.value(x.float32 - Epsilon),
        noise.value(x.float32 + Epsilon),
        "1D continuity at " & $x
      )
    for axis in 0 ..< 4:
      var
        left = [0.25'f, 0.25'f, 0.25'f, 0.25'f]
        right = left
      left[axis] -= Epsilon
      right[axis] += Epsilon
      checkClose(
        noise.value(left[0], left[1], left[2], left[3]),
        noise.value(right[0], right[1], right[2], right[3]),
        "4D continuity at tied ranks on axis " & $axis
      )
    let
      y = 0.25'f
      z = -0.125'f
      w = 0.5'f
      # The first skewed lattice coordinate wraps between -1 and 0 here.
      x = -(y + z + w) * F4 / (1 + F4)
    checkClose(
      noise.value(x - Epsilon, y, z, w),
      noise.value(x + Epsilon, y, z, w),
      "4D continuity at permutation wrapping"
    )

echo "Testing octave composition in new dimensions"
block:
  let base = initSimplex(1988)
  var noise = base
  noise.octaves = 4
  noise.frequency = 0.13'f
  noise.amplitude = 1.7'f
  noise.lacunarity = 1.9'f
  noise.gain = 0.6'f
  for point in [
    (0.125'f, -1.75'f, 3.25'f, -0.5'f),
    (-4.37'f, 9.13'f, -2.61'f, 5.19'f)
  ]:
    var
      expected1d, expected4d: float32
      frequency = noise.frequency
      amplitude = noise.amplitude
    for i in 0 ..< noise.octaves:
      expected1d += base.value(point[0] * frequency) * amplitude
      expected4d += base.value(
        point[0] * frequency,
        point[1] * frequency,
        point[2] * frequency,
        point[3] * frequency
      ) * amplitude
      frequency *= noise.lacunarity
      amplitude *= noise.gain
    checkClose(
      noise.value(point[0]),
      expected1d / noise.octaves.float32,
      "1D octave composition"
    )
    checkClose(
      noise.value(point[0], point[1], point[2], point[3]),
      expected4d / noise.octaves.float32,
      "4D octave composition"
    )

echo "Testing grid bounds, empty dimensions, and invalid input"
block:
  let
    noise = initSimplex(42)
    line = noise.grid(-0.25'f, 5)
    volume = noise.grid((0, 0, 0, 0), (2, 3, 5, 2))
    old2d = noise.grid((0, 0), (5, 7))
    old3d = noise.grid((0, 0, 0), (5, 7, 9))
  doAssert old2d.frames == 1 and old3d.frames == 1
  doAssert old2d[2, 3] == old2d[2, 3, 0, 0]
  doAssert old3d[2, 3, 4] == old3d[2, 3, 4, 0]
  for x in [-1, 5]:
    expectError(IndexDefect):
      discard line[x]
  for point in [
    (-1, 0, 0, 0), (2, 0, 0, 0),
    (0, -1, 0, 0), (0, 3, 0, 0),
    (0, 0, -1, 0), (0, 0, 5, 0),
    (0, 0, 0, -1), (0, 0, 0, 2)
  ]:
    expectError(IndexDefect):
      discard volume[point[0], point[1], point[2], point[3]]
  doAssert noise.grid(0, 0).values.len == 0
  for dimens in [(0, 2, 5, 3), (2, 0, 5, 3), (2, 3, 0, 3), (2, 3, 5, 0)]:
    let empty = noise.grid((0, 0, 0, 0), dimens)
    doAssert empty.values.len == 0
    doAssert (empty.width, empty.height, empty.depth, empty.frames) == dimens
  expectError(NoisyError):
    discard noise.grid(0, -1)
  expectError(NoisyError):
    discard noise.grid(0, int.high)
  for dimens in [
    (-1, 2, 3, 4), (1, -2, 3, 4), (1, 2, -3, 4), (1, 2, 3, -4),
    (int.high, 2, 2, 2)
  ]:
    expectError(NoisyError):
      discard noise.grid((0, 0, 0, 0), dimens)
  for octaves in [-1, 0]:
    var invalid = noise
    invalid.octaves = octaves
    expectError(NoisyError):
      discard invalid.value(0.25'f)
    expectError(NoisyError):
      discard invalid.value(0.25'f, 0.5'f, 0.75'f, 1.25'f)
    expectError(NoisyError):
      discard invalid.grid(0, 5)
    expectError(NoisyError):
      discard invalid.grid((0, 0, 0, 0), (1, 1, 5, 1))

var comparison1d, comparison4d: Comparison

echo "Testing 1D/4D scalar and SIMD parity across fixed configurations"
block:
  for seed in [-12345, 0, 1, 42, 1988]:
    for settings in [
      (1, 1.0'f, 1.0'f, 2.0'f, 0.5'f),
      (1, 0.1'f, 0.7'f, 2.0'f, 0.5'f),
      (3, 0.07'f, 1.8'f, 1.9'f, 0.65'f),
      (5, 0.013'f, 0.8'f, 2.3'f, 0.45'f),
      (3, -0.25'f, -1.5'f, 1.5'f, -0.4'f),
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
        (0.0'f, 0.0'f, 0.0'f, 0.0'f),
        (-0.125'f, -0.75'f, -0.25'f, -0.5'f),
        (-5.37'f, 9.13'f, -2.61'f, 5.19'f),
        (-16.0'f, 0.0'f, 16.0'f, -8.0'f),
        (-257.25'f, 255.125'f, -128.75'f, 256.125'f),
        (0.000001'f, -0.000001'f, 0.999999'f, -0.999999'f)
      ]:
        for width in [1, 3, 4, 5, 9, 17]:
          compareGrid(noise, start[0], width, context, comparison1d)
        for dimens in [
          (1, 1, 1, 1), (2, 3, 3, 2), (2, 3, 4, 2),
          (3, 2, 5, 4), (2, 3, 9, 3)
        ]:
          compareGrid(noise, start, dimens, context, comparison4d)

echo "Testing 1D/4D parity across reproducible random configurations"
block:
  var rng = initRand(20261118)
  for i in 0 ..< 96:
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
        (rng.rand(128.0) - 64).float32,
        (rng.rand(128.0) - 64).float32
      )
      dimens = (rng.rand(1 .. 4), rng.rand(1 .. 4), rng.rand(4 .. 11),
        rng.rand(1 .. 4))
      context = &"random case={i}, " & describe(noise, seed)
    compareGrid(noise, start[0], rng.rand(4 .. 31), context, comparison1d)
    compareGrid(noise, start, dimens, context, comparison4d)

doAssert comparison1d.simdSamples > 0 and comparison4d.simdSamples > 0
echo &"1D: {comparison1d.samples} comparisons, " &
  &"{comparison1d.simdSamples} SIMD samples, max error={comparison1d.maxError}"
echo &"4D: {comparison4d.samples} comparisons, " &
  &"{comparison4d.simdSamples} SIMD samples, max error={comparison4d.maxError}"
