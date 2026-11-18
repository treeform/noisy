import
  std/[math, os, random, strformat],
  noisy

const
  ReferenceDir = currentSourcePath().parentDir / "reference"
  NoisyReferenceDir {.strdefine.} = ""
  AbsoluteTolerance = 0.0001'f
  RelativeTolerance = 0.00001'f

when NoisyReferenceDir.len == 0:
  {.error: "Fetch references with tools/gen_references.py, then pass -d:noisyReferenceDir=/cache/path".}
static:
  doAssert fileExists(NoisyReferenceDir / "caseman/simplex23.c"),
    "Reference cache missing; run python tools/gen_references.py --cache-dir /cache/path"
{.passC: "-I" & quoteShell(NoisyReferenceDir).}
{.compile: ReferenceDir / "reference_noise.cpp".}
{.push cdecl, importc, header: ReferenceDir / "reference_noise.h".}
proc referenceNoise1(x: float32): float32
  ## Samples the unchanged Rombauts 1D kernel.
proc referenceNoise2(x, y: float32): float32
  ## Samples the unchanged Duncan 2D kernel.
proc referenceNoise3(x, y, z: float32): float32
  ## Samples the unchanged Duncan 3D kernel.
proc referenceNoise4(x, y, z, w: float32): float32
  ## Samples the unchanged Gustavson 4D kernel.
proc referencePermutation(index: int32): int32
  ## Reads the original permutation shared by the references.
{.pop.}

type
  Point = array[4, float32]
  Comparison = object
    scalarSamples, gridSamples, simdSamples: int
    maxScalarError, maxGridError: float32

proc referenceValue(noise: Simplex, point: Point, dimensions: int): float32 =
  ## Layers unchanged upstream kernels using Noisy's octave convention.
  var
    amplitude = noise.amplitude
    frequency = noise.frequency
  for octave in 0 ..< noise.octaves:
    let
      x = point[0] * frequency
      y = point[1] * frequency
      z = point[2] * frequency
      w = point[3] * frequency
      sample = case dimensions
        of 1: referenceNoise1(x)
        of 2: referenceNoise2(x, y)
        of 3: referenceNoise3(x, y, z)
        else: referenceNoise4(x, y, z, w)
    result += sample * amplitude
    amplitude *= noise.gain
    frequency *= noise.lacunarity
  result /= noise.octaves.float32

proc scalarValue(noise: Simplex, point: Point, dimensions: int): float32 =
  ## Samples the public scalar overload for the selected dimension.
  case dimensions
  of 1: noise.value(point[0])
  of 2: noise.value(point[0], point[1])
  of 3: noise.value(point[0], point[1], point[2])
  else: noise.value(point[0], point[1], point[2], point[3])

proc checkClose(
  actual, expected: float32,
  point: Point,
  dimensions: int,
  context: string,
  coordinateAllowance = 0.0'f
): float32 =
  ## Rejects nonfinite samples and differences beyond float32 rounding.
  doAssert abs(actual) < Inf and abs(expected) < Inf,
    &"{context}: {dimensions}D nonfinite sample at {point}"
  result = abs(actual - expected)
  doAssert coordinateAllowance == 0 or dimensions == 4
  let tolerance = AbsoluteTolerance +
    RelativeTolerance * max(abs(actual), abs(expected)) + coordinateAllowance
  doAssert result <= tolerance,
    &"{context}: {dimensions}D at {point}, actual={actual}, " &
      &"reference={expected}, error={result}, tolerance={tolerance}"

proc referenceRounding(point: Point): float32 =
  ## Allows for float32 coordinate resolution in the separate 4D stress test.
  const Float32Epsilon = 1.1920928955078125e-7'f
  for coordinate in point:
    result = max(result, abs(coordinate))
  result *= 8 * Float32Epsilon

proc comparePoint(
  noise: Simplex,
  point: Point,
  dimensions: int,
  comparison: var Comparison,
  largeCoordinates = false
) =
  ## Compares point sampling directly to an independent compiled reference.
  let error = checkClose(
    noise.scalarValue(point, dimensions),
    noise.referenceValue(point, dimensions),
    point,
    dimensions,
    "Scalar reference",
    if largeCoordinates: point.referenceRounding() else: 0.0'f
  )
  inc comparison.scalarSamples
  comparison.maxScalarError = max(comparison.maxScalarError, error)

proc compareGrid(
  noise: Simplex,
  start: Point,
  dimensions: int,
  comparison: var Comparison,
  largeCoordinates = false
) =
  ## Checks every vector lane and leftover edge directly against C/C++.
  let grid = case dimensions
    of 1: noise.grid(start[0], 9)
    of 2: noise.grid((start[0], start[1]), (5, 7))
    of 3: noise.grid((start[0], start[1], start[2]), (5, 5, 7))
    else: noise.grid((start[0], start[1], start[2], start[3]), (3, 2, 7, 2))
  var index = 0
  for w in 0 ..< grid.frames:
    for x in 0 ..< grid.width:
      for y in 0 ..< grid.height:
        for z in 0 ..< grid.depth:
          let
            point = [start[0] + x.float32, start[1] + y.float32,
              start[2] + z.float32, start[3] + w.float32]
            expected = noise.referenceValue(point, dimensions)
            allowance =
              if largeCoordinates: point.referenceRounding()
              else: 0.0'f
            scalarError = checkClose(
              noise.scalarValue(point, dimensions), expected,
              point, dimensions, "Scalar grid reference", allowance
            )
            gridError = checkClose(
              grid.values[index], expected,
              point, dimensions, "SIMD grid reference", allowance
            )
            vectorSample = case dimensions
              of 1: x < 8
              of 2: x < 4 and y < 4
              of 3: x < 4 and y < 4 and z < 4
              else: z < 4
          inc index
          inc comparison.scalarSamples
          inc comparison.gridSamples
          comparison.simdSamples += vectorSample.int
          comparison.maxScalarError = max(comparison.maxScalarError, scalarError)
          comparison.maxGridError = max(comparison.maxGridError, gridError)
  doAssert index == grid.values.len

var
  noise = initSimplex(1988)
  comparisons: array[4, Comparison]
  largeCoordinates: Comparison

# Match the published tables, not seed numbers from unrelated RNGs.
# Upstream kernels, constants, and normalization remain unchanged.
# Sources, pinned revisions, and licenses are in reference/README.md.
for i in 0 ..< noise.perm.len:
  noise.perm[i] = referencePermutation(i.int32)
  noise.permMod12[i] = noise.perm[i] mod 12
for i in 0 ..< 256:
  doAssert noise.perm[i] == noise.perm[i + 256]

echo "Testing C/C++ references at integer, wrapping, and tied-rank boundaries"
for dimensions in 1 .. 4:
  for code in 0 ..< 256:
    const Values = [-0.75'f, -0.25'f, 0.25'f, 0.75'f]
    var
      point: Point
      remaining = code
    for axis in 0 ..< dimensions:
      point[axis] = Values[remaining and 3]
      remaining = remaining shr 2
    noise.comparePoint(point, dimensions, comparisons[dimensions - 1])
  for boundary in [-257, -256, -255, -1, 0, 1, 255, 256, 257]:
    for delta in [-0.0001'f, 0.0'f, 0.0001'f]:
      var point: Point
      for axis in 0 ..< dimensions:
        point[axis] = boundary.float32 + delta
      if dimensions == 4 and abs(boundary) > 16:
        noise.comparePoint(point, dimensions, largeCoordinates, true)
        noise.compareGrid(point, dimensions, largeCoordinates, true)
      else:
        noise.comparePoint(point, dimensions, comparisons[dimensions - 1])
        noise.compareGrid(point, dimensions, comparisons[dimensions - 1])
  let unskew = case dimensions
    of 1: 0.0'f
    of 2: ((3 - sqrt(3.0)) / 6).float32
    of 3: (1.0 / 6.0).float32
    else: ((5 - sqrt(5.0)) / 20).float32
  for cell in [-256, -1, 0, 1, 255, 256]:
    # Invert the skew transform to sample either side of cell boundaries.
    var lattice: Point
    lattice[0] = cell.float32
    for axis in 0 ..< dimensions:
      lattice[axis] -= cell.float32 * unskew
    for axis in 0 ..< dimensions:
      for delta in [-0.0001'f, 0.0'f, 0.0001'f]:
        var point = lattice
        point[axis] += delta
        if dimensions == 4 and abs(cell) > 16:
          noise.comparePoint(point, dimensions, largeCoordinates, true)
          noise.compareGrid(point, dimensions, largeCoordinates, true)
        else:
          noise.comparePoint(point, dimensions, comparisons[dimensions - 1])
          noise.compareGrid(point, dimensions, comparisons[dimensions - 1])

echo "Testing C/C++ references across random coordinates and octave settings"
block:
  var rng = initRand(20261118)
  for settings in [
    (1, 1.0'f, 1.0'f, 2.0'f, 0.5'f),
    (1, 0.13'f, 1.7'f, 2.0'f, 0.5'f),
    (3, 0.07'f, 1.8'f, 1.9'f, 0.65'f),
    (4, -0.25'f, -1.5'f, 1.5'f, -0.4'f),
    (3, 0.0'f, 1.0'f, 2.0'f, 0.5'f)
  ]:
    noise.octaves = settings[0]
    noise.frequency = settings[1]
    noise.amplitude = settings[2]
    noise.lacunarity = settings[3]
    noise.gain = settings[4]
    # The 4D C reference promotes its skew constants to double precision.
    # Keep dense random tests local so their float32 error stays bounded.
    # The explicit tests above still cover permutation wrapping at +/-256.
    for sample in 0 ..< 1024:
      var point: Point
      for axis in 0 ..< 4:
        point[axis] = (rng.rand(32.0) - 16).float32
      for dimensions in 1 .. 4:
        noise.comparePoint(point, dimensions, comparisons[dimensions - 1])
        if sample < 16:
          noise.compareGrid(point, dimensions, comparisons[dimensions - 1])

for i, comparison in comparisons:
  doAssert comparison.scalarSamples > 5000
  doAssert comparison.gridSamples > comparison.simdSamples
  doAssert comparison.simdSamples > 0
  echo &"{i + 1}D C/C++: {comparison.scalarSamples} scalar samples, " &
    &"{comparison.gridSamples} grid samples ({comparison.simdSamples} SIMD), " &
    &"max scalar error={comparison.maxScalarError}, " &
    &"max grid error={comparison.maxGridError}"

echo "Testing 4D reference rounding at larger coordinates"
block:
  noise.octaves = 1
  noise.frequency = 1
  noise.amplitude = 1
  var rng = initRand(193982)
  for sample in 0 ..< 4096:
    var point: Point
    for axis in 0 ..< 4:
      point[axis] = (rng.rand(512.0) - 256).float32
    noise.comparePoint(point, 4, largeCoordinates, true)
    if sample < 32:
      noise.compareGrid(point, 4, largeCoordinates, true)
  noise.compareGrid(
    [208.81325'f, -240.0477'f, 170.85532'f, -242.97598'f],
    4, largeCoordinates, true
  )
  echo &"4D coordinate rounding: {largeCoordinates.scalarSamples} scalar, " &
    &"{largeCoordinates.gridSamples} grid ({largeCoordinates.simdSamples} SIMD), " &
    &"max scalar error={largeCoordinates.maxScalarError}, " &
    &"max grid error={largeCoordinates.maxGridError}"
