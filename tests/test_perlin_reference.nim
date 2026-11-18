import
  std/[random, strformat],
  noisy,
  reference/bindings

const
  AbsoluteTolerance = 0.00001'f
  RelativeTolerance = 0.00001'f

type
  Point = array[4, float32]
  Comparison = object
    samples: int
    maxError: float32

proc referenceValue(noise: Perlin, point: Point, dimensions: int): float32 =
  ## Layers unchanged C kernels using the public octave convention.
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
        of 1: referencePerlin1(x)
        of 2: referencePerlin2(x, y)
        of 3: referencePerlin3(x, y, z)
        else: referencePerlin4(x, y, z, w)
    result += sample * amplitude
    amplitude *= noise.gain
    frequency *= noise.lacunarity
  result /= noise.octaves.float32

proc scalarValue(noise: Perlin, point: Point, dimensions: int): float32 =
  ## Samples the public scalar overload for the requested dimension.
  case dimensions
  of 1: noise.value(point[0])
  of 2: noise.value(point[0], point[1])
  of 3: noise.value(point[0], point[1], point[2])
  else: noise.value(point[0], point[1], point[2], point[3])

proc compare(
  noise: Perlin,
  point: Point,
  dimensions: int,
  comparison: var Comparison
) =
  ## Rejects nonfinite outputs and discrepancies beyond float32 rounding.
  let
    actual = noise.scalarValue(point, dimensions)
    expected = noise.referenceValue(point, dimensions)
    error = abs(actual - expected)
    tolerance = AbsoluteTolerance +
      RelativeTolerance * max(abs(actual), abs(expected))
  doAssert abs(actual) < Inf and abs(expected) < Inf,
    &"Perlin {dimensions}D nonfinite result at {point}"
  doAssert error <= tolerance,
    &"Perlin {dimensions}D at {point}: actual={actual}, " &
      &"reference={expected}, error={error}, tolerance={tolerance}, " &
      &"octaves={noise.octaves}, frequency={noise.frequency}"
  inc comparison.samples
  comparison.maxError = max(comparison.maxError, error)

var
  noise = initPerlin(1988)
  comparisons: array[4, Comparison]
  rng = initRand(734911)

# Match the actual upstream table without changing any reference source.
for i in 0 ..< noise.perm.len:
  noise.perm[i] = referencePerlinPermutation(i.int32)
for i in 0 ..< 256:
  doAssert noise.perm[i] == noise.perm[i + 256]
  doAssert noise.perm[i] == referencePermutation(i.int32)

echo "Testing improved Perlin against the independent C reference"
for dimensions in 1 .. 4:
  # Exact lattice points and fractional planes expose hashing and fade errors.
  # Sample both sides of negative cells and the 256-entry permutation seam.
  for coordinate in [-1024.0'f, -512, -257, -256, -255, -2, -1, 0,
    1, 2, 255, 256, 257, 512, 1024]:
    for offset in [-0.0001'f, 0.0'f, 0.0001'f, 0.125'f, 0.5'f, 0.875'f]:
      for axis in 0 ..< dimensions:
        var point = [0.375'f, -0.625'f, 1.25'f, -2.75'f]
        point[axis] = coordinate + offset
        noise.compare(point, dimensions, comparisons[dimensions - 1])
      let point = [coordinate + offset, coordinate + offset,
        coordinate + offset, coordinate + offset]
      noise.compare(point, dimensions, comparisons[dimensions - 1])
  # Cover every hash-table entry with a nonzero fractional sample.
  for lattice in -256 ..< 256:
    let point = [lattice.float32 + 0.3125'f,
      (lattice * 7).float32 + 0.5625'f,
      (lattice * 13).float32 + 0.8125'f,
      (lattice * 17).float32 + 0.0625'f]
    noise.compare(point, dimensions, comparisons[dimensions - 1])
  for sample in 0 ..< 4000:
    var point: Point
    for axis in 0 ..< dimensions:
      point[axis] = rng.rand(-1024.0 .. 1024.0).float32
    noise.compare(point, dimensions, comparisons[dimensions - 1])

echo "Testing improved Perlin reference octave composition"
for settings in [
  (1, 0.75'f, 0.0'f, 2.0'f, 0.5'f),
  (2, 0.5'f, 0.125'f, 2.0'f, 0.5'f),
  (3, 1.25'f, 0.37'f, 1.7'f, 0.6'f),
  (4, -0.75'f, -0.25'f, 1.5'f, -0.4'f),
  (5, 0.8'f, 1.3'f, 0.75'f, 0.9'f)
]:
  (noise.octaves, noise.amplitude, noise.frequency, noise.lacunarity,
    noise.gain) = settings
  for dimensions in 1 .. 4:
    for sample in 0 ..< 1000:
      var point: Point
      for axis in 0 ..< dimensions:
        point[axis] = rng.rand(-128.0 .. 128.0).float32
      noise.compare(point, dimensions, comparisons[dimensions - 1])

for index, comparison in comparisons:
  echo &"Perlin {index + 1}D: {comparison.samples} reference comparisons, " &
    &"max error={comparison.maxError}"
