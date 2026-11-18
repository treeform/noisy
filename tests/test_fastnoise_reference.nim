import
  std/[random, strformat],
  noisy,
  reference/fastnoise

const
  AbsoluteTolerance = 0.00001'f
  RelativeTolerance = 0.00001'f
  Seeds = [0, 1, -1, 1988, -2147483647 - 1, 2147483647]

type
  Point = array[3, float32]
  Comparison = object
    samples: int
    maxError: float32

proc algorithm[T](noise: T): int32 =
  ## Maps native types to the independent wrapper's explicit algorithm IDs.
  when T is OpenSimplex2:
    0
  elif T is ValueNoise:
    1
  else:
    2

proc referenceValue[T](noise: T, point: Point, dimensions: int): float32 =
  ## Composes independent single-octave samples with Noisy's layering rules.
  var
    amplitude = noise.amplitude
    frequency = noise.frequency
    distance = 1'i32
    output = 1'i32
    jitter = 1.0'f
  when T is Cellular:
    distance = noise.distanceFunction.ord.int32
    output = noise.returnType.ord.int32
    jitter = noise.jitter
  for octave in 0 ..< noise.octaves:
    let
      x = point[0] * frequency
      y = point[1] * frequency
      z = point[2] * frequency
      sample =
        if dimensions == 2:
          referenceFastNoise2(noise.algorithm(), cast[int32](noise.seed),
            x, y, distance, output, jitter)
        else:
          referenceFastNoise3(noise.algorithm(), cast[int32](noise.seed),
            x, y, z, distance, output, jitter)
    result += sample * amplitude
    amplitude *= noise.gain
    frequency *= noise.lacunarity
  result /= noise.octaves.float32

proc compare[T](
  noise: T,
  point: Point,
  dimensions: int,
  comparison: var Comparison
) =
  ## Checks the public native sampler against the untouched C++ implementation.
  let
    actual =
      if dimensions == 2:
        noise.value(point[0], point[1])
      else:
        noise.value(point[0], point[1], point[2])
    expected = noise.referenceValue(point, dimensions)
    error = abs(actual - expected)
    tolerance = AbsoluteTolerance +
      RelativeTolerance * max(abs(actual), abs(expected))
  var context = $T & &" {dimensions}D seed={noise.seed}, point={point}, " &
    &"octaves={noise.octaves}, frequency={noise.frequency}"
  when T is Cellular:
    context.add(&" distance={noise.distanceFunction}, output={noise.returnType}, " &
      &"jitter={noise.jitter}")
  doAssert abs(actual) < Inf and abs(expected) < Inf,
    context & " produced a nonfinite sample"
  when T is Cellular:
    if noise.returnType == CellValue and noise.octaves == 1:
      doAssert actual == expected,
        context & &" selected another cell: actual={actual}, expected={expected}"
  doAssert error <= tolerance,
    context & &" actual={actual}, expected={expected}, error={error}"
  inc comparison.samples
  comparison.maxError = max(comparison.maxError, error)

proc compareBasic[T](noise: T, comparisons: var array[2, Comparison]) =
  ## Covers cell boundaries, negative coordinates, seed wrapping, and random inputs.
  var rng = initRand(48217)
  for dimensions in 2 .. 3:
    for coordinate in [-1024.0'f, -256, -2, -1, -0.5, 0,
      0.5, 1, 2, 255, 256, 1024]:
      for offset in [-0.0001'f, 0.0'f, 0.0001'f, 0.125'f, 0.5'f]:
        for axis in 0 ..< dimensions:
          var point = [0.375'f, -0.625'f, 1.25'f]
          point[axis] = coordinate + offset
          noise.compare(point, dimensions, comparisons[dimensions - 2])
        let point = [coordinate + offset, coordinate + offset,
          coordinate + offset]
        noise.compare(point, dimensions, comparisons[dimensions - 2])
    for sample in 0 ..< 1000:
      var point: Point
      for axis in 0 ..< dimensions:
        point[axis] = rng.rand(-256.0 .. 256.0).float32
      noise.compare(point, dimensions, comparisons[dimensions - 2])

proc compareOctaves[T](initial: T, comparisons: var array[2, Comparison]) =
  ## Retains one seed while layering the reference, rather than using FNL fBm.
  var
    noise = initial
    rng = initRand(341719)
  for settings in [
    (1, 0.75'f, 0.0'f, 2.0'f, 0.5'f),
    (2, 0.5'f, 0.125'f, 2.0'f, 0.5'f),
    (3, 1.25'f, 0.37'f, 1.7'f, 0.6'f),
    (4, -0.75'f, -0.25'f, 1.5'f, -0.4'f)
  ]:
    (noise.octaves, noise.amplitude, noise.frequency, noise.lacunarity,
      noise.gain) = settings
    for dimensions in 2 .. 3:
      for sample in 0 ..< 200:
        var point: Point
        for axis in 0 ..< dimensions:
          point[axis] = rng.rand(-32.0 .. 32.0).float32
        noise.compare(point, dimensions, comparisons[dimensions - 2])

var
  openComparisons, valueComparisons, cellularComparisons: array[2, Comparison]
  rng = initRand(191237)

echo "Testing OpenSimplex2, Value, and Cellular against FastNoiseLite C++"
for seed in Seeds:
  initOpenSimplex2(seed).compareBasic(openComparisons)
  initValue(seed).compareBasic(valueComparisons)
  initCellular(seed).compareBasic(cellularComparisons)
  initOpenSimplex2(seed).compareOctaves(openComparisons)
  initValue(seed).compareOctaves(valueComparisons)
  initCellular(seed).compareOctaves(cellularComparisons)

echo "Testing every Cellular distance, return mode, jitter, and exact ties"
for seed in Seeds:
  var noise = initCellular(seed)
  for distance in low(CellularDistance) .. high(CellularDistance):
    noise.distanceFunction = distance
    for output in low(CellularReturn) .. high(CellularReturn):
      noise.returnType = output
      for jitter in [0.0'f, 0.5'f, 1.0'f]:
        noise.jitter = jitter
        for dimensions in 2 .. 3:
          for point in [
            [0.0'f, 0.0'f, 0.0'f],
            [0.5'f, 0.5'f, 0.5'f],
            [-0.5'f, -0.5'f, -0.5'f],
            [1.5'f, -2.5'f, 3.5'f],
            [-256.5'f, 255.5'f, -255.5'f]
          ]:
            noise.compare(point, dimensions, cellularComparisons[dimensions - 2])
          for sample in 0 ..< 64:
            var point: Point
            for axis in 0 ..< dimensions:
              point[axis] = rng.rand(-64.0 .. 64.0).float32
            noise.compare(point, dimensions, cellularComparisons[dimensions - 2])

for index in 0 ..< 2:
  echo &"OpenSimplex2 {index + 2}D: {openComparisons[index].samples} " &
    &"reference comparisons, max error={openComparisons[index].maxError}"
  echo &"Value {index + 2}D: {valueComparisons[index].samples} " &
    &"reference comparisons, max error={valueComparisons[index].maxError}"
  echo &"Cellular {index + 2}D: {cellularComparisons[index].samples} " &
    &"reference comparisons, max error={cellularComparisons[index].maxError}"
