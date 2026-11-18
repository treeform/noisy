import
  std/[random, strformat],
  noisy/opensimplexes

type Point = array[3, float32]

proc sample(noise: OpenSimplex2, point: Point, dimensions: int): float32 =
  ## Calls the native scalar overload for the requested dimension.
  if dimensions == 2:
    noise.value(point[0], point[1])
  else:
    noise.value(point[0], point[1], point[2])

proc checkClose(actual, expected: float32, context: string) =
  ## Allows float32 rounding while rejecting nonfinite outputs.
  doAssert abs(actual) < Inf and abs(expected) < Inf, context
  let tolerance = 0.00001'f + 0.00001'f * max(abs(actual), abs(expected))
  doAssert abs(actual - expected) <= tolerance,
    context & &" actual={actual}, expected={expected}"

echo "Testing OpenSimplex2 against recorded C++ outputs"
block:
  # Generated with FastNoiseLite C++ at revision
  # 785f37a9ad76e283586a379675085f2063ae03f7.
  # Frequency is 1 and rotation uses the default orientation.
  # The unchanged reference was compiled with explicit signed wrapping.
  for recorded in [
    (0, [0.125'f, 0.25'f, 0.5'f], [0.895951331'f, -0.0898984745'f]),
    (0, [-0.75'f, 1.25'f, -2.5'f], [-0.732237697'f, -0.415868044'f]),
    (0, [255.875'f, -256.25'f, 511.5'f], [-0.870669186'f, 0.924650908'f]),
    (0, [0.0'f, 0.0'f, 0.0'f], [0.0'f, 0.0'f]),
    # Detect rounding the skew constant from an unrounded sqrt(3) literal.
    (0, [-1024.0001'f, -0.625'f, 1.25'f],
      [0.880473554'f, 0.611563087'f]),
    (-7, [0.125'f, 0.25'f, 0.5'f], [0.108032294'f, -0.212609991'f]),
    (-7, [-0.75'f, 1.25'f, -2.5'f], [-0.709360778'f, -0.476102263'f]),
    (-7, [255.875'f, -256.25'f, 511.5'f], [-0.130620778'f, 0.901816666'f]),
    (-7, [0.0'f, 0.0'f, 0.0'f], [0.0'f, 0.0'f]),
    (1337, [0.125'f, 0.25'f, 0.5'f], [0.923915744'f, 0.467908025'f]),
    (1337, [-0.75'f, 1.25'f, -2.5'f], [-0.248549402'f, 0.683826149'f]),
    (1337, [255.875'f, -256.25'f, 511.5'f], [-0.91916424'f, -0.61959368'f]),
    (1337, [0.0'f, 0.0'f, 0.0'f], [0.0'f, 0.0'f]),
    (low(int32).int, [0.125'f, 0.25'f, 0.5'f], [0.895951331'f, -0.0898984745'f]),
    (low(int32).int, [-0.75'f, 1.25'f, -2.5'f], [-0.732237697'f, -0.415868044'f]),
    (low(int32).int, [255.875'f, -256.25'f, 511.5'f],
      [-0.870669186'f, 0.924650908'f]),
    (low(int32).int, [0.0'f, 0.0'f, 0.0'f], [0.0'f, 0.0'f]),
    (2147483647, [0.125'f, 0.25'f, 0.5'f], [-0.250946581'f, 0.16036734'f]),
    (2147483647, [-0.75'f, 1.25'f, -2.5'f], [-0.709360778'f, 0.662877858'f]),
    (2147483647, [255.875'f, -256.25'f, 511.5'f],
      [0.000101377598'f, 0.289829373'f]),
    (2147483647, [0.0'f, 0.0'f, 0.0'f], [0.0'f, 0.0'f])
  ]:
    let noise = initOpenSimplex2(recorded[0])
    for dimensions in 2 .. 3:
      checkClose(noise.sample(recorded[1], dimensions),
        recorded[2][dimensions - 2],
        &"OpenSimplex2 {dimensions}D reference, seed={recorded[0]}")

echo "Testing OpenSimplex2 initialization and native scalar sampling"
block:
  static:
    doAssert not compiles(initOpenSimplex2(1).value(0.5'f))
    doAssert not compiles(initOpenSimplex2(1).value(0.5'f, 1.0'f, 2.0'f, 3.0'f))
  var rng = initRand(42187)
  for seed in [-42, 0, 7, 1988]:
    let
      noise = initOpenSimplex2(seed = seed)
      repeated = initOpenSimplex2(seed)
      other = initOpenSimplex2(seed + 1)
    doAssert noise.seed == seed
    doAssert noise.octaves == 1
    doAssert noise.amplitude == 1 and noise.frequency == 1
    doAssert noise.lacunarity == 2 and noise.gain == 0.5'f
    doAssert noise.value(0, 0) == 0
    doAssert noise.value(0, 0, 0) == 0
    for dimensions in 2 .. 3:
      var varies, differentSeed: bool
      var axisVaries: array[3, bool]
      for i in 0 ..< 128:
        var point: Point
        for axis in 0 ..< dimensions:
          point[axis] = (rng.rand(32.0) - 16).float32
        let actual = noise.sample(point, dimensions)
        doAssert abs(actual) < Inf
        doAssert actual == repeated.sample(point, dimensions)
        varies = varies or abs(actual) > 0.00001'f
        differentSeed = differentSeed or
          actual != other.sample(point, dimensions)
        for axis in 0 ..< dimensions:
          var shifted = point
          shifted[axis] += 0.25'f
          axisVaries[axis] = axisVaries[axis] or
            actual != noise.sample(shifted, dimensions)
      doAssert varies and differentSeed
      for axis in 0 ..< dimensions:
        doAssert axisVaries[axis]
    var scaled = noise
    scaled.frequency = 0.125'f
    for coordinate in -8 .. 8:
      doAssert scaled.value(coordinate, 2) ==
        scaled.value(coordinate.float32, 2.0'f)
      doAssert scaled.value(coordinate, 2, -3) ==
        scaled.value(coordinate.float32, 2.0'f, -3.0'f)

echo "Testing OpenSimplex2 continuity around negative and rounding boundaries"
block:
  const Epsilon = 1.0'f / 8192
  for seed in [-7, 0, 1988]:
    let noise = initOpenSimplex2(seed)
    for dimensions in 2 .. 3:
      for axis in 0 ..< dimensions:
        for boundary in [-256.0'f, -1.0'f, -0.5'f, 0.0'f,
          0.5'f, 1.0'f, 256.0'f]:
          var point = [0.125'f, -0.25'f, 0.5'f]
          point[axis] = boundary
          let center = noise.sample(point, dimensions)
          point[axis] = boundary - Epsilon
          let left = noise.sample(point, dimensions)
          point[axis] = boundary + Epsilon
          let right = noise.sample(point, dimensions)
          doAssert abs(left - center) < 0.01'f
          doAssert abs(right - center) < 0.01'f
          doAssert abs(left - right) < 0.01'f

echo "Testing OpenSimplex2 octave composition and invalid counts"
block:
  for dimensions in 2 .. 3:
    for parameters in [
      (1, 0.75'f, 0.25'f, 2.0'f, 0.5'f),
      (3, 1.25'f, -0.75'f, 1.5'f, -0.5'f),
      (4, -0.8'f, 1.5'f, 2.25'f, 0.35'f),
      (3, 0.0'f, 0.75'f, 2.0'f, 0.5'f),
      (3, 1.0'f, 0.0'f, 2.0'f, 0.5'f),
      (3, 1.0'f, 0.75'f, 2.0'f, 0.0'f)
    ]:
      var noise = initOpenSimplex2(1988)
      noise.octaves = parameters[0]
      noise.amplitude = parameters[1]
      noise.frequency = parameters[2]
      noise.lacunarity = parameters[3]
      noise.gain = parameters[4]
      let base = initOpenSimplex2(1988)
      for point in [[0.125'f, 0.25'f, -0.5'f], [-3.75'f, 2.125'f, 1.5'f]]:
        var
          expected = 0.0'f
          amplitude = noise.amplitude
          frequency = noise.frequency
        for octave in 0 ..< noise.octaves:
          var scaled = point
          for axis in 0 ..< dimensions:
            scaled[axis] *= frequency
          expected += base.sample(scaled, dimensions) * amplitude
          amplitude *= noise.gain
          frequency *= noise.lacunarity
        expected /= noise.octaves.float32
        checkClose(noise.sample(point, dimensions), expected,
          &"OpenSimplex2 {dimensions}D octave composition")
    for octaves in [0, -1, -3]:
      var noise = initOpenSimplex2(1988)
      noise.octaves = octaves
      var caught = false
      try:
        discard noise.sample([0.125'f, 0.25'f, 0.5'f], dimensions)
      except NoisyError:
        caught = true
      doAssert caught, &"OpenSimplex2 {dimensions}D accepted octaves={octaves}"
