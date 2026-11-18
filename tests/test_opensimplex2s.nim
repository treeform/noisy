import
  std/[random, strformat],
  noisy

type Point = array[3, float32]

proc sample(noise: OpenSimplex2S, point: Point, dimensions: int): float32 =
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

echo "Testing OpenSimplex2S against recorded C++ outputs"
block:
  # Generated with FastNoiseLite C++ at revision
  # 785f37a9ad76e283586a379675085f2063ae03f7.
  # Frequency is 1 and rotation uses the default orientation.
  # The unchanged reference was compiled with explicit signed wrapping.
  for recorded in [
    (0, [0.125'f, 0.25'f, 0.5'f],
      [0.724977612'f, 0.0815235153'f]),
    (0, [-0.75'f, 1.25'f, -2.5'f],
      [-0.610085785'f, 0.30403021'f]),
    (0, [0.75'f, 0.25'f, 0.125'f],
      [-0.0208808873'f, -0.0565822944'f]),
    (0, [-0.5'f, -0.5'f, -0.5'f],
      [-0.346485138'f, 0.0'f]),
    (0, [-1024.00012'f, -0.625'f, 255.999893'f],
      [0.692286909'f, 0.616423905'f]),
    (1, [0.125'f, 0.25'f, 0.5'f],
      [-0.178485289'f, 0.176762313'f]),
    (1, [-0.75'f, 1.25'f, -2.5'f],
      [-0.657496035'f, -0.0410478897'f]),
    (1, [0.75'f, 0.25'f, 0.125'f],
      [0.632207513'f, -0.108116336'f]),
    (1, [-0.5'f, -0.5'f, -0.5'f],
      [-0.217624471'f, 0.0'f]),
    (1, [-1024.00012'f, -0.625'f, 255.999893'f],
      [-0.111890294'f, -0.226899549'f]),
    (-7, [0.125'f, 0.25'f, 0.5'f],
      [0.105788402'f, -0.426436037'f]),
    (-7, [-0.75'f, 1.25'f, -2.5'f],
      [-0.660588026'f, -0.72275281'f]),
    (-7, [0.75'f, 0.25'f, 0.125'f],
      [-0.561867535'f, -0.166802332'f]),
    (-7, [-0.5'f, -0.5'f, -0.5'f],
      [-0.226360112'f, 0.0'f]),
    (-7, [-1024.00012'f, -0.625'f, 255.999893'f],
      [-0.429214567'f, -0.594848692'f]),
    (1988, [0.125'f, 0.25'f, 0.5'f],
      [0.447658807'f, 0.437856078'f]),
    (1988, [-0.75'f, 1.25'f, -2.5'f],
      [-0.278522432'f, 0.550321221'f]),
    (1988, [0.75'f, 0.25'f, 0.125'f],
      [-0.513240397'f, -0.13477385'f]),
    (1988, [-0.5'f, -0.5'f, -0.5'f],
      [0.282705635'f, 0.0'f]),
    (1988, [-1024.00012'f, -0.625'f, 255.999893'f],
      [-0.134435847'f, -0.314159334'f]),
    (low(int32).int, [0.125'f, 0.25'f, 0.5'f],
      [0.724977612'f, 0.0815235153'f]),
    (low(int32).int, [-0.75'f, 1.25'f, -2.5'f],
      [-0.610085785'f, 0.30403021'f]),
    (low(int32).int, [0.75'f, 0.25'f, 0.125'f],
      [-0.0208808873'f, -0.0565822944'f]),
    (low(int32).int, [-0.5'f, -0.5'f, -0.5'f],
      [-0.346485138'f, 0.0'f]),
    (low(int32).int, [-1024.00012'f, -0.625'f, 255.999893'f],
      [0.692286909'f, 0.616423905'f]),
    (2147483647, [0.125'f, 0.25'f, 0.5'f],
      [-0.0776248202'f, -0.301910043'f]),
    (2147483647, [-0.75'f, 1.25'f, -2.5'f],
      [-0.667285621'f, -0.00937497523'f]),
    (2147483647, [0.75'f, 0.25'f, 0.125'f],
      [0.419845253'f, 0.181554973'f]),
    (2147483647, [-0.5'f, -0.5'f, -0.5'f],
      [-0.216827676'f, 0.0'f]),
    (2147483647, [-1024.00012'f, -0.625'f, 255.999893'f],
      [-0.325451702'f, 0.120644413'f])
  ]:
    let noise = initOpenSimplex2S(recorded[0])
    for dimensions in 2 .. 3:
      checkClose(noise.sample(recorded[1], dimensions),
        recorded[2][dimensions - 2],
        &"OpenSimplex2S {dimensions}D reference, seed={recorded[0]}")

echo "Testing OpenSimplex2S initialization and native scalar sampling"
block:
  static:
    doAssert compiles(initOpenSimplex2S(1).grid((0, 0), (5, 5)))
    doAssert compiles(initOpenSimplex2S(1).grid((0, 0, 0), (5, 5, 5)))
    doAssert compiles(initValueCubic(1).grid((0, 0), (5, 5)))
    doAssert compiles(initValueCubic(1).grid((0, 0, 0), (5, 5, 5)))
    doAssert not compiles(initOpenSimplex2S(1).value(0.5'f))
    doAssert not compiles(initOpenSimplex2S(1).value(0.5'f, 1.0'f, 2.0'f, 3.0'f))
  var rng = initRand(42187)
  for seed in [-42, 0, 7, 1988]:
    let
      noise = initOpenSimplex2S(seed = seed)
      repeated = initOpenSimplex2S(seed)
      other = initOpenSimplex2S(seed + 1)
      original = initOpenSimplex2(seed)
    doAssert noise.seed == seed
    doAssert noise.octaves == 1
    doAssert noise.amplitude == 1 and noise.frequency == 1
    doAssert noise.lacunarity == 2 and noise.gain == 0.5'f
    for dimensions in 2 .. 3:
      var varies, differentSeed, differentAlgorithm: bool
      var axisVaries: array[3, bool]
      for i in 0 ..< 128:
        var point: Point
        for axis in 0 ..< dimensions:
          point[axis] = (rng.rand(32.0) - 16).float32
        let
          actual = noise.sample(point, dimensions)
          fast =
            if dimensions == 2:
              original.value(point[0], point[1])
            else:
              original.value(point[0], point[1], point[2])
        differentAlgorithm = differentAlgorithm or abs(actual - fast) > 0.00001'f
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
      doAssert varies and differentSeed and differentAlgorithm
      for axis in 0 ..< dimensions:
        doAssert axisVaries[axis]
    var scaled = noise
    scaled.frequency = 0.125'f
    for coordinate in -8 .. 8:
      doAssert scaled.value(coordinate, 2) ==
        scaled.value(coordinate.float32, 2.0'f)
      doAssert scaled.value(coordinate, 2, -3) ==
        scaled.value(coordinate.float32, 2.0'f, -3.0'f)

block:
  when sizeof(int) > 4:
    let
      base = initOpenSimplex2S(7)
      wrapped = initOpenSimplex2S((1'i64 shl 32).int + 7)
    doAssert base.value(0.125'f, 0.25'f) == wrapped.value(0.125'f, 0.25'f)
    doAssert base.value(0.125'f, 0.25'f, 0.5'f) ==
      wrapped.value(0.125'f, 0.25'f, 0.5'f)

echo "Testing OpenSimplex2S continuity around negative and rounding boundaries"
block:
  const Epsilon = 1.0'f / 8192
  for seed in [-7, 0, 1988]:
    let noise = initOpenSimplex2S(seed)
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

echo "Testing OpenSimplex2S octave composition and invalid counts"
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
      var noise = initOpenSimplex2S(1988)
      noise.octaves = parameters[0]
      noise.amplitude = parameters[1]
      noise.frequency = parameters[2]
      noise.lacunarity = parameters[3]
      noise.gain = parameters[4]
      let base = initOpenSimplex2S(1988)
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
          &"OpenSimplex2S {dimensions}D octave composition")
    for octaves in [0, -1, -3]:
      var noise = initOpenSimplex2S(1988)
      noise.octaves = octaves
      var caught = false
      try:
        discard noise.sample([0.125'f, 0.25'f, 0.5'f], dimensions)
      except NoisyError:
        caught = true
      doAssert caught, &"OpenSimplex2S {dimensions}D accepted octaves={octaves}"
