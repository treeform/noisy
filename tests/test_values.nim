import
  std/[random, strformat],
  noisy/values

type Point = array[3, float32]

proc sample(noise: ValueNoise, point: Point, dimensions: int): float32 =
  ## Calls the native scalar overload for the requested dimension.
  if dimensions == 2:
    noise.value(point[0], point[1])
  else:
    noise.value(point[0], point[1], point[2])

proc checkClose(actual, expected: float32, context: string) =
  ## Checks finite samples with a small float32 rounding allowance.
  doAssert abs(actual) < Inf and abs(expected) < Inf, context
  let tolerance = 0.00001'f + 0.00001'f * max(abs(actual), abs(expected))
  doAssert abs(actual - expected) <= tolerance,
    context & &" actual={actual}, expected={expected}"

echo "Testing Value noise against recorded C++ reference outputs"
block:
  # Generated from FastNoiseLite at revision
  # 785f37a9ad76e283586a379675085f2063ae03f7 using frequency 1 and no fractal.
  for recorded in [
    (0, [0.125'f, 0.25'f, 0.5'f],
      [-0.0431922786'f, -0.359443843'f]),
    (0, [-0.75'f, 1.25'f, -2.5'f],
      [0.112270102'f, -0.428141683'f]),
    (0, [255.875'f, -256.25'f, 511.5'f],
      [0.504597902'f, -0.115530431'f]),
    (0, [0'f, 0'f, 0'f],
      [0'f, 0'f]),
    (0, [-1'f, -2'f, 3'f],
      [-0.350282609'f, -0.753594995'f]),
    (1, [0.125'f, 0.25'f, 0.5'f],
      [0.269043088'f, 0.0143271983'f]),
    (1, [-0.75'f, 1.25'f, -2.5'f],
      [0.121033534'f, -0.357810259'f]),
    (1, [255.875'f, -256.25'f, 511.5'f],
      [-0.115405537'f, -0.424039662'f]),
    (1, [0'f, 0'f, 0'f],
      [0.498890638'f, 0.498890638'f]),
    (1, [-1'f, -2'f, 3'f],
      [-0.456565738'f, 0.80233556'f]),
    (-7, [0.125'f, 0.25'f, 0.5'f],
      [0.0905476734'f, -0.343570203'f]),
    (-7, [-0.75'f, 1.25'f, -2.5'f],
      [-0.119796976'f, -0.234199643'f]),
    (-7, [255.875'f, -256.25'f, 511.5'f],
      [0.340282261'f, 0.141466141'f]),
    (-7, [0'f, 0'f, 0'f],
      [0.254234493'f, 0.254234493'f]),
    (-7, [-1'f, -2'f, 3'f],
      [-0.0121755302'f, 0.422142565'f]),
    (1988, [0.125'f, 0.25'f, 0.5'f],
      [-0.0659129247'f, -0.173089385'f]),
    (1988, [-0.75'f, 1.25'f, -2.5'f],
      [-0.651719689'f, -0.597496986'f]),
    (1988, [255.875'f, -256.25'f, 511.5'f],
      [-0.36987412'f, 0.0833657905'f]),
    (1988, [0'f, 0'f, 0'f],
      [0.100688525'f, 0.100688525'f]),
    (1988, [-1'f, -2'f, 3'f],
      [-0.250574887'f, -0.110779226'f]),
    (2147483647, [0.125'f, 0.25'f, 0.5'f],
      [0.434188426'f, 0.194765583'f]),
    (2147483647, [-0.75'f, 1.25'f, -2.5'f],
      [0.508338749'f, -0.390334517'f]),
    (2147483647, [255.875'f, -256.25'f, 511.5'f],
      [-0.323415548'f, -0.0391559601'f]),
    (2147483647, [0'f, 0'f, 0'f],
      [0.498890638'f, 0.498890638'f]),
    (2147483647, [-1'f, -2'f, 3'f],
      [0.730344296'f, 0.802335501'f]),
    (low(int32).int, [0.125'f, 0.25'f, 0.5'f],
      [-0.0431922786'f, -0.359443843'f]),
    (low(int32).int, [-0.75'f, 1.25'f, -2.5'f],
      [0.112270102'f, -0.428141683'f]),
    (low(int32).int, [255.875'f, -256.25'f, 511.5'f],
      [0.504597902'f, -0.115530431'f]),
    (low(int32).int, [0'f, 0'f, 0'f],
      [0'f, 0'f]),
    (low(int32).int, [-1'f, -2'f, 3'f],
      [-0.350282609'f, -0.753594995'f])
  ]:
    let noise = initValue(recorded[0])
    for dimensions in 2 .. 3:
      checkClose(noise.sample(recorded[1], dimensions),
        recorded[2][dimensions - 2], &"{dimensions}D recorded Value reference")

echo "Testing Value initialization, overloads, and determinism"
block:
  static:
    doAssert not compiles(initValue(1).value(0.0'f))
    doAssert not compiles(initValue(1).value(0.0'f, 0.0'f, 0.0'f, 0.0'f))
  var rng = initRand(7339)
  for seed in [0, 1, -7, 1988, 2147483647, low(int32).int]:
    let
      noise = initValue(seed = seed)
      repeated = initValue(seed)
      other = initValue(seed xor 12345)
    doAssert noise.seed == seed and noise.octaves == 1
    doAssert noise.amplitude == 1 and noise.frequency == 1
    doAssert noise.lacunarity == 2 and noise.gain == 0.5'f
    for dimensions in 2 .. 3:
      var varied, differentSeed: bool
      var axisVaries: array[3, bool]
      for i in 0 ..< 128:
        var point: Point
        for axis in 0 ..< dimensions:
          point[axis] = (rng.rand(64.0) - 32).float32
        let actual = noise.sample(point, dimensions)
        doAssert actual >= -1.00001'f and actual <= 1.00001'f
        doAssert actual == repeated.sample(point, dimensions)
        varied = varied or abs(actual) > 0.00001'f
        differentSeed = differentSeed or
          actual != other.sample(point, dimensions)
        for axis in 0 ..< dimensions:
          var shifted = point
          shifted[axis] += 0.25'f
          axisVaries[axis] = axisVaries[axis] or
            actual != noise.sample(shifted, dimensions)
      doAssert varied and differentSeed
      for axis in 0 ..< dimensions:
        doAssert axisVaries[axis]
    for coordinate in -8 .. 8:
      doAssert noise.value(coordinate, 2) ==
        noise.value(coordinate.float32, 2.0'f)
      doAssert noise.value(coordinate, 2, -3) ==
        noise.value(coordinate.float32, 2.0'f, -3.0'f)
  when sizeof(int) > 4:
    let
      wrapped = initValue((1'i64 shl 32).int + 7)
      original = initValue(7)
    doAssert wrapped.value(0.125'f, 0.25'f) == original.value(0.125'f, 0.25'f)
    doAssert wrapped.value(0.125'f, 0.25'f, 0.5'f) ==
      original.value(0.125'f, 0.25'f, 0.5'f)

echo "Testing Value lattice samples and continuity"
block:
  const Epsilon = 1.0'f / 4096
  for seed in [1, 42, 1988]:
    let noise = initValue(seed)
    for dimensions in 2 .. 3:
      var latticeVaries = false
      for coordinate in [-257, -256, -1, 0, 1, 255, 256, 257]:
        let point = [coordinate.float32, -3.0'f, 7.0'f]
        latticeVaries = latticeVaries or
          abs(noise.sample(point, dimensions)) > 0.00001'f
      doAssert latticeVaries
      for axis in 0 ..< dimensions:
        for boundary in [-256.0'f, -1.0'f, 0.0'f, 1.0'f, 256.0'f]:
          var point = [0.125'f, -0.25'f, 0.5'f]
          point[axis] = boundary
          let center = noise.sample(point, dimensions)
          point[axis] = boundary - Epsilon
          let left = noise.sample(point, dimensions)
          point[axis] = boundary + Epsilon
          let right = noise.sample(point, dimensions)
          doAssert abs(left - center) < 0.003'f
          doAssert abs(right - center) < 0.003'f
          doAssert abs(left - right) < 0.003'f

echo "Testing Value octave composition and invalid octave counts"
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
      var noise = initValue(1988)
      noise.octaves = parameters[0]
      noise.amplitude = parameters[1]
      noise.frequency = parameters[2]
      noise.lacunarity = parameters[3]
      noise.gain = parameters[4]
      let base = initValue(1988)
      for point in [
        [0.125'f, 0.25'f, -0.5'f],
        [-3.75'f, 2.125'f, 1.5'f]
      ]:
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
          &"{dimensions}D Value octave composition, parameters={parameters}")
    for octaves in [0, -1, -3]:
      var noise = initValue(1988)
      noise.octaves = octaves
      var caught = false
      try:
        discard noise.sample([0.125'f, 0.25'f, 0.5'f], dimensions)
      except NoisyError:
        caught = true
      doAssert caught, &"{dimensions}D Value accepted octaves={octaves}"
