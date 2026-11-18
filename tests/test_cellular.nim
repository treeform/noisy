import noisy/cellulars

const Recorded = [
  [
    [-0.500490189'f, -0.369607627'f, -0.269942582'f,
      -0.319775105'f, -0.900334954'f, -0.769888699'f, -0.13651675'f],
    [-0.500490189'f, -0.602605462'f, -0.46701616'f,
      -0.534810781'f, -0.864410698'f, -0.894097567'f, -0.254396617'f],
    [-0.500490189'f, -0.273397386'f, 0.031386137'f,
      -0.121005654'f, -0.695216477'f, -0.625296056'f, -0.295508623'f],
    [-0.500490189'f, 0.123997211'f, 0.564369917'f,
      0.344183564'f, -0.559627295'f, -0.120826304'f, -0.281501651'f]
  ],
  [
    [-0.537805617'f, -0.633392572'f, -0.332094431'f,
      -0.482743502'f, -0.698701859'f, -0.87757045'f, -0.451108932'f],
    [-0.537805617'f, -0.865599036'f, -0.553902149'f,
      -0.709750593'f, -0.688303173'f, -0.970022023'f, -0.698718548'f],
    [-0.537805617'f, -0.381411195'f, -0.173528016'f,
      -0.277469635'f, -0.792116821'f, -0.744376838'f, -0.251530826'f],
    [-0.537805617'f, -0.247010231'f, 0.401465535'f,
      0.0772275925'f, -0.351524234'f, -0.472355366'f, -0.462712586'f]
  ]
]

proc sample(noise: Cellular, point: array[3, float32], dimensions: int): float32 =
  ## Calls the public scalar overload for the requested dimension.
  if dimensions == 2:
    noise.value(point[0], point[1])
  else:
    noise.value(point[0], point[1], point[2])

proc checkClose(actual, expected: float32) =
  ## Rejects nonfinite results and discrepancies larger than float32 rounding.
  doAssert abs(actual) < Inf and abs(expected) < Inf
  doAssert abs(actual - expected) <=
    0.00001'f + 0.00001'f * max(abs(actual), abs(expected)),
    "Cellular actual=" & $actual & ", expected=" & $expected

echo "Testing Cellular against recorded C++ reference outputs"
block:
  # Recorded from unchanged FastNoiseLite at revision
  # 785f37a9ad76e283586a379675085f2063ae03f7 with frequency 1 and seed 1988.
  for dimensions in 2 .. 3:
    var noise = initCellular(1988)
    for metric in CellularDistance:
      noise.distanceFunction = metric
      for returnType in CellularReturn:
        noise.returnType = returnType
        checkClose(
          noise.sample([0.125'f, -0.75'f, 1.5'f], dimensions),
          Recorded[dimensions - 2][metric.ord][returnType.ord]
        )

echo "Testing Cellular defaults, seeds, and scalar overloads"
block:
  static:
    doAssert not compiles(initCellular(1).grid((0, 0), (4, 4)))
    doAssert not compiles(initCellular(1).grid((0, 0, 0), (4, 4, 4)))
  for seed in [low(int32).int, -1988, 0, 1, 1988, high(int32).int]:
    let
      noise = initCellular(seed = seed)
      repeated = initCellular(seed)
      other = initCellular(seed + 1)
    doAssert noise.seed == seed and noise.octaves == 1
    doAssert noise.amplitude == 1 and noise.frequency == 1
    doAssert noise.lacunarity == 2 and noise.gain == 0.5'f
    doAssert noise.distanceFunction == EuclideanSquared
    doAssert noise.returnType == Distance and noise.jitter == 1
    for dimensions in 2 .. 3:
      var varied, seedVaries: bool
      var axesVary: array[3, bool]
      for i in -12 .. 12:
        let point = [i.float32 * 0.37'f, i.float32 * -0.125'f, 0.75'f]
        let actual = noise.sample(point, dimensions)
        doAssert abs(actual) < Inf
        doAssert actual == repeated.sample(point, dimensions)
        varied = varied or actual != noise.sample([0.0'f, 0, 0], dimensions)
        seedVaries = seedVaries or actual != other.sample(point, dimensions)
        for axis in 0 ..< dimensions:
          var shifted = point
          shifted[axis] += 0.25'f
          axesVary[axis] = axesVary[axis] or
            actual != noise.sample(shifted, dimensions)
      doAssert varied and seedVaries
      for axis in 0 ..< dimensions:
        doAssert axesVary[axis]
    for x in -8 .. 8:
      doAssert noise.value(x, -3) == noise.value(x.float32, -3.0'f)
      doAssert noise.value(x, -3, 5) ==
        noise.value(x.float32, -3.0'f, 5.0'f)
  when sizeof(int) > 4:
    let
      first = initCellular(1988)
      wrapped = initCellular(1988 + (1 shl 32))
    doAssert first.value(0.125'f, 0.75'f) ==
      wrapped.value(0.125'f, 0.75'f)

echo "Testing Cellular distances, zero-jitter ties, and continuity"
block:
  for dimensions in 2 .. 3:
    var noise = initCellular(1988)
    noise.jitter = 0
    for metric in CellularDistance:
      noise.distanceFunction = metric
      noise.returnType = Distance
      doAssert noise.sample([0.0'f, 0, 0], dimensions) == -1
      noise.returnType = Distance2Sub
      doAssert noise.sample([0.5'f, -0.5'f, 1.5'f], dimensions) == -1
      noise.returnType = Distance2Div
      doAssert noise.sample([0.5'f, -0.5'f, 1.5'f], dimensions) == 0
      noise.returnType = CellValue
      let expected = if dimensions == 2: -0.500490189'f else: -0.537805617'f
      checkClose(noise.sample([0.5'f, -0.5'f, 1.5'f], dimensions), expected)
    noise.returnType = Distance
    noise.distanceFunction = EuclideanSquared
    for jitter in [0.0'f, 0.5'f, 1.0'f]:
      noise.jitter = jitter
      for axis in 0 ..< dimensions:
        for boundary in [-256.5'f, -1.5'f, -0.5'f, 0.5'f, 1.5'f, 256.5'f]:
          var point = [0.125'f, -0.25'f, 0.75'f]
          point[axis] = boundary
          let center = noise.sample(point, dimensions)
          point[axis] = boundary - 0.0001'f
          let left = noise.sample(point, dimensions)
          point[axis] = boundary + 0.0001'f
          let right = noise.sample(point, dimensions)
          doAssert abs(left - center) < 0.001'f
          doAssert abs(right - center) < 0.001'f

echo "Testing Cellular octave composition and parameter errors"
block:
  for dimensions in 2 .. 3:
    for metric in CellularDistance:
      for returnType in CellularReturn:
        var noise = initCellular(42)
        noise.octaves = 3
        noise.amplitude = 0.75'f
        noise.frequency = -0.25'f
        noise.lacunarity = 1.5'f
        noise.gain = -0.4'f
        noise.distanceFunction = metric
        noise.returnType = returnType
        noise.jitter = 0.5'f
        var base = noise
        base.octaves = 1
        base.amplitude = 1
        base.frequency = 1
        for point in [
          [0.125'f, -0.75'f, 1.5'f],
          [-3.75'f, 2.125'f, -0.25'f]
        ]:
          var
            expected: float32
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
          checkClose(noise.sample(point, dimensions), expected)
    for octaves in [0, -1, -3]:
      var noise = initCellular(42)
      noise.octaves = octaves
      var caught = false
      try:
        discard noise.sample([0.125'f, -0.75'f, 1.5'f], dimensions)
      except NoisyError:
        caught = true
      doAssert caught
    for jitter in [-0.001'f, 1.001'f, Inf.float32, NegInf.float32, NaN.float32]:
      var noise = initCellular(42)
      noise.jitter = jitter
      var caught = false
      try:
        discard noise.sample([0.125'f, -0.75'f, 1.5'f], dimensions)
      except NoisyError:
        caught = true
      doAssert caught
