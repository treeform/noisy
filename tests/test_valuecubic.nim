import std/strformat, noisy/values

const
  Seeds = [0, 1, -1, 1988, low(int32).int, high(int32).int]
  Points = [
    [0.0'f, 0.0'f, 0.0'f],
    [0.125'f, 0.25'f, 0.5'f],
    [-0.75'f, 1.25'f, -2.5'f],
    [255.875'f, -256.25'f, 511.5'f],
    [-1024.0001'f, -0.625'f, 1.25'f]
  ]
  # Recorded from unchanged FastNoiseLite C++ at revision
  # 785f37a9ad76e283586a379675085f2063ae03f7, frequency 1, no fractal.
  Expected = [
    [
      [0.0'f, 0.0'f],
      [-0.0451375879'f, -0.136378735'f],
      [0.0681701154'f, -0.256528676'f],
      [0.319436938'f, -0.0162594225'f],
      [0.0201216005'f, -0.0539245419'f]
    ],
    [
      [0.221729174'f, 0.147819445'f],
      [0.0692128614'f, 0.00701034535'f],
      [0.0466976352'f, -0.121986054'f],
      [-0.00252707792'f, -0.0672198161'f],
      [-0.300865561'f, 0.225546524'f]
    ],
    [
      [0.221729174'f, 0.147819445'f],
      [0.247070819'f, 0.113498285'f],
      [0.311992854'f, -0.201760367'f],
      [-0.234333992'f, 0.0263744555'f],
      [-0.338119656'f, 0.0966780782'f]
    ],
    [
      [0.0447504558'f, 0.0298336372'f],
      [0.0190529004'f, -0.0270034056'f],
      [-0.420188487'f, -0.215460926'f],
      [-0.0952456668'f, 0.020236766'f],
      [-0.23294121'f, 0.260413617'f]
    ],
    [
      [0.0'f, 0.0'f],
      [-0.0451375879'f, -0.136378735'f],
      [0.0681701154'f, -0.256528676'f],
      [0.319436938'f, -0.0162594225'f],
      [0.0201216005'f, -0.0539245419'f]
    ],
    [
      [0.221729174'f, 0.147819445'f],
      [0.247070819'f, 0.113498285'f],
      [0.311992854'f, -0.201760367'f],
      [-0.234333992'f, 0.0263744555'f],
      [-0.338119656'f, 0.0966780782'f]
    ]
  ]

func sample(noise: ValueCubic, point: array[3, float32], dimensions: int): float32 =
  ## Calls the native scalar overload for the requested dimension.
  if dimensions == 2:
    noise.value(point[0], point[1])
  else:
    noise.value(point[0], point[1], point[2])

echo "Testing ValueCubic against recorded C++ outputs"
for seedIndex, seed in Seeds:
  let noise = initValueCubic(seed)
  for pointIndex, point in Points:
    for dimensions in 2 .. 3:
      let
        actual = noise.sample(point, dimensions)
        expected = Expected[seedIndex][pointIndex][dimensions - 2]
      doAssert abs(actual - expected) <= 0.00001'f,
        &"ValueCubic {dimensions}D, seed={seed}, point={point}, " &
        &"actual={actual}, expected={expected}"

echo "Testing ValueCubic settings and lattice boundaries"
block:
  let
    noise = initValueCubic(seed = 1988)
    hermite = initValue(1988)
  doAssert noise.seed == 1988 and noise.octaves == 1
  doAssert noise.amplitude == 1 and noise.frequency == 1
  doAssert noise.lacunarity == 2 and noise.gain == 0.5'f
  for x in -4 .. 4:
    doAssert noise.value(x, -2) == noise.value(x.float32, -2.0'f)
    doAssert noise.value(x, -2, 3) == noise.value(x.float32, -2.0'f, 3.0'f)
    doAssert abs(noise.value(x, -2) - hermite.value(x, -2) / 2.25'f) < 0.00001'f
    doAssert abs(noise.value(x, -2, 3) - hermite.value(x, -2, 3) / 3.375'f) <
      0.00001'f
  for dimensions in 2 .. 3:
    let center = noise.sample([-1.0'f, 0.25'f, -0.5'f], dimensions)
    for x in [-1.0001'f, -0.9999'f]:
      doAssert abs(noise.sample([x, 0.25'f, -0.5'f], dimensions) - center) <
        0.001'f
    var layered = noise
    layered.octaves = 2
    layered.amplitude = 0.75'f
    layered.frequency = -0.5'f
    layered.lacunarity = 1.5'f
    layered.gain = 0.25'f
    let
      first = noise.sample([-0.125'f, 0.375'f, -0.625'f], dimensions)
      second = noise.sample([-0.1875'f, 0.5625'f, -0.9375'f], dimensions)
      expected = (first * 0.75'f + second * 0.1875'f) / 2
      actual = layered.sample([0.25'f, -0.75'f, 1.25'f], dimensions)
    doAssert abs(actual - expected) <= 0.00001'f
    for octaves in [0, -1]:
      layered.octaves = octaves
      var caught = false
      try:
        discard layered.sample(Points[1], dimensions)
      except NoisyError:
        caught = true
      doAssert caught, &"ValueCubic {dimensions}D accepted octaves={octaves}"
