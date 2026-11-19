import std/math, noisy

const Sentinel = 123456.0'f

proc close(actual, expected: float32, tolerance = 0.00003'f) =
  doAssert abs(actual) < Inf and abs(expected) < Inf
  doAssert abs(actual - expected) <= tolerance,
    "Tile sample mismatch: " & $actual & " != " & $expected

template expectError(body: untyped) =
  block:
    var caught = false
    try:
      body
    except NoisyError:
      caught = true
    doAssert caught

func referenceCircle(coordinate, period: float64): array[2, float32] =
  # Independent formula for small coordinates used by the sampling oracle.
  let angle = coordinate * 2 * PI / period
  [(cos(angle) * period / (2 * PI)).float32,
    (sin(angle) * period / (2 * PI)).float32]

proc checkLine[T](initial: T) =
  var noise = initial
  noise.frequency = 0.37'f
  noise.amplitude = 0.8'f
  noise.gain = -0.6'f
  noise.lacunarity = 1.7'f # Non-integer lacunarity must still tile.
  for octaves in [1, 3]:
    noise.octaves = octaves
    for width in 0 .. 25:
      let
        start = -3.25
        period = 17.5
        expected = noise.tileGrid(start, width, period)
        shifted = noise.tileGrid(start + period, width, period)
        backwards = noise.tileGrid(start - 3 * period, width, period)
      doAssert expected.width == width
      doAssert (expected.height, expected.depth, expected.frames) == (1, 1, 1)
      doAssert expected.values == shifted.values and expected.values == backwards.values
      var guarded: array[40, float32]
      for value in guarded.mitems: value = Sentinel
      noise.tileGrid(start, width, period, guarded.toOpenArray(1, width + 2))
      doAssert guarded[0] == Sentinel
      for x in 0 ..< width:
        doAssert guarded[x + 1] == expected[x]
        let point = referenceCircle(start + x.float64, period)
        close(expected[x], noise.value(point[0], point[1]))
      for i in width + 1 ..< guarded.len:
        doAssert guarded[i] == Sentinel
      var exact = newSeq[float32](width)
      noise.tileGrid(start, width, period, exact)
      doAssert exact == expected.values
      when defined(nimAllocStats):
        let before = getAllocStats()
      noise.tileGrid(start + period, width, period, exact)
      when defined(nimAllocStats):
        let after = getAllocStats()
        doAssert before == after
      doAssert exact == expected.values
      if width > 0:
        for value in guarded.mitems: value = Sentinel
        expectError:
          noise.tileGrid(start, width, period, guarded.toOpenArray(1, width - 1))
        for value in guarded: doAssert value == Sentinel

  for frequency in [0'f, -0.125'f]:
    noise.frequency = frequency
    var buffer: array[19, float32]
    noise.tileGrid(-3, 19, buffer)
    doAssert @buffer == noise.tileGrid(-3, 19, 19).values
    doAssert @buffer == noise.tileGrid(16, 19).values
    noise.tileGrid(-3.0'f, 19, 19.0'f, buffer)
    doAssert @buffer == noise.tileGrid(-3, 19).values

  # Endpoint equality belongs at x=period, not the last included sample.
  noise.frequency = 0.21'f
  let loop = noise.tileGrid(0, 33, 32)
  close(loop[0], loop[32])
  doAssert loop.values[0 ..< 32] != newSeq[float32](32)
  let middle = noise.tileGrid(16, 16, 32)
  for i in 0 ..< 16: close(middle[i], loop[i + 16])

proc checkPlane[T](initial: T) =
  var noise = initial
  noise.frequency = -0.29'f
  noise.amplitude = 1.2'f
  noise.gain = 0.6'f
  noise.lacunarity = 1.7'f
  for octaves in [1, 3]:
    noise.octaves = octaves
    for height in 0 .. 25:
      let
        start = (-3.25, 2.5)
        period = (17.5, 23.25)
        shape = (3, height)
        expected = noise.tileGrid(start, shape, period)
        shiftedX = noise.tileGrid((start[0] + period[0], start[1]), shape, period)
        shiftedY = noise.tileGrid((start[0], start[1] - period[1]), shape, period)
      doAssert (expected.width, expected.height, expected.depth, expected.frames) ==
        (3, height, 1, 1)
      doAssert expected.values == shiftedX.values and expected.values == shiftedY.values
      var guarded: array[96, float32]
      for value in guarded.mitems: value = Sentinel
      let count = 3 * height
      noise.tileGrid(start, shape, period, guarded.toOpenArray(1, count + 2))
      doAssert guarded[0] == Sentinel
      for x in 0 ..< 3:
        for y in 0 ..< height:
          doAssert guarded[x * height + y + 1] == expected[x, y]
          let
            a = referenceCircle(start[0] + x.float64, period[0])
            b = referenceCircle(start[1] + y.float64, period[1])
          close(expected[x, y], noise.value(a[0], a[1], b[0], b[1]))
      for i in count + 1 ..< guarded.len: doAssert guarded[i] == Sentinel
      var exact = newSeq[float32](count)
      noise.tileGrid(start, shape, period, exact)
      doAssert exact == expected.values
      when defined(nimAllocStats):
        let before = getAllocStats()
      noise.tileGrid((start[0] + period[0], start[1]), shape, period, exact)
      when defined(nimAllocStats):
        let after = getAllocStats()
        doAssert before == after
      doAssert exact == expected.values
      if count > 0:
        for value in guarded.mitems: value = Sentinel
        expectError:
          noise.tileGrid(start, shape, period, guarded.toOpenArray(1, count - 1))
        for value in guarded: doAssert value == Sentinel

  # Both axes close, including the corner, with independent arbitrary periods.
  let whole = noise.tileGrid((0, 0), (18, 24), (17, 23))
  for x in 0 ..< 18: close(whole[x, 0], whole[x, 23])
  for y in 0 ..< 24: close(whole[0, y], whole[17, y])
  let chunk = noise.tileGrid((8, 12), (9, 11), (17, 23))
  for x in 0 ..< 9:
    for y in 0 ..< 11: close(chunk[x, y], whole[x + 8, y + 12])
  var buffer: array[17 * 23, float32]
  noise.tileGrid((-2, 3), (17, 23), buffer)
  doAssert @buffer == noise.tileGrid((-2, 3), (17, 23), (17, 23)).values
  doAssert @buffer == noise.tileGrid((15, -20), (17, 23)).values
  noise.tileGrid((-2.0'f, 3.0'f), (17, 23), buffer)
  doAssert @buffer == noise.tileGrid((-2, 3), (17, 23)).values

  # Smoothness through the wrap: samples and finite differences on both sides
  # agree with the same neighborhoods one whole period away.
  for epsilon in [0.01, 0.001]:
    for y in [2.25, 7.75]:
      let
        left = noise.tileGrid((-epsilon, y), (1, 1), (17.0, 23.0))[0, 0]
        right = noise.tileGrid((epsilon, y), (1, 1), (17.0, 23.0))[0, 0]
        before = noise.tileGrid((17 - epsilon, y), (1, 1), (17.0, 23.0))[0, 0]
        after = noise.tileGrid((17 + epsilon, y), (1, 1), (17.0, 23.0))[0, 0]
      close(left, before)
      close(right, after)
      close((right - left) / (2 * epsilon).float32,
        (after - before) / (2 * epsilon).float32, 0.001'f)

proc checkInvalid[T](initial: T) =
  var noise = initial
  var buffer: array[32, float32]
  for value in buffer.mitems: value = Sentinel
  for width in [-1, int.high]:
    expectError: noise.tileGrid(0, width, 32, buffer)
    expectError: discard noise.tileGrid(0, width, 32)
  for period in [0.0, -1.0, NaN, Inf, -Inf, 1e300]:
    expectError: noise.tileGrid(0, 8, period, buffer)
    expectError: discard noise.tileGrid(0, 8, period)
    expectError: noise.tileGrid(0, 0, period, buffer)
  for start in [NaN, Inf, -Inf]:
    expectError: noise.tileGrid(start, 8, 32, buffer)
  for octaves in [-1, 0]:
    noise.octaves = octaves
    expectError: noise.tileGrid(0, 8, buffer)
    expectError: noise.tileGrid(0, 0, buffer)
  noise = initial
  noise.frequency = Inf.float32
  expectError: noise.tileGrid(0, 8, buffer)
  noise = initial
  noise.octaves = 3
  noise.lacunarity = 1e30'f
  expectError: noise.tileGrid(0, 8, buffer)
  noise = initial
  noise.frequency = 0
  expectError: noise.tileGrid(0, 8, 1e300, buffer)
  when T is Cellular:
    noise = initial
    noise.jitter = NaN.float32
    expectError: noise.tileGrid(0, 8, buffer)
    expectError: noise.tileGrid(0, 0, buffer)
  when T is Simplex or T is Perlin:
    noise = initial
    for shape in [(-1, 3), (3, -1), (0, -1), (int.high, 2)]:
      expectError: noise.tileGrid((0, 0), shape, (32, 32), buffer)
      expectError: discard noise.tileGrid((0, 0), shape, (32, 32))
    for period in [(0.0, 32.0), (32.0, -1.0), (NaN, 32.0), (32.0, Inf)]:
      expectError: noise.tileGrid((0, 0), (3, 3), period, buffer)
      expectError: discard noise.tileGrid((0, 0), (3, 3), period)
    expectError: noise.tileGrid((Inf, 0.0), (3, 3), (32, 32), buffer)
    for shape in [(0, 17), (19, 0)]:
      noise.tileGrid((0, 0), shape, buffer)
      doAssert noise.tileGrid((0, 0), shape).values.len == 0
    noise.octaves = 0
    expectError: noise.tileGrid((0, 0), (3, 3), buffer)
    expectError: noise.tileGrid((0, 0), (0, 3), buffer)
  noise = initial
  noise.tileGrid(0, 0, buffer)
  doAssert noise.tileGrid(0, 0).values.len == 0
  for value in buffer: doAssert value == Sentinel

echo "Testing seamless tiled buffers, periods, seams, and scalar/SIMD parity"
template check(initial: untyped) =
  checkLine(initial)
  checkInvalid(initial)
check(initSimplex(1988))
check(initPerlin(42))
check(initOpenSimplex2(-7))
check(initOpenSimplex2S(1988))
check(initValue(42))
check(initValueCubic(-7))
check(initCellular(1988))
checkPlane(initSimplex(42))
checkPlane(initPerlin(1988))

block:
  var noise = initCellular(-7)
  for metric in CellularDistance:
    noise.distanceFunction = metric
    for mode in CellularReturn:
      noise.returnType = mode
      for jitter in [0'f, 1'f]:
        noise.jitter = jitter
        let loop = noise.tileGrid(-2.25, 19, 17.5)
        for x in 0 ..< 19:
          let point = referenceCircle(-2.25 + x.float64, 17.5)
          close(loop[x], noise.value(point[0], point[1]))

static:
  doAssert not compiles(initValue(42).tileGrid((0, 0), (8, 8)))
  doAssert not compiles(initCellular(42).tileGrid((0, 0), (8, 8)))
  doAssert not compiles(initSimplex(42).tileGrid((0, 0, 0), (8, 8, 8)))

when (NimMajor, NimMinor) >= (2, 4):
  static:
    let noise = initPerlin(42)
    var buffer: array[18, float32]
    noise.tileGrid((0, 0), (2, 9), (7, 8), buffer)
    for x in 0 ..< 2:
      doAssert buffer[x * 9] == buffer[x * 9 + 8]
    doAssert @buffer == noise.tileGrid((0, 0), (2, 9), (7, 8)).values
    let line = initValue(42).tileGrid(0, 9, 8)
    doAssert line[0] == line[8]
