import noisy

const Sentinel = 123456.0'f

template expectError(body: untyped) =
  block:
    var caught = false
    try:
      body
    except NoisyError:
      caught = true
    doAssert caught

proc checkFill[T, Start, Size](noise: T, start: Start, dimens: Size) =
  ## Compare all destination types, guard the slice boundaries, and reuse storage.
  let expected = noise.grid(start, dimens)
  let length = expected.values.len
  var exact = newSeq[float32](length)
  noise.grid(start, dimens, exact)
  doAssert exact == expected.values
  # Grid owns a sequence, so its existing storage can be refilled directly too.
  noise.grid(start, dimens, expected.values)
  doAssert exact == expected.values

  var buffer = newSeq[float32](length + 7)
  for value in buffer.mitems:
    value = Sentinel
  let address = buffer[0].addr
  noise.grid(start, dimens, buffer)
  doAssert buffer[0].addr == address and buffer.len == length + 7
  for i in 0 ..< length:
    doAssert buffer[i] == expected.values[i]
  for i in length ..< buffer.len:
    doAssert buffer[i] == Sentinel

  # Arrays and unaligned subviews use the same SIMD paths as sequences.
  var guarded: array[2048, float32]
  doAssert length + 3 < guarded.len
  for value in guarded.mitems:
    value = Sentinel
  noise.grid(start, dimens, guarded.toOpenArray(1, length + 2))
  doAssert guarded[0] == Sentinel
  for i in 0 ..< length:
    doAssert guarded[i + 1] == expected.values[i]
  for i in length + 1 ..< guarded.len:
    doAssert guarded[i] == Sentinel

  var changed = noise
  changed.amplitude = noise.amplitude * -0.75'f
  changed.frequency = noise.frequency * 1.3'f
  let changedExpected = changed.grid(start, dimens)
  # Warmed calls must allocate nothing when allocator instrumentation is enabled.
  when defined(nimAllocStats):
    let before = getAllocStats()
  changed.grid(start, dimens, buffer)
  when defined(nimAllocStats):
    let after = getAllocStats()
    doAssert after == before, "Buffer-taking grid allocated or freed memory"
  doAssert buffer[0].addr == address and buffer.len == length + 7
  for i in 0 ..< length:
    doAssert buffer[i] == changedExpected.values[i]
  for i in length ..< buffer.len:
    doAssert buffer[i] == Sentinel

  if length > 0:
    var short = newSeq[float32](length - 1)
    for value in short.mitems:
      value = Sentinel
    expectError:
      noise.grid(start, dimens, short)
    for value in short:
      doAssert value == Sentinel

proc checkFamily[T](initial: T) =
  var noise = initial
  noise.frequency = -0.125'f
  noise.amplitude = 0.8'f
  noise.gain = -0.6'f
  noise.lacunarity = 1.7'f
  for octaves in [1, 3]:
    noise.octaves = octaves
    # All remainders for four/eight lanes, plus empty grids and multiple batches.
    for size in 0 .. 17:
      checkFill(noise, (-1.375'f, 255.125'f), (5, size))
      checkFill(noise, (-1.375'f, 255.125'f, -0.625'f), (5, 5, size))
      when T is Simplex or T is Perlin:
        checkFill(noise, -1.375'f, size)
        checkFill(noise, (-1.375'f, 255.125'f, -0.625'f, 0.25'f), (5, 3, size, 2))
    checkFill(noise, (-1, 2), (5, 9))
    checkFill(noise, (-1, 2, -3), (5, 5, 9))
    when T is Simplex or T is Perlin:
      checkFill(noise, -1, 9)
      checkFill(noise, (-1, 2, -3, 4), (5, 3, 9, 2))

  var storage: array[27, float32]
  for value in storage.mitems:
    value = Sentinel
  # Empty dimensions must not form a pointer or enter large outer loops.
  for shape in [(0, int.high), (int.high, 0)]:
    noise.grid((0, 0), shape, storage)
    doAssert noise.grid((0, 0), shape).values.len == 0
  for shape in [(0, 3, int.high), (int.high, 0, 3), (3, int.high, 0)]:
    noise.grid((0, 0, 0), shape, storage)
    doAssert noise.grid((0, 0, 0), shape).values.len == 0
  for shape in [(-1, 3), (3, -1), (0, -1), (int.high, 2)]:
    expectError: noise.grid((0, 0), shape, storage)
    expectError: discard noise.grid((0, 0), shape)
  for shape in [(-1, 2, 3), (2, -1, 3), (2, 3, -1), (0, 2, -1), (int.high, 2, 2)]:
    expectError: noise.grid((0, 0, 0), shape, storage)
    expectError: discard noise.grid((0, 0, 0), shape)
  when T is Simplex or T is Perlin:
    for shape in [-1, int.high]:
      expectError: noise.grid(0, shape, storage)
      expectError: discard noise.grid(0, shape)
    for shape in [(-1, 2, 3, 4), (2, -1, 3, 4), (2, 3, -1, 4), (2, 3, 4, -1),
        (0, 2, 3, -1), (int.high, 2, 2, 2)]:
      expectError: noise.grid((0, 0, 0, 0), shape, storage)
      expectError: discard noise.grid((0, 0, 0, 0), shape)
    for shape in [(0, 2, 3, int.high), (int.high, 0, 3, 2),
        (2, int.high, 0, 3), (2, 3, int.high, 0)]:
      noise.grid((0, 0, 0, 0), shape, storage)
      doAssert noise.grid((0, 0, 0, 0), shape).values.len == 0

  for octaves in [-1, 0]:
    noise.octaves = octaves
    for shape in [(3, 9), (0, 9)]:
      expectError: noise.grid((0, 0), shape, storage)
      expectError: discard noise.grid((0, 0), shape)
    for shape in [(1, 3, 9), (0, 3, 9)]:
      expectError: noise.grid((0, 0, 0), shape, storage)
      expectError: discard noise.grid((0, 0, 0), shape)
    when T is Simplex or T is Perlin:
      for size in [0, 9]:
        expectError: noise.grid(0, size, storage)
        expectError: noise.grid((0, 0, 0, 0), (1, 1, size, 1), storage)
  for value in storage:
    doAssert value == Sentinel

echo "Testing grid buffers, slices, reuse, and preflight validation"
checkFamily(initSimplex(1988))
checkFamily(initPerlin(1988))
checkFamily(initOpenSimplex2(1988))
checkFamily(initOpenSimplex2S(1988))
checkFamily(initValue(1988))
checkFamily(initValueCubic(1988))
checkFamily(initCellular(1988))

block:
  var noise = initCellular(1988)
  noise.frequency = 0.075'f
  for metric in CellularDistance:
    noise.distanceFunction = metric
    for mode in CellularReturn:
      noise.returnType = mode
      checkFill(noise, (0.125'f, -0.25'f), (3, 9))
      checkFill(noise, (0.125'f, -0.25'f, 0.5'f), (3, 2, 9))
  var storage: array[27, float32]
  for value in storage.mitems:
    value = Sentinel
  for jitter in [-0.1'f, 1.1'f, NaN.float32]:
    noise.jitter = jitter
    expectError: noise.grid((0, 0), (3, 9), storage)
    expectError: noise.grid((0, 0, 0), (1, 3, 9), storage)
    expectError: noise.grid((0, 0), (0, 9), storage)
    expectError: noise.grid((0, 0, 0), (0, 3, 9), storage)
  for value in storage:
    doAssert value == Sentinel

when (NimMajor, NimMinor) >= (2, 4):
  static:
    template checkVm(initial: untyped) =
      block:
        let noise = initial
        var buffer: array[9, float32]
        noise.grid((0, 0), (1, 9), buffer)
        for i, value in buffer:
          doAssert value == noise.value(0, i)
        noise.grid((0, 0, 0), (1, 1, 9), buffer)
        for i, value in buffer:
          doAssert value == noise.value(0, 0, i)
        when typeof(noise) is Simplex or typeof(noise) is Perlin:
          noise.grid(0, 9, buffer)
          for i, value in buffer:
            doAssert value == noise.value(i)
          noise.grid((0, 0, 0, 0), (1, 1, 9, 1), buffer)
          for i, value in buffer:
            doAssert value == noise.value(0, 0, i, 0)
    checkVm(initSimplex(42))
    checkVm(initPerlin(42))
    checkVm(initOpenSimplex2(42))
    checkVm(initOpenSimplex2S(42))
    checkVm(initValue(42))
    checkVm(initValueCubic(42))
    checkVm(initCellular(42))
