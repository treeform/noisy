import
  std/strformat,
  noisy/[opensimplexes, values]

from noisy/grids import useGridNeon

template expectError(error: typedesc, body: untyped) =
  ## Checks that invalid input raises the specified exception.
  var caught = false
  try:
    body
  except error:
    caught = true
  doAssert caught

template checkSample(
  actual, expected: float32, offset, size: int, context: string
) =
  ## Allows float32 rounding in complete NEON groups and checks tails exactly.
  let
    sample = actual
    point = expected
    error = abs(sample - point)
  var tolerance = 0.0'f
  when useGridNeon:
    if offset < size div 4 * 4:
      tolerance = 0.00001'f + 0.00001'f * max(abs(sample), abs(point))
  doAssert abs(sample) < Inf and abs(point) < Inf and error <= tolerance,
    context & " grid=" & $sample & ", point=" & $point & ", error=" & $error &
      ", tolerance=" & $tolerance

proc checkGrid[T](initial: T) =
  ## Checks the scalar sampler, flat layout, and shape validation.
  var noise = initial
  noise.octaves = 3
  noise.frequency = -0.37'f
  noise.amplitude = 1.2'f
  noise.lacunarity = 1.7'f
  noise.gain = 0.6'f

  let
    start = (-5.375'f, 8.25'f, -0.125'f)
    plane: Grid = noise.grid((start[0], start[1]), (3, 5))
    volume: Grid = noise.grid(start, (3, 5, 7))
  doAssert (plane.width, plane.height, plane.depth, plane.frames) == (3, 5, 1, 1)
  doAssert (volume.width, volume.height, volume.depth, volume.frames) ==
    (3, 5, 7, 1)
  doAssert plane.values.len == 15 and volume.values.len == 105
  for x in 0 ..< 3:
    for y in 0 ..< 5:
      let
        sampleX = start[0] + x.float32
        sampleY = start[1] + y.float32
        expected = noise.value(sampleX, sampleY)
      checkSample(plane[x, y], expected, y, plane.height,
        &"{$T} plane offset=({x}, {y})")
      doAssert plane.values[x * 5 + y] == plane[x, y]
      for z in 0 ..< 7:
        let expected = noise.value(sampleX, sampleY, start[2] + z.float32)
        checkSample(volume[x, y, z], expected, z, volume.depth,
          &"{$T} volume offset=({x}, {y}, {z})")
        doAssert volume.values[(x * 5 + y) * 7 + z] == volume[x, y, z]

  doAssert noise.grid((-5, 8), (3, 5)).values ==
    noise.grid((-5.0'f, 8.0'f), (3, 5)).values
  doAssert noise.grid((-5, 8, -1), (3, 5, 7)).values ==
    noise.grid((-5.0'f, 8.0'f, -1.0'f), (3, 5, 7)).values
  expectError(IndexDefect):
    discard plane[3, 0]
  expectError(IndexDefect):
    discard volume[0, 0, -1]

  for dimens in [(0, 5), (int.high, 0)]:
    let empty = noise.grid((0, 0), dimens)
    doAssert empty.values.len == 0
    doAssert (empty.width, empty.height) == dimens
    doAssert empty.depth == 1 and empty.frames == 1
  for dimens in [(0, 2, 3), (2, 0, 3), (int.high, 2, 0)]:
    let empty = noise.grid((0, 0, 0), dimens)
    doAssert empty.values.len == 0
    doAssert (empty.width, empty.height, empty.depth) == dimens
    doAssert empty.frames == 1

  const MaxPayload = int.high div sizeof(float32)
  for dimens in [(-1, 5), (0, -1), (int.high, 5), (MaxPayload, 1)]:
    expectError(NoisyError):
      discard noise.grid((0, 0), dimens)
  for dimens in [
    (-1, 2, 3), (2, -1, 3), (0, 2, -1), (2, 3, int.high), (1, 1, MaxPayload)
  ]:
    expectError(NoisyError):
      discard noise.grid((0, 0, 0), dimens)
  for octaves in [-1, 0]:
    noise.octaves = octaves
    expectError(NoisyError):
      discard noise.grid((0, 0), (0, 0))
    expectError(NoisyError):
      discard noise.grid((0, 0, 0), (1, 2, 3))

proc compareGrid[T](noise: T, start: (float32, float32), height: int) =
  ## Compares both rows of a plane with point samples.
  let grid = noise.grid(start, (2, height))
  for x in 0 ..< grid.width:
    for y in 0 ..< height:
      let expected = noise.value(start[0] + x.float32, start[1] + y.float32)
      checkSample(grid[x, y], expected, y, height,
        &"{$T} 2D seed={noise.seed} octaves={noise.octaves} start={start} " &
          &"size={height} offset=({x}, {y})")

proc compareGrid[T](noise: T, start: (float32, float32, float32), depth: int) =
  ## Compares all depth lines of a volume with point samples.
  let grid = noise.grid(start, (2, 2, depth))
  for x in 0 ..< grid.width:
    for y in 0 ..< grid.height:
      for z in 0 ..< depth:
        let expected = noise.value(
          start[0] + x.float32,
          start[1] + y.float32,
          start[2] + z.float32
        )
        checkSample(grid[x, y, z], expected, z, depth,
          &"{$T} 3D seed={noise.seed} octaves={noise.octaves} start={start} " &
            &"size={depth} offset=({x}, {y}, {z})")

proc checkGridBatches[T](initial: T) =
  ## Exercises complete groups, all tail lengths, and octave coordinate order.
  for settings in [
    (1, 1.0'f, 1.0'f, 2.0'f, 0.5'f),
    (3, -0.37'f, 1.2'f, 1.7'f, -0.6'f),
    (2, 0.0'f, -1.5'f, 2.0'f, 0.7'f)
  ]:
    var noise = initial
    (noise.octaves, noise.frequency, noise.amplitude, noise.lacunarity,
      noise.gain) = settings
    # At the float32 integer precision limit, rounding the batch start before
    # adding the lane can produce a different coordinate from point sampling.
    for start in [(-5.37'f, 9.13'f, -2.61'f), (0.25'f, 16777215.0'f, -0.5'f)]:
      for size in 1 .. 11:
        compareGrid(noise, (start[0], start[1]), size)
        compareGrid(noise, (start[0], start[2], start[1]), size)

proc checkGridBoundaries[T](noise: T) =
  ## Places a lattice boundary in every lane of two groups and all three tails.
  const
    X = -0.125'f
    Y = 0.25'f
    Skew = 0.5'f * (1.73205077648162841796875'f - 1)
  for boundary in [-256.0'f, -1.0'f, -0.5'f, 0.0'f, 0.5'f, 255.0'f]:
    let
      y =
        when T is OpenSimplex2S: (boundary - Skew * X) / (1 + Skew)
        else: boundary
      z =
        when T is OpenSimplex2S: 2 * (X + Y) - 3 * boundary
        else: boundary
    for delta in [-0.0001'f, 0.0'f, 0.0001'f]:
      for offset in 0 .. 10:
        compareGrid(noise, (X, y + delta - offset.float32), 11)
        compareGrid(noise, (X, Y, z + delta - offset.float32), 11)

echo "Testing OpenSimplex2S and ValueCubic grids against scalar sampling"
when (NimMajor, NimMinor) >= (2, 4):
  # Nim 2.2's VM cannot evaluate the existing scalar kernels.
  static:
    let
      smooth = initOpenSimplex2S(-7)
      cubic = initValueCubic(1988)
      plane = smooth.grid((0, 0), (1, 4))
      volume = cubic.grid((0, 0, 0), (1, 1, 4))
    doAssert plane[0, 3] == smooth.value(0, 3)
    doAssert volume[0, 0, 3] == cubic.value(0, 0, 3)

checkGrid(initOpenSimplex2S(-7))
checkGrid(initValueCubic(1988))
for seed in [-7, 0, low(int32).int, high(int32).int]:
  let
    smooth = initOpenSimplex2S(seed)
    cubic = initValueCubic(seed)
  checkGridBatches(smooth)
  checkGridBatches(cubic)
  checkGridBoundaries(smooth)
  checkGridBoundaries(cubic)
