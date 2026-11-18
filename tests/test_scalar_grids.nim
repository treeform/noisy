import noisy/[opensimplexes, values]

template expectError(error: typedesc, body: untyped) =
  ## Checks that invalid input raises the specified exception.
  var caught = false
  try:
    body
  except error:
    caught = true
  doAssert caught

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
      doAssert plane[x, y] == expected
      doAssert plane.values[x * 5 + y] == expected
      for z in 0 ..< 7:
        let expected = noise.value(sampleX, sampleY, start[2] + z.float32)
        doAssert volume[x, y, z] == expected
        doAssert volume.values[(x * 5 + y) * 7 + z] == expected

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

echo "Testing OpenSimplex2S and ValueCubic grids against scalar sampling"
checkGrid(initOpenSimplex2S(-7))
checkGrid(initValueCubic(1988))
