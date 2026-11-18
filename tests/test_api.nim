import noisy

type
  Initialize = proc(seed: int): Simplex {.
    nimcall, noSideEffect, gcsafe, raises: [].}
  Sample2f = proc(simplex: Simplex, x, y: float32): float32 {.
    nimcall, noSideEffect, gcsafe, raises: [NoisyError].}
  Sample3f = proc(simplex: Simplex, x, y, z: float32): float32 {.
    nimcall, noSideEffect, gcsafe, raises: [NoisyError].}
  Sample2i = proc(simplex: Simplex, x, y: int): float32 {.
    inline, noSideEffect, gcsafe, raises: [NoisyError].}
  Sample3i = proc(simplex: Simplex, x, y, z: int): float32 {.
    inline, noSideEffect, gcsafe, raises: [NoisyError].}
  Grid2f = proc(simplex: Simplex, start: (float32, float32),
    dimens: (int, int)): Grid {.
    nimcall, noSideEffect, gcsafe, raises: [NoisyError].}
  Grid3f = proc(simplex: Simplex, start: (float32, float32, float32),
    dimens: (int, int, int)): Grid {.
    nimcall, noSideEffect, gcsafe, raises: [NoisyError].}
  Grid2i = proc(simplex: Simplex, start: (int, int),
    dimens: (int, int)): Grid {.
    inline, noSideEffect, gcsafe, raises: [NoisyError].}
  Grid3i = proc(simplex: Simplex, start: (int, int, int),
    dimens: (int, int, int)): Grid {.
    inline, noSideEffect, gcsafe, raises: [NoisyError].}
  Index = proc(g: Grid, x, y, z: int): float32 {.
    inline, noSideEffect, gcsafe, raises: [].}

proc checkLegacyApi() {.raises: [NoisyError], gcsafe.} =
  ## Exercises the API available before native 1D, 4D, and NEON support.
  let
    initialize: Initialize = initSimplex
    sample2f: Sample2f = value
    sample3f: Sample3f = value
    sample2i: Sample2i = value
    sample3i: Sample3i = value
    grid2f: Grid2f = grid
    grid3f: Grid3f = grid
    grid2i: Grid2i = grid
    grid3i: Grid3i = grid
    index: Index = `[]`
  var noise = initSimplex(seed = 1988)
  doAssert noise.perm == initialize(1988).perm
  noise.octaves = 3
  noise.frequency = 0.125'f
  noise.amplitude = 0.75'f
  noise.lacunarity = 1.5'f
  noise.gain = 0.6'f
  static:
    doAssert typeof(noise.octaves) is int
    doAssert typeof(noise.frequency) is float32
    doAssert typeof(noise.amplitude) is float32
    doAssert typeof(noise.lacunarity) is float32
    doAssert typeof(noise.gain) is float32
    doAssert typeof(noise.perm) is array[512, int32]
    doAssert typeof(noise.permMod12) is array[512, int32]
  doAssert value(simplex = noise, x = 1.25'f, y = -2.5'f) ==
    sample2f(noise, 1.25'f, -2.5'f)
  doAssert value(simplex = noise, x = 1.25'f, y = -2.5'f, z = 3.75'f) ==
    sample3f(noise, 1.25'f, -2.5'f, 3.75'f)
  doAssert sample2i(noise, 1, -2) == noise.value(1.0'f, -2.0'f)
  doAssert sample3i(noise, 1, -2, 3) == noise.value(1.0'f, -2.0'f, 3.0'f)
  doAssert value(simplex = noise, x = 1, y = -2) == sample2i(noise, 1, -2)
  doAssert value(simplex = noise, x = 1, y = -2, z = 3) ==
    sample3i(noise, 1, -2, 3)
  let
    plane = grid(simplex = noise, start = (1.25'f, -2.5'f), dimens = (5, 7))
    volume = noise.grid(start = (1.25'f, -2.5'f, 3.75'f), dimens = (5, 6, 7))
  doAssert plane.values == grid2f(noise, (1.25'f, -2.5'f), (5, 7)).values
  doAssert volume.values ==
    grid3f(noise, (1.25'f, -2.5'f, 3.75'f), (5, 6, 7)).values
  doAssert grid2i(noise, (1, -2), (5, 7)).values ==
    noise.grid((1.0'f, -2.0'f), (5, 7)).values
  doAssert grid3i(noise, (1, -2, 3), (5, 6, 7)).values ==
    noise.grid((1.0'f, -2.0'f, 3.0'f), (5, 6, 7)).values
  doAssert grid(simplex = noise, start = (1, -2), dimens = (5, 7)).values ==
    grid2i(noise, (1, -2), (5, 7)).values
  doAssert grid(simplex = noise, start = (1, -2, 3),
    dimens = (5, 6, 7)).values == grid3i(noise, (1, -2, 3), (5, 6, 7)).values
  static:
    doAssert typeof(plane.width) is int
    doAssert typeof(plane.height) is int
    doAssert typeof(plane.depth) is int
    doAssert typeof(plane.values) is seq[float32]
  doAssert (plane.width, plane.height, plane.depth) == (5, 7, 1)
  doAssert (volume.width, volume.height, volume.depth) == (5, 6, 7)
  doAssert plane[2, 3] == plane[2, 3, 0]
  doAssert plane[2, 3] == index(plane, 2, 3, 0)
  doAssert plane[2, 3] == plane.values[2 * plane.height + 3]
  doAssert volume[2, 3, 4] ==
    volume.values[(2 * volume.height + 3) * volume.depth + 4]

echo "Testing the legacy 2D/3D public API contract"
checkLegacyApi()

block:
  # Original public fields remain usable in object constructors.
  let
    original = initSimplex(42)
    copied = Simplex(
      octaves: original.octaves,
      amplitude: original.amplitude,
      frequency: original.frequency,
      lacunarity: original.lacunarity,
      gain: original.gain,
      perm: original.perm,
      permMod12: original.permMod12
    )
    plane = Grid(width: 2, height: 3, depth: 1,
      values: @[0.0'f, 1, 2, 3, 4, 5])
    volume = Grid(width: 2, height: 2, depth: 2,
      values: @[0.0'f, 1, 2, 3, 4, 5, 6, 7])
  doAssert copied.value(1.25'f, -2.5'f) == original.value(1.25'f, -2.5'f)
  doAssert copied.value(1.25'f, -2.5'f, 3.75'f) ==
    original.value(1.25'f, -2.5'f, 3.75'f)
  doAssert plane[1, 2] == 5 and plane[1, 2, 0] == 5
  doAssert volume[1, 0, 1] == 5 and volume[1, 1, 1] == 7

block:
  var noise = initSimplex(42)
  noise.octaves = 0
  for dimensions in [2, 3]:
    var caught = false
    try:
      if dimensions == 2:
        discard noise.value(0, 0)
      else:
        discard noise.value(0, 0, 0)
    except NoisyError:
      caught = true
    doAssert caught
