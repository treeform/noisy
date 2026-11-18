import
  std/[random, strformat],
  noisy

type Point = array[4, float32]

proc sample(noise: Perlin, point: Point, dimensions: int): float32 =
  ## Calls the public scalar overload for the selected dimension.
  case dimensions
  of 1: noise.value(point[0])
  of 2: noise.value(point[0], point[1])
  of 3: noise.value(point[0], point[1], point[2])
  else: noise.value(point[0], point[1], point[2], point[3])

proc checkClose(actual, expected: float32, context: string) =
  ## Allows float32 rounding while rejecting nonfinite samples.
  doAssert abs(actual) < Inf and abs(expected) < Inf, context
  let tolerance = 0.00001'f + 0.00001'f * max(abs(actual), abs(expected))
  doAssert abs(actual - expected) <= tolerance,
    context & &" actual={actual}, expected={expected}"

echo "Testing Perlin samples against recorded C reference outputs"
block:
  # Generated with Gustavson's public-domain noise1234.c at revision
  # a624f5a18395c0ebc04738191efddf1dbd2cad9d of stegu/perlin-noise.
  # Its mutable permutation was set to 0..255 repeated; kernels were unchanged.
  var noise = initPerlin(1)
  for i in 0 ..< noise.perm.len:
    noise.perm[i] = (i and 255).int32
  for recorded in [
    ([0.125'f, 0.25'f, 0.5'f, 0.75'f],
      [0.0178415831'f, 0.206157118'f, 0.207576424'f, 0.0773944035'f]),
    ([-0.75'f, 1.25'f, -2.5'f, 3.75'f],
      [-0.351673812'f, 0.286150664'f, 0.152908191'f, 0.0591864958'f]),
    ([255.875'f, -256.25'f, 511.5'f, -0.125'f],
      [-0.0442475267'f, -0.317119986'f, -0.231037319'f, -0.37563175'f]),
    ([0.0'f, 0.25'f, 0.5'f, 0.75'f],
      [0.0'f, 0.148535162'f, 0.20977734'f, 0.080736123'f]),
    ([-1.0'f, 0.25'f, 0.5'f, 0.75'f],
      [0.0'f, -0.19235304'f, -0.282445312'f, 0.097165361'f]),
    ([0.25'f, 0.25'f, 0.25'f, 0.25'f],
      [0.0129433591'f, 0.286150664'f, 0.338746816'f, 0.411872029'f]),
    ([0.75'f, 0.125'f, -0.25'f, -0.5'f],
      [-0.0696738288'f, 0.284088641'f, 0.33741951'f, -0.11714457'f])
  ]:
    for dimensions in 1 .. 4:
      checkClose(noise.sample(recorded[0], dimensions),
        recorded[1][dimensions - 1], &"{dimensions}D recorded C reference")

echo "Testing Perlin initialization, scalar overloads, and determinism"
block:
  var rng = initRand(7419)
  for seed in [-12345, 1, 42, 1988]:
    let
      noise = initPerlin(seed = seed)
      repeated = initPerlin(seed)
      other = initPerlin(seed + 1)
    doAssert noise.octaves == 1
    doAssert noise.amplitude == 1 and noise.frequency == 1
    doAssert noise.lacunarity == 2 and noise.gain == 0.5'f
    doAssert noise.perm == repeated.perm
    doAssert noise.perm == initSimplex(seed).perm
    doAssert noise.perm != other.perm
    var seen: array[256, bool]
    for i in 0 ..< 256:
      let entry = noise.perm[i]
      doAssert entry >= 0 and entry < 256
      doAssert not seen[entry]
      seen[entry] = true
      doAssert noise.perm[i + 256] == entry
    for dimensions in 1 .. 4:
      var varied, differentSeed: bool
      var axisVaries: array[4, bool]
      for i in 0 ..< 128:
        var point: Point
        for axis in 0 ..< dimensions:
          point[axis] = (rng.rand(64.0) - 32).float32
        let actual = noise.sample(point, dimensions)
        doAssert abs(actual) < Inf
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
    var scaled = noise
    scaled.frequency = 0.125'f
    for x in -8 .. 8:
      doAssert scaled.value(x) == scaled.value(x.float32)
      doAssert scaled.value(x, 2) == scaled.value(x.float32, 2.0'f)
      doAssert scaled.value(x, 2, -3) ==
        scaled.value(x.float32, 2.0'f, -3.0'f)
      doAssert scaled.value(x, 2, -3, 5) ==
        scaled.value(x.float32, 2.0'f, -3.0'f, 5.0'f)

echo "Testing Perlin lattice zeros, periodicity, and continuity"
block:
  const Epsilon = 1.0'f / 4096
  for seed in [1, 42, 1988]:
    let noise = initPerlin(seed)
    for dimensions in 1 .. 4:
      for coordinate in [-257, -256, -1, 0, 1, 255, 256, 257]:
        let point = [coordinate.float32, -3.0'f, 7.0'f, -9.0'f]
        doAssert noise.sample(point, dimensions) == 0
      for point in [
        [-3.125'f, 1.75'f, -4.5'f, 2.25'f],
        [255.875'f, -256.25'f, 511.5'f, -0.125'f]
      ]:
        let expected = noise.sample(point, dimensions)
        for axis in 0 ..< dimensions:
          for direction in [-1.0'f, 1.0'f]:
            var shifted = point
            shifted[axis] += direction * 256
            checkClose(noise.sample(shifted, dimensions), expected,
              &"{dimensions}D periodicity, seed={seed}, axis={axis}")
      for axis in 0 ..< dimensions:
        for boundary in [-256.0'f, -1.0'f, 0.0'f, 1.0'f, 256.0'f]:
          var point = [0.125'f, -0.25'f, 0.5'f, 0.75'f]
          point[axis] = boundary
          let center = noise.sample(point, dimensions)
          point[axis] = boundary - Epsilon
          let left = noise.sample(point, dimensions)
          point[axis] = boundary + Epsilon
          let right = noise.sample(point, dimensions)
          doAssert abs(left - center) < 0.003'f
          doAssert abs(right - center) < 0.003'f
          doAssert abs(left - right) < 0.003'f

echo "Testing Perlin octave composition and invalid octave counts"
block:
  for dimensions in 1 .. 4:
    for parameters in [
      (1, 0.75'f, 0.25'f, 2.0'f, 0.5'f),
      (3, 1.25'f, -0.75'f, 1.5'f, -0.5'f),
      (4, -0.8'f, 1.5'f, 2.25'f, 0.35'f),
      (3, 0.0'f, 0.75'f, 2.0'f, 0.5'f),
      (3, 1.0'f, 0.0'f, 2.0'f, 0.5'f),
      (3, 1.0'f, 0.75'f, 2.0'f, 0.0'f)
    ]:
      var noise = initPerlin(1988)
      noise.octaves = parameters[0]
      noise.amplitude = parameters[1]
      noise.frequency = parameters[2]
      noise.lacunarity = parameters[3]
      noise.gain = parameters[4]
      let base = initPerlin(1988)
      for point in [
        [0.125'f, 0.25'f, -0.5'f, 0.75'f],
        [-3.75'f, 2.125'f, 1.5'f, -0.25'f]
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
          &"{dimensions}D octave composition, parameters={parameters}")
    for octaves in [0, -1, -3]:
      var noise = initPerlin(1988)
      noise.octaves = octaves
      var caught = false
      try:
        discard noise.sample([0.125'f, 0.25'f, 0.5'f, 0.75'f], dimensions)
      except NoisyError:
        caught = true
      doAssert caught, &"{dimensions}D accepted octaves={octaves}"

echo "Testing Perlin grids against point sampling"
block:
  for seed in [1, 42, 1988]:
    var noise = initPerlin(seed)
    if seed == 42:
      # A caller can change the two halves independently.
      for i in 256 ..< noise.perm.len:
        noise.perm[i] = (noise.perm[i] + 37) and 255
    for settings in [(1, 1.0'f), (3, -0.75'f), (4, 0.125'f)]:
      noise.octaves = settings[0]
      noise.frequency = settings[1]
      noise.amplitude = 0.8'f
      noise.lacunarity = 1.75'f
      noise.gain = -0.6'f
      var starts = @[
        [-1.375'f, 2.125'f, -0.625'f, 0.375'f],
        [255.875'f, -256.25'f, 511.5'f, -256.0'f]
      ]
      when sizeof(int) == 8:
        starts.add [
          [0.125'f, 3_000_000_000.0'f, 0.625'f, -4_000_000_000.0'f],
          [-3_000_000_000.0'f, 0.125'f, 4_000_000_000.0'f, 0.875'f]
        ]
      for start in starts:
        for span in [3, 4, 9]:
          for dimensions in 1 .. 4:
            let grid = case dimensions
              of 1: noise.grid(start[0], span)
              of 2: noise.grid((start[0], start[1]), (3, span))
              of 3: noise.grid((start[0], start[1], start[2]), (3, 2, span))
              else: noise.grid(
                (start[0], start[1], start[2], start[3]), (3, 2, span, 2)
              )
            let expectedShape = case dimensions
              of 1: [span, 1, 1, 1]
              of 2: [3, span, 1, 1]
              of 3: [3, 2, span, 1]
              else: [3, 2, span, 2]
            doAssert [grid.width, grid.height, grid.depth, grid.frames] ==
              expectedShape
            var index = 0
            for w in 0 ..< grid.frames:
              for x in 0 ..< grid.width:
                for y in 0 ..< grid.height:
                  for z in 0 ..< grid.depth:
                    let point = [
                      start[0] + x.float32, start[1] + y.float32,
                      start[2] + z.float32, start[3] + w.float32
                    ]
                    let indexed = case dimensions
                      of 1: grid[x]
                      of 2: grid[x, y]
                      of 3: grid[x, y, z]
                      else: grid[x, y, z, w]
                    doAssert indexed == grid.values[index]
                    checkClose(indexed, noise.sample(point, dimensions),
                      &"Perlin {dimensions}D, seed={seed}, span={span}, point={point}")
                    inc index
            doAssert grid.values.len == index

  # Both same-cell and mixed-cell batches must retain signed and full-width
  # permutation values whenever the resulting indexes are in bounds.
  for base in [-8'i32, 300'i32]:
    var noise = initPerlin(1)
    for i in 0 ..< noise.perm.len:
      noise.perm[i] = base + int32(i mod 7)
    noise.frequency = 0.125'f
    let samples = noise.grid((97.0'f, 90.0'f, 99.0'f), (1, 1, 12))
    for z in 0 ..< samples.depth:
      checkClose(samples[0, 0, z], noise.value(97.0'f, 90.0'f, 99.0'f + z.float32),
        &"Perlin custom permutation, base={base}, z={z}")

echo "Testing Perlin grid overloads and dimensions"
block:
  var noise = initPerlin(1)
  noise.frequency = 0.125'f
  doAssert noise.grid(-1, 5).values == noise.grid(-1.0'f, 5).values
  doAssert noise.grid((-1, 2), (3, 5)).values ==
    noise.grid((-1.0'f, 2.0'f), (3, 5)).values
  doAssert noise.grid((-1, 2, -3), (2, 3, 5)).values ==
    noise.grid((-1.0'f, 2.0'f, -3.0'f), (2, 3, 5)).values
  doAssert noise.grid((-1, 2, -3, 4), (2, 3, 5, 2)).values ==
    noise.grid((-1.0'f, 2.0'f, -3.0'f, 4.0'f), (2, 3, 5, 2)).values
  doAssert noise.grid(0, 0).values.len == 0
  doAssert noise.grid((0, 0, 0, 0), (2, 3, 5, 0)).values.len == 0

  template expectGridError(expression: untyped) =
    block:
      var caught = false
      try:
        discard expression
      except NoisyError:
        caught = true
      doAssert caught

  expectGridError(noise.grid(0, -1))
  expectGridError(noise.grid((0, 0, 0, 0), (2, 3, 0, -1)))
  expectGridError(noise.grid(0, int.high))
  expectGridError(noise.grid((0, 0, 0, 0), (int.high, 3, 5, 2)))
  var invalid = noise
  invalid.octaves = 0
  expectGridError(invalid.grid(0, 0))
  expectGridError(invalid.grid((0, 0, 0, 0), (2, 3, 5, 2)))
