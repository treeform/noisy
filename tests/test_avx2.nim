import noisy

when defined(amd64) and (defined(gcc) or defined(clang) or defined(vcc)) and
    not defined(noisyNoAvx):
  import noisy/simd/runtimecheck
  let testingAvx2 = hasAvx2()
else:
  const testingAvx2 = false

when defined(noisyRequireAvx2):
  doAssert testingAvx2, "This test run requires AVX2 CPU and OS support"

template checkValue(actual, expected: float32) =
  if testingAvx2:
    doAssert actual == expected
  else:
    doAssert abs(actual - expected) < 0.0001'f

block grid_boundaries:
  # Cover empty grids, SSE2/AVX2 boundaries, every eight-lane remainder,
  # unaligned column starts, and dimensions too narrow for SSE2 tiles.
  const lengths = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
  for seed in [1, 1988, 8675309]:
    for octaves in [1, 3]:
      var s = initSimplex(seed)
      s.octaves = octaves
      s.frequency = 0.125
      s.amplitude = 0.75
      s.gain = 0.5
      s.lacunarity = 2
      for width in [0, 1, 3, 4, 5, 9]:
        for height in lengths:
          let g = s.grid((-257.25'f32, 253.5'f32), (width, height))
          doAssert (g.width, g.height, g.depth) == (width, height, 1)
          doAssert g.values.len == width * height
          for x in 0 ..< width:
            for y in 0 ..< height:
              checkValue(g[x, y], s.value(
                -257.25'f32 + x.float32, 253.5'f32 + y.float32))
        for height in [0, 1, 3, 4, 5, 9]:
          for depth in lengths:
            let g = s.grid(
              (-257.25'f32, 253.5'f32, -0.75'f32), (width, height, depth)
            )
            doAssert (g.width, g.height, g.depth) == (width, height, depth)
            doAssert g.values.len == width * height * depth
            for x in 0 ..< width:
              for y in 0 ..< height:
                for z in 0 ..< depth:
                  checkValue(g[x, y, z], s.value(
                    -257.25'f32 + x.float32, 253.5'f32 + y.float32,
                    -0.75'f32 + z.float32
                  ))

block fractional_coordinates:
  # Non-binary frequencies expose differences in coordinate scaling order.
  var s = initSimplex(12345)
  s.octaves = 4
  s.amplitude = 0.7
  s.lacunarity = 1.7
  s.gain = 0.3
  for frequency in [0'f32, 0.1'f32, -0.13'f32]:
    s.frequency = frequency
    let
      g2 = s.grid((-3.4'f32, 2.1'f32), (9, 23))
      g3 = s.grid((-3.4'f32, 2.1'f32, -0.7'f32), (5, 7, 23))
    for x in 0 ..< g2.width:
      for y in 0 ..< g2.height:
        let expected = s.value(-3.4'f32 + x.float32, 2.1'f32 + y.float32)
        if testingAvx2:
          doAssert g2[x, y] == expected
        else:
          doAssert abs(g2[x, y] - expected) < 0.00001'f32
    for x in 0 ..< g3.width:
      for y in 0 ..< g3.height:
        for z in 0 ..< g3.depth:
          let expected = s.value(
            -3.4'f32 + x.float32, 2.1'f32 + y.float32, -0.7'f32 + z.float32
          )
          if testingAvx2:
            doAssert g3[x, y, z] == expected
          else:
            doAssert abs(g3[x, y, z] - expected) < 0.00001'f32

block zero_octaves:
  var s = initSimplex(1988)
  s.octaves = 0
  for size in [1, 4, 8, 9, 16, 17]:
    var raised2d, raised3d: bool
    try:
      discard s.grid((0, 0), (size, size))
    except NoisyError:
      raised2d = true
    try:
      discard s.grid((0, 0, 0), (size, size, size))
    except NoisyError:
      raised3d = true
    doAssert raised2d and raised3d

block line_and_four_dimensions:
  for seed in [1, 1988, -42]:
    for octaves in [1, 3]:
      for frequency in [0'f, 0.125'f, -0.13'f]:
        var s = initSimplex(seed)
        s.octaves = octaves
        s.frequency = frequency
        s.gain = 0.3
        s.lacunarity = 1.7
        for size in 0 .. 17:
          let
            line = s.grid(-257.25'f, size)
            volume = s.grid((-1.25'f, 2.5'f, -0.75'f, 0.5'f), (2, 3, size, 2))
          doAssert line.values.len == size
          doAssert volume.values.len == 12 * size
          for x in 0 ..< size:
            checkValue(line[x], s.value(-257.25'f + x.float32))
          for w in 0 ..< 2:
            for x in 0 ..< 2:
              for y in 0 ..< 3:
                for z in 0 ..< size:
                  checkValue(volume[x, y, z, w], s.value(
                    -1.25'f + x.float32, 2.5'f + y.float32,
                    -0.75'f + z.float32, 0.5'f + w.float32))
