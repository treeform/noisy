import std/[random, strformat], noisy

block:
  const
    fixedSimplex = initSimplex(42)
    fixedPerlin = initPerlin(42)
  doAssert fixedSimplex.value(0.125'f, 0.25'f) ==
    initSimplex(42).value(0.125'f, 0.25'f)
  doAssert fixedPerlin.value(0.125'f, 0.25'f) ==
    initPerlin(42).value(0.125'f, 0.25'f)

proc checkClose(actual, expected: float32) =
  doAssert abs(actual - expected) <= 0.00003'f,
    &"kernel mismatch: {actual} != {expected}"

proc checkPermutation[T: Simplex | Perlin](noise, normalized: T) =
  # Complete AVX2/SSE2/NEON groups and scalar tails, in every dimension.
  for start in [-256.125'f, -0.00001'f, 255.875'f]:
    let
      line = noise.grid(start, 19)
      plane = noise.grid((start, -start), (2, 19))
      volume = noise.grid((start, -start, start), (2, 2, 19))
      hyper = noise.grid((start, -start, start, -start), (2, 2, 19, 2))
    for i in 0 ..< 19:
      let z = start + i.float32
      checkClose(line[i], normalized.value(z))
      checkClose(noise.value(z), normalized.value(z))
      for a in 0 ..< 2:
        let x = start + a.float32
        checkClose(plane[a, i], normalized.value(x, -start + i.float32))
        checkClose(noise.value(x, z), normalized.value(x, z))
        for b in 0 ..< 2:
          let y = -start + b.float32
          checkClose(volume[a, b, i], normalized.value(x, y, z))
          checkClose(noise.value(x, y, z), normalized.value(x, y, z))
          for c in 0 ..< 2:
            let w = -start + c.float32
            checkClose(hyper[a, b, i, c], normalized.value(x, y, z, w))
            checkClose(noise.value(x, y, z, w), normalized.value(x, y, z, w))

echo "Testing bounded lookups with arbitrary public permutation contents"
block:
  const Edges = [int32.low, int32.high, -513'i32, -512, -1, 0, 11, 12,
    255, 256, 511, 512, 513]
  var rng = initRand(1988)
  for trial in 0 ..< 8:
    var simplex = initSimplex(trial + 1)
    var perlin = initPerlin(trial + 1)
    for i in 0 ..< 512:
      let entry = if trial == 0: Edges[i mod Edges.len]
        else: cast[int32](rng.next().uint32)
      simplex.perm[i] = entry
      simplex.permMod12[i] = Edges[(i + trial) mod Edges.len]
      perlin.perm[i] = entry
    simplex.octaves = 3
    simplex.frequency = -0.125'f
    perlin.octaves = 3
    perlin.frequency = -0.125'f
    var normalizedSimplex = simplex
    var normalizedPerlin = perlin
    for i in 0 ..< 512:
      normalizedSimplex.perm[i] = simplex.perm[i] and 511
      normalizedPerlin.perm[i] = perlin.perm[i] and 511
      if simplex.permMod12[i] < 0 or simplex.permMod12[i] >= 12:
        normalizedSimplex.permMod12[i] = 0
    checkPermutation(simplex, normalizedSimplex)
    checkPermutation(perlin, normalizedPerlin)

echo "Testing hashed table bounds across seeds and lattice boundaries"
block:
  const Coordinates = [-65536.125'f, -256.00003'f, -256, -255.99997'f,
    -1.0000001'f, -1, -0.99999994'f, -0.000001'f, 0, 0.000001'f,
    0.99999994'f, 1, 1.0000001'f, 255.99997'f, 256, 256.00003'f, 65536.125'f]
  template exercise(initial: untyped) =
    block:
      var noise = initial
      noise.octaves = 3
      noise.frequency = 0.125'f
      for x in Coordinates:
        let
          plane = noise.grid((x, -x), (1, 9))
          volume = noise.grid((x, -x, x), (1, 1, 9))
        for i in 0 ..< 9:
          checkClose(plane[0, i], noise.value(x, -x + i.float32))
          checkClose(volume[0, 0, i], noise.value(x, -x, x + i.float32))
  for seed in [int32.low.int, -1, 0, 1, int32.high.int]:
    exercise(initOpenSimplex2(seed))
    exercise(initOpenSimplex2S(seed))
    exercise(initValue(seed))
    exercise(initValueCubic(seed))
    exercise(initCellular(seed))
