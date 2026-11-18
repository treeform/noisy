import ../lattices, avx2base, latticesavx2
include avx2flags

type
  LatticeAxis = object
    primed: array[3, Int8]
    delta: array[3, Float8]
  NearestFeatures* = object
    hash*: Int8
    first*, second*: Float8

func latticeAxis(point: Float8, prime: uint32): LatticeAxis {.avxInline.} =
  ## Matches scalar rounding before constructing the three cells of this axis.
  let
    negative = less(point, float8(0))
    rounding = select(negative, float8(-0.5'f), float8(0.5'f))
    center = truncate(point + rounding)
    stride = uints8(prime)
    primed = mm256_mullo_epi32(center, stride)
  result.primed = [mm256_sub_epi32(primed, stride), primed, mm256_add_epi32(primed, stride)]
  result.delta = [
    toFloats(mm256_add_epi32(center, int8(-1))) - point,
    toFloats(center) - point,
    toFloats(mm256_add_epi32(center, int8(1))) - point
  ]

func featureOffsets[N: static[int]](
  hash: Int8, delta: array[N, Float8], jitter: Float8
): array[N, Float8] {.avxInline.} =
  ## Gathers the jitter vectors selected by eight independent feature hashes.
  when N == 2:
    let indexes = mm256_and_si256(hash, uints8(510))
  else:
    let indexes = mm256_and_si256(hash, uints8(1020))
  for axis in 0 ..< N:
    let component =
      when N == 2: gather(RandVecs2D[axis].unsafeAddr, indexes, 4)
      else: gather(RandVecs3D[axis].unsafeAddr, indexes, 4)
    result[axis] = delta[axis] + component * jitter

func searchFeatures8*[N: static[int]](
  seed: int32, jitter: float32, point: array[N, Float8],
  squared, manhattan, trackSecond, trackHash: static[bool]
): NearestFeatures {.avxInline.} =
  ## Searches each lane's neighboring features in the scalar traversal order.
  let
    xAxis = latticeAxis(point[0], PrimeX)
    yAxis = latticeAxis(point[1], PrimeY)
    jitterScale = when N == 2: 0.43701595'f else: 0.39614353'f
    jitter8 = float8(jitter * jitterScale)
  when N == 3:
    let zAxis = latticeAxis(point[2], PrimeZ)
  result = NearestFeatures(hash: uints8(0),
    first: float8(1.0e10'f), second: float8(1.0e10'f))

  template consider(lattice: Int8, delta: array[N, Float8]) =
    let
      hash = hash8(seed, lattice)
      offset = featureOffsets(hash, delta, jitter8)
    when squared:
      let squaredXY = offset[0] * offset[0] + offset[1] * offset[1]
      let euclideanSquared =
        when N == 2: squaredXY
        else: squaredXY + offset[2] * offset[2]
    when manhattan:
      let manhattanXY = abs8(offset[0]) + abs8(offset[1])
      let manhattanDistance =
        when N == 2: manhattanXY
        else: manhattanXY + abs8(offset[2])
    let distance =
      when squared and manhattan: manhattanDistance + euclideanSquared
      elif squared: euclideanSquared
      else: manhattanDistance
    when trackSecond:
      result.second = mm256_max_ps(mm256_min_ps(result.second, distance), result.first)
    when trackHash:
      # Tied distances retain the first feature encountered, as in the scalar
      # search. Use this mask before replacing its nearest distance.
      result.hash = select(less(distance, result.first), hash, result.hash)
    result.first = mm256_min_ps(distance, result.first)

  for x in 0 .. 2:
    for y in 0 .. 2:
      let xy = mm256_xor_si256(xAxis.primed[x], yAxis.primed[y])
      when N == 2:
        consider(xy, [xAxis.delta[x], yAxis.delta[y]])
      else:
        for z in 0 .. 2:
          consider(mm256_xor_si256(xy, zAxis.primed[z]),
            [xAxis.delta[x], yAxis.delta[y], zAxis.delta[z]])
