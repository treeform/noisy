import nimsimd/neon, ../[lattices, simd], latticesneon

type
  RandomPairs {.importc: "float32x4x2_t", header: "arm_neon.h".} = object
    val: array[2, Float4]
  RandomTriples {.importc: "float32x4x3_t", header: "arm_neon.h".} = object
    val: array[3, Float4]
  LatticeAxis = object
    primed: array[3, uint32x4]
    delta: array[3, Float4]
  NearestFeatures* = object
    hash*: uint32x4
    first*, second*: Float4

func loadPairs(
  p: ptr float32, values: RandomPairs, lane: int32
): RandomPairs {.importc: "vld2q_lane_f32", header: "arm_neon.h".}
  ## Loads a feature's two jitter components into one lane.

func loadTriples(
  p: ptr float32, values: RandomTriples, lane: int32
): RandomTriples {.importc: "vld3q_lane_f32", header: "arm_neon.h".}
  ## Loads a feature's three jitter components into one lane.

func abs4(a: Float4): Float4 {.importc: "vabsq_f32", header: "arm_neon.h".}

func sqrt4*(a: Float4): Float4 {.importc: "vsqrtq_f32", header: "arm_neon.h".}

func latticeAxis(point: Float4, prime: uint32): LatticeAxis {.inline.} =
  ## Matches scalar rounding before constructing the three cells of this axis.
  let
    negative = vcltq_f32(point, float4(0))
    rounding = vbslq_f32(negative, float4(-0.5'f), float4(0.5'f))
    center = truncate(point + rounding)
    stride = vmovq_n_u32(prime)
    primed = multiply4(cast[uint32x4](center), stride)
  result.primed = [vsubq_u32(primed, stride), primed, vaddq_u32(primed, stride)]
  result.delta = [
    toFloats(addInts(center, int4(-1))) - point,
    toFloats(center) - point,
    toFloats(addInts(center, int4(1))) - point
  ]

func featureOffsets[N: static[int]](
  hash: uint32x4, delta: array[N, Float4], jitter: Float4
): array[N, Float4] {.inline.} =
  ## Gathers the jitter vectors selected by four independent feature hashes.
  when N == 2:
    let indexes = vandq_u32(hash, vmovq_n_u32(510))
    var values = RandomPairs(val: [float4(0), float4(0)])
    template loadLane(lane: static[int]) =
      values = loadPairs(RandVecs2D[vgetq_lane_u32(indexes, lane)].unsafeAddr,
        values, lane)
  else:
    let indexes = vandq_u32(hash, vmovq_n_u32(1020))
    var values = RandomTriples(val: [float4(0), float4(0), float4(0)])
    template loadLane(lane: static[int]) =
      values = loadTriples(RandVecs3D[vgetq_lane_u32(indexes, lane)].unsafeAddr,
        values, lane)
  loadLane(0)
  loadLane(1)
  loadLane(2)
  loadLane(3)
  for axis in 0 ..< N:
    result[axis] = delta[axis] + values.val[axis] * jitter

func searchFeatures4*[N: static[int]](
  seed: int32, jitter: float32, point: array[N, Float4],
  squared, manhattan, trackSecond, trackHash: static[bool]
): NearestFeatures =
  ## Searches each lane's neighboring features in the scalar traversal order.
  let
    xAxis = latticeAxis(point[0], PrimeX)
    yAxis = latticeAxis(point[1], PrimeY)
    jitterScale = when N == 2: 0.43701595'f else: 0.39614353'f
    jitter4 = float4(jitter * jitterScale)
  when N == 3:
    let zAxis = latticeAxis(point[2], PrimeZ)
  result = NearestFeatures(hash: vmovq_n_u32(0),
    first: float4(1.0e10'f), second: float4(1.0e10'f))

  template consider(lattice: uint32x4, delta: array[N, Float4]) =
    let
      hash = hash4(seed, lattice)
      offset = featureOffsets(hash, delta, jitter4)
    when squared:
      let squaredXY = offset[0] * offset[0] + offset[1] * offset[1]
      let euclideanSquared =
        when N == 2: squaredXY
        else: squaredXY + offset[2] * offset[2]
    when manhattan:
      let manhattanXY = abs4(offset[0]) + abs4(offset[1])
      let manhattanDistance =
        when N == 2: manhattanXY
        else: manhattanXY + abs4(offset[2])
    let distance =
      when squared and manhattan: manhattanDistance + euclideanSquared
      elif squared: euclideanSquared
      else: manhattanDistance
    when trackSecond:
      result.second = vmaxq_f32(vminq_f32(result.second, distance), result.first)
    when trackHash:
      # Tied distances retain the first feature encountered, as in the scalar
      # search. Use this mask before replacing its nearest distance.
      result.hash = vbslq_u32(vcltq_f32(distance, result.first), hash, result.hash)
    result.first = vminq_f32(distance, result.first)

  for x in 0 .. 2:
    for y in 0 .. 2:
      let xy = veorq_u32(xAxis.primed[x], yAxis.primed[y])
      when N == 2:
        consider(xy, [xAxis.delta[x], yAxis.delta[y]])
      else:
        for z in 0 .. 2:
          consider(veorq_u32(xy, zAxis.primed[z]),
            [xAxis.delta[x], yAxis.delta[y], zAxis.delta[z]])
