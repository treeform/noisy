import std/macros, avx2base
include avx2flags

macro unroll(index: untyped, limit: static[int], body: untyped): untyped =
  ## Expands the dimension and corner loops with constant array indexes.
  result = newStmtList()
  for i in 0 ..< limit:
    let value = newLit(i)
    result.add quote do:
      block:
        const `index` = `value`
        `body`

func permute8(perm: array[512, int32], indexes: Int8): Int8 {.avxInline.} =
  gather(perm[0].unsafeAddr, indexes, 4)

func cornerHashes[N: static[int]](
  perm: array[512, int32], lattice: array[N, array[2, Int8]]
): array[1 shl N, Int8] {.avxInline.} =
  ## Shares the nested permutation suffixes between corners.
  unroll(level, N):
    const
      axis = N - 1 - level
      count = 1 shl level
    unroll(corner, count):
      let suffix = result[corner]
      result[corner] = perm.permute8(lattice[axis][0] + suffix)
      result[corner + count] = perm.permute8(lattice[axis][1] + suffix)

func signed8(value: Float8, hash: Int8, shift: static[int]): Float8 {.avxInline.} =
  ## Negates each lane whose selected hash bit is set.
  let sign = mm256_and_si256(
    mm256_slli_epi32(hash, 31 - shift), uints8(0x80000000'u32)
  )
  value xor mm256_castsi256_ps(sign)

func gradient8[N: static[int]](
  hash: Int8, offset: array[N, Float8]
): Float8 {.avxInline.} =
  ## Computes eight Gustavson gradients using the scalar dimension's rules.
  when N == 1:
    let magnitude = toFloats(mm256_add_epi32(int8(1), hash and int8(7)))
    signed8(magnitude, hash, 3) * offset[0]
  elif N == 2:
    let
      h = hash and int8(7)
      useX = lessInt(h, uints8(4))
      u = select(useX, offset[0], offset[1])
      v = select(useX, offset[1], offset[0])
    signed8(u, hash, 0) + signed8(float8(2) * v, hash, 1)
  elif N == 3:
    let
      h = hash and int8(15)
      u = select(lessInt(h, uints8(8)), offset[0], offset[1])
      useX = mm256_cmpeq_epi32(mm256_and_si256(h, uints8(13)), uints8(12))
      v = select(lessInt(h, uints8(4)), offset[1],
        select(useX, offset[0], offset[2]))
    signed8(u, hash, 0) + signed8(v, hash, 1)
  else:
    let
      h = hash and int8(31)
      u = select(lessInt(h, uints8(24)), offset[0], offset[1])
      v = select(lessInt(h, uints8(16)), offset[1], offset[2])
      w = select(lessInt(h, uints8(8)), offset[2], offset[3])
    signed8(u, hash, 0) + signed8(v, hash, 1) + signed8(w, hash, 2)

func samplePerlinGrid8*[N: static[int]](
  perm: array[512, int32], point: array[N, Float8]
): Float8 {.avxInline.} =
  ## Samples eight grid points that vary only on the contiguous axis.
  var
    lattice: array[N, array[2, Int8]]
    offset: array[N, array[2, Float8]]
    weights: array[N, Float8]
  unroll(axis, N):
    let
      cell = mm256_floor_ps(point[axis])
      fraction = point[axis] - cell
      # Out-of-int32 float32 cells have zero low bits. CVTT returns int32.low,
      # whose low eight bits are also zero.
      lower = truncate(cell) and int8(255)
    lattice[axis] = [lower, mm256_add_epi32(lower, int8(1)) and int8(255)]
    offset[axis] = [fraction, fraction - float8(1)]
    weights[axis] = fraction * fraction * fraction *
      (fraction * (fraction * float8(6) - float8(15)) + float8(10))

  let hashes = perm.cornerHashes(lattice)

  var corners: array[1 shl N, Float8]
  unroll(corner, 1 shl N):
    var delta: array[N, Float8]
    unroll(axis, N):
      const upper = (corner shr (N - 1 - axis)) and 1
      delta[axis] = offset[axis][upper]
    corners[corner] = gradient8(hashes[corner], delta)
  unroll(level, N):
    const axis = N - 1 - level
    unroll(corner, 1 shl axis):
      let
        lower = corners[corner * 2]
        upper = corners[corner * 2 + 1]
      corners[corner] = lower + weights[axis] * (upper - lower)
  corners[0]
