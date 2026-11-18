import std/macros, nimsimd/neon, ../simd

func floor4(a: Float4): Float4 {.
    importc: "vrndmq_f32", header: "arm_neon.h".}
  ## Floors four float32 lanes without narrowing the integer range.

macro unroll(index: untyped, limit: static[int], body: untyped): untyped =
  ## Expands the dimension and corner loops with constant array indexes.
  result = newStmtList()
  for i in 0 ..< limit:
    let value = newLit(i)
    result.add quote do:
      block:
        const `index` = `value`
        `body`

func permute4(perm: array[512, int32], indexes: Int4): Int4 {.inline.} =
  ## Gathers from the original 512-entry permutation in lane order.
  cast[Int4]([
    perm[vgetq_lane_s32(indexes, 0) and 511],
    perm[vgetq_lane_s32(indexes, 1) and 511],
    perm[vgetq_lane_s32(indexes, 2) and 511],
    perm[vgetq_lane_s32(indexes, 3) and 511]
  ])

func cornerHashes[N: static[int]; Index: int | Int4](
  perm: array[512, int32], lattice: array[N, array[2, Index]]
): array[1 shl N, Index] {.inline.} =
  ## Shares the nested permutation suffixes between corners.
  unroll(level, N):
    const
      axis = N - 1 - level
      count = 1 shl level
    unroll(corner, count):
      let suffix = result[corner]
      when Index is int:
        result[corner] = perm[(lattice[axis][0] +% suffix) and 511]
        result[corner + count] = perm[(lattice[axis][1] +% suffix) and 511]
      else:
        result[corner] = perm.permute4(addInts(lattice[axis][0], suffix))
        result[corner + count] = perm.permute4(addInts(lattice[axis][1], suffix))

func signed4(value: Float4, hash: Int4, shift: static[int]): Float4 {.inline.} =
  ## Negates each lane whose selected hash bit is set.
  let sign = vandq_u32(
    vshlq_n_u32(cast[uint32x4](hash), 31 - shift), vmovq_n_u32(0x80000000'u32)
  )
  value xor cast[Float4](sign)

func gradient4[N: static[int]](
  hash: Int4, offset: array[N, Float4]
): Float4 {.inline.} =
  ## Computes four Gustavson gradients using the scalar dimension's rules.
  when N == 1:
    let magnitude = toFloats(addInts(int4(1), hash and int4(7)))
    signed4(magnitude, hash, 3) * offset[0]
  elif N == 2:
    let
      h = cast[uint32x4](hash and int4(7))
      useX = vcltq_u32(h, vmovq_n_u32(4))
      u = vbslq_f32(useX, offset[0], offset[1])
      v = vbslq_f32(useX, offset[1], offset[0])
    signed4(u, hash, 0) + signed4(float4(2) * v, hash, 1)
  elif N == 3:
    let
      h = cast[uint32x4](hash and int4(15))
      u = vbslq_f32(vcltq_u32(h, vmovq_n_u32(8)), offset[0], offset[1])
      useX = vceqq_u32(vandq_u32(h, vmovq_n_u32(13)), vmovq_n_u32(12))
      v = vbslq_f32(vcltq_u32(h, vmovq_n_u32(4)), offset[1],
        vbslq_f32(useX, offset[0], offset[2]))
    signed4(u, hash, 0) + signed4(v, hash, 1)
  else:
    let
      h = cast[uint32x4](hash and int4(31))
      u = vbslq_f32(vcltq_u32(h, vmovq_n_u32(24)), offset[0], offset[1])
      v = vbslq_f32(vcltq_u32(h, vmovq_n_u32(16)), offset[1], offset[2])
      w = vbslq_f32(vcltq_u32(h, vmovq_n_u32(8)), offset[2], offset[3])
    signed4(u, hash, 0) + signed4(v, hash, 1) + signed4(w, hash, 2)

func samplePerlinGrid4*[N: static[int]](
  perm: array[512, int32], point: array[N, Float4]
): Float4 =
  ## Samples four grid points that vary only on the contiguous axis.
  static: doAssert N in 1 .. 4
  var
    lattice: array[N, array[2, Int4]]
    offset: array[N, array[2, Float4]]
    weights: array[N, Float4]
  unroll(axis, N):
    let
      cell = floor4(point[axis])
      fraction = point[axis] - cell
      # The low eight bits are zero outside the int32 range due to float32
      # spacing. Mask positive overflow because NEON saturates to int32.high.
      fitsUpper = cast[Int4](cell < float4(2147483648.0'f))
      lower = toInts(cell) and int4(255) and fitsUpper
    lattice[axis] = [lower, addInts(lower, int4(1)) and int4(255)]
    offset[axis] = [fraction, fraction - float4(1)]
    weights[axis] = fraction * fraction * fraction *
      (fraction * (fraction * float4(6) - float4(15)) + float4(10))

  var hashes: array[1 shl N, Int4]
  when N == 1:
    hashes = perm.cornerHashes(lattice)
  else:
    # Grid batches only vary along the contiguous axis. Within one lattice
    # cell, all four lanes can share the same permutation lookups.
    const gridAxis = min(N - 1, 2)
    let
      lanes = lattice[gridAxis][0]
      sameCell = vceqq_u32(cast[uint32x4](lanes),
        cast[uint32x4](int4(vgetq_lane_s32(lanes, 0))))
    if vminvq_u32(sameCell) != 0:
      var cell: array[N, array[2, int]]
      unroll(axis, N):
        cell[axis] = [
          vgetq_lane_s32(lattice[axis][0], 0).int,
          vgetq_lane_s32(lattice[axis][1], 0).int
        ]
      let shared = perm.cornerHashes(cell)
      unroll(corner, 1 shl N):
        hashes[corner] = int4(shared[corner].int32)
    else:
      hashes = perm.cornerHashes(lattice)

  var corners: array[1 shl N, Float4]
  unroll(corner, 1 shl N):
    var delta: array[N, Float4]
    unroll(axis, N):
      const upper = (corner shr (N - 1 - axis)) and 1
      delta[axis] = offset[axis][upper]
    corners[corner] = gradient4(hashes[corner], delta)
  unroll(level, N):
    const axis = N - 1 - level
    unroll(corner, 1 shl axis):
      let
        lower = corners[corner * 2]
        upper = corners[corner * 2 + 1]
      corners[corner] = lower + weights[axis] * (upper - lower)
  corners[0]
