import nimsimd/avx2
export avx2
include avx2flags

type
  Float8* = M256
  Int8* = M256i

template float8*(a: float32): Float8 = mm256_set1_ps(a)
template int8*(a: int32): Int8 = mm256_set1_epi32(a)
template uints8*(a: uint32): Int8 = int8(cast[int32](a))
template `+`*(a, b: Float8): Float8 = mm256_add_ps(a, b)
template `-`*(a, b: Float8): Float8 = mm256_sub_ps(a, b)
template `*`*(a, b: Float8): Float8 = mm256_mul_ps(a, b)
template `/`*(a, b: Float8): Float8 = mm256_div_ps(a, b)
template `xor`*(a, b: Float8): Float8 = mm256_xor_ps(a, b)
template `and`*(a, b: Int8): Int8 = mm256_and_si256(a, b)
template `+=`*(a: var Float8, b: Float8) = a = a + b
template `+`*(a, b: Int8): Int8 = mm256_add_epi32(a, b)
template toFloats*(a: Int8): Float8 = mm256_cvtepi32_ps(a)
template truncate*(a: Float8): Int8 = mm256_cvttps_epi32(a)
template less*(a, b: Float8): Int8 =
  mm256_castps_si256(mm256_cmp_ps(a, b, 17))
template greaterEqual*(a, b: Float8): Int8 =
  mm256_castps_si256(mm256_cmp_ps(a, b, 29))
template `<`*(a, b: Float8): Int8 = less(a, b)
template lessInt*(a, b: Int8): Int8 = mm256_cmpgt_epi32(b, a)
template bitClear*(a, b: Int8): Int8 = mm256_andnot_si256(b, a)
template bitNot*(a: Int8): Int8 = mm256_xor_si256(a, int8(-1))
template select*(mask: Int8, yes, no: Float8): Float8 =
  mm256_blendv_ps(no, yes, mm256_castsi256_ps(mask))
template select*(mask: Int8, yes, no: Int8): Int8 =
  mm256_blendv_epi8(no, yes, mask)
template abs8*(a: Float8): Float8 =
  mm256_andnot_ps(mm256_set1_ps(-0.0'f), a)
template sqrt8*(a: Float8): Float8 = mm256_sqrt_ps(a)

# Typed pointers are required by the C++ intrinsic declarations.
func gather*(p: ptr int32, indexes: Int8, scale: int32): Int8
  {.importc: "_mm256_i32gather_epi32", header: "immintrin.h".}
func gather*(p: ptr float32, indexes: Int8, scale: int32): Float8
  {.importc: "_mm256_i32gather_ps", header: "immintrin.h".}

template coordinates8*(start: float32, offset: int): Float8 =
  # Preserve start + float32(index) even near the float32 precision limit.
  float8(start) + mm256_setr_ps(
    offset.float32, (offset + 1).float32, (offset + 2).float32,
    (offset + 3).float32, (offset + 4).float32, (offset + 5).float32,
    (offset + 6).float32, (offset + 7).float32)
