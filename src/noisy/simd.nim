when defined(arm64):
  import nimsimd/neon

  type
    Float4* = float32x4
    Int4* = int32x4

  func storeNeon(p: ptr float32, a: Float4) {.
      importc: "vst1q_f32", header: "arm_neon.h".}
    ## Stores lanes through the typed pointer missing in nimsimd 1.3.2.
elif defined(amd64) or defined(i386):
  import nimsimd/sse2

  type
    Float4* = M128
    Int4* = M128i
else:
  {.error: "Noisy SIMD requires ARM64 NEON or x86 SSE2.".}

func float4*(a: float32): Float4 {.inline.} =
  ## Broadcasts a float to all four lanes.
  when defined(arm64):
    vmovq_n_f32(a)
  else:
    mm_set1_ps(a)

func int4*(a: int32): Int4 {.inline.} =
  ## Broadcasts an integer to all four lanes.
  when defined(arm64):
    vmovq_n_s32(a)
  else:
    mm_set1_epi32(a)

func `and`*(a, b: Float4): Float4 {.inline.} =
  ## Computes the bitwise intersection of float lane representations.
  when defined(arm64):
    cast[Float4](vandq_u32(cast[uint32x4](a), cast[uint32x4](b)))
  else:
    mm_and_ps(a, b)

func `and`*(a, b: Int4): Int4 {.inline.} =
  ## Computes the bitwise intersection of integer lanes.
  when defined(arm64):
    cast[Int4](vandq_u32(cast[uint32x4](a), cast[uint32x4](b)))
  else:
    mm_and_si128(a, b)

func `or`*(a, b: Float4): Float4 {.inline.} =
  ## Computes the bitwise union of float lane representations.
  when defined(arm64):
    cast[Float4](vorrq_u32(cast[uint32x4](a), cast[uint32x4](b)))
  else:
    mm_or_ps(a, b)

func `xor`*(a, b: Float4): Float4 {.inline.} =
  ## Computes the bitwise exclusive union of float lane representations.
  when defined(arm64):
    cast[Float4](veorq_u32(cast[uint32x4](a), cast[uint32x4](b)))
  else:
    mm_xor_ps(a, b)

func `<`*(a, b: Float4): Float4 {.inline.} =
  ## Returns all bits set in each lane where a is less than b.
  when defined(arm64):
    cast[Float4](vcltq_f32(a, b))
  else:
    mm_cmplt_ps(a, b)

func `<=`*(a, b: Float4): Float4 {.inline.} =
  ## Returns all bits set in each lane where a is at most b.
  when defined(arm64):
    cast[Float4](vcleq_f32(a, b))
  else:
    mm_cmple_ps(a, b)

func `+`*(a, b: Float4): Float4 {.inline.} =
  ## Adds corresponding float lanes.
  when defined(arm64):
    vaddq_f32(a, b)
  else:
    mm_add_ps(a, b)

func `+=`*(a: var Float4, b: Float4) {.inline.} =
  ## Adds corresponding float lanes in place.
  a = a + b

func `-`*(a, b: Float4): Float4 {.inline.} =
  ## Subtracts corresponding float lanes.
  when defined(arm64):
    vsubq_f32(a, b)
  else:
    mm_sub_ps(a, b)

func `*`*(a, b: Float4): Float4 {.inline.} =
  ## Multiplies corresponding float lanes.
  when defined(arm64):
    vmulq_f32(a, b)
  else:
    mm_mul_ps(a, b)

func `/`*(a, b: Float4): Float4 {.inline.} =
  ## Divides corresponding float lanes using hardware division.
  when defined(arm64):
    vdivq_f32(a, b)
  else:
    mm_div_ps(a, b)

func `/=`*(a: var Float4, b: Float4) {.inline.} =
  ## Divides corresponding float lanes in place.
  a = a / b

func addInts*(a, b: Int4): Int4 {.inline.} =
  ## Adds corresponding integer lanes with wrapping arithmetic.
  when defined(arm64):
    cast[Int4](vaddq_u32(cast[uint32x4](a), cast[uint32x4](b)))
  else:
    mm_add_epi32(a, b)

func truncate*(a: Float4): Int4 {.inline.} =
  ## Converts finite float lanes in the int32 range towards zero.
  when defined(arm64):
    vcvtq_s32_f32(a)
  else:
    mm_cvttps_epi32(a)

func toInts*(a: Float4): Int4 {.inline.} =
  ## Converts integer-valued float lanes in the int32 range to integers.
  truncate(a)

func toFloats*(a: Int4): Float4 {.inline.} =
  ## Converts integer lanes to floats.
  when defined(arm64):
    vcvtq_f32_s32(a)
  else:
    mm_cvtepi32_ps(a)

func store*(p: ptr float32, a: Float4) {.inline.} =
  ## Stores four float lanes without requiring SIMD alignment.
  when defined(arm64):
    storeNeon(p, a)
  else:
    mm_storeu_ps(p, a)
