import nimsimd/neon, ../[lattices, simd]

const Gradients3DCodes = block:
  # Each component takes two bits for its magnitude and sign.
  var codes: array[64, uint8]
  for gradientIndex in 0 ..< codes.len:
    for axis in 0 ..< 3:
      let component = Gradients3D[gradientIndex * 4 + axis]
      if component != 0:
        codes[gradientIndex] = codes[gradientIndex] or (1'u8 shl (axis * 2))
      if component < 0:
        codes[gradientIndex] = codes[gradientIndex] or (2'u8 shl (axis * 2))
  codes

type GradientPair {.importc: "float32x4x2_t", header: "arm_neon.h".} = object
  val: array[2, Float4]

func multiply4*(a, b: uint32x4): uint32x4 {.
    importc: "vmulq_u32", header: "arm_neon.h".}
  ## Multiplies four unsigned lanes with wrapping.

func fastFloor4*(a: Float4): Int4 {.inline.} =
  ## Uses the preceding cell for negative integers, matching fastFloor.
  addInts(truncate(a), cast[Int4](a < float4(0)))

func hash4*(seed: int32, lattice: uint32x4): uint32x4 {.inline.} =
  ## Hashes a seed and the XOR of the primed coordinates in four lanes.
  multiply4(veorq_u32(vmovq_n_u32(cast[uint32](seed)), lattice),
    vmovq_n_u32(HashMultiplier))

func latticeValue4*(seed: int32, lattice: uint32x4): Float4 {.inline.} =
  ## Converts four wrapped hashes into signed lattice values.
  let
    hash = hash4(seed, lattice)
    square = multiply4(hash, hash)
    bits = veorq_u32(square, vshlq_n_u32(square, 19))
  toFloats(cast[Int4](bits)) * float4(1.0'f / 2147483648.0'f)

func loadPairLane(p: ptr float32, values: GradientPair, lane: int32): GradientPair {.
    importc: "vld2q_lane_f32", header: "arm_neon.h".}
  ## Loads an interleaved gradient pair into one NEON lane.

func loadBytes(p: ptr uint8): uint8x16 {.
    importc: "vld1q_u8", header: "arm_neon.h".}
  ## Loads table bytes through a typed pointer on C and C++.

func lookupBytes(table: uint8x16x4, indexes: uint8x16): uint8x16 {.
    importc: "vqtbl4q_u8", header: "arm_neon.h".}
  ## Looks up bytes across four NEON registers.

func testBits(a, b: uint32x4): uint32x4 {.
    importc: "vtstq_u32", header: "arm_neon.h".}
  ## Returns a full lane mask wherever two values share a set bit.

func gradientIndex4(seed: int32, lattice: uint32x4): uint32x4 {.inline.} =
  ## Finishes the gradient hash with the scalar sampler's signed shift.
  let hashed = hash4(seed, lattice)
  veorq_u32(hashed, cast[uint32x4](vshrq_n_s32(cast[Int4](hashed), 15)))

func gradient4*(
  seed: int32, x, y: uint32x4, dx, dy: Float4
): Float4 {.inline.} =
  ## Gathers four 2D gradient pairs and takes their dot products.
  let indexes = vandq_u32(gradientIndex4(seed, veorq_u32(x, y)), vmovq_n_u32(254))
  var gradients = GradientPair(val: [float4(0), float4(0)])
  gradients = loadPairLane(Gradients2D[vgetq_lane_u32(indexes, 0)].unsafeAddr,
    gradients, 0)
  gradients = loadPairLane(Gradients2D[vgetq_lane_u32(indexes, 1)].unsafeAddr,
    gradients, 1)
  gradients = loadPairLane(Gradients2D[vgetq_lane_u32(indexes, 2)].unsafeAddr,
    gradients, 2)
  gradients = loadPairLane(Gradients2D[vgetq_lane_u32(indexes, 3)].unsafeAddr,
    gradients, 3)
  dx * gradients.val[0] + dy * gradients.val[1]

func signedGradient(
  delta: Float4, encoded: uint32x4, axis: static[int]
): Float4 {.inline.} =
  ## Decodes one coefficient by keeping, negating, or zeroing the delta.
  let
    nonzero = testBits(encoded, vmovq_n_u32(1'u32 shl (axis * 2)))
    sign = vandq_u32(vshlq_n_u32(encoded, 30 - axis * 2), vmovq_n_u32(1'u32 shl 31))
  cast[Float4](vandq_u32(veorq_u32(cast[uint32x4](delta), sign), nonzero))

func gradient4*(
  seed: int32, x, y, z: uint32x4, dx, dy, dz: Float4
): Float4 {.inline.} =
  ## Looks up four packed 3D gradients and takes their dot products.
  let
    hash = gradientIndex4(seed, veorq_u32(veorq_u32(x, y), z))
    indexes = vandq_u32(vshrq_n_u32(hash, 2), vmovq_n_u32(63))
    table = uint8x16x4(val: [
      loadBytes(Gradients3DCodes[0].unsafeAddr),
      loadBytes(Gradients3DCodes[16].unsafeAddr),
      loadBytes(Gradients3DCodes[32].unsafeAddr),
      loadBytes(Gradients3DCodes[48].unsafeAddr)
    ])
    encoded = cast[uint32x4](lookupBytes(table, cast[uint8x16](indexes)))
  signedGradient(dx, encoded, 0) + signedGradient(dy, encoded, 1) +
    signedGradient(dz, encoded, 2)
