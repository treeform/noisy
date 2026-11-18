import ../lattices, avx2base
include avx2flags

func fastFloor8*(a: Float8): Int8 {.avxInline.} =
  ## Matches the preceding-cell convention at negative integers.
  truncate(a) + less(a, float8(0))

func hash8*(seed: int32, lattice: Int8): Int8 {.avxInline.} =
  mm256_mullo_epi32(mm256_xor_si256(int8(seed), lattice),
    uints8(HashMultiplier))

func latticeValue8*(seed: int32, lattice: Int8): Float8 {.avxInline.} =
  let
    hash = hash8(seed, lattice)
    square = mm256_mullo_epi32(hash, hash)
    bits = mm256_xor_si256(square, mm256_slli_epi32(square, 19))
  toFloats(bits) * float8(1.0'f / 2147483648.0'f)

func gradientIndex8(seed: int32, lattice: Int8): Int8 {.avxInline.} =
  let hash = hash8(seed, lattice)
  mm256_xor_si256(hash, mm256_srai_epi32(hash, 15))

func gradient8*(seed: int32, x, y: Int8, dx, dy: Float8): Float8 {.avxInline.} =
  let index = gradientIndex8(seed, mm256_xor_si256(x, y)) and int8(254)
  dx * gather(Gradients2D[0].unsafeAddr, index, 4) +
    dy * gather(Gradients2D[1].unsafeAddr, index, 4)

func gradient8*(
  seed: int32, x, y, z: Int8, dx, dy, dz: Float8
): Float8 {.avxInline.} =
  let index = gradientIndex8(seed, mm256_xor_si256(mm256_xor_si256(x, y), z)) and
    int8(252)
  dx * gather(Gradients3D[0].unsafeAddr, index, 4) +
    dy * gather(Gradients3D[1].unsafeAddr, index, 4) +
    dz * gather(Gradients3D[2].unsafeAddr, index, 4)
