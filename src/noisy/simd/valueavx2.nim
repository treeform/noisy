import ../lattices, avx2base, latticesavx2
include avx2flags

func adjacent2(
  cell: Int8, stride: uint32
): tuple[first, second: Int8] {.avxInline.} =
  ## Primes the two neighboring lattice points in each lane.
  let
    step = uints8(stride)
    first = mm256_mullo_epi32(cell, step)
  (first, mm256_add_epi32(first, step))

func lerp8(a, b, t: Float8): Float8 {.avxInline.} =
  ## Interpolates along one axis using the scalar arithmetic order.
  a + t * (b - a)

func hermite8(t: Float8): Float8 {.avxInline.} =
  ## Applies the same cubic Hermite fade as the scalar sampler.
  t * t * (float8(3) - float8(2) * t)

func plane8(
  seed: int32, x0, x1, y0, y1: Int8, xs, ys: Float8
): Float8 {.avxInline.} =
  ## Interpolates four corner hashes along x, then y.
  let
    row0 = lerp8(
      latticeValue8(seed, mm256_xor_si256(x0, y0)),
      latticeValue8(seed, mm256_xor_si256(x1, y0)), xs
    )
    row1 = lerp8(
      latticeValue8(seed, mm256_xor_si256(x0, y1)),
      latticeValue8(seed, mm256_xor_si256(x1, y1)), xs
    )
  lerp8(row0, row1, ys)

func sampleValue8*(seed: int32, x, y: Float8): Float8 {.avxInline.} =
  ## Samples eight 2D Value noise points in parallel.
  let
    ix = fastFloor8(x)
    iy = fastFloor8(y)
    (x0, x1) = adjacent2(ix, PrimeX)
    (y0, y1) = adjacent2(iy, PrimeY)
    xs = hermite8(x - toFloats(ix))
    ys = hermite8(y - toFloats(iy))
  plane8(seed, x0, x1, y0, y1, xs, ys)

func sampleValue8*(seed: int32, x, y, z: Float8): Float8 {.avxInline.} =
  ## Samples eight 3D Value noise points in parallel.
  let
    ix = fastFloor8(x)
    iy = fastFloor8(y)
    iz = fastFloor8(z)
    (x0, x1) = adjacent2(ix, PrimeX)
    (y0, y1) = adjacent2(iy, PrimeY)
    (z0, z1) = adjacent2(iz, PrimeZ)
    xs = hermite8(x - toFloats(ix))
    ys = hermite8(y - toFloats(iy))
    zs = hermite8(z - toFloats(iz))
    plane0 = plane8(seed, x0, x1, mm256_xor_si256(y0, z0),
      mm256_xor_si256(y1, z0), xs, ys)
    plane1 = plane8(seed, x0, x1, mm256_xor_si256(y0, z1),
      mm256_xor_si256(y1, z1), xs, ys)
  lerp8(plane0, plane1, zs)
