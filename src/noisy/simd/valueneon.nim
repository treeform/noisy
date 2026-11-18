import nimsimd/neon, ../[lattices, simd], latticesneon

func adjacent2(
  cell: Int4, stride: uint32
): tuple[first, second: uint32x4] {.inline.} =
  ## Primes the two neighboring lattice points in each lane.
  let
    step = vmovq_n_u32(stride)
    first = multiply4(cast[uint32x4](cell), step)
  (first, vaddq_u32(first, step))

func lerp4(a, b, t: Float4): Float4 {.inline.} =
  ## Interpolates along one axis using the scalar arithmetic order.
  a + t * (b - a)

func hermite4(t: Float4): Float4 {.inline.} =
  ## Applies the same cubic Hermite fade as the scalar sampler.
  t * t * (float4(3) - float4(2) * t)

func plane4(
  seed: int32, x0, x1, y0, y1: uint32x4, xs, ys: Float4
): Float4 {.inline.} =
  ## Interpolates four corner hashes along x, then y.
  let
    row0 = lerp4(
      latticeValue4(seed, veorq_u32(x0, y0)),
      latticeValue4(seed, veorq_u32(x1, y0)), xs
    )
    row1 = lerp4(
      latticeValue4(seed, veorq_u32(x0, y1)),
      latticeValue4(seed, veorq_u32(x1, y1)), xs
    )
  lerp4(row0, row1, ys)

func sampleValue4*(seed: int32, x, y: Float4): Float4 {.inline.} =
  ## Samples four 2D Value noise points in parallel.
  let
    ix = fastFloor4(x)
    iy = fastFloor4(y)
    (x0, x1) = adjacent2(ix, PrimeX)
    (y0, y1) = adjacent2(iy, PrimeY)
    xs = hermite4(x - toFloats(ix))
    ys = hermite4(y - toFloats(iy))
  plane4(seed, x0, x1, y0, y1, xs, ys)

func sampleValue4*(seed: int32, x, y, z: Float4): Float4 {.inline.} =
  ## Samples four 3D Value noise points in parallel.
  let
    ix = fastFloor4(x)
    iy = fastFloor4(y)
    iz = fastFloor4(z)
    (x0, x1) = adjacent2(ix, PrimeX)
    (y0, y1) = adjacent2(iy, PrimeY)
    (z0, z1) = adjacent2(iz, PrimeZ)
    xs = hermite4(x - toFloats(ix))
    ys = hermite4(y - toFloats(iy))
    zs = hermite4(z - toFloats(iz))
    plane0 = plane4(seed, x0, x1, veorq_u32(y0, z0), veorq_u32(y1, z0), xs, ys)
    plane1 = plane4(seed, x0, x1, veorq_u32(y0, z1), veorq_u32(y1, z1), xs, ys)
  lerp4(plane0, plane1, zs)
