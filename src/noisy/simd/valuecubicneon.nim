import nimsimd/neon, ../[lattices, simd], latticesneon

func adjacent4(cell: Int4, stride: uint32): array[4, uint32x4] {.inline.} =
  ## Computes the primed coordinates of the four interpolation points.
  let
    step = vmovq_n_u32(stride)
    center = multiply4(cast[uint32x4](cell), step)
    next = vaddq_u32(center, step)
  [vsubq_u32(center, step), center, next, vaddq_u32(next, step)]

func cubic4(a, b, c, d, t, t2, t3: Float4): Float4 {.inline.} =
  ## Uses the scalar interpolant's arithmetic order across four grid cells.
  let
    difference = a - b
    p = (d - c) - difference
  t3 * p + t2 * (difference - p) + t * (c - a) + b

func plane4(
  seed: int32,
  columns, rows: array[4, uint32x4],
  layer: uint32x4,
  x, y: Float4
): Float4 {.inline.} =
  ## Interpolates four independent sets of lattice values along x then y.
  let
    x2 = x * x
    x3 = x2 * x
    y2 = y * y
  var interpolated: array[4, Float4]
  for row in 0 ..< 4:
    let lattice = veorq_u32(rows[row], layer)
    interpolated[row] = cubic4(
      latticeValue4(seed, veorq_u32(columns[0], lattice)),
      latticeValue4(seed, veorq_u32(columns[1], lattice)),
      latticeValue4(seed, veorq_u32(columns[2], lattice)),
      latticeValue4(seed, veorq_u32(columns[3], lattice)),
      x, x2, x3
    )
  cubic4(
    interpolated[0], interpolated[1], interpolated[2], interpolated[3],
    y, y2, y2 * y
  )

func sampleCubic4*(seed: int32, x, y: Float4): Float4 =
  ## Evaluates four 2D cubic Value samples with NEON.
  let
    ix = fastFloor4(x)
    iy = fastFloor4(y)
  plane4(
    seed, adjacent4(ix, PrimeX), adjacent4(iy, PrimeY), vmovq_n_u32(0),
    x - toFloats(ix), y - toFloats(iy)
  ) * float4(1.0'f / (1.5'f * 1.5'f))

func sampleCubic4*(seed: int32, x, y, z: Float4): Float4 =
  ## Evaluates four 3D cubic Value samples with NEON.
  let
    ix = fastFloor4(x)
    iy = fastFloor4(y)
    iz = fastFloor4(z)
    columns = adjacent4(ix, PrimeX)
    rows = adjacent4(iy, PrimeY)
    layers = adjacent4(iz, PrimeZ)
    xs = x - toFloats(ix)
    ys = y - toFloats(iy)
    zs = z - toFloats(iz)
    z2 = zs * zs
  var interpolated: array[4, Float4]
  for layer in 0 ..< 4:
    interpolated[layer] = plane4(seed, columns, rows, layers[layer], xs, ys)
  cubic4(
    interpolated[0], interpolated[1], interpolated[2], interpolated[3],
    zs, z2, z2 * zs
  ) * float4(1.0'f / (1.5'f * 1.5'f * 1.5'f))
