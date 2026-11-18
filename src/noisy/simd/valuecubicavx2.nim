import ../lattices, avx2base, latticesavx2
include avx2flags

func adjacent4(cell: Int8, stride: uint32): array[4, Int8] {.avxInline.} =
  ## Computes the primed coordinates of the four interpolation points.
  let
    step = uints8(stride)
    center = mm256_mullo_epi32(cell, step)
    next = mm256_add_epi32(center, step)
  [mm256_sub_epi32(center, step), center, next, mm256_add_epi32(next, step)]

func cubic8(a, b, c, d, t, t2, t3: Float8): Float8 {.avxInline.} =
  ## Uses the scalar interpolant's arithmetic order across eight grid cells.
  let
    difference = a - b
    p = (d - c) - difference
  t3 * p + t2 * (difference - p) + t * (c - a) + b

func plane8(
  seed: int32,
  columns, rows: array[4, Int8],
  layer: Int8,
  x, y: Float8
): Float8 {.avxInline.} =
  ## Interpolates eight independent sets of lattice values along x then y.
  let
    x2 = x * x
    x3 = x2 * x
    y2 = y * y
  var interpolated: array[4, Float8]
  for row in 0 ..< 4:
    let lattice = mm256_xor_si256(rows[row], layer)
    interpolated[row] = cubic8(
      latticeValue8(seed, mm256_xor_si256(columns[0], lattice)),
      latticeValue8(seed, mm256_xor_si256(columns[1], lattice)),
      latticeValue8(seed, mm256_xor_si256(columns[2], lattice)),
      latticeValue8(seed, mm256_xor_si256(columns[3], lattice)),
      x, x2, x3
    )
  cubic8(
    interpolated[0], interpolated[1], interpolated[2], interpolated[3],
    y, y2, y2 * y
  )

func sampleCubic8*(seed: int32, x, y: Float8): Float8 {.avxInline.} =
  ## Evaluates eight 2D cubic Value samples with AVX2.
  let
    ix = fastFloor8(x)
    iy = fastFloor8(y)
  plane8(
    seed, adjacent4(ix, PrimeX), adjacent4(iy, PrimeY), uints8(0),
    x - toFloats(ix), y - toFloats(iy)
  ) * float8(1.0'f / (1.5'f * 1.5'f))

func sampleCubic8*(seed: int32, x, y, z: Float8): Float8 {.avxInline.} =
  ## Evaluates eight 3D cubic Value samples with AVX2.
  let
    ix = fastFloor8(x)
    iy = fastFloor8(y)
    iz = fastFloor8(z)
    columns = adjacent4(ix, PrimeX)
    rows = adjacent4(iy, PrimeY)
    layers = adjacent4(iz, PrimeZ)
    xs = x - toFloats(ix)
    ys = y - toFloats(iy)
    zs = z - toFloats(iz)
    z2 = zs * zs
  var interpolated: array[4, Float8]
  for layer in 0 ..< 4:
    interpolated[layer] = plane8(seed, columns, rows, layers[layer], xs, ys)
  cubic8(
    interpolated[0], interpolated[1], interpolated[2], interpolated[3],
    zs, z2, z2 * zs
  ) * float8(1.0'f / (1.5'f * 1.5'f * 1.5'f))
