import common
from simplexcommon import Grid

export Grid

const useGridNeon* = defined(arm64) and not defined(noisyNoNeon) and
  not defined(js) and not defined(nimscript)

when useGridNeon:
  import simd

  func gridCoordinates(start: float32, offset: int): Float4 {.inline.} =
    ## Adds each integer offset independently, as scalar grid sampling does.
    float4(start) + cast[Float4]([
      offset.float32,
      (offset + 1).float32,
      (offset + 2).float32,
      (offset + 3).float32
    ])

type GridSampler = enum
  sampleScalar

when defined(release):
  {.push checks: off.}

func valueIndex*(g: Grid, x, y, z: int): int {.inline.} =
  ## Computes the flat index with z as the contiguous axis.
  z + y * g.depth + x * g.height * g.depth

func valueIndex*(g: Grid, x, y, z, w: int): int {.inline.} =
  ## Stores fourth-axis slices consecutively, retaining the 3D layout.
  w * g.width * g.height * g.depth + g.valueIndex(x, y, z)

func `[]`*(g: Grid, x, y: int, z = 0): float32 {.inline.} =
  ## Returns the noise value at (x, y) or (x, y, z).
  if x < 0 or x >= g.width:
    raise newException(IndexDefect, "Index x out of bounds")
  if y < 0 or y >= g.height:
    raise newException(IndexDefect, "Index y out of bounds")
  if z < 0 or z >= g.depth:
    raise newException(IndexDefect, "Index z out of bounds")
  g.values[g.valueIndex(x, y, z)]

func `[]`*(g: Grid, x: int): float32 {.inline.} =
  ## Returns the noise value at x in a 1D grid.
  g[x, 0, 0]

func `[]`*(g: Grid, x, y, z, w: int): float32 {.inline.} =
  ## Returns a 4D sample, with w selecting a whole 3D slice.
  if x < 0 or x >= g.width:
    raise newException(IndexDefect, "Index x out of bounds")
  if y < 0 or y >= g.height:
    raise newException(IndexDefect, "Index y out of bounds")
  if z < 0 or z >= g.depth:
    raise newException(IndexDefect, "Index z out of bounds")
  if w < 0 or w >= g.frames:
    raise newException(IndexDefect, "Index w out of bounds")
  g.values[g.valueIndex(x, y, z, w)]

when defined(release):
  {.pop.}

func gridLength*(dimens: openArray[int]): int =
  ## Validates dimensions and computes a representable storage size.
  # Leave room for sequence headers and allocator rounding near int.high.
  const MaxSamples = (int.high - 64 * 1024) div sizeof(float32)
  var empty = false
  for size in dimens:
    if size < 0:
      raise newException(NoisyError, "Grid dimensions must be >= 0")
    if size == 0:
      empty = true
  if empty:
    return 0
  result = 1
  for size in dimens:
    if result > MaxSamples div size:
      raise newException(NoisyError, "Grid dimensions are too large")
    result *= size

func fillGridImpl[T, Batch](
  noise: T, start: (float32, float32), dimens: (int, int), batch: static[Batch]
): Grid =
  ## Fills a plane using the supplied batch sampler and scalar tails.
  mixin value
  if noise.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  let
    (width, height) = dimens
    length = gridLength([width, height])
  result = Grid(
    width: width, height: height, depth: 1, frames: 1,
    values: newSeq[float32](length)
  )
  if length == 0:
    return
  let scalarStart =
    when Batch is GridSampler: 0
    else: height - height mod 4
  var index = 0
  for x in 0 ..< width:
    let sampleX = start[0] + x.float32
    when Batch isnot GridSampler:
      let batchX = float4(sampleX)
      for y in countup(0, scalarStart - 4, 4):
        store(
          result.values[index].unsafeAddr,
          batch(noise, batchX, gridCoordinates(start[1], y))
        )
        index += 4
    for y in scalarStart ..< height:
      result.values[index] = noise.value(sampleX, start[1] + y.float32)
      inc index

func fillGridImpl[T, Batch](
  noise: T,
  start: (float32, float32, float32),
  dimens: (int, int, int),
  batch: static[Batch]
): Grid =
  ## Fills a volume using the supplied batch sampler and scalar tails.
  mixin value
  if noise.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  let
    (width, height, depth) = dimens
    length = gridLength([width, height, depth])
  result = Grid(
    width: width, height: height, depth: depth, frames: 1,
    values: newSeq[float32](length)
  )
  if length == 0:
    return
  let scalarStart =
    when Batch is GridSampler: 0
    else: depth - depth mod 4
  var index = 0
  for x in 0 ..< width:
    let sampleX = start[0] + x.float32
    when Batch isnot GridSampler:
      let batchX = float4(sampleX)
    for y in 0 ..< height:
      let sampleY = start[1] + y.float32
      when Batch isnot GridSampler:
        let batchY = float4(sampleY)
        for z in countup(0, scalarStart - 4, 4):
          store(
            result.values[index].unsafeAddr,
            batch(noise, batchX, batchY, gridCoordinates(start[2], z))
          )
          index += 4
      for z in scalarStart ..< depth:
        result.values[index] = noise.value(sampleX, sampleY, start[2] + z.float32)
        inc index

func scalarGrid*[T](
  noise: T, start: (float32, float32), dimens: (int, int)
): Grid {.inline.} =
  ## Fills a plane with scalar samples at unit intervals.
  noise.fillGridImpl(start, dimens, sampleScalar)

func scalarGrid*[T](
  noise: T, start: (float32, float32, float32), dimens: (int, int, int)
): Grid {.inline.} =
  ## Fills a volume with scalar samples at unit intervals.
  noise.fillGridImpl(start, dimens, sampleScalar)

when useGridNeon:
  type
    GridBatch2[T] = proc(noise: T, x, y: Float4): Float4 {.noSideEffect.}
    GridBatch3[T] = proc(noise: T, x, y, z: Float4): Float4 {.noSideEffect.}

  func fillGrid*[T](
    noise: T,
    start: (float32, float32),
    dimens: (int, int),
    batch: static[GridBatch2[T]]
  ): Grid {.inline.} =
    ## Generates a 2D grid using four contiguous y samples per batch.
    when nimvm:
      noise.scalarGrid(start, dimens)
    else:
      noise.fillGridImpl(start, dimens, batch)

  func fillGrid*[T](
    noise: T,
    start: (float32, float32, float32),
    dimens: (int, int, int),
    batch: static[GridBatch3[T]]
  ): Grid {.inline.} =
    ## Generates a 3D grid using four contiguous z samples per batch.
    when nimvm:
      noise.scalarGrid(start, dimens)
    else:
      noise.fillGridImpl(start, dimens, batch)
