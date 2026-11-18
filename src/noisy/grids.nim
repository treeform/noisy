import common
from simplexcommon import Grid

export Grid

const useGridNeon* = defined(arm64) and not defined(noisyNoNeon) and
  not defined(js) and not defined(nimscript)

const useGridAvx2* = defined(amd64) and
  (defined(gcc) or defined(clang) or defined(vcc)) and
  not defined(noisyNoAvx) and not defined(js) and not defined(nimscript)

when useGridAvx2:
  import simd/[gridavx2, runtimecheck]
  export AvxKind
else:
  # Keep the generic filler independent of platform-specific vector types.
  type
    AvxKind = enum avxDisabled
    AvxSettings = object

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
  noise: T, start: float32, width: int,
  batch: static[Batch], avxSettings: AvxSettings = default(AvxSettings)
): Grid =
  ## Fills a line using the supplied batch sampler and scalar tails.
  mixin value
  if noise.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  let length = gridLength([width])
  result = Grid(
    width: width, height: 1, depth: 1, frames: 1,
    values: newSeq[float32](length)
  )
  let scalarStart =
    when Batch is GridSampler: 0
    elif Batch is AvxKind: width - width mod 8
    else: width - width mod 4
  when Batch isnot GridSampler and Batch isnot AvxKind:
    for x in countup(0, scalarStart - 4, 4):
      store(result.values[x].unsafeAddr, batch(noise, gridCoordinates(start, x)))
  when Batch is AvxKind:
    if scalarStart > 0:
      fillRowAvx2(AvxBatch(kind: batch, settings: avxSettings),
        [start, 0'f, 0'f, 0'f], 1, width,
        result.values[0].unsafeAddr)
  for x in scalarStart ..< width:
    result.values[x] = noise.value(start + x.float32)

func fillGridImpl[T, Batch](
  noise: T, start: (float32, float32), dimens: (int, int),
  batch: static[Batch], avxSettings: AvxSettings = default(AvxSettings)
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
    elif Batch is AvxKind: height - height mod 8
    else: height - height mod 4
  var index = 0
  for x in 0 ..< width:
    let sampleX = start[0] + x.float32
    when Batch isnot GridSampler and Batch isnot AvxKind:
      let batchX = float4(sampleX)
      for y in countup(0, scalarStart - 4, 4):
        store(
          result.values[index].unsafeAddr,
          batch(noise, batchX, gridCoordinates(start[1], y))
        )
        index += 4
    when Batch is AvxKind:
      if scalarStart > 0:
        fillRowAvx2(AvxBatch(kind: batch, settings: avxSettings),
          [sampleX, start[1], 0'f, 0'f], 2, height,
          result.values[index].unsafeAddr)
      index += scalarStart
    for y in scalarStart ..< height:
      result.values[index] = noise.value(sampleX, start[1] + y.float32)
      inc index

func fillGridImpl[T, Batch](
  noise: T,
  start: (float32, float32, float32),
  dimens: (int, int, int),
  batch: static[Batch], avxSettings: AvxSettings = default(AvxSettings)
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
    elif Batch is AvxKind: depth - depth mod 8
    else: depth - depth mod 4
  var index = 0
  for x in 0 ..< width:
    let sampleX = start[0] + x.float32
    when Batch isnot GridSampler and Batch isnot AvxKind:
      let batchX = float4(sampleX)
    for y in 0 ..< height:
      let sampleY = start[1] + y.float32
      when Batch isnot GridSampler and Batch isnot AvxKind:
        let batchY = float4(sampleY)
        for z in countup(0, scalarStart - 4, 4):
          store(
            result.values[index].unsafeAddr,
            batch(noise, batchX, batchY, gridCoordinates(start[2], z))
          )
          index += 4
      when Batch is AvxKind:
        if scalarStart > 0:
          fillRowAvx2(AvxBatch(kind: batch, settings: avxSettings),
            [sampleX, sampleY, start[2], 0'f], 3, depth,
            result.values[index].unsafeAddr)
        index += scalarStart
      for z in scalarStart ..< depth:
        result.values[index] = noise.value(sampleX, sampleY, start[2] + z.float32)
        inc index

func fillGridImpl[T, Batch](
  noise: T,
  start: (float32, float32, float32, float32),
  dimens: (int, int, int, int),
  batch: static[Batch], avxSettings: AvxSettings = default(AvxSettings)
): Grid =
  ## Fills each volume slice using the supplied batch sampler and scalar tails.
  mixin value
  if noise.octaves <= 0:
    raise newException(NoisyError, "Octaves must be > 0")
  let
    (width, height, depth, frames) = dimens
    length = gridLength([width, height, depth, frames])
  result = Grid(
    width: width, height: height, depth: depth, frames: frames,
    values: newSeq[float32](length)
  )
  if length == 0:
    return
  let scalarStart =
    when Batch is GridSampler: 0
    elif Batch is AvxKind: depth - depth mod 8
    else: depth - depth mod 4
  var index = 0
  for w in 0 ..< frames:
    let sampleW = start[3] + w.float32
    when Batch isnot GridSampler and Batch isnot AvxKind:
      let batchW = float4(sampleW)
    for x in 0 ..< width:
      let sampleX = start[0] + x.float32
      when Batch isnot GridSampler and Batch isnot AvxKind:
        let batchX = float4(sampleX)
      for y in 0 ..< height:
        let sampleY = start[1] + y.float32
        when Batch isnot GridSampler and Batch isnot AvxKind:
          let batchY = float4(sampleY)
          for z in countup(0, scalarStart - 4, 4):
            store(
              result.values[index].unsafeAddr,
              batch(noise, batchX, batchY, gridCoordinates(start[2], z), batchW)
            )
            index += 4
        when Batch is AvxKind:
          if scalarStart > 0:
            fillRowAvx2(AvxBatch(kind: batch, settings: avxSettings),
              [sampleX, sampleY, start[2], sampleW], 4, depth,
              result.values[index].unsafeAddr)
          index += scalarStart
        for z in scalarStart ..< depth:
          result.values[index] = noise.value(
            sampleX, sampleY, start[2] + z.float32, sampleW
          )
          inc index

func scalarGrid*[T](noise: T, start: float32, width: int): Grid {.inline.} =
  ## Fills a line with scalar samples at unit intervals.
  noise.fillGridImpl(start, width, sampleScalar)

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

func scalarGrid*[T](
  noise: T,
  start: (float32, float32, float32, float32),
  dimens: (int, int, int, int)
): Grid {.inline.} =
  ## Fills volume slices with scalar samples at unit intervals.
  noise.fillGridImpl(start, dimens, sampleScalar)

when useGridNeon:
  type
    GridBatch1[T] = proc(noise: T, x: Float4): Float4 {.noSideEffect.}
    GridBatch2[T] = proc(noise: T, x, y: Float4): Float4 {.noSideEffect.}
    GridBatch3[T] = proc(noise: T, x, y, z: Float4): Float4 {.noSideEffect.}
    GridBatch4[T] = proc(noise: T, x, y, z, w: Float4): Float4 {.noSideEffect.}

  func fillGrid*[T](
    noise: T, start: float32, width: int, batch: static[GridBatch1[T]]
  ): Grid {.inline.} =
    ## Generates a 1D grid using four contiguous samples per batch.
    when nimvm:
      noise.scalarGrid(start, width)
    else:
      noise.fillGridImpl(start, width, batch)

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

  func fillGrid*[T](
    noise: T,
    start: (float32, float32, float32, float32),
    dimens: (int, int, int, int),
    batch: static[GridBatch4[T]]
  ): Grid {.inline.} =
    ## Generates a 4D grid using four contiguous z samples per batch.
    when nimvm:
      noise.scalarGrid(start, dimens)
    else:
      noise.fillGridImpl(start, dimens, batch)

when useGridAvx2:
  func avxGrid*[T, Start, Size](
    noise: T, start: Start, dimens: Size, kind: static[AvxKind]
  ): Grid =
    ## Dispatches through a scalar ABI; unsupported machines use point samples.
    when nimvm:
      result = noise.scalarGrid(start, dimens)
    else:
      if not hasAvx2():
        return noise.scalarGrid(start, dimens)
      var settings = AvxSettings(
        octaves: noise.octaves, amplitude: noise.amplitude,
        frequency: noise.frequency, lacunarity: noise.lacunarity, gain: noise.gain)
      when kind == avxPerlin:
        settings.perm = noise.perm.unsafeAddr
      else:
        settings.seed = cast[int32](noise.seed)
      when kind == avxCellular:
        settings.jitter = noise.jitter
        settings.distanceFunction = ord(noise.distanceFunction)
        settings.returnType = ord(noise.returnType)
      result = noise.fillGridImpl(start, dimens, kind, settings)
