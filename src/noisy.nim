import
  std/[math, random],
  noisy/[perlins, opensimplexes, values, cellulars, grids, simd, simplexcommon]

export Simplex, Grid, NoisyError, `[]`
export perlins except NoisyError
export opensimplexes except NoisyError
export values except NoisyError
export cellulars except NoisyError

const allowAvx2 = defined(amd64) and
  (defined(gcc) or defined(clang) or defined(vcc)) and not defined(noisyNoAvx)

when allowAvx2:
  import noisy/simd/[avx2, runtimecheck]

when defined(release):
  {.push checks: off.}

template failOctaves() =
  raise newException(NoisyError, "Octaves must be > 0")

func floor(a: Float4): Float4 {.inline.} =
  const one = [1.float32, 1, 1, 1]
  let tmp = toFloats(truncate(a))
  tmp - ((a < tmp) and cast[Float4](one))

func blend(a, b, mask: Float4): Float4 {.inline.} =
  ((a xor b) and mask) xor a

func initSimplex*(seed: int): Simplex =
  ## Initializes simplex noise with a seeded, repeated permutation table.
  result.octaves = 1
  result.amplitude = 1
  result.frequency = 1
  result.lacunarity = 2
  result.gain = 0.5

  var
    permutation: array[256, int32]
    rng = initRand(seed)
  for i in 0 ..< permutation.len:
    permutation[i] = i.int32
  shuffle(rng, permutation)
  # Duplicate the permutation so masked lattice coordinates wrap consistently.
  for i in 0 ..< result.perm.len:
    result.perm[i] = permutation[i and 255]
    result.permMod12[i] = result.perm[i] mod 12

func `[]=`(g: var Grid, x, y, z: int, value: float32) {.inline.} =
  g.values[g.valueIndex(x, y, z)] = value

func dot(g: array[3, float32], x, y: float32): float32 {.inline.} =
  g[0] * x + g[1] * y

func dot(g: array[3, float32], x, y, z: float32): float32 {.inline.} =
  g[0] * x + g[1] * y + g[2] * z

func noise(simplex: Simplex, x, y: float32): float32 =
  let
    s = (x + y) * F2
    i = floor(x + s)
    j = floor(y + s)
    t = (i + j) * G2
    x0 = x - (i - t)
    y0 = y - (j - t)
    gt = (x0 > y0).int32
    i1 = gt
    j1 = (not gt) and 1
    x1 = x0 - i1.float32 + G2
    y1 = y0 - j1.float32 + G2
    x2 = x0 - 1.float32 + 2.float32 * G2
    y2 = y0 - 1.float32 + 2.float32 * G2
    ii = i.int32 and 255
    jj = j.int32 and 255
    t0 = 0.5.float32 - x0 * x0 - y0 * y0
    t1 = 0.5.float32 - x1 * x1 - y1 * y1
    t2 = 0.5.float32 - x2 * x2 - y2 * y2

  var n0, n1, n2: float32
  if t0 > 0:
    n0 = t0 * t0 * t0 * t0 * dot(
      grad3[simplex.permMod12[ii + simplex.perm[jj]]], x0, y0
    )

  if t1 > 0:
    n1 = t1 * t1 * t1 * t1 * dot(
      grad3[simplex.permMod12[ii + i1 + simplex.perm[jj + j1]]], x1, y1
    )

  if t2 > 0:
    n2 = t2 * t2 * t2 * t2 * dot(
      grad3[simplex.permMod12[ii + 1 + simplex.perm[jj + 1]]], x2, y2
    )

  70.float32 * (n0 + n1 + n2)

func noise(simplex: Simplex, x, y, z: float32): float32 =
  let
    s = (x + y + z) * F3
    i = floor(x + s)
    j = floor(y + s)
    k = floor(z + s)
    t = (i + j + k) * G3
    x0 = x - (i - t)
    y0 = y - (j - t)
    z0 = z - (k - t)
    i1 = (x0 >= y0 and x0 >= z0).int32
    j1 = (x0 < y0 and y0 >= z0).int32
    k1 = (x0 < z0 and y0 < z0).int32
    i2 = (x0 >= y0 or x0 >= z0).int32
    j2 = (x0 < y0 or y0 >= z0).int32
    k2 = (x0 < z0 or y0 < z0).int32
    x1 = x0 - i1.float32 + G3
    y1 = y0 - j1.float32 + G3
    z1 = z0 - k1.float32 + G3
    x2 = x0 - i2.float32 + 2.float32 * G3
    y2 = y0 - j2.float32 + 2.float32 * G3
    z2 = z0 - k2.float32 + 2.float32 * G3
    x3 = x0 - 1.float32 + 3.float32 * G3
    y3 = y0 - 1.float32 + 3.float32 * G3
    z3 = z0 - 1.float32 + 3.float32 * G3
    ii = i.int32 and 255
    jj = j.int32 and 255
    kk = k.int32 and 255
    t0 = 0.6.float32 - x0 * x0 - y0 * y0 - z0 * z0
    t1 = 0.6.float32 - x1 * x1 - y1 * y1 - z1 * z1
    t2 = 0.6.float32 - x2 * x2 - y2 * y2 - z2 * z2
    t3 = 0.6.float32 - x3 * x3 - y3 * y3 - z3 * z3

  var n0, n1, n2, n3: float32
  if t0 >= 0:
    let gi0 = simplex.permMod12[
      ii + simplex.perm[jj + simplex.perm[kk]]
    ]
    n0 = t0 * t0 * t0 * t0 * dot(grad3[gi0], x0, y0, z0)

  if t1 >= 0:
    let gi1 = simplex.permMod12[
      ii + i1 + simplex.perm[jj + j1 + simplex.perm[kk + k1]]
    ]
    n1 = t1 * t1 * t1 * t1 * dot(grad3[gi1], x1, y1, z1)

  if t2 >= 0:
    let gi2 = simplex.permMod12[
      ii + i2 + simplex.perm[jj + j2 + simplex.perm[kk + k2]]
    ]
    n2 = t2 * t2 * t2 * t2 * dot(grad3[gi2], x2, y2, z2)

  if t3 >= 0:
    let gi3 = simplex.permMod12[
      ii + 1 + simplex.perm[jj + 1 + simplex.perm[kk + 1]]
    ]
    n3 = t3 * t3 * t3 * t3 * dot(grad3[gi3], x3, y3, z3)

  32.float32 * (n0 + n1 + n2 + n3)

func value*(simplex: Simplex, x, y: float32): float32 =
  ## Generates the 2D noise value at (x, y) based on the Simplex parameters.

  if simplex.octaves == 0:
    failOctaves()

  var
    total: float32
    amplitude = simplex.amplitude
    frequency = simplex.frequency

  for _ in 0 ..< simplex.octaves:
    total += simplex.noise(x * frequency, y * frequency) * amplitude
    amplitude *= simplex.gain
    frequency *= simplex.lacunarity

  total / simplex.octaves.float32

func value*(simplex: Simplex, x, y, z: float32): float32 =
  ## Generates the 3D noise value at (x, y, z) based on the Simplex parameters.

  if simplex.octaves == 0:
    raise newException(NoisyError, "Octaves must be > 0")

  var
    total: float32
    amplitude = simplex.amplitude
    frequency = simplex.frequency

  for _ in 0 ..< simplex.octaves:
    total += simplex.noise(
      x * frequency, y * frequency, z * frequency
    ) * amplitude
    amplitude *= simplex.gain
    frequency *= simplex.lacunarity

  total / simplex.octaves.float32

func value*(simplex: Simplex, x, y: int): float32 {.inline.} =
  ## Helper for working with ints.
  simplex.value(x.float32, y.float32)

func value*(simplex: Simplex, x, y, z: int): float32 {.inline.} =
  ## Helper for working with ints
  simplex.value(x.float32, y.float32, z.float32)

func coordinates(start: float32, offset: int): Float4 =
  ## Forms world coordinates in the same addition order as scalar sampling.
  var values: array[4, float32]
  for i in 0 ..< 4:
    values[i] = start + (offset + i).float32
  cast[Float4](values)

func column4(simplex: Simplex, x: float32, y: Float4): Float4 =
  ## Evaluates four 2D samples at already scaled coordinates.
  let
    F2 = float4(F2)
    G2 = float4(G2)
    x = float4(x)
    vec0 = float4(0)
    vec1 = float4(1)
    vec2 = float4(2)
    vec0dot5 = float4(0.5)
    vec255i = int4(255)
    vec1i = int4(1)

  let
    s = (x + y) * F2
    i = floor(x + s)
    j = floor(y + s)
    t = (i + j) * G2
    x0 = x - (i - t)
    y0 = y - (j - t)
    gt = x0 > y0
    i1 = blend(vec0, vec1, gt)
    j1 = blend(vec1, vec0, gt)
    x1 = x0 - i1 + G2
    y1 = y0 - j1 + G2
    x2 = x0 - vec1 + vec2 * G2
    y2 = y0 - vec1 + vec2 * G2
    ii = toInts(i) and vec255i
    jj = toInts(j) and vec255i
    t0 = vec0dot5 - x0 * x0 - y0 * y0
    t1 = vec0dot5 - x1 * x1 - y1 * y1
    t2 = vec0dot5 - x2 * x2 - y2 * y2
    t0gt = t0 > vec0
    t1gt = t1 > vec0
    t2gt = t2 > vec0

  # Set up the gradient vectors
  let
    jji = cast[array[4, int32]](jj)
    jjj1 = cast[array[4, int32]](addInts(jj, toInts(j1)))
    jjplus1 = cast[array[4, int32]](addInts(jj, vec1i))
    jjperm = cast[Int4]([
      simplex.perm[jji[0]],
      simplex.perm[jji[1]],
      simplex.perm[jji[2]],
      simplex.perm[jji[3]]
    ])
    jjj1perm = cast[Int4]([
      simplex.perm[jjj1[0]],
      simplex.perm[jjj1[1]],
      simplex.perm[jjj1[2]],
      simplex.perm[jjj1[3]]
    ])
    jjplus1perm = cast[Int4]([
      simplex.perm[jjplus1[0]],
      simplex.perm[jjplus1[1]],
      simplex.perm[jjplus1[2]],
      simplex.perm[jjplus1[3]]
    ])
    iiplusjjperm = cast[array[4, int32]](addInts(ii, jjperm))
    iii1plusjjj1perm = cast[array[4, int32]](addInts(
      addInts(ii, toInts(i1)), jjj1perm)
    )
    iiplus1plusjjplus1perm = cast[array[4, int32]](
      addInts(addInts(ii, vec1i), jjplus1perm)
    )
    gx0 = cast[Float4]([
      grad3[simplex.permMod12[iiplusjjperm[0]]][0],
      grad3[simplex.permMod12[iiplusjjperm[1]]][0],
      grad3[simplex.permMod12[iiplusjjperm[2]]][0],
      grad3[simplex.permMod12[iiplusjjperm[3]]][0],
    ])
    gy0 = cast[Float4]([
      grad3[simplex.permMod12[iiplusjjperm[0]]][1],
      grad3[simplex.permMod12[iiplusjjperm[1]]][1],
      grad3[simplex.permMod12[iiplusjjperm[2]]][1],
      grad3[simplex.permMod12[iiplusjjperm[3]]][1],
    ])
    gx1 = cast[Float4]([
      grad3[simplex.permMod12[iii1plusjjj1perm[0]]][0],
      grad3[simplex.permMod12[iii1plusjjj1perm[1]]][0],
      grad3[simplex.permMod12[iii1plusjjj1perm[2]]][0],
      grad3[simplex.permMod12[iii1plusjjj1perm[3]]][0]
    ])
    gy1 = cast[Float4]([
      grad3[simplex.permMod12[iii1plusjjj1perm[0]]][1],
      grad3[simplex.permMod12[iii1plusjjj1perm[1]]][1],
      grad3[simplex.permMod12[iii1plusjjj1perm[2]]][1],
      grad3[simplex.permMod12[iii1plusjjj1perm[3]]][1]
    ])
    gx2 = cast[Float4]([
      grad3[simplex.permMod12[iiplus1plusjjplus1perm[0]]][0],
      grad3[simplex.permMod12[iiplus1plusjjplus1perm[1]]][0],
      grad3[simplex.permMod12[iiplus1plusjjplus1perm[2]]][0],
      grad3[simplex.permMod12[iiplus1plusjjplus1perm[3]]][0]
    ])
    gy2 = cast[Float4]([
      grad3[simplex.permMod12[iiplus1plusjjplus1perm[0]]][1],
      grad3[simplex.permMod12[iiplus1plusjjplus1perm[1]]][1],
      grad3[simplex.permMod12[iiplus1plusjjplus1perm[2]]][1],
      grad3[simplex.permMod12[iiplus1plusjjplus1perm[3]]][1]
    ])

  let
    n0 = blend(vec0, vec1, t0gt) * t0 * t0 * t0 * t0 * (gx0 * x0 + gy0 * y0)
    n1 = blend(vec0, vec1, t1gt) * t1 * t1 * t1 * t1 * (gx1 * x1 + gy1 * y1)
    n2 = blend(vec0, vec1, t2gt) * t2 * t2 * t2 * t2 * (gx2 * x2 + gy2 * y2)

  float4(70) * (n0 + n1 + n2)

func row4(
  simplex: Simplex, x: array[4, float32], y: Float4
): array[4, Float4] =
  ## Evaluates four columns using the supplied scaled coordinates.
  for i in 0 ..< 4:
    result[i] = simplex.column4(x[i], y)

func grid4(
  simplex: Simplex, x, y: float32, offsetX, offsetY: int
): array[4, Float4] =
  ## Generates a 4x4 2D noise grid based on the Simplex parameters.
  ## Applies integer grid offsets before scaling coordinates for each octave.
  ## Uses SIMD instructions.

  if simplex.octaves == 0:
    failOctaves()

  let
    xs = coordinates(x, offsetX)
    ys = coordinates(y, offsetY)
  var
    totals: array[4, Float4]
    amplitude = float4(simplex.amplitude)
    gain = float4(simplex.gain)
    frequency = simplex.frequency

  for _ in 0 ..< simplex.octaves:
    let
      scale = float4(frequency)
      rows = simplex.row4(cast[array[4, float32]](xs * scale), ys * scale)
    for i in 0 ..< 4:
      totals[i] += rows[i] * amplitude

    amplitude = amplitude * gain
    frequency *= simplex.lacunarity

  if simplex.octaves > 1:
    let octaves = float4(simplex.octaves.float32)
    for i in 0 ..< 4:
      totals[i] /= octaves

  totals

func layer4(simplex: Simplex, x, y: float32, z: Float4): Float4 =
  ## Evaluates four 3D samples at already scaled coordinates.
  let
    F3 = float4(F3)
    G3 = float4(G3)
    x = float4(x)
    y = float4(y)
    vec0 = float4(0)
    vec1 = float4(1)
    vec2 = float4(2)
    vec3 = float4(3)
    vec0dot6 = float4(0.6)
    vec255i = int4(255)
    vec1i = int4(1)

  let
    s = (x + y + z) * F3
    i = floor(x + s)
    j = floor(y + s)
    k = floor(z + s)
    t = (i + j + k) * G3
    x0 = x - (i - t)
    y0 = y - (j - t)
    z0 = z - (k - t)
    i1 = blend(vec0, vec1, (x0 >= y0 and x0 >= z0))
    j1 = blend(vec0, vec1, (x0 < y0 and y0 >= z0))
    k1 = blend(vec0, vec1, (x0 < z0 and y0 < z0))
    i2 = blend(vec0, vec1, (x0 >= y0 or x0 >= z0))
    j2 = blend(vec0, vec1, (x0 < y0 or y0 >= z0))
    k2 = blend(vec0, vec1, (x0 < z0 or y0 < z0))
    x1 = x0 - i1 + G3
    y1 = y0 - j1 + G3
    z1 = z0 - k1 + G3
    x2 = x0 - i2 + vec2 * G3
    y2 = y0 - j2 + vec2 * G3
    z2 = z0 - k2 + vec2 * G3
    x3 = x0 - vec1 + vec3 * G3
    y3 = y0 - vec1 + vec3 * G3
    z3 = z0 - vec1 + vec3 * G3
    ii = toInts(i) and vec255i
    jj = toInts(j) and vec255i
    kk = toInts(k) and vec255i
    t0 = vec0dot6 - x0 * x0 - y0 * y0 - z0 * z0
    t1 = vec0dot6 - x1 * x1 - y1 * y1 - z1 * z1
    t2 = vec0dot6 - x2 * x2 - y2 * y2 - z2 * z2
    t3 = vec0dot6 - x3 * x3 - y3 * y3 - z3 * z3
    t0gt = t0 > vec0
    t1gt = t1 > vec0
    t2gt = t2 > vec0
    t3gt = t3 > vec0

  # Set up the gradient vectors
  let
    kki = cast[array[4, int32]](kk)
    kkk1 = cast[array[4, int32]](addInts(kk, toInts(k1)))
    kkk2 = cast[array[4, int32]](addInts(kk, toInts(k2)))
    kkplus1 = cast[array[4, int32]](addInts(kk, vec1i))
    kkperm = cast[Int4]([
      simplex.perm[kki[0]],
      simplex.perm[kki[1]],
      simplex.perm[kki[2]],
      simplex.perm[kki[3]]
    ])
    kkk1perm = cast[Int4]([
      simplex.perm[kkk1[0]],
      simplex.perm[kkk1[1]],
      simplex.perm[kkk1[2]],
      simplex.perm[kkk1[3]]
    ])
    kkk2perm = cast[Int4]([
      simplex.perm[kkk2[0]],
      simplex.perm[kkk2[1]],
      simplex.perm[kkk2[2]],
      simplex.perm[kkk2[3]]
    ])
    kkplus1perm = cast[Int4]([
      simplex.perm[kkplus1[0]],
      simplex.perm[kkplus1[1]],
      simplex.perm[kkplus1[2]],
      simplex.perm[kkplus1[3]]
    ])
    jjpluskkperm = cast[array[4, int32]](addInts(jj, kkperm))
    jjj1pluskkk1perm = cast[array[4, int32]](
      addInts(addInts(jj, toInts(j1)), kkk1perm)
    )
    jjj2pluskkk2perm = cast[array[4, int32]](
      addInts(addInts(jj, toInts(j2)), kkk2perm)
    )
    jjplus1pluskkplus1perm = cast[array[4, int32]](
      addInts(addInts(jj, vec1i), kkplus1perm)
    )
    jjpluskpermperm = cast[Int4]([
      simplex.perm[jjpluskkperm[0]],
      simplex.perm[jjpluskkperm[1]],
      simplex.perm[jjpluskkperm[2]],
      simplex.perm[jjpluskkperm[3]]
    ])
    jjj1pluskkk1permperm = cast[Int4]([
      simplex.perm[jjj1pluskkk1perm[0]],
      simplex.perm[jjj1pluskkk1perm[1]],
      simplex.perm[jjj1pluskkk1perm[2]],
      simplex.perm[jjj1pluskkk1perm[3]]
    ])
    jjj2pluskkk2permperm = cast[Int4]([
      simplex.perm[jjj2pluskkk2perm[0]],
      simplex.perm[jjj2pluskkk2perm[1]],
      simplex.perm[jjj2pluskkk2perm[2]],
      simplex.perm[jjj2pluskkk2perm[3]],
    ])
    jjplus1pluskkplus1permperm = cast[Int4]([
      simplex.perm[jjplus1pluskkplus1perm[0]],
      simplex.perm[jjplus1pluskkplus1perm[1]],
      simplex.perm[jjplus1pluskkplus1perm[2]],
      simplex.perm[jjplus1pluskkplus1perm[3]]
    ])
    iiplusjjpluskpermperm = cast[array[4, int32]](
      addInts(ii, jjpluskpermperm)
    )
    iii1plusjjj1pluskkk1permperm = cast[array[4, int32]](
      addInts(addInts(ii, toInts(i1)), jjj1pluskkk1permperm)
    )
    iii2plusjjj2pluskkk2permperm = cast[array[4, int32]](
      addInts(addInts(ii, toInts(i2)), jjj2pluskkk2permperm)
    )
    iiplus1plusjjplus1pluskkplus1permperm = cast[array[4, int32]](
      addInts(addInts(ii, vec1i), jjplus1pluskkplus1permperm)
    )
    gx0 = cast[Float4]([
      grad3[simplex.permMod12[iiplusjjpluskpermperm[0]]][0],
      grad3[simplex.permMod12[iiplusjjpluskpermperm[1]]][0],
      grad3[simplex.permMod12[iiplusjjpluskpermperm[2]]][0],
      grad3[simplex.permMod12[iiplusjjpluskpermperm[3]]][0]
    ])
    gy0 = cast[Float4]([
      grad3[simplex.permMod12[iiplusjjpluskpermperm[0]]][1],
      grad3[simplex.permMod12[iiplusjjpluskpermperm[1]]][1],
      grad3[simplex.permMod12[iiplusjjpluskpermperm[2]]][1],
      grad3[simplex.permMod12[iiplusjjpluskpermperm[3]]][1]
    ])
    gz0 = cast[Float4]([
      grad3[simplex.permMod12[iiplusjjpluskpermperm[0]]][2],
      grad3[simplex.permMod12[iiplusjjpluskpermperm[1]]][2],
      grad3[simplex.permMod12[iiplusjjpluskpermperm[2]]][2],
      grad3[simplex.permMod12[iiplusjjpluskpermperm[3]]][2]
    ])
    gx1 = cast[Float4]([
      grad3[simplex.permMod12[iii1plusjjj1pluskkk1permperm[0]]][0],
      grad3[simplex.permMod12[iii1plusjjj1pluskkk1permperm[1]]][0],
      grad3[simplex.permMod12[iii1plusjjj1pluskkk1permperm[2]]][0],
      grad3[simplex.permMod12[iii1plusjjj1pluskkk1permperm[3]]][0]
    ])
    gy1 = cast[Float4]([
      grad3[simplex.permMod12[iii1plusjjj1pluskkk1permperm[0]]][1],
      grad3[simplex.permMod12[iii1plusjjj1pluskkk1permperm[1]]][1],
      grad3[simplex.permMod12[iii1plusjjj1pluskkk1permperm[2]]][1],
      grad3[simplex.permMod12[iii1plusjjj1pluskkk1permperm[3]]][1]
    ])
    gz1 = cast[Float4]([
      grad3[simplex.permMod12[iii1plusjjj1pluskkk1permperm[0]]][2],
      grad3[simplex.permMod12[iii1plusjjj1pluskkk1permperm[1]]][2],
      grad3[simplex.permMod12[iii1plusjjj1pluskkk1permperm[2]]][2],
      grad3[simplex.permMod12[iii1plusjjj1pluskkk1permperm[3]]][2]
    ])
    gx2 = cast[Float4]([
      grad3[simplex.permMod12[iii2plusjjj2pluskkk2permperm[0]]][0],
      grad3[simplex.permMod12[iii2plusjjj2pluskkk2permperm[1]]][0],
      grad3[simplex.permMod12[iii2plusjjj2pluskkk2permperm[2]]][0],
      grad3[simplex.permMod12[iii2plusjjj2pluskkk2permperm[3]]][0]
    ])
    gy2 = cast[Float4]([
      grad3[simplex.permMod12[iii2plusjjj2pluskkk2permperm[0]]][1],
      grad3[simplex.permMod12[iii2plusjjj2pluskkk2permperm[1]]][1],
      grad3[simplex.permMod12[iii2plusjjj2pluskkk2permperm[2]]][1],
      grad3[simplex.permMod12[iii2plusjjj2pluskkk2permperm[3]]][1]
    ])
    gz2 = cast[Float4]([
      grad3[simplex.permMod12[iii2plusjjj2pluskkk2permperm[0]]][2],
      grad3[simplex.permMod12[iii2plusjjj2pluskkk2permperm[1]]][2],
      grad3[simplex.permMod12[iii2plusjjj2pluskkk2permperm[2]]][2],
      grad3[simplex.permMod12[iii2plusjjj2pluskkk2permperm[3]]][2]
    ])
    gx3 = cast[Float4]([
      grad3[simplex.permMod12[iiplus1plusjjplus1pluskkplus1permperm[0]]][0],
      grad3[simplex.permMod12[iiplus1plusjjplus1pluskkplus1permperm[1]]][0],
      grad3[simplex.permMod12[iiplus1plusjjplus1pluskkplus1permperm[2]]][0],
      grad3[simplex.permMod12[iiplus1plusjjplus1pluskkplus1permperm[3]]][0]
    ])
    gy3 = cast[Float4]([
      grad3[simplex.permMod12[iiplus1plusjjplus1pluskkplus1permperm[0]]][1],
      grad3[simplex.permMod12[iiplus1plusjjplus1pluskkplus1permperm[1]]][1],
      grad3[simplex.permMod12[iiplus1plusjjplus1pluskkplus1permperm[2]]][1],
      grad3[simplex.permMod12[iiplus1plusjjplus1pluskkplus1permperm[3]]][1]
    ])
    gz3 = cast[Float4]([
      grad3[simplex.permMod12[iiplus1plusjjplus1pluskkplus1permperm[0]]][2],
      grad3[simplex.permMod12[iiplus1plusjjplus1pluskkplus1permperm[1]]][2],
      grad3[simplex.permMod12[iiplus1plusjjplus1pluskkplus1permperm[2]]][2],
      grad3[simplex.permMod12[iiplus1plusjjplus1pluskkplus1permperm[3]]][2]
    ])

  let
    n0 = blend(vec0, vec1, t0gt) *
      t0 * t0 * t0 * t0 * (gx0 * x0 + gy0 * y0 + gz0 * z0)
    n1 = blend(vec0, vec1, t1gt) *
      t1 * t1 * t1 * t1 * (gx1 * x1 + gy1 * y1 + gz1 * z1)
    n2 = blend(vec0, vec1, t2gt) *
      t2 * t2 * t2 * t2 * (gx2 * x2 + gy2 * y2 + gz2 * z2)
    n3 = blend(vec0, vec1, t3gt) *
      t3 * t3 * t3 * t3 * (gx3 * x3 + gy3 * y3 + gz3 * z3)

  float4(32) * (n0 + n1 + n2 + n3)

func column4(
  simplex: Simplex, x: float32, y: array[4, float32], z: Float4
): array[4, Float4] =
  ## Evaluates four layers using the supplied scaled coordinates.
  for i in 0 ..< 4:
    result[i] = simplex.layer4(x, y[i], z)

func row4(
  simplex: Simplex, x, y: array[4, float32], z: Float4
): array[4, array[4, Float4]] =
  ## Evaluates four columns using the supplied scaled coordinates.
  for i in 0 ..< 4:
    result[i] = simplex.column4(x[i], y, z)

func grid4(
  simplex: Simplex, x, y, z: float32, offsetX, offsetY, offsetZ: int
): array[4, array[4, Float4]] =
  ## Generates a 4x4x4 3D noise grid based on the Simplex parameters.
  ## Applies integer grid offsets before scaling coordinates for each octave.
  ## Uses SIMD instructions.

  if simplex.octaves == 0:
    failOctaves()

  let
    xs = coordinates(x, offsetX)
    ys = coordinates(y, offsetY)
    zs = coordinates(z, offsetZ)
  var
    totals: array[4, array[4, Float4]]
    amplitude = float4(simplex.amplitude)
    gain = float4(simplex.gain)
    frequency = simplex.frequency

  for _ in 0 ..< simplex.octaves:
    let
      scale = float4(frequency)
      rows = simplex.row4(
        cast[array[4, float32]](xs * scale),
        cast[array[4, float32]](ys * scale),
        zs * scale
      )
    for i in 0 ..< 4:
      for j in 0 ..< 4:
        totals[i][j] += rows[i][j] * amplitude

    amplitude = amplitude * gain
    frequency *= simplex.lacunarity

  if simplex.octaves > 1:
    let octaves = float4(simplex.octaves.float32)
    for i in 0 ..< 4:
      for j in 0 ..< 4:
        totals[i][j] /= octaves

  totals

func grid*(
  simplex: Simplex, start: (float32, float32), dimens: (int, int)
): Grid =
  ## Beginning at position start, generate a grid of 2D noise based on the
  ## Simplex parameters. The width and height of the grid is set by the dimens
  ## parameter.

  let
    (x, y) = start
    (width, height) = dimens

  result = Grid()
  result.width = width
  result.height = height
  result.depth = 1
  result.frames = 1
  result.values.setLen(width * height)

  var widthDone, heightDone: int
  template fillSse2OrNeon() =
    for i in countup(0, width - 4, 4):
      for j in countup(0, height - 4, 4):
        let grid4 = simplex.grid4(x, y, i, j)
        for a in 0 ..< 4:
          store(
            result.values[result.valueIndex(a + i, j, 0)].unsafeAddr,
            grid4[a]
          )
        widthDone = i + 4
        heightDone = j + 4

  when allowAvx2:
    if hasAvx2() and width > 0 and height >= 8:
      simplex.fillGridAvx2(start, result)
      widthDone = width
      heightDone = height - height mod 8
    else:
      fillSse2OrNeon()
  else:
    fillSse2OrNeon()

  # Fill any width leftover
  for i in widthDone ..< width:
    for j in 0 ..< heightDone:
      result[i, j, 0] = simplex.value(x + i.float32, y + j.float32)

  # Fill any height leftover
  for i in 0 ..< width:
    for j in heightDone ..< height:
      result[i, j, 0] = simplex.value(x + i.float32, y + j.float32)

func grid*(
  simplex: Simplex, start: (int, int), dimens: (int, int)
): Grid {.inline.} =
  ## Helper for working with ints.
  simplex.grid((start[0].float32, start[1].float32), dimens)

func grid*(
  simplex: Simplex, start: (float32, float32, float32), dimens: (int, int, int)
): Grid =
  ## Beginning at position start, generate a grid of 3D noise based on the
  ## Simplex parameters. The width, depth, and height of the grid is set by
  ## the dimens parameter.

  let
    (x, y, z) = start
    (width, height, depth) = dimens

  result = Grid()
  result.width = width
  result.height = height
  result.depth = depth
  result.frames = 1
  result.values.setLen(width * height * depth)

  var widthDone, heightDone, depthDone: int
  template fillSse2OrNeon() =
    for i in countup(0, width - 4, 4):
      for j in countup(0, height - 4, 4):
        for k in countup(0, depth - 4, 4):
          let grid4 = simplex.grid4(x, y, z, i, j, k)
          for a in 0 ..< 4:
            for b in 0 ..< 4:
              store(
                result.values[result.valueIndex(a + i, b + j, k)].unsafeAddr,
                grid4[a][b]
              )
          widthDone = i + 4
          heightDone = j + 4
          depthDone = k + 4

  when allowAvx2:
    if hasAvx2() and width > 0 and height > 0 and depth >= 8:
      simplex.fillGridAvx2(start, result)
      widthDone = width
      heightDone = height
      depthDone = depth - depth mod 8
    else:
      fillSse2OrNeon()
  else:
    fillSse2OrNeon()

  # Fill any width leftover
  for i in widthDone ..< width:
    for j in 0 ..< heightDone:
      for k in 0 ..< depthDone:
        result[i, j, k] = simplex.value(
          x + i.float32, y + j.float32, z + k.float32
        )

  # Fill any height leftover
  for i in 0 ..< width:
    for j in heightDone ..< height:
      for k in 0 ..< depthDone:
        result[i, j, k] = simplex.value(
          x + i.float32, y + j.float32, z + k.float32
        )

  # Fill any depth leftover
  for i in 0 ..< width:
    for j in 0 ..< height:
      for k in depthDone ..< depth:
        result[i, j, k] = simplex.value(
          x + i.float32, y + j.float32, z + k.float32
        )

func grid*(
  simplex: Simplex, start: (int, int, int), dimens: (int, int, int)
): Grid {.inline.} =
  ## Helper for working with ints.
  simplex.grid((start[0].float32, start[1].float32, start[2].float32), dimens)

func gradient1(hash: int32): float32 {.inline.} =
  ## Selects one of the sixteen signed 1D gradients.
  result = (1 + (hash and 7)).float32
  if (hash and 8) != 0:
    result = -result

func noise(simplex: Simplex, x: float32): float32 =
  ## Evaluates native 1D simplex noise from its two lattice endpoints.
  let
    cell = floor(x)
    index = cell.int32 and 255
    x0 = x - cell
    x1 = x0 - 1.0'f
    t0 = 1.0'f - x0 * x0
    t1 = 1.0'f - x1 * x1
    t0Squared = t0 * t0
    t1Squared = t1 * t1
    n0 = t0Squared * t0Squared * (gradient1(simplex.perm[index]) * x0)
    n1 = t1Squared * t1Squared * (gradient1(simplex.perm[index + 1]) * x1)
  # The unscaled maximum is 8 * (3 / 4)^4 = 2.53125.
  0.395'f * (n0 + n1)

func value*(simplex: Simplex, x: float32): float32 =
  ## Generates the 1D noise value at x based on the Simplex parameters.
  if simplex.octaves <= 0:
    failOctaves()
  var
    total: float32
    amplitude = simplex.amplitude
    frequency = simplex.frequency
  for _ in 0 ..< simplex.octaves:
    total += simplex.noise(x * frequency) * amplitude
    amplitude *= simplex.gain
    frequency *= simplex.lacunarity
  total / simplex.octaves.float32

func value*(simplex: Simplex, x: int): float32 {.inline.} =
  ## Samples native 1D simplex noise at an integer coordinate.
  simplex.value(x.float32)

func noise1x4(simplex: Simplex, x: Float4): Float4 =
  ## Evaluates four native 1D samples at already scaled coordinates.
  let
    cell = floor(x)
    indices = cast[array[4, int32]](toInts(cell) and int4(255))
    x0 = x - cell
    x1 = x0 - float4(1)
    t0 = float4(1) - x0 * x0
    t1 = float4(1) - x1 * x1
    t0Squared = t0 * t0
    t1Squared = t1 * t1
  var gradients0, gradients1: array[4, float32]
  for i in 0 ..< 4:
    gradients0[i] = gradient1(simplex.perm[indices[i]])
    gradients1[i] = gradient1(simplex.perm[indices[i] + 1])
  let
    n0 = t0Squared * t0Squared * (cast[Float4](gradients0) * x0)
    n1 = t1Squared * t1Squared * (cast[Float4](gradients1) * x1)
  float4(0.395'f) * (n0 + n1)

func values1x4(simplex: Simplex, x: Float4): Float4 =
  ## Evaluates four world coordinates with the configured octave layering.
  if simplex.octaves <= 0:
    failOctaves()
  var
    total = float4(0)
    amplitude = float4(simplex.amplitude)
    frequency = simplex.frequency
  let gain = float4(simplex.gain)
  for _ in 0 ..< simplex.octaves:
    total += simplex.noise1x4(x * float4(frequency)) * amplitude
    amplitude = amplitude * gain
    frequency *= simplex.lacunarity
  total / float4(simplex.octaves.float32)

func noise(simplex: Simplex, x, y, z, w: float32): float32 =
  ## Evaluates native 4D simplex noise at five simplex corners.
  let
    point = [x, y, z, w]
    skew = (x + y + z + w) * F4
    cells = [floor(x + skew), floor(y + skew),
      floor(z + skew), floor(w + skew)]
    unskew = (cells[0] + cells[1] + cells[2] + cells[3]) * G4
  var
    origin: array[4, float32]
    ranks: array[4, int]
  for axis in 0 ..< 4:
    origin[axis] = point[axis] - (cells[axis] - unskew)
  for i in 0 ..< 3:
    for j in i + 1 ..< 4:
      if origin[i] > origin[j]:
        inc ranks[i]
      else:
        inc ranks[j]

  for corner in 0 .. 4:
    var
      coordinates: array[4, int32]
      distances: array[4, float32]
      falloff = 0.5'f
    for axis in 0 ..< 4:
      let offset = (ranks[axis] >= 4 - corner).int32
      coordinates[axis] = (cells[axis].int32 and 255) + offset
      distances[axis] = origin[axis] - offset.float32 + corner.float32 * G4
      falloff -= distances[axis] * distances[axis]
    if falloff > 0:
      let gradient = Grad4[simplex.perm[
        coordinates[0] + simplex.perm[
          coordinates[1] + simplex.perm[
            coordinates[2] + simplex.perm[coordinates[3]]
          ]
        ]
      ] and 31]
      var projection: float32
      for axis in 0 ..< 4:
        projection += gradient[axis] * distances[axis]
      result += falloff * falloff * falloff * falloff * projection
  result *= 62

func noise4x4(simplex: Simplex, x, y, z, w: Float4): Float4 =
  ## Evaluates four native 4D simplex samples using SIMD arithmetic.
  let
    point = [x, y, z, w]
    zero = float4(0)
    one = float4(1)
    skew = (x + y + z + w) * float4(F4)
    cells = [floor(x + skew), floor(y + skew),
      floor(z + skew), floor(w + skew)]
    unskew = (cells[0] + cells[1] + cells[2] + cells[3]) * float4(G4)
  var
    origin, ranks: array[4, Float4]
    indices: array[4, Int4]
    total = zero
  for axis in 0 ..< 4:
    origin[axis] = point[axis] - (cells[axis] - unskew)
    ranks[axis] = zero
    indices[axis] = toInts(cells[axis]) and int4(255)
  for i in 0 ..< 3:
    for j in i + 1 ..< 4:
      let greater = origin[i] > origin[j]
      ranks[i] += blend(zero, one, greater)
      ranks[j] += blend(one, zero, greater)

  for corner in 0 .. 4:
    let shift = float4(corner.float32 * G4)
    var
      coordinates: array[4, array[4, int32]]
      distances: array[4, Float4]
      gradients: array[4, array[4, float32]]
      falloff = float4(0.5'f)
      projection = zero
    for axis in 0 ..< 4:
      let offset = blend(zero, one, ranks[axis] >= float4((4 - corner).float32))
      coordinates[axis] = cast[array[4, int32]](addInts(
        indices[axis], toInts(offset)
      ))
      distances[axis] = origin[axis] - offset + shift
      falloff = falloff - distances[axis] * distances[axis]
    for lane in 0 ..< 4:
      let gradient = Grad4[simplex.perm[
        coordinates[0][lane] + simplex.perm[
          coordinates[1][lane] + simplex.perm[
            coordinates[2][lane] + simplex.perm[coordinates[3][lane]]
          ]
        ]
      ] and 31]
      for axis in 0 ..< 4:
        gradients[axis][lane] = gradient[axis]
    for axis in 0 ..< 4:
      projection += cast[Float4](gradients[axis]) * distances[axis]
    total += blend(
      zero,
      falloff * falloff * falloff * falloff * projection,
      falloff > zero
    )
  float4(62) * total

func value*(simplex: Simplex, x, y, z, w: float32): float32 =
  ## Generates the 4D noise value using the configured octave parameters.
  if simplex.octaves <= 0:
    failOctaves()
  var
    total: float32
    amplitude = simplex.amplitude
    frequency = simplex.frequency
  for _ in 0 ..< simplex.octaves:
    total += simplex.noise(
      x * frequency,
      y * frequency,
      z * frequency,
      w * frequency
    ) * amplitude
    amplitude *= simplex.gain
    frequency *= simplex.lacunarity
  total / simplex.octaves.float32

func value*(simplex: Simplex, x, y, z, w: int): float32 {.inline.} =
  ## Samples 4D noise at integer coordinates.
  simplex.value(x.float32, y.float32, z.float32, w.float32)

func values4x4(simplex: Simplex, x, y, z, w: Float4): Float4 =
  ## Samples four world positions with the configured 4D noise octaves.
  if simplex.octaves <= 0:
    failOctaves()
  var
    total = float4(0)
    amplitude = float4(simplex.amplitude)
    frequency = simplex.frequency
  let gain = float4(simplex.gain)
  for _ in 0 ..< simplex.octaves:
    let scale = float4(frequency)
    total += simplex.noise4x4(
      x * scale,
      y * scale,
      z * scale,
      w * scale
    ) * amplitude
    amplitude = amplitude * gain
    frequency *= simplex.lacunarity
  total / float4(simplex.octaves.float32)

func grid*(simplex: Simplex, start: float32, width: int): Grid =
  ## Generates 1D noise at unit intervals using SIMD for complete groups of samples.
  if simplex.octaves <= 0:
    failOctaves()
  let length = gridLength([width])
  result = Grid(
    width: width,
    height: 1,
    depth: 1,
    frames: 1,
    values: newSeq[float32](length)
  )
  var widthDone: int
  when allowAvx2:
    if hasAvx2() and width >= 8:
      simplex.fillGridAvx2(start, result)
      widthDone = width - width mod 8
  for x in countup(widthDone, width - 4, 4):
    store(
      result.values[x].unsafeAddr,
      simplex.values1x4(coordinates(start, x))
    )
  for x in width div 4 * 4 ..< width:
    result.values[x] = simplex.value(start + x.float32)

func grid*(simplex: Simplex, start, width: int): Grid {.inline.} =
  ## Generates a 1D grid beginning at an integer coordinate.
  simplex.grid(start.float32, width)

func grid*(
  simplex: Simplex,
  start: (float32, float32, float32, float32),
  dimens: (int, int, int, int)
): Grid =
  ## Generates a 4D grid at unit intervals, with SIMD along the z axis.
  ## The fourth dimension is frames; each frame contains one 3D slice.
  if simplex.octaves <= 0:
    failOctaves()
  let
    (width, height, depth, frames) = dimens
    length = gridLength([width, height, depth, frames])
  result = Grid(
    width: width,
    height: height,
    depth: depth,
    frames: frames,
    values: newSeq[float32](length)
  )
  if length == 0:
    return
  var depthDone: int
  when allowAvx2:
    if hasAvx2() and depth >= 8:
      simplex.fillGridAvx2(start, result)
      depthDone = depth - depth mod 8
  for w in 0 ..< frames:
    let
      worldW = start[3] + w.float32
      ws = float4(worldW)
    for x in 0 ..< width:
      let
        worldX = start[0] + x.float32
        xs = float4(worldX)
      for y in 0 ..< height:
        let
          worldY = start[1] + y.float32
          ys = float4(worldY)
        for z in countup(depthDone, depth - 4, 4):
          store(
            result.values[result.valueIndex(x, y, z, w)].unsafeAddr,
            simplex.values4x4(xs, ys, coordinates(start[2], z), ws)
          )
        for z in depth div 4 * 4 ..< depth:
          result.values[result.valueIndex(x, y, z, w)] = simplex.value(
            worldX,
            worldY,
            start[2] + z.float32,
            worldW
          )

func grid*(
  simplex: Simplex,
  start: (int, int, int, int),
  dimens: (int, int, int, int)
): Grid {.inline.} =
  ## Generates a 4D grid beginning at integer coordinates.
  simplex.grid(
    (start[0].float32, start[1].float32, start[2].float32, start[3].float32),
    dimens
  )

when defined(release):
  {.pop.}
