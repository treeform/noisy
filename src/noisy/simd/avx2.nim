import nimsimd/avx2, ../simplexcommon

when defined(gcc) or defined(clang):
  # Preserve the scalar arithmetic order; AVX2 does not require FMA support.
  {.localPassC: "-mavx2 -ffp-contract=off".}
elif defined(vcc):
  {.localPassC: "/arch:AVX2 /fp:strict".}

include ../kernelchecks

# Templates keep vector values within this module and avoid a vector ABI.
template m256(a: float32): M256 = mm256_set1_ps(a)
template m256i(a: int32): M256i = mm256_set1_epi32(a)
template `and`(a, b: M256): M256 = mm256_and_ps(a, b)
template `and`(a, b: M256i): M256i = mm256_and_si256(a, b)
template `or`(a, b: M256): M256 = mm256_or_ps(a, b)
template `>`(a, b: M256): M256 = mm256_cmp_ps(a, b, 30) # _CMP_GT_OQ
template `>=`(a, b: M256): M256 = mm256_cmp_ps(a, b, 29) # _CMP_GE_OQ
template `<`(a, b: M256): M256 = mm256_cmp_ps(a, b, 17) # _CMP_LT_OQ
template `+`(a, b: M256): M256 = mm256_add_ps(a, b)
template `-`(a, b: M256): M256 = mm256_sub_ps(a, b)
template `*`(a, b: M256): M256 = mm256_mul_ps(a, b)
template `/`(a, b: M256): M256 = mm256_div_ps(a, b)
template `+`(a, b: M256i): M256i = mm256_add_epi32(a, b)
template `+=`(a: var M256, b: M256) = a = a + b

template blend(a, b, mask: M256): M256 =
  mm256_blendv_ps(a, b, mask)

# Typed pointers also work with the C++ backend (nimsimd's gather bindings use
# void pointers, which C++ cannot implicitly convert to the intrinsic types).
func gather(p: ptr int32, indices: M256i, scale: int32): M256i
  {.importc: "_mm256_i32gather_epi32", header: "immintrin.h".}
func gather(p: ptr float32, indices: M256i, scale: int32): M256
  {.importc: "_mm256_i32gather_ps", header: "immintrin.h".}

template lookupPerm(simplex: Simplex, indices: M256i): M256i =
  gather(simplex.perm[0].unsafeAddr, indices and m256i(511), 4) and m256i(511)

template gradients(simplex: Simplex, indices: M256i): M256i =
  block:
    let values = gather(simplex.permMod12[0].unsafeAddr, indices and m256i(511), 4)
    values and mm256_cmpgt_epi32(values, m256i(-1)) and
      mm256_cmpgt_epi32(m256i(12), values)

func gradientAxis(axis: int): array[12, float32] =
  for i in 0 ..< grad3.len:
    result[i] = grad3[i][axis]

const
  gradX = gradientAxis(0)
  gradY = gradientAxis(1)
  gradZ = gradientAxis(2)

template dot(indices: M256i, x, y: M256): M256 =
  gather(cast[ptr float32](gradX.unsafeAddr), indices, 4) * x +
    gather(cast[ptr float32](gradY.unsafeAddr), indices, 4) * y

template dot(indices: M256i, x, y, z: M256): M256 =
  dot(indices, x, y) +
    gather(cast[ptr float32](gradZ.unsafeAddr), indices, 4) * z

template noise8(simplex: Simplex, xIn, yIn: M256): M256 =
  block:
    let x = xIn
    let y = yIn
    let
      F2 = m256(F2)
      G2 = m256(G2)
      vec0 = m256(0)
      vec1 = m256(1)
      vec2 = m256(2)
      vec0dot5 = m256(0.5)
      vec255i = m256i(255)
      vec1i = m256i(1)

    let
      s = (x + y) * F2
      i = mm256_floor_ps(x + s)
      j = mm256_floor_ps(y + s)
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
      ii = mm256_cvttps_epi32(i) and vec255i
      jj = mm256_cvttps_epi32(j) and vec255i
      t0 = vec0dot5 - x0 * x0 - y0 * y0
      t1 = vec0dot5 - x1 * x1 - y1 * y1
      t2 = vec0dot5 - x2 * x2 - y2 * y2

    let
      gi0 = simplex.gradients(ii + simplex.lookupPerm(jj))
      gi1 = simplex.gradients(ii + mm256_cvttps_epi32(i1) +
        simplex.lookupPerm(jj + mm256_cvttps_epi32(j1)))
      gi2 = simplex.gradients(ii + vec1i + simplex.lookupPerm(jj + vec1i))
      n0 = blend(vec0, vec1, t0 > vec0) * t0 * t0 * t0 * t0 * dot(gi0, x0, y0)
      n1 = blend(vec0, vec1, t1 > vec0) * t1 * t1 * t1 * t1 * dot(gi1, x1, y1)
      n2 = blend(vec0, vec1, t2 > vec0) * t2 * t2 * t2 * t2 * dot(gi2, x2, y2)

    m256(70) * (n0 + n1 + n2)

template noise8(simplex: Simplex, xIn, yIn, zIn: M256): M256 =
  block:
    let x = xIn
    let y = yIn
    let z = zIn
    let
      F3 = m256(F3)
      G3 = m256(G3)
      vec0 = m256(0)
      vec1 = m256(1)
      vec2 = m256(2)
      vec3 = m256(3)
      vec0dot6 = m256(0.6)
      vec255i = m256i(255)
      vec1i = m256i(1)

    let
      s = (x + y + z) * F3
      i = mm256_floor_ps(x + s)
      j = mm256_floor_ps(y + s)
      k = mm256_floor_ps(z + s)
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
      ii = mm256_cvttps_epi32(i) and vec255i
      jj = mm256_cvttps_epi32(j) and vec255i
      kk = mm256_cvttps_epi32(k) and vec255i
      t0 = vec0dot6 - x0 * x0 - y0 * y0 - z0 * z0
      t1 = vec0dot6 - x1 * x1 - y1 * y1 - z1 * z1
      t2 = vec0dot6 - x2 * x2 - y2 * y2 - z2 * z2
      t3 = vec0dot6 - x3 * x3 - y3 * y3 - z3 * z3

    let
      gi0 = simplex.gradients(ii + simplex.lookupPerm(jj + simplex.lookupPerm(kk)))
      gi1 = simplex.gradients(ii + mm256_cvttps_epi32(i1) +
        simplex.lookupPerm(jj + mm256_cvttps_epi32(j1) +
          simplex.lookupPerm(kk + mm256_cvttps_epi32(k1))))
      gi2 = simplex.gradients(ii + mm256_cvttps_epi32(i2) +
        simplex.lookupPerm(jj + mm256_cvttps_epi32(j2) +
          simplex.lookupPerm(kk + mm256_cvttps_epi32(k2))))
      gi3 = simplex.gradients(ii + vec1i +
        simplex.lookupPerm(jj + vec1i + simplex.lookupPerm(kk + vec1i)))
      n0 = blend(vec0, vec1, t0 > vec0) *
        t0 * t0 * t0 * t0 * dot(gi0, x0, y0, z0)
      n1 = blend(vec0, vec1, t1 > vec0) *
        t1 * t1 * t1 * t1 * dot(gi1, x1, y1, z1)
      n2 = blend(vec0, vec1, t2 > vec0) *
        t2 * t2 * t2 * t2 * dot(gi2, x2, y2, z2)
      n3 = blend(vec0, vec1, t3 > vec0) *
        t3 * t3 * t3 * t3 * dot(gi3, x3, y3, z3)

    m256(32) * (n0 + n1 + n2 + n3)

template coordinates(start: float32, index: int): M256 =
  # Add the start before scaling, in the same order as scalar value().
  m256(start) + mm256_setr_ps(
    index.float32, (index + 1).float32, (index + 2).float32,
    (index + 3).float32, (index + 4).float32, (index + 5).float32,
    (index + 6).float32, (index + 7).float32
  )

template gradient1(hash: M256i): M256 =
  block:
    let
      magnitude = mm256_cvtepi32_ps((hash and m256i(7)) + m256i(1))
      sign = mm256_slli_epi32(hash and m256i(8), 28)
    mm256_xor_ps(magnitude, mm256_castsi256_ps(sign))

template noise8(simplex: Simplex, xIn: M256): M256 =
  block:
    let
      x = xIn
      cell = mm256_floor_ps(x)
      index = mm256_cvttps_epi32(cell) and m256i(255)
      x0 = x - cell
      x1 = x0 - m256(1)
      t0 = m256(1) - x0 * x0
      t1 = m256(1) - x1 * x1
      t0Squared = t0 * t0
      t1Squared = t1 * t1
      hash0 = simplex.lookupPerm(index)
      hash1 = simplex.lookupPerm(index + m256i(1))
      n0 = t0Squared * t0Squared * (gradient1(hash0) * x0)
      n1 = t1Squared * t1Squared * (gradient1(hash1) * x1)
    m256(0.395'f) * (n0 + n1)

const grad4Axes = block:
  var axes: array[4, array[32, float32]]
  for axis in 0 ..< 4:
    for i in 0 ..< Grad4.len:
      axes[axis][i] = Grad4[i][axis]
  axes

template noise8(simplex: Simplex, xIn, yIn, zIn, wIn: M256): M256 =
  block:
    let
      point = [xIn, yIn, zIn, wIn]
      zero = m256(0)
      one = m256(1)
      skew = (point[0] + point[1] + point[2] + point[3]) * m256(F4)
      cells = [mm256_floor_ps(point[0] + skew),
        mm256_floor_ps(point[1] + skew), mm256_floor_ps(point[2] + skew),
        mm256_floor_ps(point[3] + skew)]
      unskew = (cells[0] + cells[1] + cells[2] + cells[3]) * m256(G4)
    var
      origin, ranks: array[4, M256]
      indices: array[4, M256i]
      total = zero
    for axis in 0 ..< 4:
      origin[axis] = point[axis] - (cells[axis] - unskew)
      ranks[axis] = zero
      indices[axis] = mm256_cvttps_epi32(cells[axis]) and m256i(255)
    for i in 0 ..< 3:
      for j in i + 1 ..< 4:
        let greater = origin[i] > origin[j]
        ranks[i] += blend(zero, one, greater)
        ranks[j] += blend(one, zero, greater)
    for corner in 0 .. 4:
      let shift = m256(corner.float32 * G4)
      var
        coordinates: array[4, M256i]
        distances: array[4, M256]
        falloff = m256(0.5'f)
        projection = zero
      for axis in 0 ..< 4:
        let offset = blend(zero, one, ranks[axis] >= m256((4 - corner).float32))
        coordinates[axis] = indices[axis] + mm256_cvttps_epi32(offset)
        distances[axis] = origin[axis] - offset + shift
        falloff = falloff - distances[axis] * distances[axis]
      let gradient = simplex.lookupPerm(coordinates[0] +
        simplex.lookupPerm(coordinates[1] +
          simplex.lookupPerm(coordinates[2] +
            simplex.lookupPerm(coordinates[3])))) and m256i(31)
      for axis in 0 ..< 4:
        projection += gather(cast[ptr float32](grad4Axes[axis].unsafeAddr),
          gradient, 4) * distances[axis]
      total += blend(zero, falloff * falloff * falloff * falloff * projection,
        falloff > zero)
    m256(62) * total

func fillGridAvx2*(simplex: Simplex, start: float32, grid: Grid) =
  ## Fills complete groups of eight along X; the caller fills the remainder.
  for x in countup(0, grid.width - 8, 8):
    let vx = coordinates(start, x)
    var
      total = mm256_setzero_ps()
      amplitude = m256(simplex.amplitude)
      frequency = simplex.frequency
    for octave in 0 ..< simplex.octaves:
      total += noise8(simplex, vx * m256(frequency)) * amplitude
      amplitude = amplitude * m256(simplex.gain)
      frequency *= simplex.lacunarity
    total = total / m256(simplex.octaves.float32)
    mm256_storeu_ps(grid.values[x].unsafeAddr, total)

func fillGridAvx2*(simplex: Simplex, start: (float32, float32), grid: Grid) =
  ## Fills complete groups of eight along Y; the caller fills the scalar tail.
  if simplex.octaves == 0:
    raise newException(NoisyError, "Octaves must be > 0")

  for x in 0 ..< grid.width:
    let vx = m256(start[0] + x.float32)
    for y in countup(0, grid.height - 8, 8):
      let vy = coordinates(start[1], y)
      var
        total = mm256_setzero_ps()
        amplitude = m256(simplex.amplitude)
        frequency = simplex.frequency
      for octave in 0 ..< simplex.octaves:
        let f = m256(frequency)
        total = total + noise8(simplex, vx * f, vy * f) * amplitude
        amplitude = amplitude * m256(simplex.gain)
        frequency *= simplex.lacunarity
      if simplex.octaves > 1:
        total = total / m256(simplex.octaves.float32)
      mm256_storeu_ps(grid.values[x * grid.height + y].unsafeAddr, total)

func fillGridAvx2*(
  simplex: Simplex, start: (float32, float32, float32), grid: Grid
) =
  ## Fills complete groups of eight along Z; the caller fills the scalar tail.
  if simplex.octaves == 0:
    raise newException(NoisyError, "Octaves must be > 0")

  for x in 0 ..< grid.width:
    let vx = m256(start[0] + x.float32)
    for y in 0 ..< grid.height:
      let vy = m256(start[1] + y.float32)
      for z in countup(0, grid.depth - 8, 8):
        let vz = coordinates(start[2], z)
        var
          total = mm256_setzero_ps()
          amplitude = m256(simplex.amplitude)
          frequency = simplex.frequency
        for octave in 0 ..< simplex.octaves:
          let f = m256(frequency)
          total = total + noise8(simplex, vx * f, vy * f, vz * f) * amplitude
          amplitude = amplitude * m256(simplex.gain)
          frequency *= simplex.lacunarity
        if simplex.octaves > 1:
          total = total / m256(simplex.octaves.float32)
        mm256_storeu_ps(
          grid.values[(x * grid.height + y) * grid.depth + z].unsafeAddr,
          total
        )

func fillGridAvx2*(
  simplex: Simplex, start: (float32, float32, float32, float32), grid: Grid
) =
  ## Fills eight Z samples at a time within each fourth-axis slice.
  for w in 0 ..< grid.frames:
    let vw = m256(start[3] + w.float32)
    for x in 0 ..< grid.width:
      let vx = m256(start[0] + x.float32)
      for y in 0 ..< grid.height:
        let vy = m256(start[1] + y.float32)
        for z in countup(0, grid.depth - 8, 8):
          let vz = coordinates(start[2], z)
          var
            total = mm256_setzero_ps()
            amplitude = m256(simplex.amplitude)
            frequency = simplex.frequency
          for octave in 0 ..< simplex.octaves:
            let f = m256(frequency)
            total += noise8(simplex, vx * f, vy * f, vz * f, vw * f) * amplitude
            amplitude = amplitude * m256(simplex.gain)
            frequency *= simplex.lacunarity
          total = total / m256(simplex.octaves.float32)
          let index = ((w * grid.width + x) * grid.height + y) * grid.depth + z
          mm256_storeu_ps(grid.values[index].unsafeAddr, total)

{.pop.}
