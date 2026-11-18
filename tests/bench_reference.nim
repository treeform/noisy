import
  std/[algorithm, json, math, monotimes, strformat, times],
  noisy,
  reference/bindings

when not defined(release):
  {.error: "Run benchmarks with -d:release.".}

const
  NoisyBenchmarkRounds {.intdefine.} = 11
  NoisyBenchmarkBatches {.intdefine.} = 4
  NoisyBenchmarkOutput {.strdefine.} = ""
  Start = [-1.375'f, 2.125'f, -0.625'f, 0.375'f]
  CompilerBackend =
    when defined(cpp): "C++"
    else: "C"
  Shapes = [
    [262144, 1, 1, 1],
    [512, 512, 1, 1],
    [64, 64, 64, 1],
    [32, 32, 32, 8]
  ]

when defined(amd64) and (defined(gcc) or defined(clang) or defined(vcc)) and
    not defined(noisyNoAvx):
  import noisy/simd/runtimecheck
  let SimdName = if hasAvx2(): "AVX2" else: "SSE2"
else:
  const SimdName = when defined(arm64): "NEON" else: "SSE2"

type
  Backend = enum
    Reference, Scalar, Vector
  Shape = array[4, int]

var checksum: float64
var results = newJArray()

proc consume(values: seq[float32]) =
  ## Consumes every output outside the timer to prevent dead-code removal.
  var total: float64
  for value in values:
    total += value.float64
  checksum += total

proc generate(
  dimensions: static[int],
  backend: Backend,
  noise: Simplex,
  shape: Shape
): seq[float32] =
  ## Allocates and fills the same output layout through each public path.
  case backend
  of Reference:
    result = newSeq[float32](shape[0] * shape[1] * shape[2] * shape[3])
    referenceGrid(
      result[0].addr,
      dimensions.int32,
      shape[0].int32, shape[1].int32, shape[2].int32, shape[3].int32,
      Start[0].unsafeAddr,
      noise.octaves.int32,
      noise.amplitude, noise.frequency, noise.lacunarity, noise.gain
    )
  of Scalar:
    result = newSeq[float32](shape[0] * shape[1] * shape[2] * shape[3])
    var index = 0
    when dimensions == 1:
      for x in 0 ..< shape[0]:
        result[index] = noise.value(Start[0] + x.float32)
        inc index
    elif dimensions == 2:
      for x in 0 ..< shape[0]:
        let worldX = Start[0] + x.float32
        for y in 0 ..< shape[1]:
          result[index] = noise.value(worldX, Start[1] + y.float32)
          inc index
    elif dimensions == 3:
      for x in 0 ..< shape[0]:
        let worldX = Start[0] + x.float32
        for y in 0 ..< shape[1]:
          let worldY = Start[1] + y.float32
          for z in 0 ..< shape[2]:
            result[index] = noise.value(worldX, worldY, Start[2] + z.float32)
            inc index
    else:
      for w in 0 ..< shape[3]:
        let worldW = Start[3] + w.float32
        for x in 0 ..< shape[0]:
          let worldX = Start[0] + x.float32
          for y in 0 ..< shape[1]:
            let worldY = Start[1] + y.float32
            for z in 0 ..< shape[2]:
              result[index] = noise.value(
                worldX, worldY, Start[2] + z.float32, worldW
              )
              inc index
  of Vector:
    when dimensions == 1:
      result = noise.grid(Start[0], shape[0]).values
    elif dimensions == 2:
      result = noise.grid((Start[0], Start[1]), (shape[0], shape[1])).values
    elif dimensions == 3:
      result = noise.grid(
        (Start[0], Start[1], Start[2]), (shape[0], shape[1], shape[2])
      ).values
    else:
      result = noise.grid(
        (Start[0], Start[1], Start[2], Start[3]),
        (shape[0], shape[1], shape[2], shape[3])
      ).values

proc validate(dimensions: static[int], noise: Simplex, shape: Shape) =
  ## Checks every timed input against the reference before benchmarking.
  let expected = generate(dimensions, Reference, noise, shape)
  for backend in [Scalar, Vector]:
    let actual = generate(dimensions, backend, noise, shape)
    doAssert actual.len == expected.len
    for i in 0 ..< actual.len:
      let
        error = abs(actual[i] - expected[i])
        tolerance = 0.0001'f + 0.00001'f *
          max(abs(actual[i]), abs(expected[i]))
      doAssert abs(actual[i]) < Inf and abs(expected[i]) < Inf and
        error <= tolerance,
        &"{dimensions}D {noise.octaves} octaves {backend} sample {i}: " &
          &"actual={actual[i]}, reference={expected[i]}, error={error}"
    consume(actual)
  consume(expected)

proc measure(
  dimensions: static[int],
  backend: Backend,
  noise: Simplex,
  shape: Shape
): float64 =
  ## Times generation, excluding checksums and output buffer disposal.
  for batch in 0 ..< NoisyBenchmarkBatches:
    let started = getMonoTime()
    let values = generate(dimensions, backend, noise, shape)
    let elapsed = getMonoTime() - started
    result += elapsed.inNanoseconds.float64 / 1_000_000.0
    consume(values)
  result /= NoisyBenchmarkBatches.float64

proc runCase(dimensions: static[int], noise: Simplex, shape: Shape) =
  ## Warms the paths, rotates measurement order, and reports sorted medians.
  validate(dimensions, noise, shape)
  for warmup in 0 ..< 2:
    for backend in Backend:
      consume(generate(dimensions, backend, noise, shape))
  var timings: array[Backend, seq[float64]]
  for round in 0 ..< NoisyBenchmarkRounds:
    for offset in 0 ..< 3:
      let backend = Backend((round + offset) mod 3)
      timings[backend].add measure(dimensions, backend, noise, shape)
  var rawTimings = newJObject()
  for backend in Backend:
    rawTimings[$backend] = %timings[backend]
  var medians: array[Backend, float64]
  for backend in Backend:
    timings[backend].sort()
    let middle = timings[backend].len div 2
    medians[backend] =
      if timings[backend].len mod 2 == 0:
        (timings[backend][middle - 1] + timings[backend][middle]) / 2
      else:
        timings[backend][middle]
  let
    reference = medians[Reference]
    scalar = medians[Scalar]
    vector = medians[Vector]
  results.add %*{
    "dimensions": dimensions,
    "octaves": noise.octaves,
    "shape": shape,
    "round_means_ms": rawTimings,
    "reference_ms": reference,
    "scalar_ms": scalar,
    "simd_ms": vector,
    "simd_vs_reference": reference / vector,
    "simd_vs_scalar": scalar / vector
  }
  echo &"| {dimensions}D | {noise.octaves} | {reference:.3f} | " &
    &"{scalar:.3f} | {vector:.3f} | {reference / scalar:.2f}x | " &
    &"{reference / vector:.2f}x | {scalar / vector:.2f}x |"
  echo &"  ms min..max: reference {timings[Reference][0]:.3f}.." &
    &"{timings[Reference][^1]:.3f}; scalar {timings[Scalar][0]:.3f}.." &
    &"{timings[Scalar][^1]:.3f}; {SimdName} {timings[Vector][0]:.3f}.." &
    &"{timings[Vector][^1]:.3f}"

static:
  doAssert NoisyBenchmarkRounds >= 3
  doAssert NoisyBenchmarkBatches >= 1

var noise = initSimplex(1988)
for i in 0 ..< noise.perm.len:
  noise.perm[i] = referencePermutation(i.int32)
  noise.permMod12[i] = noise.perm[i] mod 12
noise.frequency = 0.075'f
noise.amplitude = 1
noise.lacunarity = 2
noise.gain = 0.5

echo &"Nim {NimVersion}, {hostOS}/{hostCPU}, {CompilerBackend}, " &
  &"release, {SimdName}; 262144 samples per grid"
echo &"{NoisyBenchmarkRounds} rounds, {NoisyBenchmarkBatches} grids per round; " &
  "times are medians of round means, in milliseconds per grid"
echo &"Start={Start}; frequency={noise.frequency}; " &
  &"amplitude={noise.amplitude}; lacunarity={noise.lacunarity}; gain={noise.gain}"
echo &"Shapes (x,y,z,w): {Shapes}"
echo "Includes output allocation and writes; validation/checksums excluded."
echo "Reference and scalar paths fill buffers; references use one batch call."
echo "Speedup = baseline time / measured time; above 1 means faster."
echo "| Dim | Octaves | Reference ms | Scalar ms | " & SimdName &
  " ms | Scalar vs ref | " & SimdName & " vs ref | " & SimdName & " vs scalar |"
echo "| --- | --- | --- | --- | --- | --- | --- | --- |"
for octaves in [1, 3]:
  noise.octaves = octaves
  runCase(1, noise, Shapes[0])
  runCase(2, noise, Shapes[1])
  runCase(3, noise, Shapes[2])
  runCase(4, noise, Shapes[3])
echo &"Consumed checksum: {checksum:.6f}"
if NoisyBenchmarkOutput.len > 0:
  writeFile(NoisyBenchmarkOutput, pretty(%*{
    "nim_version": NimVersion,
    "os": hostOS,
    "cpu": hostCPU,
    "compiler_backend": CompilerBackend,
    "simd_backend": SimdName,
    "rounds": NoisyBenchmarkRounds,
    "batches": NoisyBenchmarkBatches,
    "samples_per_grid": 262144,
    "start": Start,
    "frequency": noise.frequency,
    "amplitude": noise.amplitude,
    "lacunarity": noise.lacunarity,
    "gain": noise.gain,
    "checksum": checksum,
    "cases": results
  }) & "\n")
