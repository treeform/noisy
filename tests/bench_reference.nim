import
  std/[algorithm, json, math, monotimes, strformat, strutils, times],
  noisy,
  noisy/grids,
  reference/[bindings, fastnoise]

when not defined(release):
  {.error: "Run benchmarks with -d:release.".}

const
  NoisyBenchmarkRounds {.intdefine.} = 11
  NoisyBenchmarkBatches {.intdefine.} = 4
  NoisyBenchmarkOutput {.strdefine.} = ""
  Seed = 1988
  Frequency = 0.075'f
  Amplitude = 1.0'f
  Lacunarity = 2.0'f
  Gain = 0.5'f
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

template supportsVector(noise: typed): bool =
  when useGridAvx2:
    noise is Simplex or hasAvx2()
  else:
    noise is Simplex or useGridNeon

proc consume(values: seq[float32]) =
  ## Consumes every output outside the timer to prevent dead-code removal.
  var total: float64
  for value in values:
    total += value.float64
  checksum += total

proc generate(
  dimensions: static[int],
  backend: Backend,
  noise: auto,
  shape: Shape
): seq[float32] =
  ## Allocates and fills the same output layout through each public path.
  case backend
  of Reference:
    result = newSeq[float32](shape[0] * shape[1] * shape[2] * shape[3])
    when noise is Simplex or noise is Perlin:
      template fill(reference: untyped) =
        reference(
          result[0].addr,
          dimensions.int32,
          shape[0].int32, shape[1].int32, shape[2].int32, shape[3].int32,
          Start[0].unsafeAddr,
          noise.octaves.int32,
          noise.amplitude, noise.frequency, noise.lacunarity, noise.gain
        )
      when noise is Simplex:
        fill(referenceGrid)
      else:
        fill(referencePerlinGrid)
    else:
      const algorithm =
        when noise is OpenSimplex2: 0'i32
        elif noise is ValueNoise: 1'i32
        elif noise is Cellular: 2'i32
        elif noise is OpenSimplex2S: 3'i32
        elif noise is ValueCubic: 4'i32
        else: {.error: "Missing FastNoiseLite reference mapping".}
      referenceFastNoiseGrid(
        result[0].addr,
        algorithm, cast[int32](noise.seed), dimensions.int32,
        shape[0].int32, shape[1].int32, shape[2].int32,
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
    if not supportsVector(noise):
      raise newException(AssertionDefect, "No SIMD grid for this noise")
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

proc validate(
  dimensions: static[int],
  name: string,
  noise: auto,
  shape: Shape,
  backends: openArray[Backend]
) =
  ## Checks every timed input against the reference before benchmarking.
  let expected = generate(dimensions, Reference, noise, shape)
  const absoluteTolerance = when noise is Simplex: 0.0001'f else: 0.00001'f
  for backend in backends:
    let actual = generate(dimensions, backend, noise, shape)
    doAssert actual.len == expected.len
    for i in 0 ..< actual.len:
      let
        error = abs(actual[i] - expected[i])
        tolerance = absoluteTolerance + 0.00001'f *
          max(abs(actual[i]), abs(expected[i]))
      doAssert abs(actual[i]) < Inf and abs(expected[i]) < Inf and
        error <= tolerance,
        &"{name} {dimensions}D {noise.octaves} octaves {backend} sample {i}: " &
          &"actual={actual[i]}, reference={expected[i]}, error={error}"
    consume(actual)
  consume(expected)

proc measure(
  dimensions: static[int],
  backend: Backend,
  noise: auto,
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

proc runCase(dimensions: static[int], name: string, noise: auto) =
  ## Warms the paths, rotates measurement order, and reports sorted medians.
  const shape = Shapes[dimensions - 1]
  let backends =
    if supportsVector(noise): @[Reference, Scalar, Vector]
    else: @[Reference, Scalar]
  validate(dimensions, name, noise, shape, backends.toOpenArray(1, backends.high))
  for warmup in 0 ..< 2:
    for backend in backends:
      consume(generate(dimensions, backend, noise, shape))
  var timings: array[Backend, seq[float64]]
  for round in 0 ..< NoisyBenchmarkRounds:
    for offset in 0 ..< backends.len:
      let backend = backends[(round + offset) mod backends.len]
      timings[backend].add measure(dimensions, backend, noise, shape)
  var rawTimings = newJObject()
  var medians: array[Backend, float64]
  var ranges: seq[string]
  for backend in backends:
    rawTimings[$backend] = %timings[backend]
    timings[backend].sort()
    let middle = timings[backend].len div 2
    medians[backend] =
      if timings[backend].len mod 2 == 0:
        (timings[backend][middle - 1] + timings[backend][middle]) / 2
      else:
        timings[backend][middle]
    ranges.add &"{backend} {timings[backend][0]:.3f}..{timings[backend][^1]:.3f}"
  let
    reference = medians[Reference]
    scalar = medians[Scalar]
  var row = %*{
    "noise": name,
    "dimensions": dimensions,
    "octaves": noise.octaves,
    "shape": shape,
    "round_means_ms": rawTimings,
    "reference_ms": reference,
    "scalar_ms": scalar
  }
  let line = &"| {name} | {dimensions}D | {noise.octaves} | " &
    &"{reference:.3f} | {scalar:.3f} | "
  if supportsVector(noise):
    let vector = medians[Vector]
    row["simd_ms"] = %vector
    row["simd_vs_reference"] = %(reference / vector)
    row["simd_vs_scalar"] = %(scalar / vector)
    echo line & &"{vector:.3f} | {reference / scalar:.2f}x | " &
      &"{reference / vector:.2f}x | {scalar / vector:.2f}x |"
  else:
    row["simd_ms"] = newJNull()
    row["simd_vs_reference"] = newJNull()
    row["simd_vs_scalar"] = newJNull()
    echo line & &"— | {reference / scalar:.2f}x | — | — |"
  echo "  ms min..max: " & ranges.join(", ")
  results.add row

proc runNoise(name: string, noise: auto, octaves: int) =
  ## Runs each dimension supported by this noise with the shared settings.
  var configured = noise
  configured.octaves = octaves
  configured.frequency = Frequency
  configured.amplitude = Amplitude
  configured.lacunarity = Lacunarity
  configured.gain = Gain
  when noise is Simplex or noise is Perlin:
    runCase(1, name, configured)
  runCase(2, name, configured)
  runCase(3, name, configured)
  when noise is Simplex or noise is Perlin:
    runCase(4, name, configured)

static:
  doAssert NoisyBenchmarkRounds >= 3
  doAssert NoisyBenchmarkBatches >= 1

var
  simplex = initSimplex(Seed)
  perlin = initPerlin(Seed)
for i in 0 ..< simplex.perm.len:
  simplex.perm[i] = referencePermutation(i.int32)
  simplex.permMod12[i] = simplex.perm[i] mod 12
for i in 0 ..< perlin.perm.len:
  perlin.perm[i] = referencePerlinPermutation(i.int32)

echo &"Nim {NimVersion}, {hostOS}/{hostCPU}, {CompilerBackend}, " &
  &"release, {SimdName}; 262144 samples per grid"
echo &"{NoisyBenchmarkRounds} rounds, {NoisyBenchmarkBatches} grids per round; " &
  "times are medians of round means, in milliseconds per grid"
echo &"Start={Start}; frequency={Frequency}; " &
  &"amplitude={Amplitude}; lacunarity={Lacunarity}; gain={Gain}"
echo &"Shapes (x,y,z,w): {Shapes}"
echo "Includes output allocation and writes; validation/checksums excluded."
echo "Reference and scalar paths fill buffers; references use one batch call."
echo "Speedup = baseline time / measured time; above 1 means faster."
echo "| Noise | Dim | Octaves | Reference ms | Scalar ms | " & SimdName &
  " ms | Scalar vs ref | " & SimdName & " vs ref | " & SimdName & " vs scalar |"
echo "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"
for octaves in [1, 3]:
  runNoise("Simplex", simplex, octaves)
  runNoise("Improved Perlin", perlin, octaves)
  runNoise("OpenSimplex2", initOpenSimplex2(Seed), octaves)
  runNoise("OpenSimplex2S", initOpenSimplex2S(Seed), octaves)
  runNoise("Value", initValue(Seed), octaves)
  runNoise("ValueCubic", initValueCubic(Seed), octaves)
  runNoise("Cellular", initCellular(Seed), octaves)
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
    "seed": Seed,
    "frequency": Frequency,
    "amplitude": Amplitude,
    "lacunarity": Lacunarity,
    "gain": Gain,
    "checksum": checksum,
    "cases": results
  }) & "\n")
