import std/os

const
  ReferenceDir = currentSourcePath().parentDir
  NoisyReferenceDir {.strdefine.} = ""

when NoisyReferenceDir.len == 0:
  {.error: "Fetch references with tools/gen_references.py, then pass " &
    "-d:noisyReferenceDir=/cache/path".}
static:
  doAssert fileExists(NoisyReferenceDir / "fastnoise/FastNoiseLite.h"),
    "FastNoiseLite reference missing; rerun tools/gen_references.py"

# FastNoiseLite intentionally wraps signed hashes. Define that behavior only
# for this unchanged upstream translation unit, leaving other builds alone.
{.compile(ReferenceDir / "fastnoise.cpp",
  "-fwrapv -I" & quoteShell(NoisyReferenceDir)).}
{.push cdecl, importc, header: ReferenceDir / "fastnoise.h".}
proc referenceFastNoise2*(
  algorithm, seed: int32,
  x, y: float32,
  distance, output: int32,
  jitter: float32
): float32
  ## Samples the unchanged public FastNoiseLite 2D API.
proc referenceFastNoise3*(
  algorithm, seed: int32,
  x, y, z: float32,
  distance, output: int32,
  jitter: float32
): float32
  ## Samples the unchanged public FastNoiseLite 3D API.
proc referenceFastNoiseGrid*(
  output: ptr float32,
  noiseType, seed, dimensions, width, height, depth: int32,
  start: ptr float32,
  octaves: int32,
  amplitude, frequency, lacunarity, gain: float32
)
  ## Fills a grid through the public C++ API with default Cellular settings.
{.pop.}
