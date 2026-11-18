import std/os

const
  ReferenceDir = currentSourcePath().parentDir
  NoisyReferenceDir* {.strdefine.} = ""

when NoisyReferenceDir.len == 0:
  {.error: "Fetch references with tools/gen_references.py, then pass " &
    "-d:noisyReferenceDir=/cache/path".}
static:
  doAssert fileExists(NoisyReferenceDir / "stegu/noise1234.c"),
    "Perlin reference missing; rerun tools/gen_references.py with --cache-dir"
  doAssert fileExists(NoisyReferenceDir / "caseman/simplex23.c"),
    "Reference cache missing; run tools/gen_references.py with --cache-dir"

{.passC: "-I" & quoteShell(NoisyReferenceDir).}
{.compile: ReferenceDir / "reference_noise.cpp".}
{.push cdecl, importc, header: ReferenceDir / "reference_noise.h".}
proc referenceNoise1*(x: float32): float32
  ## Samples the unchanged Rombauts 1D kernel.
proc referenceNoise2*(x, y: float32): float32
  ## Samples the unchanged Duncan 2D kernel.
proc referenceNoise3*(x, y, z: float32): float32
  ## Samples the unchanged Duncan 3D kernel.
proc referenceNoise4*(x, y, z, w: float32): float32
  ## Samples the unchanged Gustavson 4D kernel.
proc referencePerlin1*(x: float32): float32
  ## Samples Gustavson's unchanged improved Perlin 1D kernel.
proc referencePerlin2*(x, y: float32): float32
  ## Samples Gustavson's unchanged improved Perlin 2D kernel.
proc referencePerlin3*(x, y, z: float32): float32
  ## Samples Gustavson's unchanged improved Perlin 3D kernel.
proc referencePerlin4*(x, y, z, w: float32): float32
  ## Samples Gustavson's unchanged improved Perlin 4D kernel.
proc referencePerlinPermutation*(index: int32): int32
  ## Reads the published permutation from the Perlin reference itself.
proc referencePermutation*(index: int32): int32
  ## Reads the original permutation shared by the references.
proc referenceGrid*(
  output: ptr float32,
  dimensions, width, height, depth, frames: int32,
  start: ptr float32,
  octaves: int32,
  amplitude, frequency, lacunarity, gain: float32
)
  ## Fills caller-owned storage using one C++ call and Noisy's octave rules.
proc referencePerlinGrid*(
  output: ptr float32,
  dimensions, width, height, depth, frames: int32,
  start: ptr float32,
  octaves: int32,
  amplitude, frequency, lacunarity, gain: float32
)
  ## Fills a grid using the C reference's kernels and Noisy's octave rules.
{.pop.}
