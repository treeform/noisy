# Independent noise references

These tests compare Noisy's scalar and SIMD output directly against established
C and C++ implementations. This directory contains our wrapper and source
metadata only. Third-party source code is downloaded into a verified cache and
is not committed to this repository. By default, the cache is the ignored
`tests/reference/.cache/` directory inside the checkout; `--cache-dir` also
accepts an external location.

| Noise | Dimensions | Upstream source | Pinned revision | License |
| --- | --- | --- | --- | --- |
| Simplex | 1D | [SRombauts/SimplexNoise](https://github.com/SRombauts/SimplexNoise/blob/175add2d2a9ba99a43320b66bb87bb21b879b0f5/src/SimplexNoise.cpp) | `175add2d2a9ba99a43320b66bb87bb21b879b0f5` | [MIT](https://github.com/SRombauts/SimplexNoise/blob/175add2d2a9ba99a43320b66bb87bb21b879b0f5/LICENSE.txt) |
| Simplex | 2D, 3D | [caseman/noise](https://github.com/caseman/noise/blob/bb32991ab97e90882d0e46e578060717c5b90dc5/_simplex.c) | `bb32991ab97e90882d0e46e578060717c5b90dc5` | [MIT](https://github.com/caseman/noise/blob/bb32991ab97e90882d0e46e578060717c5b90dc5/LICENSE.txt) |
| Simplex | 4D | [stegu/perlin-noise](https://github.com/stegu/perlin-noise/blob/a624f5a18395c0ebc04738191efddf1dbd2cad9d/src/simplexnoise1234.c) | `a624f5a18395c0ebc04738191efddf1dbd2cad9d` | Public domain declaration in the source and header |
| Improved Perlin | 1D–4D | [stegu/perlin-noise](https://github.com/stegu/perlin-noise/blob/a624f5a18395c0ebc04738191efddf1dbd2cad9d/src/noise1234.c) | `a624f5a18395c0ebc04738191efddf1dbd2cad9d` | Public domain declaration in the source and header |
| OpenSimplex2, OpenSimplex2S, Value, ValueCubic, Cellular | 2D, 3D | [FastNoiseLite](https://github.com/Auburn/FastNoiseLite/blob/785f37a9ad76e283586a379675085f2063ae03f7/Cpp/FastNoiseLite.h) | `785f37a9ad76e283586a379675085f2063ae03f7` | [MIT](https://github.com/Auburn/FastNoiseLite/blob/785f37a9ad76e283586a379675085f2063ae03f7/LICENSE) |

## Reproduce the comparison

Install Noisy's dependencies with Nimby as usual. The reference check additionally
needs Python 3.8 or newer for downloading, a C++ compiler, and network access for
the initial fetch. From the repository root on macOS or Linux:

```sh
python tools/gen_references.py
nim r -f -d:noisyReferenceDir=tests/reference/.cache tests/test_reference.nim
nim cpp -r -f -d:release -d:noisyReferenceDir=tests/reference/.cache tests/test_reference.nim
```

Use `python3` if that is the name of your Python 3 executable. Quote the whole
`-d:noisyReferenceDir=...` argument when the path contains spaces. The fetch
script also prints the exact command for its resolved cache path.

The fetch verifies all SHA-256 hashes from `sources.json`, including already
cached files. `--offline` verifies a cache without making network requests and
fails if any pinned file is missing or altered. `-f` forces recompilation so
changes to external include files cannot leave stale compiled reference code.
To include these comparisons in the complete suite, pass the same definition to
`nim r -f tests/tests.nim`. Without the definition, the normal test suite runs
offline and reports that external comparisons were skipped. CI explicitly
fetches and runs the C and C++ reference comparisons on all three platforms.

The test prints comparison counts and maximum absolute errors for each Simplex
dimension. The other tests and their settings are described in the
[Improved Perlin](../../docs/perlin.md) and [scalar noise](../../docs/scalar-noises.md)
guides.

## Timing

Run this from the repository root to download the pinned sources and benchmark
all supported noise families and dimensions with one and three octaves.

```sh
python tools/gen_references.py --benchmark
```

On x86, automatic dispatch selects AVX2 or SSE2 for Simplex. When it selects
AVX2, the runner also benchmarks a separate `-d:noisyNoAvx` SSE2 build. On
ARM64, Simplex, OpenSimplex2S, and ValueCubic use NEON. The table shows a dash
for a noise without the given SIMD backend. Its reference and scalar medians
come from the automatic-dispatch run. Git records the checked-out source
revision in the summary.

All generated files are ignored under `tests/reference/.cache/`. The `results/`
subdirectory contains `auto.json`, `sse2.json` (when both x86 backends run),
compiler and run logs, `summary.json` with commands and source revision, and
`table.md`. The JSON includes the individual round means and their medians.
The logs include the C/C++ compiler commands. Each run overwrites the output
files for the modes it runs.

Useful variants:

```sh
# Reuse verified downloads without network access.
python tools/gen_references.py --offline --benchmark
# Short smoke test; use the default settings for reported measurements.
python tools/gen_references.py --benchmark --rounds 3 --batches 1
# Run one configuration directly, optionally adding -d:noisyNoAvx.
nim r -f -d:release -d:noisyReferenceDir=tests/reference/.cache tests/bench_reference.nim
```

`bench_reference.nim` uses Nim's standard library. It checks every output against
the reference before timing. Each path allocates and fills a buffer of 262,144
samples with matched coordinates, seeds or permutations, and octave settings.
Shapes are 262144 in 1D, 512×512 in 2D, 64×64×64 in 3D, and 32×32×32×8 in 4D.
Cellular uses the default squared Euclidean distance, `Distance` return mode,
and jitter of 1.

Each reference grid uses one call into the C++ bridge. The bridge applies
Noisy's octave composition and uses FastNoiseLite's public `GetNoise()` API
for its five noise types, with one configuration per grid. Noisy scalar
measurements fill a buffer with `value()` calls. SIMD measurements call the
public `grid()` API. The timer
includes allocation and output writes. Validation, checksums, and freeing
the output buffer happen outside the timer. A printed checksum consumes
the output and prevents dead-code elimination.

After two warmups per path, the runner rotates measurement order across 11 rounds,
each averaging four grids per path. It reports the sorted median of those round
means and their min/max range. `-d:noisyBenchmarkRounds=N` and
`-d:noisyBenchmarkBatches=N` change these settings. A speedup is baseline time
divided by the measured path's time, so values above 1 mean faster. These are
batch grid throughput measurements; arbitrary point queries can behave differently.
Compiler optimizations of the C/C++ references remain enabled, and the 4D
reference retains its double-precision intermediates.

### AMD Ryzen 9 5900X Simplex results

Recorded on Windows x86-64 with Nim 2.2.10, MinGW-w64 GCC/G++ 11.1.0,
the Nim C backend, and `-d:release`. These use the default 11 rounds of four
grids per path, with 262,144 samples per grid. Both SIMD configurations passed
the full reference validation before timing. The full suite also passes with
automatic AVX2 and forced SSE2.

| Dimensions | Octaves | C/C++ reference ms | Noisy scalar ms | Noisy SSE2 ms | Noisy AVX2 ms | AVX2 vs reference |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1D | 1 | 1.299 | 8.681 | 3.421 | 1.157 | 1.12x |
| 2D | 1 | 4.414 | 18.551 | 2.965 | 2.054 | 2.15x |
| 3D | 1 | 7.324 | 24.837 | 4.144 | 3.316 | 2.21x |
| 4D | 1 | 13.017 | 24.217 | 7.542 | 6.937 | 1.88x |
| 1D | 3 | 3.267 | 11.069 | 8.436 | 1.518 | 2.15x |
| 2D | 3 | 12.618 | 41.304 | 7.810 | 4.778 | 2.64x |
| 3D | 3 | 22.614 | 62.805 | 11.496 | 8.403 | 2.69x |
| 4D | 3 | 58.749 | 68.583 | 19.028 | 18.640 | 3.15x |

AVX2 was 1.12x to 3.15x faster than these selected C/C++ references. Its gain
over SSE2 varied substantially: 1D was about 3x to 5.6x faster, while the 4D
three-octave case was roughly tied (a measured advantage of about 2%).
Noisy scalar was slower than the
references in every case on this machine. These results compare the pinned,
matching Simplex implementations above, not all noise libraries or Perlin.

[Raw round means and metadata](results/ryzen-5900x.json) include both runs'
reference and scalar baselines, which help assess variation between runs.
The source commit is recorded there; subsequent benchmark/documentation changes
do not change the measured noise kernels.

### Apple M4 Pro Simplex results

Recorded with Nim 2.4.0, Apple Clang 21.0.0, the Nim C backend, `-d:release`, and
native ARM64 NEON. Times below are milliseconds per grid from the first complete
run; a second run confirmed the same faster/slower conclusions. Timing varies
with hardware, compiler, and system load.

| Dimensions | Octaves | C/C++ reference ms | Noisy scalar ms | Noisy NEON ms | Scalar vs reference | NEON vs reference |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1D | 1 | 1.001 | 0.595 | 0.377 | 1.68x | 2.65x |
| 2D | 1 | 1.312 | 1.800 | 0.968 | 0.73x | 1.36x |
| 3D | 1 | 2.388 | 2.940 | 2.115 | 0.81x | 1.13x |
| 4D | 1 | 7.175 | 6.577 | 3.059 | 1.09x | 2.35x |
| 1D | 3 | 2.910 | 1.462 | 0.913 | 1.99x | 3.19x |
| 2D | 3 | 3.800 | 4.670 | 2.637 | 0.81x | 1.44x |
| 3D | 3 | 9.471 | 8.476 | 6.126 | 1.12x | 1.55x |
| 4D | 3 | 31.176 | 32.003 | 8.807 | 0.97x | 3.54x |

NEON was faster than each selected reference in these workloads. Scalar results
were mixed: notably, scalar 2D was slower at both octave counts. These numbers
compare the specific implementations and workloads above, not every noise library.

## Comparable algorithms and inputs

Simplex implementations have different gradient and normalization conventions.
The selected kernels match Noisy's conventions without changing either side:
1D uses scale 0.395, 2D and 3D use the 12-entry gradient table with scales 70 and
32, and 4D uses support 0.5 and scale 62. Other kernels from these libraries are
not interchangeable. For example, Gustavson's C 1D scale is 0.25 and its 2D/3D
gradient selection differs from Noisy's.

All three sources contain the same original Perlin permutation table. The tests
install that table into Noisy instead of modifying the references. Seed numbers
from different libraries need not produce matching permutations. The reference
functions evaluate one octave; tests compose those samples using Noisy's
amplitude, frequency, and octave rules.

The core comparisons allow `1e-4 + 1e-5 * max(abs(actual), abs(expected))` and
cover random coordinates, octave settings, tied ranks, negative cells, every
SIMD lane, and scalar tails. Dense random inputs lie in [-16, 16]; explicit
boundary cases also exercise permutation wrapping near +/-256.

Gustavson's 4D C kernel uses double-precision skew constants, whereas Noisy uses
float32 arithmetic. Their rounding differences grow with coordinate magnitude.
Larger 4D boundary grids and a separate random stress check near +/-256 therefore
add `8 * float32_epsilon * max(abs(coordinate))` to the tolerance and report their
errors separately. This is a conservative regression allowance for that tested
range, not a guarantee for arbitrary coordinates. No reference constants are
changed and the normal scalar/SIMD parity tests retain their existing tolerance.

## Recorded validation

With Nim 2.4.0 and Apple Clang, the C debug and C++ release checks passed on
native ARM64 NEON and x86-64 SSE2 via Rosetta. The suite made 130,114 direct
reference comparisons, including 23,768 SIMD samples. Maximum absolute errors
across those builds were:

| Comparison | Maximum error |
| --- | ---: |
| 1D core, scalar and grid | 0.000000418 |
| 2D core, scalar and grid | 0.000000597 |
| 3D core, scalar and grid | 0.00000120 |
| 4D core, scalar and grid | 0.00000567 |
| Separate larger-coordinate 4D cases | 0.000338 |

These results validate the selected conventions and tested input ranges. They
do not assert identical seed mappings, bitwise equality, or arbitrary-coordinate
error bounds.

## Source integrity and extraction

`sources.json` records exact upstream URLs and SHA-256 hashes. The fetch script
keeps the complete original sources and licenses in the external cache. Only
Casey Duncan's Python extension needs extraction for standalone compilation:

- Copy the prefix before `#include "Python.h"` (lines 1-4).
- Append the bytes from `// 2D simplex skew factors` up to, but not including,
  `static inline float` followed by `fbm_noise3` (lines 10-128).
- Verify the resulting `caseman/simplex23.c` against its separately pinned hash.

This preserves the 2D/3D kernels verbatim while excluding Python bindings and
unrelated functions. Our `reference_noise.cpp` places the fetched sources in
separate namespaces and exports a C ABI. It does not alter gradient tables,
permutation tables, numeric constants, normalization, or kernel arithmetic.
