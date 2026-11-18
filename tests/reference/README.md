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

Install Noisy's dependencies with Nimby as usual. Reference checks need a C++
compiler. The initial fetch uses curl and requires network access. Run these
commands from the repository root.

```sh
nim r tools/gen_references.nim
nim r -f -d:noisyReferenceDir=tests/reference/.cache tests/test_reference.nim
nim cpp -r -f -d:release -d:noisyReferenceDir=tests/reference/.cache tests/test_reference.nim
```

Quote the whole `-d:noisyReferenceDir=...` argument when the path contains spaces.
The Nim runner also prints the command for its resolved cache path.

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
[Improved Perlin](../../docs/perlin.md) and [FastNoiseLite-compatible noise](../../docs/scalar-noises.md)
guides.

## Timing

Run this from the repository root to download the pinned sources and benchmark
all supported noise families and dimensions with one and three octaves.

```sh
nim r tools/gen_references.nim --benchmark
```

On x86, automatic dispatch selects AVX2 for every noise family when supported.
The fallback uses SSE2 for Simplex and scalar sampling for the other families.
When AVX2 is available, the runner also benchmarks a separate `-d:noisyNoAvx`
fallback build. On ARM64, every noise family uses NEON. The table shows a dash for a noise without
the given SIMD backend. Reference and scalar medians come from the
automatic-dispatch run. Git records the checked-out source
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
nim r tools/gen_references.nim --offline --benchmark
# Short smoke test; use the default settings for reported measurements.
nim r tools/gen_references.nim --benchmark --rounds 3 --batches 1
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
public `grid()` API. The timer includes allocation and output writes.
Validation, checksums, and freeing
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

### Scalar kernel optimization results

Recorded on 2026-11-18 at commit `3bd36ec` with a clean working tree, on the
same Ryzen 9 5900X, Windows, Nim 2.2.10, and GCC/G++ 11.1.0 setup above.
Both measurements use release builds, ORC, **threads enabled**, and the default
goto exception backend. No global danger/checks-off flag is used. Each median
covers 11 rounds of four grids of 262,144 samples, including allocation and
output writes. All timed outputs passed reference validation.

The before column is the earlier recorded `0de468a` run; the after column
uses the audited kernels with checks and exception bookkeeping removed.
These are separate runs on the same machine, rather than interleaved revisions.
Scalar speedups range from **1.17x to 2.56x** across the 36 workloads.
The fixed lookup tables were already constants; these gains do not come from
changing their storage. See [kernel safety](../../docs/kernels.md).

| Workload (dimensions / octaves) | Before ms | After ms | Speedup |
| --- | ---: | ---: | ---: |
| Simplex 1D / 1 | 8.299 | 5.722 | 1.45x |
| Simplex 2D / 1 | 17.066 | 9.122 | 1.87x |
| Simplex 3D / 1 | 23.466 | 10.430 | 2.25x |
| Simplex 4D / 1 | 23.420 | 16.898 | 1.39x |
| Improved Perlin 1D / 1 | 10.249 | 6.409 | 1.60x |
| Improved Perlin 2D / 1 | 13.623 | 7.143 | 1.91x |
| Improved Perlin 3D / 1 | 19.021 | 13.502 | 1.41x |
| Improved Perlin 4D / 1 | 29.743 | 22.778 | 1.31x |
| OpenSimplex2 2D / 1 | 11.323 | 7.160 | 1.58x |
| OpenSimplex2 3D / 1 | 13.828 | 9.558 | 1.45x |
| OpenSimplex2S 2D / 1 | 11.992 | 7.030 | 1.71x |
| OpenSimplex2S 3D / 1 | 15.817 | 10.393 | 1.52x |
| Value 2D / 1 | 10.904 | 6.304 | 1.73x |
| Value 3D / 1 | 12.084 | 7.293 | 1.66x |
| ValueCubic 2D / 1 | 13.805 | 8.930 | 1.55x |
| ValueCubic 3D / 1 | 24.508 | 16.544 | 1.48x |
| Cellular 2D / 1 | 24.346 | 17.672 | 1.38x |
| Cellular 3D / 1 | 60.782 | 49.255 | 1.23x |
| Simplex 1D / 3 | 10.624 | 8.092 | 1.31x |
| Simplex 2D / 3 | 39.420 | 16.100 | 2.45x |
| Simplex 3D / 3 | 60.320 | 23.549 | 2.56x |
| Simplex 4D / 3 | 67.656 | 56.241 | 1.20x |
| Improved Perlin 1D / 3 | 19.728 | 9.002 | 2.19x |
| Improved Perlin 2D / 3 | 32.726 | 17.108 | 1.91x |
| Improved Perlin 3D / 3 | 49.626 | 33.972 | 1.46x |
| Improved Perlin 4D / 3 | 85.769 | 64.702 | 1.33x |
| OpenSimplex2 2D / 3 | 22.684 | 12.243 | 1.85x |
| OpenSimplex2 3D / 3 | 31.704 | 23.901 | 1.33x |
| OpenSimplex2S 2D / 3 | 24.493 | 13.767 | 1.78x |
| OpenSimplex2S 3D / 3 | 39.576 | 25.730 | 1.54x |
| Value 2D / 3 | 20.667 | 8.152 | 2.54x |
| Value 3D / 3 | 23.756 | 12.511 | 1.90x |
| ValueCubic 2D / 3 | 29.123 | 17.305 | 1.68x |
| ValueCubic 3D / 3 | 60.847 | 40.942 | 1.49x |
| Cellular 2D / 3 | 76.598 | 62.147 | 1.23x |
| Cellular 3D / 3 | 196.239 | 167.742 | 1.17x |

[Raw results](results/ryzen-5900x-kernel-optimization.json) include every round
and current C/C++ reference, scalar, and AVX2 measurements, plus the SSE2/scalar
fallback run. The prior results remain in
[the original record](results/ryzen-5900x-all-noises-avx2.json).

### AMD Ryzen 9 5900X expanded AVX2 results

Recorded on 2026-11-18 at commit `0de468a` with a clean working tree, Windows
x86-64, Nim 2.2.10, MinGW-w64 GCC/G++ 11.1.0, the Nim C backend, and `-d:release`.
These are medians from 11 rounds of four grids per path, with 262,144 samples
per grid. Every timed output was checked against the matching C/C++ reference.
Times include allocation and output writes; lower is better.

| Noise | Dimensions | Octaves | C/C++ reference ms | Noisy scalar ms | Noisy AVX2 ms | AVX2 vs reference | AVX2 vs scalar |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Improved Perlin | 1D | 1 | 1.323 | 10.249 | 1.024 | 1.29x | 10.01x |
| Improved Perlin | 2D | 1 | 2.482 | 13.623 | 1.606 | 1.55x | 8.48x |
| Improved Perlin | 3D | 1 | 3.914 | 19.021 | 2.767 | 1.41x | 6.87x |
| Improved Perlin | 4D | 1 | 8.247 | 29.743 | 4.743 | 1.74x | 6.27x |
| OpenSimplex2 | 2D | 1 | 2.982 | 11.323 | 1.513 | 1.97x | 7.48x |
| OpenSimplex2 | 3D | 1 | 6.263 | 13.828 | 2.499 | 2.51x | 5.53x |
| OpenSimplex2S | 2D | 1 | 3.392 | 11.992 | 1.869 | 1.81x | 6.41x |
| OpenSimplex2S | 3D | 1 | 7.044 | 15.817 | 3.976 | 1.77x | 3.98x |
| Value | 2D | 1 | 2.345 | 10.904 | 0.971 | 2.41x | 11.22x |
| Value | 3D | 1 | 4.098 | 12.084 | 1.170 | 3.50x | 10.33x |
| ValueCubic | 2D | 1 | 5.498 | 13.805 | 1.504 | 3.65x | 9.18x |
| ValueCubic | 3D | 1 | 16.534 | 24.508 | 3.682 | 4.49x | 6.66x |
| Cellular | 2D | 1 | 5.100 | 24.346 | 2.574 | 1.98x | 9.46x |
| Cellular | 3D | 1 | 15.661 | 60.782 | 8.252 | 1.90x | 7.37x |
| Improved Perlin | 1D | 3 | 3.295 | 19.728 | 1.445 | 2.28x | 13.66x |
| Improved Perlin | 2D | 3 | 9.870 | 32.726 | 3.373 | 2.93x | 9.70x |
| Improved Perlin | 3D | 3 | 13.557 | 49.626 | 6.357 | 2.13x | 7.81x |
| Improved Perlin | 4D | 3 | 29.931 | 85.769 | 11.975 | 2.50x | 7.16x |
| OpenSimplex2 | 2D | 3 | 8.576 | 22.684 | 3.109 | 2.76x | 7.30x |
| OpenSimplex2 | 3D | 3 | 21.053 | 31.704 | 5.781 | 3.64x | 5.48x |
| OpenSimplex2S | 2D | 3 | 9.867 | 24.493 | 4.127 | 2.39x | 5.93x |
| OpenSimplex2S | 3D | 3 | 23.671 | 39.576 | 10.146 | 2.33x | 3.90x |
| Value | 2D | 3 | 6.416 | 20.667 | 1.366 | 4.70x | 15.13x |
| Value | 3D | 3 | 10.880 | 23.756 | 1.912 | 5.69x | 12.43x |
| ValueCubic | 2D | 3 | 14.640 | 29.123 | 2.802 | 5.22x | 10.39x |
| ValueCubic | 3D | 3 | 47.472 | 60.847 | 8.707 | 5.45x | 6.99x |
| Cellular | 2D | 3 | 14.849 | 76.598 | 6.143 | 2.42x | 12.47x |
| Cellular | 3D | 3 | 46.516 | 196.239 | 23.136 | 2.01x | 8.48x |

For these newly accelerated workloads, AVX2 was 3.90x to 15.13x faster than
Noisy scalar and 1.29x to 5.69x faster than the selected matching references.
Cellular uses squared Euclidean distance, `Distance`, and jitter of 1; its
other modes are covered by correctness tests, not this timing table.

[Raw round means and metadata](results/ryzen-5900x-all-noises-avx2.json) also
include Simplex and a separate `-d:noisyNoAvx` run. Only Simplex has an SSE2
backend; the remaining families fall back to scalar on that build. These
numbers describe this machine and these grid workloads, rather than arbitrary
point queries or every noise library.

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
