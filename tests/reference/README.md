# Independent simplex references

These tests compare Noisy's scalar and SIMD output directly against established
C and C++ implementations. This directory contains our wrapper and source
metadata only. Third-party source code is downloaded into an external cache and
is not committed to this repository.

| Dimensions | Upstream source | Pinned revision | License |
| --- | --- | --- | --- |
| 1D | [SRombauts/SimplexNoise](https://github.com/SRombauts/SimplexNoise/blob/175add2d2a9ba99a43320b66bb87bb21b879b0f5/src/SimplexNoise.cpp) | `175add2d2a9ba99a43320b66bb87bb21b879b0f5` | [MIT](https://github.com/SRombauts/SimplexNoise/blob/175add2d2a9ba99a43320b66bb87bb21b879b0f5/LICENSE.txt) |
| 2D, 3D | [caseman/noise](https://github.com/caseman/noise/blob/bb32991ab97e90882d0e46e578060717c5b90dc5/_simplex.c) | `bb32991ab97e90882d0e46e578060717c5b90dc5` | [MIT](https://github.com/caseman/noise/blob/bb32991ab97e90882d0e46e578060717c5b90dc5/LICENSE.txt) |
| 4D | [stegu/perlin-noise](https://github.com/stegu/perlin-noise/blob/a624f5a18395c0ebc04738191efddf1dbd2cad9d/src/simplexnoise1234.c) | `a624f5a18395c0ebc04738191efddf1dbd2cad9d` | Public domain declaration in the source and header |

## Reproduce the comparison

Install Noisy's dependencies with Nimby as usual. The reference check additionally
needs Python 3.8 or newer for downloading, a C++ compiler, and network access for
the initial fetch. From the repository root on macOS or Linux:

```sh
python3 tools/gen_references.py --cache-dir /tmp/noisy-references
nim r -f -d:noisyReferenceDir=/tmp/noisy-references tests/test_reference.nim
nim cpp -r -f -d:release -d:noisyReferenceDir=/tmp/noisy-references tests/test_reference.nim
```

On Windows, use `python` and an external cache such as `C:/Temp/noisy-references`.
Quote the whole `-d:noisyReferenceDir=...` argument when the path contains spaces.
The fetch script also prints the exact command for its resolved cache path.

The fetch verifies all SHA-256 hashes from `sources.json`, including already
cached files. A valid cache can be reused offline. `-f` forces recompilation so
changes to external include files cannot leave stale compiled reference code.
To include these comparisons in the complete suite, pass the same definition to
`nim r -f tests/tests.nim`. Without the definition, the normal test suite runs
offline and reports that external comparisons were skipped. CI explicitly
fetches and runs the C and C++ reference comparisons on all three platforms.

The test prints comparison counts and maximum absolute errors for each dimension.
It is a numerical correctness comparison, not a performance timing benchmark.

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
