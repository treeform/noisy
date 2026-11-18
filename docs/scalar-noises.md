# OpenSimplex2, Value, and Cellular noise

Noisy includes native scalar 2D and 3D implementations of these three families,
matching the C++ implementation in
[FastNoiseLite](https://github.com/Auburn/FastNoiseLite/tree/785f37a9ad76e283586a379675085f2063ae03f7).
These types have integer and float32 `value()` overloads. They have no `grid()`
or SIMD implementations, and no 1D/4D overloads.

| Type | Constructor | Algorithm |
| --- | --- | --- |
| `OpenSimplex2` | `initOpenSimplex2(seed)` | FastNoiseLite OpenSimplex2 with its default coordinate transforms |
| `ValueNoise` | `initValue(seed)` | Hashed lattice values with cubic Hermite interpolation |
| `Cellular` | `initCellular(seed)` | Nearest-feature noise with configurable distance and return modes |

```nim
import noisy

var terrain = initOpenSimplex2(1988)
terrain.frequency = 0.01
terrain.octaves = 3
let height = terrain.value(12.5'f, -7.25'f)

var variation = initValue(42)
variation.frequency = 0.1
let density = variation.value(12.5'f, -7.25'f, 3.0'f)

var cells = initCellular(42)
cells.distanceFunction = Euclidean
cells.returnType = Distance2Sub
cells.jitter = 0.75
let edge = cells.value(0.25'f, 0.75'f)
```

The standalone modules `noisy/opensimplexes`, `noisy/values`, and
`noisy/cellulars` expose the respective types without importing the Simplex
SIMD backend. The ports retain the FastNoiseLite license and attribution in
[`LICENSE.fastnoiselite.txt`](../LICENSE.fastnoiselite.txt).

## Settings and seed compatibility

All three types expose `seed` and the existing Noisy octave controls:

| Setting | Default |
| --- | ---: |
| `octaves` | 1 |
| `amplitude` | 1 |
| `frequency` | 1 |
| `lacunarity` | 2 |
| `gain` | 0.5 |

`seed` accepts an integer and uses its low 32 bits for hashing. Zero and negative
seeds are supported. The implementation explicitly wraps hash arithmetic so
hash results agree between Nim debug and release builds. The same settings and seed match the
selected FastNoiseLite single-octave kernel. Distinct seeds are not guaranteed
to produce distinct maps: for example, Value noise preserves FastNoiseLite's
collision between seed 0 and the minimum signed 32-bit seed.

Octaves follow Noisy's convention: keep the same seed, scale coordinates by
frequency, multiply each sample by amplitude, then divide the sum by the octave
count. Frequency and amplitude are multiplied by lacunarity and gain after each
layer. FastNoiseLite's own fractal API has different seed progression and
normalization, so the reference tests compose its single-octave samples using
Noisy's rules. Noisy's default frequency is 1 rather than FastNoiseLite's 0.01.

Nonpositive octave counts raise `NoisyError`. Use finite coordinates and
settings, keeping scaled and transformed lattice coordinates comfortably
inside the signed 32-bit range. Float32 loses fractional resolution at large
coordinates. OpenSimplex2 uses FastNoiseLite's default 2D skew and 3D rotation;
this API does not add OpenSimplex2S or alternate plane rotations. Value selects
the Hermite-interpolated algorithm, not the separate ValueCubic mode.

## Cellular settings

`distanceFunction` accepts `Euclidean`, `EuclideanSquared`, `Manhattan`, or
`Hybrid` (Manhattan plus squared Euclidean distance). The default is
`EuclideanSquared`. `jitter` defaults to 1 and accepts values from 0 through 1;
invalid values raise `NoisyError`. Zero places feature points at lattice sites.

`returnType` defaults to `Distance`. Let F1 and F2 be the nearest and second
nearest feature distances under the selected metric. Returned samples follow
the reference's signed conventions:

| Return type | Single-octave result |
| --- | --- |
| `CellValue` | Signed hash value of the nearest feature |
| `Distance` | F1 - 1 |
| `Distance2` | F2 - 1 |
| `Distance2Add` | (F1 + F2) / 2 - 1 |
| `Distance2Sub` | F2 - F1 - 1 |
| `Distance2Mul` | F1 * F2 / 2 - 1 |
| `Distance2Div` | F1 / F2 - 1 |

Cell values can jump at Voronoi boundaries. Distance modes and custom octave
settings are not promised to stay in [-1, 1]. These are noise samples, not an
API that returns feature positions or an exact distance to a Voronoi edge.

## Reproduce the comparison

The upstream C++ header and license are pinned by commit and SHA-256 in
`tests/reference/sources.json` and downloaded outside this repository.
Only the native Nim implementations, their necessary tables, attribution, and
our test bridge are committed.

```sh
python3 tools/gen_references.py --cache-dir /tmp/noisy-references
nim r -f -d:noisyReferenceDir=/tmp/noisy-references tests/test_fastnoise_reference.nim
nim cpp -r -f -d:release -d:noisyReferenceDir=/tmp/noisy-references tests/test_fastnoise_reference.nim
```

The unchanged C++ reference is compiled with `-fwrapv` for its signed hash
arithmetic; this flag applies only to that reference translation unit. Tests
compare the public APIs with matched seeds, coordinates, and settings, including
negative cells, lattice boundaries, octave combinations, and all Cellular
metrics and return modes. The normal offline suite contains recorded outputs
from the compiled reference as well as behavior tests:

```sh
nim check tests/tests.nim
nim r tests/tests.nim
```


## Recorded validation

With Nim 2.4.0 and Apple Clang, 141,912 comparisons per build passed in native
ARM64 C debug and C++ release, and x86-64 C++ release under Rosetta. The maximum
errors across these builds were:

| Noise | Comparisons per build | Maximum absolute error |
| --- | ---: | ---: |
| OpenSimplex2 | 24,120 | 6.56e-7 |
| Value | 24,120 | 0 |
| Cellular | 93,672 | 4.77e-7 |

Tests use `1e-5 + 1e-5 * max(abs(actual), abs(expected))` for float32 rounding.
Single-octave `CellValue` samples additionally require exact equality, so
selecting a different feature cell cannot pass through a numeric tolerance.
The ordinary offline suite includes 158 recorded C++ outputs. The new noise
tests and reference suite also pass with Nim 1.4.8 on ARM64. Existing API tests
and the complete Nim 2.4 offline suite pass on ARM64 and x86-64.
