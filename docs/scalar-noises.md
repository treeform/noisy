# FastNoiseLite-compatible noise families

Noisy includes native scalar 2D and 3D implementations of these five types,
matching the C++ implementation in
[FastNoiseLite](https://github.com/Auburn/FastNoiseLite/tree/785f37a9ad76e283586a379675085f2063ae03f7).
Use the integer or float32 `value()` overloads for point sampling.

| Type | Constructor | Algorithm |
| --- | --- | --- |
| `OpenSimplex2` | `initOpenSimplex2(seed)` | FastNoiseLite OpenSimplex2 with its default coordinate transforms |
| `OpenSimplex2S` | `initOpenSimplex2S(seed)` | FastNoiseLite OpenSimplex2S with its default coordinate transforms |
| `ValueNoise` | `initValue(seed)` | Hashed lattice values with cubic Hermite interpolation |
| `ValueCubic` | `initValueCubic(seed)` | Cubic interpolation over four lattice samples per axis |
| `Cellular` | `initCellular(seed)` | Nearest-feature noise with configurable distance and return modes |

```nim
import noisy

var terrain = initOpenSimplex2(1988)
terrain.frequency = 0.01
terrain.octaves = 3
let height = terrain.value(12.5'f, -7.25'f)

var smooth = initOpenSimplex2S(1988)
smooth.frequency = 0.01
let smoothHeight = smooth.value(12.5'f, -7.25'f)

var variation = initValue(42)
variation.frequency = 0.1
let density = variation.value(12.5'f, -7.25'f, 3.0'f)

let cubic = initValueCubic(42)
let cubicDensity = cubic.value(1.25'f, -0.75'f, 0.5'f)

var cells = initCellular(42)
cells.distanceFunction = Euclidean
cells.returnType = Distance2Sub
cells.jitter = 0.75
let edge = cells.value(0.25'f, 0.75'f)
```

The standalone modules `noisy/opensimplexes`, `noisy/values`, and
`noisy/cellulars` expose their respective types. The ports retain the
FastNoiseLite license and attribution in
[`LICENSE.fastnoiselite.txt`](../LICENSE.fastnoiselite.txt).

## Grids

All five types support `grid(start, dimens)` in 2D and 3D. Pass integer or
float32 coordinates in `start` and integer sizes in `dimens`. Samples are one
unit apart before frequency scaling. Index the result with zero-based offsets
from `start`, using the same layout as a Simplex grid.

On ARM64, the grids automatically use NEON for groups of four contiguous
samples. Remaining samples and grids on other architectures use scalar
sampling. Build with `-d:noisyNoNeon` to select scalar grids for these types on
ARM64.

```nim
import noisy

let surface = initOpenSimplex2S(1988).grid((-2.5'f, 4.0'f), (16, 24))
let volume = initValueCubic(42).grid((0, 0, 0), (16, 24, 8))
let height = surface[2, 3]
let density = volume[2, 3, 4]
```

A zero dimension returns an empty grid. Negative dimensions, sizes that
overflow storage, and nonpositive octave counts raise `NoisyError`.

## Settings and seed compatibility

All five types expose `seed` and share these octave defaults.

| Setting | Default |
| --- | ---: |
| `octaves` | 1 |
| `amplitude` | 1 |
| `frequency` | 1 |
| `lacunarity` | 2 |
| `gain` | 0.5 |

`seed` accepts an integer and uses its low 32 bits for hashing. Zero and negative
seeds are supported. The implementation explicitly wraps hash arithmetic so
hash results agree between Nim debug and release builds. The same settings and
seed match the selected FastNoiseLite single-octave kernel. Value noise
preserves FastNoiseLite's collision between seed 0 and the minimum signed
32-bit seed.

Octaves follow Noisy's convention: keep the same seed, scale coordinates by
frequency, multiply each sample by amplitude, then divide the sum by the octave
count. Frequency and amplitude are multiplied by lacunarity and gain after each
layer. FastNoiseLite's own fractal API has different seed progression and
normalization, so the reference tests compose its single-octave samples using
Noisy's rules. Noisy's default frequency is 1 rather than FastNoiseLite's 0.01.

Nonpositive octave counts raise `NoisyError`. Use finite coordinates and
settings, keeping scaled and transformed lattice coordinates comfortably
inside the signed 32-bit range. Float32 loses fractional resolution at large
coordinates. Both OpenSimplex2 and OpenSimplex2S use FastNoiseLite's default
2D skew and 3D rotation.

`OpenSimplex2S` samples a modified triangular lattice in 2D and two offset
cube lattices in 3D. `ValueNoise` uses Hermite interpolation between two
lattice samples per axis. `ValueCubic` uses four samples per axis with the
reference's amplitude scale of `1 / 1.5^dimensions`.

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
nim r tools/gen_references.nim --cache-dir /tmp/noisy-references
nim r -f -d:noisyReferenceDir=/tmp/noisy-references tests/test_fastnoise_reference.nim
nim cpp -r -f -d:release -d:noisyReferenceDir=/tmp/noisy-references tests/test_fastnoise_reference.nim
```

The unchanged C++ reference is compiled with `-fwrapv` for its signed hash
arithmetic. This flag applies only to the reference translation unit. Tests
compare the public APIs with matched seeds, coordinates, and settings, including
negative cells, lattice boundaries, octave combinations, and all Cellular
metrics and return modes. The normal offline suite contains recorded outputs
from the compiled reference as well as behavior tests:

```sh
nim check tests/tests.nim
nim r tests/tests.nim
```

## Recorded validation

With Nim 2.4.0 and Apple Clang, 190,152 comparisons per build passed in native
ARM64 C debug and C++ release, and x86-64 C++ release under Rosetta.

| Noise | Comparisons per build | Maximum absolute error |
| --- | ---: | ---: |
| OpenSimplex2 | 24,120 | 6.56e-7 |
| OpenSimplex2S | 24,120 | 4.18e-7 |
| Value | 24,120 | 0 |
| ValueCubic | 24,120 | 0 |
| Cellular | 93,672 | 4.77e-7 |

Tests use `1e-5 + 1e-5 * max(abs(actual), abs(expected))` for float32 rounding.
Single-octave `CellValue` samples additionally require exact equality, so
selecting a different feature cell cannot pass through a numeric tolerance.
The offline suite includes 278 recorded C++ outputs.
