<img src="docs/noise.png" alt="2D simplex noise" width="480">

# Noisy - SIMD-accelerated simplex noise.

`nimby install https://github.com/treeform/noisy`

![Github Actions](https://github.com/treeform/noisy/workflows/Github%20Actions/badge.svg)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/treeform/noisy)
![GitHub Repo stars](https://img.shields.io/github/stars/treeform/noisy)
![GitHub](https://img.shields.io/github/license/treeform/noisy)
![GitHub issues](https://img.shields.io/github/issues/treeform/noisy)

[API reference](https://treeform.github.io/noisy)

## About

Noisy generates seeded simplex noise in one, two, three, and four dimensions.
It supports octave layering through `octaves`, `amplitude`, `frequency`,
`lacunarity`, and `gain`.

Point sampling with `value()` uses scalar arithmetic. Grid generation uses SSE2
for groups of samples, with scalar handling for leftover edges. Noisy depends on
Nimsimd and supports both Nim's C and C++ backends.

Originally created by Ryan Oldenburg in
[guzba/noisy](https://github.com/guzba/noisy).

> **AI disclaimer: Portions of this library and its tests were developed with AI assistance.**

## Examples

```nim
import noisy

var simplex = initSimplex(1988)
simplex.frequency = 0.1

let line = simplex.grid(-4.5'f, 16)
echo line[3]

let surface = simplex.grid((0, 0), (16, 16))
echo surface[3, 4]

let volume = simplex.grid((0, 0, 0), (8, 8, 8))
echo volume[1, 2, 3]

let animation = simplex.grid((0, 0, 0, 0), (8, 8, 8, 4))
echo animation[1, 2, 3, 0]
```

Grid coordinates advance in unit intervals before applying `frequency`.
The fourth grid dimension is `frames`. Each frame stores one 3D slice, with z
as its fastest-changing coordinate. A 1D grid has height, depth, and frames set
to one. All sampling methods accept integer or float32 coordinates.

See [examples/basic.nim](examples/basic.nim) and
[examples/dimensions.nim](examples/dimensions.nim) for runnable examples.

## Development

Install dependencies from the workspace directory containing your checkout:

```sh
cd ..
nimby create
nimby install noisy/noisy.nimble
cd noisy
```

Run the full test suite through its single entry point:

```sh
nim r tests/tests.nim
```

The suite covers independent reference values, continuity, grid bounds, octave
composition, and scalar/SIMD comparisons across seeds and parameters.

SSE2 requires an x86 target. On Apple Silicon, compile for x86_64 and run through
Rosetta:

```sh
nim r --cpu:amd64 --passC:"-arch x86_64" --passL:"-arch x86_64" tests/tests.nim
```

## Benchmarks and image generation

Optional development tools have separate dependencies. From the workspace
parent, install `benchy` for benchmarks or `pixie` for image generation with
`nimby install benchy` or `nimby install pixie`.

Run these commands from the repository root:

```sh
nim r -d:release tests/bench_noisy.nim
nim r tools/gen_noise.nim
```

The generator writes the preview to `docs/noise.png`. On Apple Silicon, use the
same x86_64 compiler flags shown above for these commands.

## Documentation

The docs workflow builds the [API reference](https://treeform.github.io/noisy)
from `src/noisy.nim`. To build it locally:

```sh
nim doc --index:on --project --git.url:https://github.com/treeform/noisy \
  --git.commit:master --out:.gh-pages src/noisy.nim
```

## License

Noisy is distributed under the [MIT license](LICENSE).
