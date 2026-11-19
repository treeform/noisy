# Seamless grids

Import `noisy` and use `tileGrid()` with the same origin, dimensions, and optional
destination buffer style as `grid()`. The optional `period` argument goes before
the buffer. Without it, each period defaults to the corresponding dimension
(one for an empty axis).

Run the complete example from the repository root with
`nim r -d:release --path:src examples/tiles.nim`.

```nim
var noise = initPerlin(42)
noise.frequency = 0.05
var samples = newSeq[float32](128 * 128)
noise.tileGrid((0, 0), (128, 128), samples)
noise.tileGrid((128, 0), (128, 128), samples) # Same phase, same tile.
let owned = noise.tileGrid((0, 0), (128, 128))

var strip = newSeq[float32](64)
initValue(42).tileGrid(0, 64, strip)
```

1D loops are supported for all seven families using their existing 2D sampler.
2D tiles are supported for Simplex and Perlin using their existing 4D sampler.
Unsupported combinations have no overload. This mapping requires six noise
coordinates for a fully periodic 3D volume; Noisy has no 6D sampler. It does not
add lattice-wrapping variants of individual algorithms.

## Mapping and periods

For a world coordinate `x` and a positive world-space period `p`, the mapped
coordinates are `(r * cos(theta), r * sin(theta))`, where `r = p / (2*pi)` and
`theta = 2*pi * (x modulo p) / p`. Two such circles supply the four coordinates
of a 2D tile. This is the standard
[4D torus mapping](https://github.com/KdotJPG/OpenSimplex2#4d-noise-torus-mapping).
Phase is reduced before trigonometry using float64; the noise samplers receive
float32 coordinates as usual.

The radius keeps the local sampling scale in world units. Frequency scales the
mapped coordinates, and the existing octave, amplitude, lacunarity, and gain
settings apply unchanged. The phase repeats before octave scaling, so negative
or zero frequency and non-integer lacunarity still repeat. This is a different
slice of noise from the ordinary 1D/2D `grid()` output.

The origin is a phase in world space; it does not reset the period. For example:

```nim
let a = noise.tileGrid((0, 0), (64, 128), (128, 128))
let b = noise.tileGrid((64, 0), (64, 128), (128, 128))
```

These are adjacent halves of the same periodic 128 x 128 field. Calling the
second form without the explicit period would choose a new 64 x 128 period.
The 2D layout remains `x * height + y`.

A grid samples the half-open interval: `start` through `start + size - 1`.
There is no duplicated endpoint. To check a period-128 seam directly, request
129 samples with an explicit period of 128 and compare positions 0 and 128.
Floating-point rounding limits equality for origins that cannot represent the
period shift or neighboring sample positions accurately.

## Storage, validation, and execution

The mutable-buffer form accepts sequences, fixed arrays, and `toOpenArray`
slices. It requires at least the product of the dimensions, writes only that
prefix, and does not resize or allocate storage. The returning form allocates
a `Grid` and calls the same filler. Zero dimensions produce no writes, but
explicit parameters are still validated.

Before writing, tiled sampling rejects invalid octave/jitter settings, negative
or overflowing sizes, undersized buffers, non-finite origins/settings, and
non-positive or non-finite periods. Circle radii must fit in finite float32.
At every octave, `sum(period / (2*pi)) * abs(frequency)` must be at most
268,435,456, leaving margin for simplex skew and neighboring int32 cells.
Periods and frequency/lacunarity that exceed this conservative domain raise
`NoisyError`. Extremely large amplitudes/gains can still overflow float32 noise
values, just as with ordinary sampling.

Complete groups of eight use the existing AVX2 noise kernels after circle/torus
coordinates are computed. The AVX2 boundary carries only scalar arrays and
pointers; runtime CPU/OS checks remain in baseline code. Smaller groups and
machines without AVX2 use scalar sampling. Tiled sampling does not currently
use the regular-grid SSE2/NEON optimizations. No temporary grid or coordinate
sequence is allocated. Cellular preserves the continuity characteristics of
its selected return mode: for example, `CellValue` has internal cell boundaries.

`tests/test_tiles.nim` verifies positive/negative period shifts, both seams and
their corner, independent circle/torus sample calculations, chunk stitching,
buffer guards, invalid-input behavior, all Cellular modes, and scalar/AVX2
agreement. Allocation counters are checked with `-d:nimAllocStats` in CI.
