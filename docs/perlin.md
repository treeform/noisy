# Improved Perlin noise

Noisy provides native scalar 1D, 2D, 3D, and 4D Improved Perlin noise. It follows
Stefan Gustavson's public-domain C
[`noise1234.c`](https://github.com/stegu/perlin-noise/blob/a624f5a18395c0ebc04738191efddf1dbd2cad9d/src/noise1234.c)
implementation of the 2002 algorithm, using quintic interpolation. The 1D, 2D,
and 4D functions are Gustavson's dimensional extensions. This is separate from
both Simplex noise and the older cubic-interpolation Perlin algorithm.

```nim
import noisy

var noise = initPerlin(1988)
noise.frequency = 0.01
noise.octaves = 3

let
  line = noise.value(12.5'f)
  plane = noise.value(12.5'f, -7.25'f)
  volume = noise.value(12.5'f, -7.25'f, 3.0'f)
  animation = noise.value(12.5'f, -7.25'f, 3.0'f, 0.5'f)
```

`import noisy/perlins` also exposes this API without importing the Simplex SIMD
backend. Perlin has no `grid()` overloads or explicit SIMD implementation.

| Setting | Default | Meaning |
| --- | ---: | --- |
| `octaves` | 1 | Number of layers, must be positive |
| `amplitude` | 1 | First layer's amplitude |
| `frequency` | 1 | First layer's coordinate multiplier |
| `lacunarity` | 2 | Frequency multiplier after each layer |
| `gain` | 0.5 | Amplitude multiplier after each layer |

`value()` accepts `float32` coordinates or integer coordinates in each dimension
and returns `float32`. Invalid octave counts raise `NoisyError`. Octave samples
are summed and divided by the octave count, matching Simplex's convention.
Use finite coordinates and settings whose scaled coordinates fit the platform's
integer range. The normal float32 loss of fractional precision at large
coordinates still applies.

Initialization shuffles a 256-entry permutation and repeats it into the public
`perm: array[512, int32]`, using the same seed mapping as `initSimplex` in the
same Nim runtime. Seed numbers from other libraries do not imply matching
permutations. As with Simplex, Nim 1.4.8's random generator rejects seed zero.
The underlying single-octave lattice repeats every 256 units along each axis;
frequency scaling changes that period in world coordinates.

The reference's scale factors are retained: 0.188 (1D), 0.507 (2D), 0.936 (3D),
and 0.87 (4D). In particular, the 3D result includes Gustavson's 0.936 multiplier
and is not the unscaled output of Ken Perlin's original 3D listing. These scales
provide the reference's signed output convention, not a strict range guarantee
for arbitrary octave parameters.

## Reproduce the reference comparison

Upstream sources remain outside this repository. The existing downloader pins
both files by commit and SHA-256; the bridge compiles them unchanged.

```sh
python3 tools/gen_references.py --cache-dir /tmp/noisy-references
nim r -f -d:noisyReferenceDir=/tmp/noisy-references tests/test_perlin_reference.nim
nim cpp -r -f -d:release -d:noisyReferenceDir=/tmp/noisy-references tests/test_perlin_reference.nim
```

The comparison installs the reference permutation into Noisy and tests random
coordinates, negative cells, integer planes, permutation boundaries, and
multiple octave configurations in all four dimensions. It allows
`1e-5 + 1e-5 * max(abs(actual), abs(expected))` for float32 rounding and prints
comparison counts and maximum errors. The existing combined reference test and
CI also run these checks. Ordinary offline tests run with `nim r tests/tests.nim`.

With Nim 2.4.0 and Apple Clang on ARM64, all 39,308 comparisons passed in C debug
and C++ release builds with zero observed error. The offline suite also contains
28 C-generated values using an identity permutation, so it checks numerical
results without downloading reference sources. New Perlin and legacy API tests
also pass with Nim 1.4.8.
