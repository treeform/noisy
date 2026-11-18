# Sampling kernels

Scalar sampling uses allocation-free internal functions with a `raises: []`
contract (the `noiseKernel` pragma). Public `value()` and `grid()` entry points
still reject invalid octave settings, cellular jitter, and grid dimensions with
`NoisyError`. Grid indexing still checks its coordinates.

In release builds, `kernelchecks.nim` disables checks within the audited sampling
sections. Debug builds retain them. `-d:noisyCheckKernels` restores checks and
normal exception propagation even with `-d:release`; CI runs this configuration.
The surrounding application does not need `--threads:off`, `-d:danger`, or a
different exception backend.

## Exception code generation

Nim's `raises: []` checks catchable exceptions, but does not exclude defects such
as `IndexDefect`. With Nim 2.2's goto exception backend, ordinary calls can
therefore still generate exception-flag reads even with an empty raises list.
On Windows/MinGW, these reads can require thread-local storage lookups in every
sample and octave.

For audited kernels only, the release `noiseKernel` pragma combines `raises: []`
with Nim's [quirky pragma](https://nim-lang.org/docs/manual_experimental.html#exceptions-quirky-pragma).
The latter suppresses exception bookkeeping; it does **not** prove safety.
These functions contain arithmetic, bounded array operations, and non-throwing
math/SIMD calls, with no allocation or callbacks. Debug and checked release builds
omit `quirky`, so defects remain detectable. The compiler continues to verify
the empty catchable-exception effect in every configuration.

Validation remains outside these kernels. Public scalar octave loops spell out
floating arithmetic instead of using generic compound-assignment procedures,
which also attract conservative exception checks in Nim 2.2. Generated release
C was checked for exception-flag accesses and bounds/overflow/range handlers in
the scalar sampling functions. This removes per-sample TLS accesses from those
functions with threads enabled and goto exceptions unchanged.

## Bounds invariants

- FastNoiseLite gradient tables are constant arrays of 256 floats. Their hash
  indexes are masked with 254 in 2D and 252 in 3D; adding components with bitwise
  OR gives indexes at most 255 and 254, respectively.
- Cellular random-vector tables contain 512 and 1024 floats. Hash masks of 510
  and 1020 plus component indexes 0..1 or 0..2 stay inside those arrays. The
  dimension is statically restricted to 2 or 3. Its neighborhood loops visit
  exactly 9 or 27 cells.
- Perlin restricts dimensions to 1..4. Corner arrays contain `2^N` values. Each
  interpolation pass halves the active count, so both child indexes remain in
  the previous pass's range. Cubic Value uses fixed four-element axis arrays;
  OpenSimplex2S uses fixed three-element axis arrays.
- Simplex's 4D rank loops only address four axes and five corners. Gradient
  indexes are masked to 0..31. Its 2D/3D gradients use the guarded 0..11 lookup
  described below. SIMD lane loops use their fixed vector widths.
- Grid allocation and dimension multiplication remain validated. Vector stores
  cover complete groups of eight or four samples; scalar code fills the tail.
  No speculative vector load or store is used beyond the allocated output.

Coordinates and octave settings must keep the scaled/rotated lattice coordinates
finite and representable by the sampler's integer conversion. FastNoiseLite and
Simplex use int32 cells; scalar Perlin uses native `int`. This optimization does
not extend the numeric domain or define output for NaN, infinity, or integer
conversion overflow.

## Seeded tables

Gradient and random-vector tables are already `const`. Simplex and Perlin
permutation arrays are generated from each instance's seed and remain public
for compatibility. A fixed instance can be constructed with `const noise =
initPerlin(42)` or `const noise = initSimplex(42)`; this precomputes initialization,
but does not guarantee that runtime-indexed lookups disappear.

Because callers can edit permutation arrays, every scalar and SIMD permutation
address is now masked to 0..511. Perlin preserves full signed stored values and
uses wrapping addition before masking addresses. Simplex keeps the low nine
bits of intermediate permutations; invalid `permMod12` entries select gradient
zero. Valid lookups retain their previous output, including independently edited
table halves. Formerly out-of-bounds custom lookups now wrap within the array.

`test_kernel_bounds.nim` exercises extreme signed permutation values, random
table contents, lattice boundaries, seeds, all dimensions, complete vector
batches, and scalar tails. The existing suite also checks parameter errors,
all cellular modes, octave composition, and independent C/C++ reference values.

```sh
nim r tests/tests.nim
nim r -d:release -d:noisyCheckKernels tests/tests.nim
nim r -d:release -d:noisyCheckKernels -d:noisyNoAvx -d:noisyNoNeon tests/tests.nim
```

Measured results and raw timings are recorded in the
[four-column scalar comparison](../tests/reference/README.md#scalar-kernel-optimization-results).
