# Include only around audited sampling code, then restore with {.pop.}.
# Fixed-size loops and masked table indexes make bounds checks redundant.
# Nim's raises: [] excludes Defects; goto codegen otherwise still checks TLS.
# quirky removes that bookkeeping only for allocation-free, audited kernels.
# raises: [] still verifies their CatchableError effects. Debug/check builds
# keep all checks and ordinary exception propagation so mistakes are detectable.
when defined(release) and not defined(noisyCheckKernels):
  {.push checks: off.}
  {.pragma: noiseKernel, raises: [], quirky.}
else:
  {.push checks: on.}
  {.pragma: noiseKernel, raises: [].}
