{.push warning[UnusedImport]: off.}
import test_api, test_noisy, test_dimensions, test_simd, test_perlin
when defined(noisyReferenceDir):
  import test_reference
else:
  echo "External C/C++ comparisons skipped; see tests/reference/README.md"
{.pop.}

when defined(arm64):
  echo "All noisy tests passed (ARM64 NEON)"
else:
  echo "All noisy tests passed (x86 SSE2)"
