{.push warning[UnusedImport]: off.}
import test_noisy, test_dimensions, test_simd
{.pop.}

when defined(arm64):
  echo "All noisy tests passed (ARM64 NEON)"
else:
  echo "All noisy tests passed (x86 SSE2)"
