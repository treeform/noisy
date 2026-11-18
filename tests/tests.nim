{.push warning[UnusedImport]: off.}
import test_api, test_noisy, test_dimensions, test_simd, test_perlin, test_avx2
import test_opensimplex2, test_opensimplex2s
import test_values, test_valuecubic, test_cellular
when defined(noisyReferenceDir):
  import test_reference
else:
  echo "External C/C++ comparisons skipped; see tests/reference/README.md"
{.pop.}

when defined(arm64):
  echo "All noisy tests passed (ARM64 NEON)"
else:
  when defined(amd64) and (defined(gcc) or defined(clang) or defined(vcc)) and
      not defined(noisyNoAvx):
    import noisy/simd/runtimecheck
    echo "All noisy tests passed (x86 ", (if hasAvx2(): "AVX2" else: "SSE2"), ")"
  else:
    echo "All noisy tests passed (x86 SSE2)"
