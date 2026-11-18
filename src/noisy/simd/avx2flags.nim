# Include in every translation unit that operates on 256-bit vectors.
# Baseline dispatch and scalar code must not inherit these compiler options.
when defined(gcc) or defined(clang):
  {.localPassC: "-mavx2 -ffp-contract=off".}
elif defined(vcc):
  {.localPassC: "/arch:AVX2 /fp:strict".}

# These helpers own no managed memory. Avoid ARC/ORC generating vector copy
# hooks in nimsimd's baseline translation unit, which has a different AVX ABI.
# Force inlining also avoids MinGW's vector return storage alignment problem.
when defined(gcc) or defined(clang):
  {.pragma: avxInline, inline, nodestroy,
    codegenDecl: "static inline __attribute__((always_inline)) $# $#$#".}
elif defined(vcc):
  {.pragma: avxInline, inline, nodestroy,
    codegenDecl: "static __forceinline $# $#$#".}
