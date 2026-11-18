## AVX2 needs both CPU support and OS support for saving XMM/YMM registers.
## Keep this module compiled for the baseline ISA, separate from avx2.nim.

proc cpuid(leaf, subleaf: int32): array[4, int32] =
  when defined(vcc):
    proc cpuidEx(info: ptr int32, leaf, subleaf: int32)
      {.cdecl, importc: "__cpuidex", header: "intrin.h".}
    cpuidEx(result[0].addr, leaf, subleaf)
  else:
    var eax, ebx, ecx, edx: int32
    asm """
      cpuid
      :"=a"(`eax`), "=b"(`ebx`), "=c"(`ecx`), "=d"(`edx`)
      :"a"(`leaf`), "c"(`subleaf`)"""
    result = [eax, ebx, ecx, edx]

proc avxStateEnabled(): bool =
  when defined(vcc):
    proc xgetbv(index: uint32): uint64
      {.cdecl, importc: "_xgetbv", header: "immintrin.h".}
    result = (xgetbv(0) and 6) == 6
  else:
    var eax, edx: uint32
    asm """
      xgetbv
      :"=a"(`eax`), "=d"(`edx`)
      :"c"(0)"""
    result = (eax and 6) == 6

proc checkAvx2(): bool =
  if cpuid(0, 0)[0] < 7:
    return false
  # XSAVE, OSXSAVE and AVX must be present before executing XGETBV.
  const avxFlags = (1'i32 shl 26) or (1'i32 shl 27) or (1'i32 shl 28)
  if (cpuid(1, 0)[2] and avxFlags) != avxFlags:
    return false
  if not avxStateEnabled():
    return false
  (cpuid(7, 0)[1] and (1'i32 shl 5)) != 0

let cpuHasAvx2 = checkAvx2()

func hasAvx2*(): bool {.inline.} =
  # This immutable capability is initialized once, before any grid is generated.
  {.cast(noSideEffect).}:
    result = cpuHasAvx2
