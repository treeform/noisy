#include "reference_noise.h"
#include <math.h>
#include "srombauts/SimplexNoise.cpp"

// Keep independent upstream tables and functions in separate namespaces.
namespace CasemanReference {
#include "caseman/_noise.h"
#include "caseman/simplex23.c"
}

#undef F2
#undef G2
#undef F3
#undef G3
#undef dot3
#undef ASSIGN

namespace SteguReference {
#include "stegu/simplexnoise1234.c"
}

float referenceNoise1(float x) {
  return SimplexNoise::noise(x);
}

float referenceNoise2(float x, float y) {
  return CasemanReference::noise2(x, y);
}

float referenceNoise3(float x, float y, float z) {
  return CasemanReference::noise3(x, y, z);
}

float referenceNoise4(float x, float y, float z, float w) {
  return SteguReference::snoise4(x, y, z, w);
}

int32_t referencePermutation(int32_t index) {
  return CasemanReference::PERM[index & 255];
}
