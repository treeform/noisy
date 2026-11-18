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

// Dimension dispatch is outside the sampling loop, so the compiler can inline
// the unchanged upstream kernel without a per-sample foreign-function call.
template<int Dimensions>
static float layeredReference(float x, float y, float z, float w,
                              int32_t octaves, float amplitude, float frequency,
                              float lacunarity, float gain) {
  float total = 0.0f;
  for (int32_t octave = 0; octave < octaves; ++octave) {
    float sample;
    if (Dimensions == 1) {
      sample = SimplexNoise::noise(x * frequency);
    } else if (Dimensions == 2) {
      sample = CasemanReference::noise2(x * frequency, y * frequency);
    } else if (Dimensions == 3) {
      sample = CasemanReference::noise3(
          x * frequency, y * frequency, z * frequency);
    } else {
      sample = SteguReference::snoise4(
          x * frequency, y * frequency, z * frequency, w * frequency);
    }
    total += sample * amplitude;
    amplitude *= gain;
    frequency *= lacunarity;
  }
  return total / static_cast<float>(octaves);
}

template<int Dimensions>
static void fillReference(float *output, int32_t width, int32_t height,
                          int32_t depth, int32_t frames, const float *start,
                          int32_t octaves, float amplitude, float frequency,
                          float lacunarity, float gain) {
  const int32_t axisHeight = Dimensions >= 2 ? height : 1;
  const int32_t axisDepth = Dimensions >= 3 ? depth : 1;
  const int32_t axisFrames = Dimensions >= 4 ? frames : 1;
  float *next = output;
  for (int32_t w = 0; w < axisFrames; ++w) {
    const float worldW = start[3] + static_cast<float>(w);
    for (int32_t x = 0; x < width; ++x) {
      const float worldX = start[0] + static_cast<float>(x);
      for (int32_t y = 0; y < axisHeight; ++y) {
        const float worldY = start[1] + static_cast<float>(y);
        for (int32_t z = 0; z < axisDepth; ++z) {
          const float worldZ = start[2] + static_cast<float>(z);
          *next++ = layeredReference<Dimensions>(
              worldX, worldY, worldZ, worldW, octaves, amplitude, frequency,
              lacunarity, gain);
        }
      }
    }
  }
}

void referenceGrid(float *output, int32_t dimensions, int32_t width,
                   int32_t height, int32_t depth, int32_t frames,
                   const float *start, int32_t octaves, float amplitude,
                   float frequency, float lacunarity, float gain) {
  switch (dimensions) {
    case 1:
      fillReference<1>(output, width, height, depth, frames, start, octaves,
                       amplitude, frequency, lacunarity, gain);
      break;
    case 2:
      fillReference<2>(output, width, height, depth, frames, start, octaves,
                       amplitude, frequency, lacunarity, gain);
      break;
    case 3:
      fillReference<3>(output, width, height, depth, frames, start, octaves,
                       amplitude, frequency, lacunarity, gain);
      break;
    case 4:
      fillReference<4>(output, width, height, depth, frames, start, octaves,
                       amplitude, frequency, lacunarity, gain);
      break;
  }
}
