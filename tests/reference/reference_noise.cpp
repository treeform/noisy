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

// The Perlin reference has a different integer-boundary floor macro.
#undef FASTFLOOR
namespace SteguPerlinReference {
#include "stegu/noise1234.c"
}
#undef FADE
#undef FASTFLOOR
#undef LERP

float referencePerlin1(float x) {
  return SteguPerlinReference::noise1(x);
}

float referencePerlin2(float x, float y) {
  return SteguPerlinReference::noise2(x, y);
}

float referencePerlin3(float x, float y, float z) {
  return SteguPerlinReference::noise3(x, y, z);
}

float referencePerlin4(float x, float y, float z, float w) {
  return SteguPerlinReference::noise4(x, y, z, w);
}

int32_t referencePerlinPermutation(int32_t index) {
  return SteguPerlinReference::perm[index & 255];
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

// Kernel and dimension dispatch are outside the sampling loop, so the compiler
// can inline the upstream code without a per-sample foreign-function call.
template<bool Perlin, int Dimensions>
static float layeredReference(float x, float y, float z, float w,
                              int32_t octaves, float amplitude, float frequency,
                              float lacunarity, float gain) {
  float total = 0.0f;
  for (int32_t octave = 0; octave < octaves; ++octave) {
    const float sampleX = x * frequency;
    const float sampleY = y * frequency;
    const float sampleZ = z * frequency;
    const float sampleW = w * frequency;
    float sample;
    if (Dimensions == 1) {
      sample = Perlin ? SteguPerlinReference::noise1(sampleX)
                      : SimplexNoise::noise(sampleX);
    } else if (Dimensions == 2) {
      sample = Perlin ? SteguPerlinReference::noise2(sampleX, sampleY)
                      : CasemanReference::noise2(sampleX, sampleY);
    } else if (Dimensions == 3) {
      sample = Perlin ? SteguPerlinReference::noise3(sampleX, sampleY, sampleZ)
                      : CasemanReference::noise3(sampleX, sampleY, sampleZ);
    } else {
      sample = Perlin
          ? SteguPerlinReference::noise4(sampleX, sampleY, sampleZ, sampleW)
          : SteguReference::snoise4(sampleX, sampleY, sampleZ, sampleW);
    }
    total += sample * amplitude;
    amplitude *= gain;
    frequency *= lacunarity;
  }
  return total / static_cast<float>(octaves);
}

template<bool Perlin, int Dimensions>
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
          *next++ = layeredReference<Perlin, Dimensions>(
              worldX, worldY, worldZ, worldW, octaves, amplitude, frequency,
              lacunarity, gain);
        }
      }
    }
  }
}

template<bool Perlin>
static void dispatchReference(float *output, int32_t dimensions, int32_t width,
                              int32_t height, int32_t depth, int32_t frames,
                              const float *start, int32_t octaves,
                              float amplitude, float frequency,
                              float lacunarity, float gain) {
  switch (dimensions) {
    case 1:
      fillReference<Perlin, 1>(output, width, height, depth, frames, start,
                               octaves, amplitude, frequency, lacunarity, gain);
      break;
    case 2:
      fillReference<Perlin, 2>(output, width, height, depth, frames, start,
                               octaves, amplitude, frequency, lacunarity, gain);
      break;
    case 3:
      fillReference<Perlin, 3>(output, width, height, depth, frames, start,
                               octaves, amplitude, frequency, lacunarity, gain);
      break;
    case 4:
      fillReference<Perlin, 4>(output, width, height, depth, frames, start,
                               octaves, amplitude, frequency, lacunarity, gain);
      break;
  }
}

void referenceGrid(float *output, int32_t dimensions, int32_t width,
                   int32_t height, int32_t depth, int32_t frames,
                   const float *start, int32_t octaves, float amplitude,
                   float frequency, float lacunarity, float gain) {
  dispatchReference<false>(output, dimensions, width, height, depth, frames,
                           start, octaves, amplitude, frequency, lacunarity,
                           gain);
}

void referencePerlinGrid(float *output, int32_t dimensions, int32_t width,
                         int32_t height, int32_t depth, int32_t frames,
                         const float *start, int32_t octaves, float amplitude,
                         float frequency, float lacunarity, float gain) {
  dispatchReference<true>(output, dimensions, width, height, depth, frames,
                          start, octaves, amplitude, frequency, lacunarity,
                          gain);
}
