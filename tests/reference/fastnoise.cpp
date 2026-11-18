#include "fastnoise.h"
#include "fastnoise/FastNoiseLite.h"

static FastNoiseLite configuredNoise(int32_t algorithm, int32_t seed,
                                    int32_t distance, int32_t output,
                                    float jitter) {
  FastNoiseLite noise(seed);
  switch (algorithm) {
    case 0: noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2); break;
    case 1: noise.SetNoiseType(FastNoiseLite::NoiseType_Value); break;
    case 2: noise.SetNoiseType(FastNoiseLite::NoiseType_Cellular); break;
    case 3: noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2S); break;
    case 4: noise.SetNoiseType(FastNoiseLite::NoiseType_ValueCubic); break;
  }
  noise.SetFrequency(1.0f);
  noise.SetFractalType(FastNoiseLite::FractalType_None);
  noise.SetRotationType3D(FastNoiseLite::RotationType3D_None);
  noise.SetCellularDistanceFunction(
      static_cast<FastNoiseLite::CellularDistanceFunction>(distance));
  noise.SetCellularReturnType(
      static_cast<FastNoiseLite::CellularReturnType>(output));
  noise.SetCellularJitter(jitter);
  return noise;
}

float referenceFastNoise2(int32_t algorithm, int32_t seed, float x, float y,
                         int32_t distance, int32_t output, float jitter) {
  return configuredNoise(algorithm, seed, distance, output, jitter)
      .GetNoise(x, y);
}

float referenceFastNoise3(int32_t algorithm, int32_t seed,
                         float x, float y, float z,
                         int32_t distance, int32_t output, float jitter) {
  return configuredNoise(algorithm, seed, distance, output, jitter)
      .GetNoise(x, y, z);
}

template<int Dimensions>
static void fillFastNoise(float *output, const FastNoiseLite &noise,
                          int32_t width, int32_t height, int32_t depth,
                          const float *start, int32_t octaves, float amplitude,
                          float frequency, float lacunarity, float gain) {
  const int32_t axisDepth = Dimensions == 3 ? depth : 1;
  for (int32_t x = 0; x < width; ++x) {
    const float worldX = start[0] + static_cast<float>(x);
    for (int32_t y = 0; y < height; ++y) {
      const float worldY = start[1] + static_cast<float>(y);
      for (int32_t z = 0; z < axisDepth; ++z) {
        const float worldZ = start[2] + static_cast<float>(z);
        float total = 0.0f;
        float layerAmplitude = amplitude;
        float layerFrequency = frequency;
        for (int32_t octave = 0; octave < octaves; ++octave) {
          const float sample = Dimensions == 2
              ? noise.GetNoise(worldX * layerFrequency, worldY * layerFrequency)
              : noise.GetNoise(worldX * layerFrequency, worldY * layerFrequency,
                               worldZ * layerFrequency);
          total += sample * layerAmplitude;
          layerAmplitude *= gain;
          layerFrequency *= lacunarity;
        }
        *output++ = total / static_cast<float>(octaves);
      }
    }
  }
}

void referenceFastNoiseGrid(float *output, int32_t algorithm, int32_t seed,
                            int32_t dimensions, int32_t width, int32_t height,
                            int32_t depth, const float *start, int32_t octaves,
                            float amplitude, float frequency,
                            float lacunarity, float gain) {
  const auto noise = configuredNoise(algorithm, seed,
      FastNoiseLite::CellularDistanceFunction_EuclideanSq,
      FastNoiseLite::CellularReturnType_Distance, 1.0f);
  if (dimensions == 2) {
    fillFastNoise<2>(output, noise, width, height, depth, start,
                     octaves, amplitude, frequency, lacunarity, gain);
  } else if (dimensions == 3) {
    fillFastNoise<3>(output, noise, width, height, depth, start,
                     octaves, amplitude, frequency, lacunarity, gain);
  }
}
