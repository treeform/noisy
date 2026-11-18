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
