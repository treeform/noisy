#ifndef NOISY_FASTNOISE_REFERENCE_H
#define NOISY_FASTNOISE_REFERENCE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Algorithm IDs are 0 OpenSimplex2, 1 Value, 2 Cellular, 3 OpenSimplex2S,
 * and 4 ValueCubic. Cellular settings use FastNoiseLite's enum ordering.
 * Samples use frequency 1 and no fractal. */
float referenceFastNoise2(int32_t algorithm, int32_t seed, float x, float y,
                         int32_t distance, int32_t output, float jitter);
float referenceFastNoise3(int32_t algorithm, int32_t seed,
                         float x, float y, float z,
                         int32_t distance, int32_t output, float jitter);

/* Fills the output in x/y/z order using Noisy's octave rules. Cellular uses
 * squared Euclidean distance, Distance return mode, and jitter 1. The caller
 * supplies three start coordinates, positive octaves, and depth 1 in 2D. */
void referenceFastNoiseGrid(float *output, int32_t algorithm, int32_t seed,
                            int32_t dimensions, int32_t width, int32_t height,
                            int32_t depth, const float *start, int32_t octaves,
                            float amplitude, float frequency,
                            float lacunarity, float gain);

#ifdef __cplusplus
}
#endif

#endif
