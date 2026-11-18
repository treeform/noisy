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

#ifdef __cplusplus
}
#endif

#endif
