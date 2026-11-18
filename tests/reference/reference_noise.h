#ifndef NOISY_REFERENCE_NOISE_H
#define NOISY_REFERENCE_NOISE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Samples the independent upstream kernels with their original permutation. */
float referenceNoise1(float x);
float referenceNoise2(float x, float y);
float referenceNoise3(float x, float y, float z);
float referenceNoise4(float x, float y, float z, float w);

/* Returns one entry of the original, repeating 256-entry permutation. */
int32_t referencePermutation(int32_t index);

#ifdef __cplusplus
}
#endif

#endif
