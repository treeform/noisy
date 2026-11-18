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

/* Samples Gustavson's unchanged improved Perlin kernels. */
float referencePerlin1(float x);
float referencePerlin2(float x, float y);
float referencePerlin3(float x, float y, float z);
float referencePerlin4(float x, float y, float z, float w);
int32_t referencePerlinPermutation(int32_t index);

/* Returns one entry of the original, repeating 256-entry permutation. */
int32_t referencePermutation(int32_t index);

/* Fills caller-owned storage in Noisy's w/x/y/z order, with one bridge call.
 * The caller supplies dimensions 1..4, positive octaves, nonnegative sizes,
 * four start coordinates, and width*height*depth*frames output floats.
 * Unused shape extents must be 1. */
void referenceGrid(float *output, int32_t dimensions, int32_t width,
                   int32_t height, int32_t depth, int32_t frames,
                   const float *start, int32_t octaves, float amplitude,
                   float frequency, float lacunarity, float gain);

#ifdef __cplusplus
}
#endif

#endif
