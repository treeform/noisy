import std/math, common

export NoisyError

const
  F2* = (0.5 * (sqrt(3.0) - 1)).float32
  G2* = ((3 - sqrt(3.0)) / 6).float32
  F3* = (1.0 / 3.0).float32
  G3* = (1.0 / 6.0).float32
  F4* = ((sqrt(5.0) - 1) / 4).float32
  G4* = ((5 - sqrt(5.0)) / 20).float32
  # Gradients follow Gustavson's public-domain simplexnoise1234 reference.
  Grad4* = block:
    var gradients: array[32, array[4, float32]]
    for i in 0 ..< 32:
      let
        first = if i < 24: 0 else: 1
        second = if i < 16: 1 else: 2
        third = if i < 8: 2 else: 3
      gradients[i][first] = if (i and 1) == 0: 1.0'f else: -1.0'f
      gradients[i][second] = if (i and 2) == 0: 1.0'f else: -1.0'f
      gradients[i][third] = if (i and 4) == 0: 1.0'f else: -1.0'f
    gradients

  grad3* = [
    [1.float32, 1, 0], [-1.float32, 1, 0], [1.float32, -1, 0],
    [-1.float32, -1, 0], [1.float32, 0, 1], [-1.float32, 0, 1],
    [1.float32, 0, -1], [-1.float32, 0, -1], [0.float32, 1, 1],
    [0.float32, -1, 1], [0.float32, 1, -1], [0.float32, -1, -1]
  ]

type
  Simplex* = object
    octaves*: int
    amplitude*, frequency*, lacunarity*, gain*: float32
    perm*: array[512, int32]
    permMod12*: array[512, int32]

  Grid* = ref object
    width*, height*, depth*: int
    frames*: int ## Number of fourth-axis slices; one for 1D through 3D grids.
    values*: seq[float32]
