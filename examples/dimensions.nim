import noisy

block:
  var simplex = initSimplex(1988)
  simplex.frequency = 0.1
  simplex.octaves = 3

  # Sample a native 1D field, or generate a line of samples at unit intervals.
  echo "1D point: ", simplex.value(0.25'f)
  let line = simplex.grid(-4.5'f, 9)
  echo "1D grid: ", line.values
  doAssert abs(line[3] - simplex.value(-1.5'f)) < 0.0001'f

  # Animate a 3D volume by stepping through the fourth coordinate.
  echo "4D point: ", simplex.value(0.25'f, -1.0'f, 2.0'f, 0.5'f)
  let volume = simplex.grid(
    (-1.0'f, -2.0'f, -3.0'f, 0.5'f),
    (2, 3, 5, 4)
  )
  echo "Frames: ", volume.frames
  echo "4D grid sample: ", volume[1, 2, 3, 1]
  doAssert abs(volume[1, 2, 3, 1] - simplex.value(0'f, 0'f, 0'f, 1.5'f)) <
    0.0001'f
