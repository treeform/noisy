# Run from the repository root: nim r -d:release --path:src examples/tiles.nim
import noisy

var noise = initSimplex(1988)
noise.frequency = 0.025
noise.octaves = 3
noise.lacunarity = 1.7

# Both forms sample the same seamless 128 x 128 field.
let texture = noise.tileGrid((0, 0), (128, 128))
var buffer = newSeq[float32](128 * 128)
noise.tileGrid((128, 0), (128, 128), buffer)
doAssert buffer == texture.values

# Request the endpoint explicitly when checking the seam.
let seam = noise.tileGrid((0, 0), (129, 129), (128, 128))
for i in 0 ..< 129:
  doAssert abs(seam[0, i] - seam[128, i]) < 0.00003'f
  doAssert abs(seam[i, 0] - seam[i, 128]) < 0.00003'f

# Separate chunks retain the full field's period.
let rightHalf = noise.tileGrid((64, 0), (64, 128), (128, 128))
for x in 0 ..< 64:
  for y in 0 ..< 128:
    doAssert abs(rightHalf[x, y] - texture[x + 64, y]) < 0.00003'f

# All noise families also support 1D loops through their 2D samplers.
let loop = initValue(42).tileGrid(0, 65, 64)
doAssert abs(loop[0] - loop[64]) < 0.00003'f
echo "Seamless tiles, reusable buffers, and adjacent chunks verified."
