import
  std/os,
  noisy, pixie

block:
  var simplex = initSimplex(1988)
  simplex.frequency = 0.1
  let
    image = newImage(480, 480)
    grid = simplex.grid((0, 0), (image.width, image.height))
  for x in 0 ..< image.width:
    for y in 0 ..< image.height:
      let value = uint8(clamp((grid[x, y] + 1) * 127.5'f, 0.0'f, 255.0'f))
      image[x, y] = rgba(value, value, value, 255)
  createDir("docs")
  image.writeFile("docs/noise.png")
