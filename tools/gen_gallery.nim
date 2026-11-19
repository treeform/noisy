## Generates the README gallery from actual Noisy grid samples.
## Run from the repository root; optionally pass a .ttf or .otf label font.
import std/os, noisy, pixie

const
  Columns = 2
  Rows = 4
  TileWidth = 384
  TileHeight = 256
  Margin = 24
  Gap = 24
  LabelHeight = 40
  HeaderHeight = 88
  FooterHeight = 48
  Seed = 1988
  Frequency = 0.025'f
  Width = 2 * Margin + Columns * TileWidth + (Columns - 1) * Gap
  Height = HeaderHeight + Rows * (LabelHeight + TileHeight) +
    (Rows - 1) * Gap + FooterHeight

proc fontPath(): string =
  if paramCount() > 0:
    return paramStr(1)
  # Fonts are used only for labels; none are copied into the repository.
  for path in [
    getEnv("WINDIR", "C:/Windows") / "Fonts/arial.ttf",
    "/System/Library/Fonts/Supplemental/Arial.ttf",
    "/Library/Fonts/Arial.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf"
  ]:
    if fileExists(path):
      return path
  quit("No label font found. Pass a .ttf/.otf path after tools/gen_gallery.nim.")

let
  gallery = newImage(Width, Height)
  font = readFont(fontPath())
gallery.fill(rgba(248, 249, 251, 255))

proc label(text: string, x, y: int, size: float32, muted = false) =
  font.size = size
  font.paint.color = if muted: color(0.36, 0.40, 0.46)
    else: color(0.10, 0.14, 0.20)
  gallery.fillText(font, text, translate(vec2(x.float32, y.float32)))

label("NOISY / SEVEN NOISE FAMILIES", Margin, 14, 27)
label("2D samples  /  seed 1988  /  frequency 0.025  /  one octave",
  Margin, 49, 19, true)

proc tile[T](index: int, title: string, initial: T) =
  var noise = initial
  noise.frequency = Frequency
  noise.octaves = 1
  let
    x0 = Margin + (index mod Columns) * (TileWidth + Gap)
    y0 = HeaderHeight + (index div Columns) * (LabelHeight + TileHeight + Gap)
    samples = noise.grid((0, 0), (TileWidth, TileHeight))
  label(title, x0, y0, 25)
  # Use one fixed transfer function for every tile, without normalizing
  # individual images: -1 = black, 0 = mid-gray, +1 = white.
  for x in 0 ..< TileWidth:
    for y in 0 ..< TileHeight:
      let value = uint8(clamp((samples[x, y] + 1) * 127.5'f, 0.0'f, 255.0'f))
      gallery[x0 + x, y0 + LabelHeight + y] = rgba(value, value, value, 255)

tile(0, "Simplex", initSimplex(Seed))
tile(1, "Improved Perlin", initPerlin(Seed))
tile(2, "OpenSimplex2", initOpenSimplex2(Seed))
tile(3, "OpenSimplex2S", initOpenSimplex2S(Seed))
tile(4, "Value", initValue(Seed))
tile(5, "ValueCubic", initValueCubic(Seed))
tile(6, "Cellular / distance", initCellular(Seed))
var cells = initCellular(Seed)
cells.returnType = CellValue
tile(7, "Cellular / cell value", cells)

label("Fixed grayscale: -1 black / 0 gray / +1 white",
  Margin, Height - 32, 19, true)
createDir("docs")
gallery.writeFile("docs/noise-gallery.png")
echo "Generated docs/noise-gallery.png (", Width, " x ", Height, ")"
