import
  std/[json, os, osproc, streams, strformat, strutils, tempfiles, times],
  sha256

const
  Repo = currentSourcePath().parentDir.parentDir
  DefaultCache = Repo / "tests/reference/.cache"
  ManifestPath = Repo / "tests/reference/sources.json"
  Usage = """
Fetch verified C/C++ noise references and optionally benchmark Noisy.

Usage: nim r tools/gen_references.nim [options]
  --cache-dir PATH  Write to PATH (default: tests/reference/.cache)
  --offline         Verify existing cached sources without network access
  --benchmark       Benchmark all noise families against the pinned references
  --rounds N        Use N benchmark rounds (default: 11, minimum: 3)
  --batches N       Use N grids per benchmark round (default: 4, minimum: 1)
  --help, -h        Show this help message
"""

type Options = object
  cache: string
  offline, benchmark: bool
  rounds, batches: int

proc parseOptions(args: openArray[string]): Options =
  result = Options(cache: DefaultCache, rounds: 11, batches: 4)
  var i = 0
  while i < args.len:
    let arg = args[i].split('=', 1)
    inc i
    case arg[0]
    of "--cache-dir", "--rounds", "--batches":
      let value =
        if arg.len == 2:
          arg[1]
        elif i < args.len and not args[i].startsWith("--"):
          inc i
          args[i - 1]
        else:
          raise newException(ValueError, "missing argument for " & arg[0])
      case arg[0]
      of "--cache-dir":
        result.cache = value
      of "--rounds":
        result.rounds = parseInt(value)
      of "--batches":
        result.batches = parseInt(value)
      else:
        discard
    else:
      case args[i - 1]
      of "--offline":
        result.offline = true
      of "--benchmark":
        result.benchmark = true
      of "--help", "-h":
        stdout.write Usage
        quit(QuitSuccess)
      else:
        raise newException(ValueError, "unknown argument: " & args[i - 1])
  if result.rounds < 3 or result.batches < 1:
    raise newException(ValueError, "--rounds must be at least 3 and " &
      "--batches must be at least 1")

proc resolvePath(path: string): string =
  ## Resolves existing ancestors before appending a new cache directory name.
  if fileExists(path) or dirExists(path):
    return expandFilename(path)
  let parent = path.parentDir
  if parent.len == 0 or parent == path:
    return path
  resolvePath(parent) / path.lastPathPart

proc within(path, directory: string): bool =
  let relative = relativePath(path, directory, '/')
  relative != ".." and not relative.startsWith("../") and not relative.isAbsolute()

proc runCommand(command: string, args: openArray[string], logPath = ""): string =
  ## Captures combined output or streams it to a log while the command runs.
  let log = if logPath.len > 0: open(logPath, fmWrite) else: nil
  defer:
    if log != nil:
      log.close()
  let process = startProcess(command, workingDir = Repo, args = args,
    options = {poUsePath, poStdErrToStdOut})
  defer: process.close()
  let output = process.outputStream()
  if log == nil:
    result = output.readAll()
  else:
    var buffer: array[8192, char]
    while true:
      let count = output.readData(addr buffer[0], buffer.len)
      if count == 0:
        break
      if log.writeBuffer(addr buffer[0], count) != count:
        raise newException(IOError, "could not write " & logPath)
  let status = process.waitForExit()
  if status != 0:
    let detail =
      if logPath.len > 0: "see " & logPath
      else: result.strip()
    raise newException(IOError, &"{command} exited with status {status}: {detail}")

proc fetch(url: string): string =
  ## Uses the platform curl TLS setup and reads the response without text conversion.
  let curl = findExe("curl")
  if curl.len == 0:
    raise newException(IOError, "curl is required to fetch reference sources")
  let (file, path) = createTempFile("noisy-reference-", ".download")
  file.close()
  defer: removeFile(path)
  discard runCommand(curl, [
    "--disable", "--fail", "--location", "--silent", "--show-error",
    "--max-time", "30", "--output", path, "--", url
  ])
  readFile(path)

proc writeVerified(path, data, expected: string) =
  ## Checks exact bytes before atomically replacing the file in the cache.
  if sha256(data) != expected:
    raise newException(ValueError, "SHA-256 mismatch: " & path.extractFilename)
  createDir(path.parentDir)
  let (file, temporary) = createTempFile("reference-", ".tmp", path.parentDir)
  defer: removeFile(temporary)
  # Reopen in binary mode because tempfiles uses a text stream on Windows.
  file.close()
  writeFile(temporary, data)
  moveFile(temporary, path)

proc populateCache(cache: string, manifest: JsonNode, offline: bool) =
  for source in manifest["sources"]:
    let
      relative = source["path"].getStr()
      target = cache / relative
      expected = source["sha256"].getStr()
    if fileExists(target) and sha256(readFile(target)) == expected:
      continue
    if offline:
      raise newException(ValueError, "missing or invalid cached source: " & target)
    echo "Fetching ", relative
    writeVerified(target, fetch(source["url"].getStr()), expected)

  # Keep the original file in the cache and preserve the exact extracted bytes.
  let
    original = readFile(cache / "caseman/_simplex.c")
    prefixEnd = original.find("#include \"Python.h\"")
    start = original.find("// 2D simplex skew factors")
    stop = original.find("static inline float\nfbm_noise3", max(0, start))
    extracted = manifest["extracted"]
  if prefixEnd < 0 or start < 0 or stop < 0:
    raise newException(ValueError, "cannot find the pinned 2D and 3D simplex kernels")
  writeVerified(cache / extracted["path"].getStr(),
    original[0 ..< prefixEnd] & original[start ..< stop],
    extracted["sha256"].getStr())

proc markdownRow(cells: openArray[string]): string =
  "| " & cells.join(" | ") & " |\n"

proc timingTable(reports: openArray[JsonNode]): string =
  var
    header = @["Noise", "Dim", "Octaves", "C/C++ reference ms", "Noisy scalar ms"]
    alignment = @["---", "---", "---:", "---:", "---:"]
  for report in reports:
    header.add "Noisy " & report["simd_backend"].getStr() & " ms"
    alignment.add "---:"
  result = markdownRow(header) & markdownRow(alignment)
  let cases = reports[0]["cases"]
  for i in 0 ..< cases.len:
    let sample = cases[i]
    var row = @[
      sample["noise"].getStr(),
      $sample["dimensions"].getInt() & "D",
      $sample["octaves"].getInt(),
      sample["reference_ms"].getFloat().formatFloat(ffDecimal, 3),
      sample["scalar_ms"].getFloat().formatFloat(ffDecimal, 3)
    ]
    for report in reports:
      let timing = report["cases"][i]["simd_ms"]
      row.add:
        if timing.kind == JNull: "—"
        else: timing.getFloat().formatFloat(ffDecimal, 3)
    result.add markdownRow(row)

proc benchmark(cache: string, options: Options, manifest: string) =
  let nim = findExe("nim")
  if nim.len == 0:
    raise newException(IOError, "nim is not on PATH")
  let output = cache / "results"
  createDir(output)
  var
    reports: seq[JsonNode]
    commands: seq[seq[string]]
  for mode in ["auto", "sse2"]:
    if mode == "sse2" and reports[0]["simd_backend"].getStr() != "AVX2":
      break
    let
      reportPath = output / (mode & ".json")
      logPath = output / (mode & ".log")
      exeSuffix = when defined(windows): ".exe" else: ""
    var args = @[
      "c", "-r", "-f", "-d:release", "--listCmd",
      "--nimcache:" & (cache / "nimcache" / mode),
      "--out:" & (cache / ("bench_" & mode & exeSuffix)),
      "-d:noisyReferenceDir=" & cache,
      "-d:noisyBenchmarkRounds=" & $options.rounds,
      "-d:noisyBenchmarkBatches=" & $options.batches,
      "-d:noisyBenchmarkOutput=" & reportPath
    ]
    if mode == "sse2":
      args.add "-d:noisyNoAvx"
    args.add "tests/bench_reference.nim"
    commands.add @[nim] & args
    echo &"Benchmarking {mode}; log: {logPath}"
    stdout.flushFile()
    discard runCommand(nim, args, logPath)
    let data = parseFile(reportPath)
    reports.add data
    echo "Completed ", data["simd_backend"].getStr()

  let report = %*{
    "recorded_utc": $now().utc,
    "git_revision": runCommand("git", ["rev-parse", "HEAD"]).strip(),
    "working_tree_dirty": runCommand("git", ["status", "--porcelain"]).strip().len > 0,
    "platform": (when defined(posix): runCommand("uname", ["-srm"]).strip()
                 else: getEnv("OS", hostOS)),
    "processor": getEnv("PROCESSOR_IDENTIFIER", hostCPU),
    "source_manifest_sha256": sha256(manifest),
    "commands": commands,
    "reports": reports
  }
  writeFile(output / "summary.json", report.pretty() & "\n")
  let table = timingTable(reports)
  writeFile(output / "table.md", table)
  stdout.write table
  echo "Raw round timings, commands, and logs: ", output

proc main() =
  let
    options = parseOptions(commandLineParams())
    cache = resolvePath(options.cache.expandTilde().absolutePath().normalizedPath())
    repoRoot = resolvePath(Repo)
  if cache.within(repoRoot) and not cache.within(resolvePath(DefaultCache)):
    raise newException(ValueError, "use tests/reference/.cache or a cache " &
      "outside the repository")
  let manifest = readFile(ManifestPath)
  populateCache(cache, parseJson(manifest), options.offline)
  echo "Verified references: ", cache
  echo "Run from the Noisy repository:"
  echo quoteShellCommand([
    "nim", "r", "-f", "-d:noisyReferenceDir=" & cache, "tests/test_reference.nim"
  ])
  if options.benchmark:
    benchmark(cache, options, manifest)

when isMainModule:
  try:
    main()
  except CatchableError as error:
    quit("Reference workflow failed: " & error.msg, QuitFailure)
