"""Fetch verified C/C++ noise references and optionally benchmark Noisy."""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from urllib.request import urlopen


REPO = Path(__file__).resolve().parents[1]
DEFAULT_CACHE = REPO / "tests/reference/.cache"


def benchmark(cache, rounds, batches):
    """Compare automatic SIMD and the fallback with the same pinned references."""
    nim = shutil.which("nim")
    if nim is None:
        raise ValueError("nim is not on PATH; install Nim and dependencies with Nimby")
    output = cache / "results"
    output.mkdir(parents=True, exist_ok=True)
    modes = ["auto"]
    if platform.machine().lower() in ("amd64", "x86_64", "i386", "i686"):
        modes.append("sse2")
    reports, commands = [], []
    for mode in modes:
        report = output / f"{mode}.json"
        command = [
            nim, "c", "-r", "-f", "-d:release", "--listCmd",
            f"--nimcache:{cache / 'nimcache' / mode}",
            f"--out:{cache / ('bench_' + mode + ('.exe' if sys.platform == 'win32' else ''))}",
            f"-d:noisyReferenceDir={cache}",
            f"-d:noisyBenchmarkRounds={rounds}",
            f"-d:noisyBenchmarkBatches={batches}",
            f"-d:noisyBenchmarkOutput={report}",
        ]
        if mode == "sse2":
            command.append("-d:noisyNoAvx")
        command.append("tests/bench_reference.nim")
        commands.append(command)
        print(f"Benchmarking {mode}; log: {output / (mode + '.log')}", flush=True)
        with (output / f"{mode}.log").open("w", encoding="utf-8") as log:
            subprocess.run(command, cwd=REPO, stdout=log,
                           stderr=subprocess.STDOUT, check=True)
        data = json.loads(report.read_text())
        reports.append(data)
        print(f"Completed {data['simd_backend']}", flush=True)
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
    status = subprocess.check_output(
        ["git", "status", "--porcelain"], cwd=REPO, text=True).strip()
    manifest = REPO / "tests/reference/sources.json"
    summary = {
        "recorded_utc": datetime.now(timezone.utc).isoformat(),
        "git_revision": revision,
        "working_tree_dirty": bool(status),
        "platform": platform.platform(),
        "processor": platform.processor(),
        "source_manifest_sha256": digest(manifest.read_bytes()),
        "commands": commands,
        "reports": reports,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    header = ["Dim", "Octaves", "Reference ms", "Scalar ms"]
    header += [f"{report['simd_backend']} ms" for report in reports]
    lines = ["| " + " | ".join(header) + " |",
             "| " + " | ".join(["---"] + ["---:"] * (len(header) - 1)) + " |"]
    for i, case in enumerate(reports[0]["cases"]):
        row = [f"{case['dimensions']}D", str(case["octaves"]),
               f"{case['reference_ms']:.3f}", f"{case['scalar_ms']:.3f}"]
        row += [f"{report['cases'][i]['simd_ms']:.3f}" for report in reports]
        lines.append("| " + " | ".join(row) + " |")
    table = "\n".join(lines)
    (output / "table.md").write_text(table + "\n", encoding="utf-8")
    print(table)
    print(f"Raw round timings, commands, and logs: {output}")


def digest(data):
    """Return the digest used by the checked-in source manifest."""
    return hashlib.sha256(data).hexdigest()


def write_verified(path, data, expected):
    """Check downloaded or extracted bytes before replacing a cache file."""
    if digest(data) != expected:
        raise ValueError(f"SHA-256 mismatch: {path.name}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        stream.write(data)
    try:
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def main():
    """Populate an ignored cache without committing third-party sources."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cache-dir", type=Path,
        default=DEFAULT_CACHE,
        help="Cache directory (default: ignored tests/reference/.cache)",
    )
    parser.add_argument("--offline", action="store_true",
                        help="Verify cached files without downloading anything")
    parser.add_argument("--benchmark", action="store_true",
                        help="Build and run the reference benchmark with automatic SIMD and SSE2")
    parser.add_argument("--rounds", type=int, default=11)
    parser.add_argument("--batches", type=int, default=4)
    args = parser.parse_args()
    if args.rounds < 3 or args.batches < 1:
        parser.error("--rounds must be at least 3 and --batches at least 1")
    cache = args.cache_dir.expanduser().resolve()
    try:
        cache.relative_to(REPO)
    except ValueError:
        pass
    else:
        try:
            cache.relative_to(DEFAULT_CACHE)
        except ValueError:
            parser.error("use tests/reference/.cache or a cache outside the repository")

    manifest = json.loads((REPO / "tests/reference/sources.json").read_text())
    for source in manifest["sources"]:
        target = cache / source["path"]
        if target.is_file() and digest(target.read_bytes()) == source["sha256"]:
            continue
        if args.offline:
            raise ValueError(f"missing or invalid cached source: {target}")
        print(f"Fetching {source['path']}")
        with urlopen(source["url"], timeout=30) as response:
            data = response.read()
        write_verified(target, data, source["sha256"])

    # Keep the complete original for verification, and extract only the two
    # standalone kernels. The copied bytes retain their original formatting.
    original = (cache / "caseman/_simplex.c").read_bytes()
    prefix = original[:original.index(b'#include "Python.h"')]
    start = original.index(b"// 2D simplex skew factors")
    end = original.index(b"static inline float\nfbm_noise3", start)
    extracted = manifest["extracted"]
    write_verified(cache / extracted["path"], prefix + original[start:end],
                   extracted["sha256"])
    print(f"Verified references: {cache}")
    print('Run from the Noisy repository:')
    print(f'nim r -f "-d:noisyReferenceDir={cache}" tests/test_reference.nim')
    if args.benchmark:
        benchmark(cache, args.rounds, args.batches)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        sys.exit(f"Reference workflow failed: {error}")
