"""Fetch pinned C/C++ noise references into a verified external cache."""

import argparse
import hashlib
import json
from pathlib import Path
import sys
import tempfile
from urllib.request import urlopen


REPO = Path(__file__).resolve().parents[1]


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
    """Populate the cache without storing third-party code in the repository."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cache-dir", type=Path,
        default=Path(tempfile.gettempdir()) / "noisy-references",
        help="External directory for the verified reference source files",
    )
    args = parser.parse_args()
    cache = args.cache_dir.expanduser().resolve()
    try:
        cache.relative_to(REPO)
    except ValueError:
        pass
    else:
        parser.error("the reference cache must be outside the repository")

    manifest = json.loads((REPO / "tests/reference/sources.json").read_text())
    for source in manifest["sources"]:
        target = cache / source["path"]
        if target.is_file() and digest(target.read_bytes()) == source["sha256"]:
            continue
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


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        sys.exit(f"Reference fetch failed: {error}")
