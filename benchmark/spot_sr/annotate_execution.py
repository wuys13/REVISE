"""Record which formal iStar runs used the qz CPU thread limit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


ORIGINAL_THREADS = {
    ("part1", 50),
    ("part2", 50),
    ("part2", 200),
    ("part3", 50),
    ("part3", 100),
}
ISTAR_REVISION = "3cb0e5352a86df337f41c5841d0808da0003457b"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared", type=Path, required=True)
    args = parser.parse_args()
    for path in sorted(args.prepared.glob("part*/spot_*/istar_run.json")):
        part = path.parent.parent.name
        size = int(path.parent.name.removeprefix("spot_"))
        data = json.loads(path.read_text())
        data["upstream_revision"] = ISTAR_REVISION
        data["blas_openmp_threads"] = (
            "unrestricted" if (part, size) in ORIGINAL_THREADS else 1
        )
        data["execution_note"] = (
            "Original iStar process started before the thread-limit patch"
            if (part, size) in ORIGINAL_THREADS
            else "Local thread-limit patch applied before this process started"
        )
        path.write_text(json.dumps(data, indent=2) + "\n")
        print(part, size, data["blas_openmp_threads"])


if __name__ == "__main__":
    main()
