#!/usr/bin/env python3
"""Source-checkout entrypoint for one REVISE benchmark family."""

import sys
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from revise.benchmark.cli import main  # noqa: E402


if __name__ == "__main__":
    main()
