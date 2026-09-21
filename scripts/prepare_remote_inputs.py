"""Build the fixed, portable input bundle for remote reconstruction work.

The bundle contains copies, never symlinks.  It does not upload anything.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import stat
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


REVISE_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_ROOT = REVISE_ROOT.parent / "REVISE_Analysis_Agent"
DEFAULT_OUTPUT = REVISE_ROOT / "upload" / "remote-inputs"
EXPECTED_TOTAL_BYTES = 3_114_093_428
CHUNK_SIZE = 1024 * 1024


@dataclass(frozen=True)
class BundleInput:
    """One fixed source and its portable path inside the bundle."""

    repo: str
    relative_path: Path
    purpose: str

    @property
    def target(self) -> Path:
        return Path(self.repo) / self.relative_path


BUNDLE_INPUTS = (
    BundleInput(
        "REVISE",
        Path("raw_data/Real_application/P2CRC_Xenium.h5ad"),
        "P2CRC Xenium spatial input",
    ),
    BundleInput(
        "REVISE",
        Path("raw_data/Real_application/P1CRC_HD.h5ad"),
        "P1CRC Visium HD spatial input",
    ),
    BundleInput(
        "REVISE",
        Path("raw_data/Real_application/P2CRC_Visium.h5ad"),
        "P2CRC Visium spatial input",
    ),
    BundleInput(
        "REVISE",
        Path("raw_data/Real_application/adata_sc_all_reanno.h5ad"),
        "annotated single-cell reference input",
    ),
    BundleInput(
        "REVISE",
        Path("output/mini-acceptance/20260920/assembly/original_spatial.h5ad"),
        "historical labels and coordinates, zero-gene comparison baseline (not full reconstructed SVC)",
    ),
    BundleInput(
        "REVISE_Analysis_Agent",
        Path("resources/h.all.v2025.1.Hs.symbols.gmt"),
        "Hallmark gene-set resource",
    ),
)


@dataclass(frozen=True)
class PreparedInput:
    item: BundleInput
    source: Path
    size: int
    sha256: str


def file_sha256(path: Path) -> str:
    """Return the SHA-256 digest of one regular file."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(CHUNK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def _source_for(item: BundleInput, revise_root: Path, analysis_root: Path) -> Path:
    root = revise_root if item.repo == "REVISE" else analysis_root
    return root / item.relative_path


def _preflight(
    revise_root: Path, analysis_root: Path, expected_total_size: int | None
) -> tuple[PreparedInput, ...]:
    """Validate every source before creating a destination directory."""

    prepared: list[PreparedInput] = []
    problems: list[str] = []
    for item in BUNDLE_INPUTS:
        source = _source_for(item, revise_root, analysis_root)
        if source.is_symlink():
            problems.append(f"source is a symlink: {source}")
            continue
        try:
            mode = source.lstat().st_mode
        except FileNotFoundError:
            problems.append(f"missing source: {source}")
            continue
        if not stat.S_ISREG(mode):
            problems.append(f"source is not a regular file: {source}")
            continue
        prepared.append(
            PreparedInput(item, source, source.stat().st_size, file_sha256(source))
        )

    if problems:
        raise FileNotFoundError("Remote-input bundle preflight failed:\n" + "\n".join(problems))

    total_size = sum(record.size for record in prepared)
    if expected_total_size is not None and total_size != expected_total_size:
        raise ValueError(
            "Source byte total changed: "
            f"expected {expected_total_size}, found {total_size}"
        )
    return tuple(prepared)


def copy_regular_file(source: Path, destination: Path) -> None:
    """Copy bytes without preserving a link relationship to the source."""

    shutil.copyfile(source, destination)


def _write_readme(destination: Path) -> None:
    (destination / "README.md").write_text(
        "# Remote reconstruction inputs\n\n"
        "This directory is an upload-ready copy of the fixed inputs from two "
        "sibling repositories: `REVISE/` and `REVISE_Analysis_Agent/`. "
        "It does not upload data.\n\n"
        "Verify copied payloads from this directory with:\n\n"
        "    sha256sum -c SHA256SUMS\n\n"
        "When merging this bundle into another checkout, copy a target only when "
        "the destination is missing or its SHA-256 is the same. Stop for any "
        "different existing destination; do not overwrite it.\n",
        encoding="utf-8",
    )


def _write_metadata(destination: Path, prepared: Sequence[PreparedInput]) -> None:
    records = [
        {
            "purpose": record.item.purpose,
            "source": str(record.source.absolute()),
            "repo_relative_target": record.item.target.as_posix(),
            "size": record.size,
            "sha256": record.sha256,
        }
        for record in prepared
    ]
    manifest = {
        "purpose": (
            "Portable bundle for the three REVISE reconstruction routes, the "
            "single-cell analysis reference, the historical-label and coordinate "
            "zero-gene comparison baseline, and the Hallmark resource. Upload is "
            "out of scope."
        ),
        "inputs": records,
        "total_size": sum(record.size for record in prepared),
    }
    (destination / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (destination / "SHA256SUMS").write_text(
        "".join(
            f"{record.sha256}  {record.item.target.as_posix()}\n"
            for record in prepared
        ),
        encoding="utf-8",
    )


def build_bundle(
    output: Path = DEFAULT_OUTPUT,
    *,
    revise_root: Path = REVISE_ROOT,
    analysis_root: Path = ANALYSIS_ROOT,
    expected_total_size: int | None = EXPECTED_TOTAL_BYTES,
) -> dict[str, object]:
    """Copy and verify the six fixed inputs into a new bundle directory.

    A pre-existing destination is never touched, including a destination symlink.
    On failure, only the destination created by this invocation is removed.
    """

    destination = Path(output).expanduser().absolute()
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"Refusing existing destination: {destination}")

    prepared = _preflight(
        Path(revise_root), Path(analysis_root), expected_total_size
    )
    created_destination = False
    try:
        destination.mkdir(parents=True)
        created_destination = True
        for record in prepared:
            target = destination / record.item.target
            target.parent.mkdir(parents=True, exist_ok=True)
            copy_regular_file(record.source, target)
            source_after_copy = file_sha256(record.source)
            target_digest = file_sha256(target)
            if source_after_copy != record.sha256:
                raise RuntimeError(f"Source changed during copy: {record.source}")
            if target_digest != record.sha256:
                raise RuntimeError(f"Copied file hash mismatch: {target}")
        _write_readme(destination)
        _write_metadata(destination, prepared)
    except BaseException:
        if created_destination:
            shutil.rmtree(destination)
        raise

    return {
        "output": str(destination),
        "file_count": len(prepared),
        "total_size": sum(record.size for record in prepared),
        "manifest": str(destination / "manifest.json"),
    }


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="new bundle directory (default: REVISE/upload/remote-inputs)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    result = build_bundle(args.output)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
