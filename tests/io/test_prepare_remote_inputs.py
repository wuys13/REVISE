"""Focused behavior checks for the remote-input bundle builder."""

import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import pytest

from scripts import prepare_remote_inputs as bundle


def _write_fixture_sources(revise_root: Path, analysis_root: Path) -> dict[Path, bytes]:
    payloads: dict[Path, bytes] = {}
    for number, item in enumerate(bundle.BUNDLE_INPUTS, start=1):
        root = revise_root if item.repo == "REVISE" else analysis_root
        path = root / item.relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = f"fixture payload {number}\n".encode()
        path.write_bytes(payload)
        payloads[item.target] = payload
    return payloads


@pytest.fixture
def fixture_roots(tmp_path):
    revise_root = tmp_path / "REVISE"
    analysis_root = tmp_path / "REVISE_Analysis_Agent"
    payloads = _write_fixture_sources(revise_root, analysis_root)
    return revise_root, analysis_root, payloads


def _build(output: Path, revise_root: Path, analysis_root: Path):
    return bundle.build_bundle(
        output,
        revise_root=revise_root,
        analysis_root=analysis_root,
        expected_total_size=None,
    )


def _assert_sources_unchanged(revise_root, analysis_root, payloads):
    for target, payload in payloads.items():
        source_root = revise_root if target.parts[0] == "REVISE" else analysis_root
        assert (source_root / Path(*target.parts[1:])).read_bytes() == payload


def test_missing_source_fails_before_creating_destination(fixture_roots, tmp_path):
    revise_root, analysis_root, _ = fixture_roots
    (revise_root / bundle.BUNDLE_INPUTS[0].relative_path).unlink()
    output = tmp_path / "bundle"

    with pytest.raises(FileNotFoundError, match="missing source"):
        _build(output, revise_root, analysis_root)

    assert not output.exists()


def test_existing_directory_and_symlinks_are_refused(fixture_roots, tmp_path):
    revise_root, analysis_root, _ = fixture_roots
    existing = tmp_path / "existing"
    existing.mkdir()
    with pytest.raises(FileExistsError):
        _build(existing, revise_root, analysis_root)

    target = tmp_path / "target"
    target.mkdir()
    link = tmp_path / "bundle-link"
    link.symlink_to(target, target_is_directory=True)
    with pytest.raises(FileExistsError):
        _build(link, revise_root, analysis_root)
    assert link.is_symlink()

    dangling_link = tmp_path / "dangling-bundle-link"
    dangling_link.symlink_to(tmp_path / "missing-target", target_is_directory=True)
    with pytest.raises(FileExistsError):
        _build(dangling_link, revise_root, analysis_root)
    assert dangling_link.is_symlink()


def test_copy_failure_removes_only_new_destination(fixture_roots, tmp_path, monkeypatch):
    revise_root, analysis_root, payloads = fixture_roots
    output = tmp_path / "bundle"

    def fail_copy(source, destination):
        raise OSError("injected copy failure")

    monkeypatch.setattr(bundle, "copy_regular_file", fail_copy)
    with pytest.raises(OSError, match="injected copy failure"):
        _build(output, revise_root, analysis_root)
    assert not output.exists()
    _assert_sources_unchanged(revise_root, analysis_root, payloads)


def test_corrupt_copy_fails_hash_verification_and_cleans_up(
    fixture_roots, tmp_path, monkeypatch
):
    revise_root, analysis_root, payloads = fixture_roots
    output = tmp_path / "bundle"

    def corrupt_copy(source, destination):
        destination.write_bytes(b"corrupt")

    monkeypatch.setattr(bundle, "copy_regular_file", corrupt_copy)
    with pytest.raises(RuntimeError, match="Copied file hash mismatch"):
        _build(output, revise_root, analysis_root)
    assert not output.exists()
    _assert_sources_unchanged(revise_root, analysis_root, payloads)


def test_source_mutation_during_copy_is_detected(fixture_roots, tmp_path, monkeypatch):
    revise_root, analysis_root, _ = fixture_roots
    output = tmp_path / "bundle"

    def mutate_source_after_copy(source, destination):
        shutil.copyfile(source, destination)
        source.write_bytes(b"changed during copy")

    monkeypatch.setattr(bundle, "copy_regular_file", mutate_source_after_copy)
    with pytest.raises(RuntimeError, match="Source changed during copy"):
        _build(output, revise_root, analysis_root)
    assert not output.exists()


def test_normal_result_contains_six_real_copies_and_portable_checksums(
    fixture_roots, tmp_path
):
    revise_root, analysis_root, payloads = fixture_roots
    output = tmp_path / "bundle"

    result = _build(output, revise_root, analysis_root)

    assert result["file_count"] == 6
    assert result["total_size"] == sum(len(payload) for payload in payloads.values())
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["total_size"] == result["total_size"]
    assert {entry["repo_relative_target"] for entry in manifest["inputs"]} == {
        path.as_posix() for path in payloads
    }
    for target, payload in payloads.items():
        copied = output / target
        assert copied.read_bytes() == payload
        assert not copied.is_symlink()
        entry = next(
            row for row in manifest["inputs"] if row["repo_relative_target"] == target.as_posix()
        )
        source_root = revise_root if target.parts[0] == "REVISE" else analysis_root
        assert entry["source"] == str(source_root / Path(*target.parts[1:]))
        assert entry["sha256"] == hashlib.sha256(payload).hexdigest()
        assert entry["size"] == len(payload)

    sums = (output / "SHA256SUMS").read_text()
    assert all(f"  {target.as_posix()}" in sums for target in payloads)
    if shutil.which("sha256sum"):
        command = ["sha256sum", "-c", "SHA256SUMS"]
    else:
        command = ["shasum", "-a", "256", "-c", "SHA256SUMS"]
    verified = subprocess.run(command, cwd=output, check=False, capture_output=True, text=True)
    assert verified.returncode == 0, verified.stdout + verified.stderr
    readme = (output / "README.md").read_text()
    assert "sha256sum -c SHA256SUMS" in readme
    assert "missing or its SHA-256 is the same" in readme
