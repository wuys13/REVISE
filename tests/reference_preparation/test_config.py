from pathlib import Path
import tempfile
import unittest

import yaml

from revise.reference_preparation.config import (
    ConfigError,
    PairedConfig,
    atomic_write_json,
    load_preparation_config,
    load_reference_config,
    write_reference_config,
)


class ConfigTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()

    def write(self, name, payload):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(yaml.safe_dump(payload), encoding="utf-8")
        return path

    def paired(self, **changes):
        return {
            "schema_version": 1, "mode": "paired", "source": "../data/all.h5ad",
            "pair_column": "pair_id", "pair_key": " 001 ", "output_dir": "../output",
            **changes,
        }

    def screen(self, candidates):
        # The core must not parse, compile, or rewrite this host-owned file.
        (self.root / "original.yaml").write_text("host-owned opaque contents: [", encoding="utf-8")
        return self.write("prepare.yaml", {
            "schema_version": 1, "mode": "screen",
            "reconstruction_config": "original.yaml", "candidates": candidates,
            "output_dir": "output",
        })

    def test_paired_paths_are_relative_to_declaring_file_and_key_is_exact(self):
        path = self.write("configs/prepare.yaml", self.paired())
        config = load_preparation_config(path)
        self.assertIsInstance(config, PairedConfig)
        self.assertEqual(config.source, self.root / "data/all.h5ad")
        self.assertEqual(config.output_dir, self.root / "output")
        self.assertEqual(config.pair_key, " 001 ")
        whitespace_key = load_preparation_config(self.write("space.yaml", self.paired(pair_key=" ")))
        self.assertEqual(whitespace_key.pair_key, " ")

    def test_directory_discovers_only_immediate_h5ad_files_in_stable_order(self):
        directory = self.root / "candidates"
        directory.mkdir()
        for name in ("z.h5ad", "a.h5ad", "b.H5AD", "notes.txt"):
            (directory / name).touch()
        (directory / "nested").mkdir()
        (directory / "nested/hidden.h5ad").touch()
        config = load_preparation_config(self.screen({"directory": "candidates"}))
        self.assertEqual([c.candidate_id for c in config.candidates], ["a.h5ad", "b.H5AD", "z.h5ad"])
        self.assertEqual(config.reconstruction_config, self.root / "original.yaml")

    def test_list_paths_resolve_against_list_and_missing_file_is_retained(self):
        self.write("lists/references.yaml", {
            "schema_version": 1,
            "candidates": [{"id": "external", "path": "../data/missing.h5ad"}],
        })
        config = load_preparation_config(self.screen({"list": "lists/references.yaml"}))
        self.assertEqual(config.candidates[0].path, self.root / "data/missing.h5ad")

    def test_duplicate_ids_or_resolved_paths_are_rejected(self):
        cases = [
            [{"id": "A", "path": "a.h5ad"}, {"id": "A", "path": "b.h5ad"}],
            [{"id": "A", "path": "a.h5ad"}, {"id": "B", "path": "./a.h5ad"}],
        ]
        for entries in cases:
            with self.subTest(entries=entries):
                self.write("candidates.yaml", {"schema_version": 1, "candidates": entries})
                with self.assertRaisesRegex(ConfigError, "unique"):
                    load_preparation_config(self.screen({"list": "candidates.yaml"}))

    def test_symlink_aliases_cannot_score_one_file_twice(self):
        directory = self.root / "candidates"
        directory.mkdir()
        (directory / "a.h5ad").touch()
        (directory / "alias.h5ad").symlink_to(directory / "a.h5ad")
        with self.assertRaisesRegex(ConfigError, "paths must be unique"):
            load_preparation_config(self.screen({"directory": "candidates"}))

    def test_rejects_ambiguous_or_unknown_preparation_fields(self):
        cases = [
            self.paired(pair_key=1),
            self.paired(schema_version=True),
            self.paired(schema_version=2),
            self.paired(mode="automatic"),
            self.paired(patient_id="inferred"),
            self.paired(pair_column=""),
        ]
        for document in cases:
            with self.subTest(document=document), self.assertRaises(ConfigError):
                load_preparation_config(self.write("bad.yaml", document))
        with self.assertRaisesRegex(ConfigError, "exactly one"):
            load_preparation_config(self.screen({"directory": "d", "list": "l"}))

    def test_duplicate_yaml_keys_fail_instead_of_silently_overriding(self):
        path = self.root / "duplicate.yaml"
        path.write_text("schema_version: 1\nmode: paired\nmode: screen\n", encoding="utf-8")
        with self.assertRaisesRegex(ConfigError, "Duplicate YAML key"):
            load_preparation_config(path)

    def test_empty_candidate_sources_are_valid_requests(self):
        (self.root / "empty").mkdir()
        self.assertEqual(load_preparation_config(self.screen({"directory": "empty"})).candidates, ())
        self.write("empty.yaml", {"schema_version": 1, "candidates": []})
        self.assertEqual(load_preparation_config(self.screen({"list": "empty.yaml"})).candidates, ())

    def test_resolved_reference_round_trip_is_small_and_relocatable(self):
        data = self.root / "data"
        data.mkdir()
        reference = data / "chosen.h5ad"
        reference.touch()
        output = self.root / "run"
        output.mkdir()
        path = output / "reference.yaml"
        write_reference_config(path, reference)
        document = yaml.safe_load(path.read_text(encoding="utf-8"))
        self.assertEqual(document, {
            "schema_version": 1, "reference": {"path": "../data/chosen.h5ad", "format": "h5ad"},
        })
        self.assertEqual(load_reference_config(path).path, reference)

    def test_reference_config_rejects_filters_original_params_and_missing_files(self):
        reference = self.root / "chosen.h5ad"
        reference.touch()
        base = {"schema_version": 1, "reference": {"path": "chosen.h5ad", "format": "h5ad"}}
        cases = [
            {**base, "mode": "screen"},
            {**base, "reference": {**base["reference"], "filter_column": "old"}},
            {**base, "reference": {"path": "absent.h5ad", "format": "h5ad"}},
            {**base, "reference": {"path": "chosen.h5ad", "format": "csv"}},
        ]
        for document in cases:
            with self.subTest(document=document), self.assertRaises(ConfigError):
                load_reference_config(self.write("bad-reference.yaml", document))

    def test_output_writers_never_overwrite_and_clean_temporary_files(self):
        path = self.root / "report.json"
        path.write_text("keep manual result", encoding="utf-8")
        with self.assertRaises(FileExistsError):
            atomic_write_json(path, {"new": "result"})
        self.assertEqual(path.read_text(encoding="utf-8"), "keep manual result")
        self.assertEqual(list(self.root.iterdir()), [path])
        with self.assertRaises(ValueError):
            atomic_write_json(self.root / "invalid.json", {"score": float("nan")})
        self.assertFalse((self.root / "invalid.json").exists())


if __name__ == "__main__":
    unittest.main()
