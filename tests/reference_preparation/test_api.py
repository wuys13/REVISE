import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest import mock

import anndata as ad
import numpy as np
import pandas as pd
import yaml

from revise.reference_preparation import (
    NoUsableReferenceError,
    load_reference_config,
    prepare_reference,
)
from revise.reference_preparation.ga_contract import GAResponse, GlobalAnchoringResult


class ApiTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()

    def write_yaml(self, name, document):
        path = self.root / name
        path.write_text(yaml.safe_dump(document), encoding="utf-8")
        return path

    def reference(self, name):
        path = self.root / name
        data = ad.AnnData(
            X=np.array([[1., 2.], [3., 4.], [5., 6.]]),
            obs=pd.DataFrame({"pair_id": ["A", "B", "A"]}, index=["c1", "c2", "c3"]),
            var=pd.DataFrame(index=["g1", "g2"]),
        )
        data.write_h5ad(path)
        return path

    def paired_request(self, **changes):
        return self.write_yaml("prepare.yaml", {
            "schema_version": 1, "mode": "paired", "source": "all.h5ad",
            "pair_column": "pair_id", "pair_key": "A", "output_dir": "run",
            **changes,
        })

    def screen_request(self, entries):
        original = self.write_yaml("original.yaml", {
            "inputs": {"reference": {"path": "old.h5ad", "filter_column": "Patient", "filter_value": "old"}},
            "host_method": "unchanged",
        })
        self.write_yaml("candidates.yaml", {"schema_version": 1, "candidates": entries})
        request = self.write_yaml("prepare.yaml", {
            "schema_version": 1, "mode": "screen",
            "reconstruction_config": "original.yaml",
            "candidates": {"list": "candidates.yaml"}, "output_dir": "run",
        })
        return request, original

    @staticmethod
    def response(probability=.9):
        return GAResponse(
            GlobalAnchoringResult(
                np.array([[probability, 1-probability], [probability, 1-probability]]),
                ["u1", "u2"], ["type1", "type2"], {},
            ),
            ["u1", "u2"],
        )

    def test_paired_end_to_end_and_copied_result_load_without_original_source(self):
        source = self.reference("all.h5ad")
        source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
        runner = mock.Mock(side_effect=AssertionError("paired must not call GA"))
        result = prepare_reference(self.paired_request(), ga_runner=runner)
        runner.assert_not_called()
        self.assertEqual(load_reference_config(result.reference_config_path).path, result.reference_path)
        selected = ad.read_h5ad(result.reference_path)
        self.assertEqual(list(selected.obs_names), ["c1", "c3"])
        np.testing.assert_array_equal(selected.X, [[1, 2], [5, 6]])
        self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), source_hash)
        report = json.loads(result.report_path.read_text(encoding="utf-8"))
        self.assertEqual(report["selected_cells"], 2)
        copied = self.root / "relocated"
        shutil.copytree(result.reference_path.parent, copied)
        source.unlink()
        self.assertEqual(load_reference_config(copied / "reference.yaml").path, copied / "reference.h5ad")

    def test_partial_failure_selects_success_without_modifying_base_or_candidates(self):
        good = self.reference("good.h5ad")
        bad = self.reference("bad.h5ad")
        request, original = self.screen_request([
            {"id": "bad", "path": bad.name}, {"id": "good", "path": good.name},
        ])
        before = {p: p.read_bytes() for p in (good, bad, original, request)}
        calls = []

        def runner(reference, reconstruction_config):
            calls.append((reference, reconstruction_config))
            if reference == bad:
                raise ValueError("missing annotation")
            return self.response()

        result = prepare_reference(request, ga_runner=runner)
        self.assertEqual(calls, [(bad, original), (good, original)])
        self.assertEqual(result.reference_path, good)
        self.assertFalse((result.report_path.parent / "reference.h5ad").exists())
        self.assertEqual(load_reference_config(result.reference_config_path).path, good)
        self.assertEqual({p: p.read_bytes() for p in before}, before)
        report = json.loads(result.report_path.read_text(encoding="utf-8"))
        self.assertEqual(report["selected_candidate_id"], "good")
        self.assertEqual(report["candidates"][0]["status"], "failed")
        self.assertEqual(set(report["candidates"][1]["scores"]), {
            "max_median", "certainty_median", "margin_median",
        })
        resolved = yaml.safe_load(result.reference_config_path.read_text(encoding="utf-8"))
        self.assertEqual(set(resolved["reference"]), {"path", "format"})

    def test_missing_ga_fails_before_outputs_even_for_empty_candidates(self):
        request, _ = self.screen_request([])
        with self.assertRaisesRegex(ValueError, "ga_runner"):
            prepare_reference(request)
        self.assertFalse((self.root / "run").exists())

    def test_all_failures_leave_report_and_no_reference_configuration(self):
        request, _ = self.screen_request([{"id": "missing", "path": "absent.h5ad"}])
        runner = mock.Mock(side_effect=AssertionError("absent candidate must not reach GA"))
        with self.assertRaises(NoUsableReferenceError) as caught:
            prepare_reference(request, ga_runner=runner)
        runner.assert_not_called()
        report = json.loads(caught.exception.report_path.read_text(encoding="utf-8"))
        self.assertEqual(report["status"], "no_reference")
        self.assertEqual(report["candidates"][0]["status"], "failed")
        self.assertFalse((self.root / "run/reference.yaml").exists())

    def test_empty_candidates_are_inspectable_no_reference_result(self):
        request, _ = self.screen_request([])
        runner = mock.Mock()
        with self.assertRaises(NoUsableReferenceError) as caught:
            prepare_reference(request, ga_runner=runner)
        runner.assert_not_called()
        report = json.loads(caught.exception.report_path.read_text(encoding="utf-8"))
        self.assertEqual(report["candidates"], [])
        self.assertFalse((self.root / "run/reference.yaml").exists())

    def test_existing_output_directory_cannot_be_overwritten(self):
        self.reference("all.h5ad")
        request = self.paired_request()
        result = prepare_reference(request)
        before = {p: p.read_bytes() for p in result.report_path.parent.iterdir()}
        with self.assertRaises(FileExistsError):
            prepare_reference(request)
        self.assertEqual({p: p.read_bytes() for p in before}, before)

    def test_paired_failure_is_reported_without_a_resolved_config(self):
        self.reference("all.h5ad")
        with self.assertRaises(ValueError):
            prepare_reference(self.paired_request(pair_key="missing"))
        self.assertFalse((self.root / "run/reference.yaml").exists())
        self.assertFalse((self.root / "run/reference.h5ad").exists())
        report = json.loads((self.root / "run/report.json").read_text(encoding="utf-8"))
        self.assertEqual(report["status"], "failed")

    def test_report_write_failure_never_publishes_reference_config(self):
        self.reference("all.h5ad")
        request = self.paired_request()
        with mock.patch("revise.reference_preparation.api.atomic_write_json", side_effect=OSError("disk full")):
            with self.assertRaisesRegex(OSError, "disk full"):
                prepare_reference(request)
        self.assertFalse((self.root / "run/reference.yaml").exists())

    def test_final_config_write_failure_updates_the_owned_report(self):
        self.reference("all.h5ad")
        request = self.paired_request()
        with mock.patch("revise.reference_preparation.api.write_reference_config", side_effect=PermissionError("publish denied")):
            with self.assertRaisesRegex(PermissionError, "publish denied"):
                prepare_reference(request)
        self.assertFalse((self.root / "run/reference.yaml").exists())
        report = json.loads((self.root / "run/report.json").read_text(encoding="utf-8"))
        self.assertEqual(report["status"], "failed")
        self.assertIn("publish denied", report["failure_reason"])


if __name__ == "__main__":
    unittest.main()
