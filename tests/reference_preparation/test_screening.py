from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np

from revise.reference_preparation.config import Candidate
from revise.reference_preparation.ga_contract import GAResponse, GlobalAnchoringResult
from revise.reference_preparation.screening import screen_candidates


class ScreeningTest(unittest.TestCase):
    def test_invalid_ga_isolated_and_stronger_success_is_selected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidates = [Candidate(name, root / f"{name}.h5ad") for name in ("weak", "invalid", "strong")]
            for candidate in candidates:
                candidate.path.touch()
            attempted = []

            def runner(path, config):
                attempted.append(path.stem)
                probability = {"weak": .6, "invalid": .9, "strong": .95}[path.stem]
                result = GlobalAnchoringResult(
                    np.array([[probability, 1-probability]]), ["unit"], ["A", "B"], {},
                )
                expected = ["wrong-order"] if path.stem == "invalid" else ["unit"]
                return GAResponse(result, expected)

            report = screen_candidates(candidates, root / "base.yaml", runner)
            self.assertEqual(attempted, ["weak", "invalid", "strong"])
            self.assertEqual(report["selected_candidate_id"], "strong")
            self.assertEqual(report["candidates"][0]["ranks"]["max_median"], 2)
            self.assertEqual(report["candidates"][1]["status"], "failed")
            self.assertIn("input order", report["candidates"][1]["failure_reason"])
            self.assertAlmostEqual(report["candidates"][2]["scores"]["max_median"], .95)

    def test_attempts_every_candidate_and_ranks_ties_by_candidate_id(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            alpha = directory / "alpha.h5ad"
            beta = directory / "beta.h5ad"
            zeta = directory / "zeta.h5ad"
            alpha.touch()
            beta.touch()
            zeta.touch()
            reconstruction = directory / "reconstruction.yaml"
            reconstruction.touch()
            calls: list[Path] = []

            def runner(path: Path, config: Path) -> GAResponse:
                self.assertEqual(config, reconstruction)
                calls.append(path)
                if path == beta:
                    raise RuntimeError("host unavailable")
                return GAResponse(
                    result=GlobalAnchoringResult(
                        distribution=np.array([[0.9, 0.1], [0.8, 0.2]]),
                        st_unit_ids=["s1", "s2"],
                        cell_type_labels=["A", "B"],
                        metadata={},
                    ),
                    expected_st_unit_ids=["s1", "s2"],
                )

            report = screen_candidates(
                [
                    Candidate("z-missing", directory / "missing.h5ad"),
                    Candidate("beta", beta),
                    Candidate("zeta", zeta),
                    Candidate("alpha", alpha),
                ],
                reconstruction,
                runner,
            )

            self.assertEqual(calls, [beta, zeta, alpha])
            self.assertEqual(report["status"], "selected")
            self.assertEqual(report["selected_candidate_id"], "alpha")
            self.assertEqual(report["selected_reference_path"], str(alpha))
            self.assertEqual(
                [row["status"] for row in report["candidates"]],
                ["failed", "failed", "success", "success"],
            )
            self.assertEqual(report["candidates"][3]["ranks"], {
                "max_median": 1,
                "certainty_median": 1,
                "margin_median": 1,
            })
            self.assertEqual(report["candidates"][2]["ranks"], {
                "max_median": 2,
                "certainty_median": 2,
                "margin_median": 2,
            })
            self.assertIn("FileNotFoundError", report["candidates"][0]["failure_reason"])
            self.assertIn("RuntimeError", report["candidates"][1]["failure_reason"])

    def test_no_successful_candidate_returns_no_reference_report(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            report = screen_candidates(
                [Candidate("missing", directory / "missing.h5ad")],
                directory / "reconstruction.yaml",
                lambda *_: self.fail("runner must not receive a missing file"),
            )

            self.assertEqual(report["status"], "no_reference")
            self.assertIsNone(report["selected_candidate_id"])
            self.assertIsNone(report["selected_reference_path"])
            self.assertTrue(
                all(value is None for value in report["candidates"][0]["ranks"].values())
            )


if __name__ == "__main__":
    unittest.main()
