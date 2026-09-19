from __future__ import annotations

import unittest

import numpy as np

from revise.reference_preparation.ga_contract import (
    GAResponse,
    GlobalAnchoringResult,
    validate_global_anchoring_result,
)


class GlobalAnchoringContractTest(unittest.TestCase):
    def test_preserves_original_invalid_probability_and_shape_cases(self):
        invalid = [
            np.array([.5, .5]),
            np.array([[.5, .3, .2]]),
            np.array([[1.1, -.1]]),
            np.array([[8., 2.]]),
            np.array([[float("nan"), .5]]),
            np.array([[float("inf"), .5]]),
        ]
        for matrix in invalid:
            with self.subTest(matrix=matrix), self.assertRaises(ValueError):
                validate_global_anchoring_result(
                    GlobalAnchoringResult(matrix, ["u1"], ["A", "B"], {}), ["u1"],
                )
        with self.assertRaisesRegex(ValueError, "unique"):
            validate_global_anchoring_result(
                GlobalAnchoringResult(np.array([[.5, .5]]), ["u1"], ["A", "A"], {}), ["u1"],
            )
        with self.assertRaisesRegex(ValueError, "unique"):
            validate_global_anchoring_result(
                GlobalAnchoringResult(np.ones((2, 2)) / 2, ["u1", "u1"], ["A", "B"], {}),
                ["u1", "u2"],
            )

    def test_validates_separate_expected_axis_and_copies_metadata(self):
        result = GlobalAnchoringResult(
            distribution=np.array([[0.8, 0.2], [0.3, 0.7]]),
            st_unit_ids=["u1", "u2"],
            cell_type_labels=["A", "B"],
            metadata={"engine": "test"},
        )
        response = GAResponse(result=result, expected_st_unit_ids=["u1", "u2"])

        validated = validate_global_anchoring_result(
            response.result, response.expected_st_unit_ids
        )

        self.assertEqual(validated.metadata["normalization"], "adapter_declared_probability")
        self.assertNotIn("normalization", result.metadata)
        self.assertIsInstance(validated.distribution, np.ndarray)

    def test_rejects_invalid_result_and_empty_or_mismatched_axes(self):
        valid = GlobalAnchoringResult(
            distribution=np.array([[0.8, 0.2]]),
            st_unit_ids=["u1"],
            cell_type_labels=["A", "B"],
            metadata={},
        )
        cases = [
            (object(), ["u1"], TypeError, "GlobalAnchoringResult"),
            (
                GlobalAnchoringResult(
                    np.array([[0.8, 0.2]]), [], ["A", "B"], {}
                ),
                ["u1"],
                ValueError,
                "nonempty",
            ),
            (valid, [""], ValueError, "nonempty string"),
            (valid, ["other"], ValueError, "input order"),
            (
                GlobalAnchoringResult(
                    np.array([[0.8, 0.2]]), ["u1"], ["A", ""], {}
                ),
                ["u1"],
                ValueError,
                "nonempty string",
            ),
            (
                GlobalAnchoringResult(
                    np.array([[0.8, 0.2]]), ["u1"], ["A"], {}
                ),
                ["u1"],
                ValueError,
                "at least two",
            ),
        ]
        for result, expected, error_type, message in cases:
            with self.subTest(message=message):
                with self.assertRaisesRegex(error_type, message):
                    validate_global_anchoring_result(result, expected)

    def test_retains_atlas_probability_validation_tolerance(self):
        within_tolerance = GlobalAnchoringResult(
            distribution=np.array([[0.800004, 0.2]]),
            st_unit_ids=["u1"],
            cell_type_labels=["A", "B"],
            metadata={},
        )
        outside_tolerance = GlobalAnchoringResult(
            distribution=np.array([[0.8001, 0.2]]),
            st_unit_ids=["u1"],
            cell_type_labels=["A", "B"],
            metadata={},
        )

        validate_global_anchoring_result(within_tolerance, ["u1"])
        with self.assertRaisesRegex(ValueError, "row-normalized"):
            validate_global_anchoring_result(outside_tolerance, ["u1"])


if __name__ == "__main__":
    unittest.main()
