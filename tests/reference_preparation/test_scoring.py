from __future__ import annotations

import unittest

import numpy as np

from revise.reference_preparation.scoring import SCORE_METHODS, row_metrics, summarize_scores


class ScoringTest(unittest.TestCase):
    def test_preserves_original_atlas_scoring_example(self):
        matrix = np.array([[0.9, 0.1], [0.6, 0.4]])
        metrics = row_metrics(matrix)
        scores = summarize_scores(metrics)
        np.testing.assert_allclose(metrics["max_confidence"], [0.9, 0.6])
        np.testing.assert_allclose(metrics["top2_margin"], [0.8, 0.2])
        self.assertAlmostEqual(scores["max_median"], 0.75)
        self.assertAlmostEqual(scores["margin_median"], 0.5)
        expected_certainty = [
            1 + (0.9 * np.log(0.9) + 0.1 * np.log(0.1)) / np.log(2),
            1 + (0.6 * np.log(0.6) + 0.4 * np.log(0.4)) / np.log(2),
        ]
        self.assertAlmostEqual(scores["certainty_median"], float(np.median(expected_certainty)))

    def test_zero_probabilities_and_uniform_distribution_have_defined_scores(self):
        with np.errstate(divide="raise", invalid="raise"):
            metrics = row_metrics(np.array([[1., 0., 0.], [1/3, 1/3, 1/3]]))
        np.testing.assert_allclose(metrics["normalized_certainty"], [1., 0.], atol=1e-15)
        scores = summarize_scores(metrics)
        self.assertAlmostEqual(scores["max_median"], 2/3)
        self.assertAlmostEqual(scores["certainty_median"], .5)
        self.assertAlmostEqual(scores["margin_median"], .5)

    def test_row_metrics_and_median_summary_match_ranking_contract(self):
        matrix = np.array([[0.8, 0.2], [0.5, 0.5], [0.1, 0.9]])

        metrics = row_metrics(matrix)

        np.testing.assert_allclose(metrics["max_confidence"], [0.8, 0.5, 0.9])
        np.testing.assert_allclose(metrics["top2_margin"], [0.6, 0.0, 0.8])
        self.assertEqual(set(SCORE_METHODS), set(summarize_scores(metrics)))
        self.assertEqual(summarize_scores(metrics)["max_median"], 0.8)


if __name__ == "__main__":
    unittest.main()
