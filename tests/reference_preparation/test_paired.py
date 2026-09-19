from __future__ import annotations

from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd

from revise.reference_preparation.paired import extract_pair


class PairedExtractionTest(unittest.TestCase):
    def _write_source(self, directory: Path) -> Path:
        source = ad.AnnData(
            X=np.array([[1, 2], [3, 4], [5, 6]], dtype=float),
            obs=pd.DataFrame(
                {"pair": [" A", "A", "A"]}, index=["c1", "c2", "c3"]
            ),
            var=pd.DataFrame({"feature": ["x", "y"]}, index=["g1", "g2"]),
            layers={"counts": np.array([[1, 2], [3, 4], [5, 6]])},
            obsm={"embedding": np.array([[1, 0], [0, 1], [1, 1]])},
        )
        source.raw = source.copy()
        path = directory / "source.h5ad"
        source.write_h5ad(path)
        return path

    def test_extracts_exact_pair_and_preserves_anndata_content(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = self._write_source(directory)
            output = directory / "reference.h5ad"

            stats = extract_pair(source, "pair", "A", output)
            selected = ad.read_h5ad(output)

            self.assertEqual(list(selected.obs_names), ["c2", "c3"])
            self.assertEqual(list(selected.var_names), ["g1", "g2"])
            np.testing.assert_array_equal(selected.layers["counts"], [[3, 4], [5, 6]])
            np.testing.assert_array_equal(selected.obsm["embedding"], [[0, 1], [1, 1]])
            self.assertIsNotNone(selected.raw)
            self.assertEqual(stats, {
                "source_path": str(source),
                "pair_column": "pair",
                "pair_key": "A",
                "source_cells": 3,
                "selected_cells": 2,
                "genes": 2,
            })

    def test_missing_or_zero_match_does_not_publish_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = self._write_source(directory)
            missing_output = directory / "missing.h5ad"
            zero_output = directory / "zero.h5ad"

            with self.assertRaisesRegex(KeyError, "pair column"):
                extract_pair(source, "missing", "A", missing_output)
            with self.assertRaisesRegex(ValueError, "exactly match"):
                extract_pair(source, "pair", "missing", zero_output)

            self.assertFalse(missing_output.exists())
            self.assertFalse(zero_output.exists())

    def test_serialization_failure_cleans_up_temporary_and_final_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = self._write_source(directory)
            output = directory / "reference.h5ad"

            with patch.object(ad.AnnData, "write_h5ad", side_effect=OSError("write failed")):
                with self.assertRaisesRegex(OSError, "write failed"):
                    extract_pair(source, "pair", "A", output)

            self.assertFalse(output.exists())
            self.assertEqual(list(directory.glob(".reference.h5ad.*")), [])


if __name__ == "__main__":
    unittest.main()
