from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd
from anndata import AnnData
import yaml

import revise.preprocess.sim2real_pseudospot as pseudospot
from revise.preprocess.sim2real_pseudospot import aggregation
from revise.preprocess.sim2real_pseudospot import cli
from revise.preprocess.sim2real_pseudospot import regions
from revise.preprocess.sim2real_pseudospot import workflow


ROOT = Path(__file__).resolve().parents[2]


class Sim2RealPseudospotContractTest(unittest.TestCase):
    def test_is_a_packaged_cli_capability(self):
        module = ROOT / "revise/preprocess/sim2real_pseudospot"
        pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")

        self.assertTrue(module.is_dir())
        self.assertTrue((module / "workflow.py").is_file())
        self.assertTrue((module / "regions.py").is_file())
        self.assertTrue((module / "aggregation.py").is_file())
        self.assertTrue((module / "cli.py").is_file())
        self.assertIn(
            'revise-prepare-sim2real-pseudospots = '
            '"revise.preprocess.sim2real_pseudospot.cli:main"',
            pyproject,
        )

    def test_package_exports_the_two_workflow_operations_and_cli_stages(self):
        self.assertTrue(hasattr(pseudospot, "propose_regions"))
        self.assertTrue(hasattr(pseudospot, "build_real_pseudospots"))
        self.assertTrue(hasattr(cli, "build_parser"))

        parser = cli.build_parser()
        proposed = parser.parse_args(
            ["propose", "--config", "config.yaml", "--sample", "P1CRC"]
        )
        built = parser.parse_args(
            [
                "build",
                "--config",
                "config.yaml",
                "--sample",
                "P1CRC",
                "--confirmation",
                "confirmed_regions.yaml",
            ]
        )

        self.assertEqual(proposed.stage, "propose")
        self.assertEqual(built.stage, "build")

    def test_real_aggregation_sums_cells_and_keeps_empty_spots_only_in_distribution(self):
        self.assertTrue(hasattr(aggregation, "aggregate_real_cells"))

        cells = AnnData(
            X=np.array([[1, 2], [3, 4], [5, 6], [7, 8]], dtype=float),
            obs=pd.DataFrame(
                {"cell_id": ["c1", "c2", "c3", "c4"]},
                index=["c1", "c2", "c3", "c4"],
            ),
            var=pd.DataFrame(index=["g1", "g2"]),
        )
        cells.obsm["spatial"] = np.array(
            [[0.0, 0.0], [1.0, 0.0], [11.0, 0.0], [0.0, 11.0]]
        )

        spots, distribution = aggregation.aggregate_real_cells(cells, spot_size=10)

        self.assertEqual(spots.obs_names.tolist(), ["SPOT_0_0", "SPOT_0_1", "SPOT_1_0"])
        np.testing.assert_array_equal(
            spots.X.toarray(), np.array([[4, 6], [7, 8], [5, 6]], dtype=float)
        )
        self.assertEqual(spots.uns["all_cells_in_spot"]["SPOT_1_1"], None)
        self.assertEqual(spots.uns["all_cells_in_spot"]["SPOT_0_0"], ["c1", "c2"])
        self.assertEqual(distribution["count"].to_dict(), {0: 1, 1: 2, 2: 1})
        self.assertEqual(list(spots.layers.keys()), [])

    def test_config_resolves_paths_relative_to_its_file_and_rejects_invalid_sizes(self):
        self.assertTrue(hasattr(workflow, "load_config"))

        with self.subTest("relative paths"):
            with tempfile.TemporaryDirectory() as temporary_directory:
                config_dir = Path(temporary_directory)
                config_path = config_dir / "config.yaml"
                config_path.write_text(
                    yaml.safe_dump(
                        {
                            "schema_version": 1,
                            "reference": {"path": "sc.h5ad", "patient_column": "Patient"},
                            "samples": {
                                "P1CRC": {
                                    "xenium_path": "P1.h5ad",
                                    "output_dir": "out/P1",
                                }
                            },
                            "preprocessing": {
                                "transcript_counts_min": 60,
                                "gene_min_cells": 100,
                                "label_key": "Level1",
                                "unknown_label": "Unknown",
                                "min_cells_per_type": 60,
                            },
                            "template_root": "templates",
                            "proposal": {
                                "base_width": 2740,
                                "base_height": 1370,
                                "scales": [0.75, 1.0, 1.25],
                                "step": 250,
                                "min_cells": 1000,
                                "max_iou": 0.25,
                            },
                            "spot_sizes": [50, 100, 150, 200],
                            "seed": 42,
                        }
                    ),
                    encoding="utf-8",
                )

                config = workflow.load_config(config_path)

                self.assertEqual(
                    config.samples["P1CRC"].xenium_path,
                    (config_dir / "P1.h5ad").resolve(),
                )
                self.assertEqual(
                    config.samples["P1CRC"].output_dir,
                    (config_dir / "out/P1").resolve(),
                )

                invalid = yaml.safe_load(config_path.read_text(encoding="utf-8"))
                invalid["spot_sizes"] = [50, 20]
                config_path.write_text(yaml.safe_dump(invalid), encoding="utf-8")

                with self.assertRaisesRegex(ValueError, "spot_sizes"):
                    workflow.load_config(config_path)

    def test_region_candidates_filter_labels_before_iou_suppression(self):
        self.assertTrue(hasattr(regions, "filter_region_cells"))
        self.assertTrue(hasattr(regions, "Candidate"))
        self.assertTrue(hasattr(regions, "suppress_overlaps"))

        cells = AnnData(
            X=np.ones((4, 1)),
            obs=pd.DataFrame(
                {"Level1": ["Tumor", "Tumor", "Unknown", "Immune"]},
                index=["c1", "c2", "c3", "c4"],
            ),
            var=pd.DataFrame(index=["g1"]),
        )
        cells.obsm["spatial"] = np.array(
            [[1.0, 1.0], [2.0, 2.0], [3.0, 3.0], [4.0, 4.0]]
        )

        retained = regions.filter_region_cells(
            cells,
            bounds=(0.0, 10.0, 0.0, 10.0),
            label_key="Level1",
            unknown_label="Unknown",
            min_cells_per_type=1,
        )

        self.assertEqual(retained.obs_names.tolist(), ["c1", "c2"])

        candidates = [
            regions.Candidate("normal_core-1", "normal_core", (0, 10, 0, 10), 0.9),
            regions.Candidate("normal_core-2", "normal_core", (1, 11, 0, 10), 0.8),
            regions.Candidate("normal_core-3", "normal_core", (20, 30, 0, 10), 0.7),
        ]
        selected = regions.suppress_overlaps(candidates, max_iou=0.25, limit=3)

        self.assertEqual([candidate.candidate_id for candidate in selected], [
            "normal_core-1",
            "normal_core-3",
        ])

    def test_candidate_search_returns_three_deterministic_nonoverlapping_proposals(self):
        cells = AnnData(
            X=np.ones((15, 1)),
            obs=pd.DataFrame({"Level1": ["Tumor"] * 15}, index=[f"c{index}" for index in range(15)]),
            var=pd.DataFrame(index=["g1"]),
        )
        cells.obsm["spatial"] = np.array(
            [[5.0 + 20.0 * group, 5.0] for group in range(3) for _ in range(5)]
        )
        proposal = workflow.ProposalConfig(
            base_width=20,
            base_height=20,
            scales=(1.0,),
            step=20,
            min_cells=2,
            max_iou=0.1,
        )

        first = regions.propose_candidates(
            cells,
            role="tumor_core",
            template_composition={"Tumor": 1.0},
            config=proposal,
            label_key="Level1",
            unknown_label="Unknown",
            min_cells_per_type=1,
        )
        second = regions.propose_candidates(
            cells,
            role="tumor_core",
            template_composition={"Tumor": 1.0},
            config=proposal,
            label_key="Level1",
            unknown_label="Unknown",
            min_cells_per_type=1,
        )

        self.assertEqual(len(first), 3)
        self.assertEqual(
            [(candidate.candidate_id, candidate.bounds) for candidate in first],
            [(candidate.candidate_id, candidate.bounds) for candidate in second],
        )

    def test_build_requires_a_confirmed_region_manifest_before_reading_data(self):
        self.assertTrue(hasattr(workflow, "build_real_pseudospots"))

        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            config_path = root / "config.yaml"
            config_path.write_text(
                yaml.safe_dump(
                    {
                        "schema_version": 1,
                        "reference": {"path": "sc.h5ad", "patient_column": "Patient"},
                        "samples": {
                            "P1CRC": {
                                "xenium_path": "P1.h5ad",
                                "output_dir": "out/P1",
                            }
                        },
                        "preprocessing": {
                            "transcript_counts_min": 60,
                            "gene_min_cells": 100,
                            "label_key": "Level1",
                            "unknown_label": "Unknown",
                            "min_cells_per_type": 60,
                        },
                        "template_root": "templates",
                        "proposal": {
                            "base_width": 2740,
                            "base_height": 1370,
                            "scales": [0.75, 1.0, 1.25],
                            "step": 250,
                            "min_cells": 1000,
                            "max_iou": 0.25,
                        },
                        "spot_sizes": [50, 100, 150, 200],
                        "seed": 42,
                    }
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(FileNotFoundError, "confirmation"):
                workflow.build_real_pseudospots(
                    config_path, "P1CRC", root / "missing-confirmed-regions.yaml"
                )

    def test_build_writes_only_real_standard_output_contract(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            labels = []
            coordinates = []
            cell_ids = []
            for index, label in enumerate(["Tumor", "Immune", "Fibroblast"]):
                for cell in range(61):
                    labels.append(label)
                    coordinates.append([index * 100 + cell % 10, cell // 10])
                    cell_ids.append(f"{label}-{cell}")
            annotated = AnnData(
                X=np.ones((len(cell_ids), 2)),
                obs=pd.DataFrame(
                    {"cell_id": cell_ids, "Level1": labels}, index=cell_ids
                ),
                var=pd.DataFrame(index=["g1", "g2"]),
            )
            annotated.obsm["spatial"] = np.asarray(coordinates, dtype=float)
            region_dir = root / "out/P1/region_selection"
            region_dir.mkdir(parents=True)
            annotated_path = region_dir / "annotated_xenium.h5ad"
            annotated.write_h5ad(annotated_path)

            reference = AnnData(
                X=np.ones((9, 2)),
                obs=pd.DataFrame(
                    {
                        "Patient": ["P1CRC"] * 9,
                        "Level1": ["Tumor"] * 3 + ["Immune"] * 3 + ["Fibroblast"] * 3,
                    },
                    index=[f"ref-{index}" for index in range(9)],
                ),
                var=pd.DataFrame(index=["g1", "g2"]),
            )
            reference.write_h5ad(root / "sc.h5ad")
            proposal = {
                "schema_version": 1,
                "sample": "P1CRC",
                "annotated_xenium": str(annotated_path),
                "regions": {
                    "leading_edge": [{"candidate_id": "leading_edge-1", "bounds": [-1, 10, -1, 7]}],
                    "normal_core": [{"candidate_id": "normal_core-1", "bounds": [99, 110, -1, 7]}],
                    "tumor_core": [{"candidate_id": "tumor_core-1", "bounds": [199, 210, -1, 7]}],
                },
            }
            (region_dir / "proposal.yaml").write_text(yaml.safe_dump(proposal), encoding="utf-8")
            confirmation_path = region_dir / "confirmed_regions.yaml"
            confirmation_path.write_text(
                yaml.safe_dump(
                    {
                        "schema_version": 1,
                        "sample": "P1CRC",
                        "regions": {
                            "leading_edge": {"candidate_id": "leading_edge-1"},
                            "normal_core": {"candidate_id": "normal_core-1"},
                            "tumor_core": {"candidate_id": "tumor_core-1"},
                        },
                    }
                ),
                encoding="utf-8",
            )
            config_path = root / "config.yaml"
            config_path.write_text(
                yaml.safe_dump(
                    {
                        "schema_version": 1,
                        "reference": {"path": "sc.h5ad", "patient_column": "Patient"},
                        "samples": {"P1CRC": {"xenium_path": "unused.h5ad", "output_dir": "out/P1"}},
                        "preprocessing": {
                            "transcript_counts_min": 60,
                            "gene_min_cells": 100,
                            "label_key": "Level1",
                            "unknown_label": "Unknown",
                            "min_cells_per_type": 1,
                        },
                        "template_root": "templates",
                        "proposal": {
                            "base_width": 2740,
                            "base_height": 1370,
                            "scales": [0.75, 1.0, 1.25],
                            "step": 250,
                            "min_cells": 1000,
                            "max_iou": 0.25,
                        },
                        "spot_sizes": [50, 100, 150, 200],
                        "seed": 42,
                    }
                ),
                encoding="utf-8",
            )

            result = workflow.build_real_pseudospots(config_path, "P1CRC", confirmation_path)

            self.assertTrue(result.output_dir.is_dir())
            for part in ("part1", "part2", "part3"):
                part_dir = result.output_dir / part
                self.assertTrue((part_dir / "selected_xenium.h5ad").is_file())
                self.assertTrue((part_dir / "real_sc_ref_part.h5ad").is_file())
                self.assertTrue((part_dir / "cut.png").is_file())
                for size in (50, 100, 150, 200):
                    self.assertTrue((part_dir / f"spot_{size}" / "xenium_spot.h5ad").is_file())
                    self.assertTrue((part_dir / f"spot_{size}" / "cell_num_distribution.csv").is_file())
            self.assertEqual(list(result.output_dir.rglob("*simulated*")), [])
