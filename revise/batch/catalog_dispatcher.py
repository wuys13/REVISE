"""Catalog Dispatcher for Batch Multi-Sample Routing in REVISE."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Dict, List
import yaml

class CatalogBatchDispatcher:
    """Resolves ST Bank sample manifests and single-cell reference catalogs
    to dispatch standardized reconstruction pipelines.
    
    Routing Principles:
      - iST (Imaging-based, e.g. Xenium, CosMx) -> sc-SVC (cluster / within_cluster)
      - sST (Spot-based, e.g. Visium 55um) -> sc-SVC-sr (super-resolution virtual single-cell restoration)
      - HST (High-res bin, e.g. VisiumHD, Stereo-seq) -> sp-SVC (continuous 2um grid super-resolution)
    """

    DEFAULT_ROUTING = {
        "iST": "sc-SVC",
        "sST": "sc-SVC-sr",
        "HST": "sp-SVC",
    }

    def __init__(self, spec_path: str | Path) -> None:
        self.spec_path = Path(spec_path).resolve()
        self.base_dir = self.spec_path.parent.parent
        with open(self.spec_path, "r", encoding="utf-8") as f:
            self.spec = yaml.safe_load(f)

        self.manifest_path = self.base_dir / self.spec.get("manifest_path", "data/input/samples_manifest.csv")
        self.catalog_path = self.base_dir / self.spec.get("catalog_path", "data/input/reference_catalog.yaml")
        self.st_bank_root = self.base_dir / self.spec.get("st_bank_root", "data/input/st_bank")
        self.sc_ref_root = self.base_dir / self.spec.get("sc_ref_root", "data/input/sc_ref")
        self.output_root = self.base_dir / self.spec.get("output_root", "data/output")
        self.routing_policy = self.spec.get("routing_policy", self.DEFAULT_ROUTING)

    def load_catalog(self) -> Dict[str, Any]:
        if not self.catalog_path.exists():
            raise FileNotFoundError(f"Reference catalog not found: {self.catalog_path}")
        with open(self.catalog_path, "r", encoding="utf-8") as f:
            return yaml.safe_load(f)

    def load_manifest(self) -> List[Dict[str, str]]:
        if not self.manifest_path.exists():
            raise FileNotFoundError(f"Samples manifest not found: {self.manifest_path}")
        with open(self.manifest_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            return list(reader)

    def plan_dispatches(self) -> List[Dict[str, Any]]:
        catalog = self.load_catalog()
        samples = self.load_manifest()
        plans = []
        for s in samples:
            cid = s.get("cancer_type", "UNKNOWN")
            mod = s.get("tech_class", "iST")
            sid = s.get("sample_id", "sample")
            ref_entry = catalog.get(cid, {})
            engine = self.routing_policy.get(mod, "sc-SVC")
            target_out = self.output_root / cid / sid

            plans.append({
                "sample_id": sid,
                "cancer_type": cid,
                "tech_class": mod,
                "engine": engine,
                "reference_h5ad": ref_entry.get("reference_h5ad", ""),
                "target_out": str(target_out),
            })
        return plans
