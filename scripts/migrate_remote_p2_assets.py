"""Safely migrate remote P2CRC Xenium assets to Schema 2 Batch standard output and clean up old directories."""

import os
from pathlib import Path
import shutil
import json
import yaml

def migrate(repo_root: Path):
    print("Starting migration and cleanup...")
    
    # Paths
    old_p2_run = repo_root / "results/p2-xenium-impact-20260922/iter-001"
    old_p2_xen = old_p2_run / "P2CRC_Xenium"
    analysis_repo = repo_root.parent / "REVISE_Analysis_Agent"
    old_analysis_run = analysis_repo / "output/p2-xenium-impact-20260922/iter-001"
    old_rendered = analysis_repo / "output/rendered_notebooks"
    
    # Destination
    output_sample = repo_root / "output/CRC/P2_CRC_Xenium"
    output_sample.mkdir(parents=True, exist_ok=True)
    revise_meta = output_sample / ".revise"
    revise_meta.mkdir(parents=True, exist_ok=True)
    analysis_dir = output_sample / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    methods_dir = output_sample / "methods"
    methods_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Migrate primary whole-sample reconstruction (mean mode)
    if (old_p2_xen / "mean/SVC.h5ad").exists() and not (output_sample / "SVC.h5ad").exists():
        print("Migrating mean SVC.h5ad to primary SVC.h5ad...")
        shutil.move(str(old_p2_xen / "mean/SVC.h5ad"), str(output_sample / "SVC.h5ad"))
    
    if (old_p2_xen / "mean/raw.h5ad").exists() and not (output_sample / "raw.h5ad").exists():
        print("Migrating mean raw.h5ad to raw.h5ad...")
        shutil.move(str(old_p2_xen / "mean/raw.h5ad"), str(output_sample / "raw.h5ad"))

    # Update sample.yaml
    sample_yaml_dest = output_sample / "sample.yaml"
    if (old_p2_xen / "mean/sample.yaml").exists() and not sample_yaml_dest.exists():
        with open(old_p2_xen / "mean/sample.yaml", "r", encoding="utf-8") as f:
            sample_doc = yaml.safe_load(f)
        sample_doc["sample_id"] = "P2_CRC_Xenium"
        sample_doc["metadata"] = {
            "patient_id": "P2",
            "cancer_type": "CRC",
            "assay": "Xenium",
            "modality": "iST",
            "technology": "in_situ_sequencing",
        }
        with open(sample_yaml_dest, "w", encoding="utf-8") as f:
            yaml.dump(sample_doc, f, sort_keys=False)
        print("Created standard sample.yaml with updated sample_id and metadata")

    # Move receipt to reconstruction.json
    if (old_p2_run / "receipt.json").exists() and not (output_sample / "reconstruction.json").exists():
        shutil.move(str(old_p2_run / "receipt.json"), str(output_sample / "reconstruction.json"))
        print("Migrated receipt.json -> reconstruction.json")

    # Move log to .revise/reconstruction.log
    if (old_p2_run / "reconstruction.log").exists() and not (revise_meta / "reconstruction.log").exists():
        shutil.move(str(old_p2_run / "reconstruction.log"), str(revise_meta / "reconstruction.log"))
        print("Migrated reconstruction.log -> .revise/reconstruction.log")

    # 2. Preserve methods (random, within_cluster, mean remaining carriers)
    for m in ["random", "within_cluster", "mean"]:
        src_m = old_p2_xen / m
        dst_m = methods_dir / m
        if src_m.exists() and not dst_m.exists():
            shutil.move(str(src_m), str(dst_m))
            print(f"Moved method directory {m} -> methods/{m}")

    # 3. Migrate Analysis Agent outputs
    # 3.1 reconstruction_impact
    impact_src = old_analysis_run / "reconstruction_impact"
    impact_dst = analysis_dir / "reconstruction_impact"
    if impact_src.exists() and not impact_dst.exists():
        shutil.move(str(impact_src), str(impact_dst))
        print("Migrated reconstruction_impact -> analysis/reconstruction_impact")

    # 3.2 HTML reports from rendered_notebooks
    aspect_map = {
        "00_master_reconstruction_synthesis.html": ("master_synthesis", "00_master_reconstruction_synthesis.html"),
        "01_reconstruction_methods_benchmark.html": ("methods_benchmark", "01_methods_benchmark.html"),
        "02_spatial_autocorrelation_and_smoothness.html": ("spatial_continuity", "02_spatial_smoothness_moran.html"),
        "03_microenvironment_diversity_and_gain.html": ("microenvironment_diversity", "03_microenvironment_diversity.html"),
        "04_pathway_and_spatial_communication.html": ("pathway_communication", "04_pathway_and_spatial_cci.html"),
    }
    for old_name, (aspect, new_name) in aspect_map.items():
        src_html = old_rendered / old_name
        dst_aspect_dir = analysis_dir / aspect
        dst_aspect_dir.mkdir(parents=True, exist_ok=True)
        dst_html = dst_aspect_dir / new_name
        if src_html.exists() and not dst_html.exists():
            shutil.move(str(src_html), str(dst_html))
            print(f"Migrated report {old_name} -> analysis/{aspect}/{new_name}")

    # 4. Generate analysis.json summary
    analysis_index = {
        "schema_version": 1,
        "sample_id": "P2_CRC_Xenium",
        "aspects": {
            "reconstruction_impact": {
                "status": "succeeded" if (analysis_dir / "reconstruction_impact/result.json").exists() else "pending",
                "result": "reconstruction_impact/result.json" if (analysis_dir / "reconstruction_impact/result.json").exists() else None,
                "report": "reconstruction_impact/report.html" if (analysis_dir / "reconstruction_impact/report.html").exists() else None,
            },
            "master_synthesis": {
                "status": "succeeded" if (analysis_dir / "master_synthesis/00_master_reconstruction_synthesis.html").exists() else "pending",
                "report": "master_synthesis/00_master_reconstruction_synthesis.html" if (analysis_dir / "master_synthesis/00_master_reconstruction_synthesis.html").exists() else None,
            },
            "methods_benchmark": {
                "status": "succeeded" if (analysis_dir / "methods_benchmark/01_methods_benchmark.html").exists() else "pending",
                "report": "methods_benchmark/01_methods_benchmark.html" if (analysis_dir / "methods_benchmark/01_methods_benchmark.html").exists() else None,
            },
            "spatial_continuity": {
                "status": "succeeded" if (analysis_dir / "spatial_continuity/02_spatial_smoothness_moran.html").exists() else "pending",
                "report": "spatial_continuity/02_spatial_smoothness_moran.html" if (analysis_dir / "spatial_continuity/02_spatial_smoothness_moran.html").exists() else None,
            },
            "microenvironment_diversity": {
                "status": "succeeded" if (analysis_dir / "microenvironment_diversity/03_microenvironment_diversity.html").exists() else "pending",
                "report": "microenvironment_diversity/03_microenvironment_diversity.html" if (analysis_dir / "microenvironment_diversity/03_microenvironment_diversity.html").exists() else None,
            },
            "pathway_communication": {
                "status": "succeeded" if (analysis_dir / "pathway_communication/04_pathway_and_spatial_cci.html").exists() else "pending",
                "report": "pathway_communication/04_pathway_and_spatial_cci.html" if (analysis_dir / "pathway_communication/04_pathway_and_spatial_cci.html").exists() else None,
            },
        }
    }
    with open(analysis_dir / "analysis.json", "w", encoding="utf-8") as f:
        json.dump(analysis_index, f, indent=2)
    print("Created analysis/analysis.json")

    # 5. Safe cleanup of old directories
    old_results_dir = repo_root / "results/p2-xenium-impact-20260922"
    if old_results_dir.exists():
        shutil.rmtree(old_results_dir)
        print(f"Cleaned up {old_results_dir}")

    old_analysis_impact_dir = analysis_repo / "output/p2-xenium-impact-20260922"
    if old_analysis_impact_dir.exists():
        shutil.rmtree(old_analysis_impact_dir)
        print(f"Cleaned up {old_analysis_impact_dir}")

    if old_rendered.exists() and not any(old_rendered.iterdir()):
        shutil.rmtree(old_rendered)
        print(f"Cleaned up empty {old_rendered}")

    print("Migration and cleanup finished successfully!")

if __name__ == "__main__":
    import sys
    root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(".").resolve()
    migrate(root)
