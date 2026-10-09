"""Setup standard Batch (Schema 2) input structure and metadata for REVISE."""

from pathlib import Path
import yaml

def setup_batch(revise_root: Path, raw_data_dir: Path):
    # 1. Root batch.yaml
    root_batch = {
        "schema_version": 2,
        "input_root": "ST",
        "output_root": "output",
        "modality": "iST",
        "coordinates": {
            "unit": "pixel",
            "microns_per_coordinate": 0.2125,
        },
        "algorithm": {
            "ot_method": "tacco",
        },
        "preprocessing": {
            "spatial": {
                "min_transcript_counts": 60,
                "min_cell_counts": 100,
            },
            "reference": {
                "min_transcript_counts": None,
                "min_cell_counts": 100,
            },
        },
        "global_anchoring": {
            "broad_column": "Level1",
        },
        "local_refinement": {
            "subtype_column": "Level2",
            "alpha": 0.2,
            "resolutions": [0.6, 0.7, 0.8],
        },
        "output": {
            "ist_mapping": "within_cluster",
        },
        "execution": {
            "seed": 42,
        },
        "analysis": {},
    }
    
    with (revise_root / "batch.yaml").open("w", encoding="utf-8") as f:
        yaml.dump(root_batch, f, sort_keys=False)
    print(f"Created {revise_root / 'batch.yaml'}")

    # 2. SC directory
    sc_crc = revise_root / "SC" / "CRC"
    sc_crc.mkdir(parents=True, exist_ok=True)
    sc_ref_target = raw_data_dir / "adata_sc_all_reanno.h5ad"
    sc_ref_link = sc_crc / "adata_sc_all_reanno.h5ad"
    if sc_ref_target.exists() and not sc_ref_link.exists():
        sc_ref_link.symlink_to(sc_ref_target.resolve())
        print(f"Linked {sc_ref_link} -> {sc_ref_target}")

    # 3. ST/CRC/batch.yaml
    st_crc = revise_root / "ST" / "CRC"
    st_crc.mkdir(parents=True, exist_ok=True)
    crc_batch = {
        "inputs": {
            "reference": {
                "path": "../../SC/CRC/adata_sc_all_reanno.h5ad",
                "filter_column": "Patient",
                "filter_value": "P2CRC",
            }
        }
    }
    with (st_crc / "batch.yaml").open("w", encoding="utf-8") as f:
        yaml.dump(crc_batch, f, sort_keys=False)
    print(f"Created {st_crc / 'batch.yaml'}")

    # 4. Samples:
    # 4.1 P2_CRC_Xenium (iST)
    sample_xen = st_crc / "P2_CRC_Xenium"
    sample_xen.mkdir(parents=True, exist_ok=True)
    xen_data = raw_data_dir / "P2CRC_Xenium.h5ad"
    xen_link = sample_xen / "spatial.h5ad"
    if xen_data.exists() and not xen_link.exists():
        xen_link.symlink_to(xen_data.resolve())
        print(f"Linked {xen_link} -> {xen_data}")
    xen_config = {
        "metadata": {
            "patient_id": "P2",
            "cancer_type": "CRC",
            "assay": "Xenium",
            "modality": "iST",
            "technology": "in_situ_sequencing",
        },
        "modality": "iST",
        "coordinates": {
            "unit": "pixel",
            "microns_per_coordinate": 0.2125,
        },
        "output": {
            "ist_mapping": "within_cluster",
        },
    }
    with (sample_xen / "batch.yaml").open("w", encoding="utf-8") as f:
        yaml.dump(xen_config, f, sort_keys=False)
    print(f"Created {sample_xen / 'batch.yaml'}")

    # 4.2 P2_CRC_VisiumHD (hST)
    sample_hd = st_crc / "P2_CRC_VisiumHD"
    sample_hd.mkdir(parents=True, exist_ok=True)
    hd_data = raw_data_dir / "P2CRC_HD.h5ad"
    hd_link = sample_hd / "spatial.h5ad"
    if hd_data.exists() and not hd_link.exists():
        hd_link.symlink_to(hd_data.resolve())
        print(f"Linked {hd_link} -> {hd_data}")
    hd_config = {
        "metadata": {
            "patient_id": "P2",
            "cancer_type": "CRC",
            "assay": "VisiumHD",
            "modality": "hST",
            "technology": "high_definition_spatial_transcriptomics",
        },
        "modality": "hST",
    }
    with (sample_hd / "batch.yaml").open("w", encoding="utf-8") as f:
        yaml.dump(hd_config, f, sort_keys=False)
    print(f"Created {sample_hd / 'batch.yaml'}")

    # 4.3 P2_CRC_Visium (sST)
    sample_visium = st_crc / "P2_CRC_Visium"
    sample_visium.mkdir(parents=True, exist_ok=True)
    visium_data = raw_data_dir / "P2CRC_Visium.h5ad"
    visium_link = sample_visium / "spatial.h5ad"
    if visium_data.exists() and not visium_link.exists():
        visium_link.symlink_to(visium_data.resolve())
        print(f"Linked {visium_link} -> {visium_data}")
    visium_config = {
        "metadata": {
            "patient_id": "P2",
            "cancer_type": "CRC",
            "assay": "Visium",
            "modality": "sST",
            "technology": "spatial_transcriptomics",
        },
        "modality": "sST",
    }
    with (sample_visium / "batch.yaml").open("w", encoding="utf-8") as f:
        yaml.dump(visium_config, f, sort_keys=False)
    print(f"Created {sample_visium / 'batch.yaml'}")

if __name__ == "__main__":
    import sys
    revise = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    raw = Path(sys.argv[2]) if len(sys.argv) > 2 else revise / "raw_data/Real_application"
    setup_batch(revise, raw)
