from pathlib import Path
import hashlib

from revise.benchmark.cli import (
    _read_benchmark_request,
    _read_benchmark_request_with_metadata,
)
from revise.config.authority import ENGINE_DEFAULTS
from revise.utils.provenance import hash_jsonable


ROUTES = (
    "segmentation",
    "bin2cell",
    "batch_effect",
    "spot_size",
    "gene_panel",
    "gene_dropout",
)


def test_benchmark_template_inventory_bytes_cardinality_and_imputation_contract():
    repo_dir = Path("configs/benchmark")
    package_dir = Path("revise/benchmark/templates")
    assert {path.name for path in repo_dir.glob("*.yaml")} == {
        f"{route}.yaml" for route in ROUTES
    }
    assert {path.name for path in package_dir.glob("*.yaml")} == {
        f"{route}.yaml" for route in ROUTES
    }

    requests = {}
    for route in ROUTES:
        name = f"{route}.yaml"
        assert (repo_dir / name).read_bytes() == (package_dir / name).read_bytes()
        requests[route] = _read_benchmark_request(str(repo_dir / name))
        assert requests[route]["route"] == route

    assert len(requests["segmentation"]["cases"]["segmentation_methods"]) == 4
    assert len(requests["bin2cell"]["cases"]["segmentation_methods"]) == 1
    expected_spot_sizes = [50, 100, 150, 200]
    assert requests["batch_effect"]["cases"]["spot_sizes"] == expected_spot_sizes
    assert len(requests["batch_effect"]["cases"]["spot_sizes"]) * len(
        requests["batch_effect"]["cases"]["batches"]
    ) == 16
    assert requests["spot_size"]["cases"]["spot_sizes"] == expected_spot_sizes
    assert requests["gene_panel"]["cases"] == {}
    assert requests["gene_dropout"]["cases"] == {}
    for request in requests.values():
        assert "data_root" not in request["io"]
        assert "sample_name" not in request["io"]

    expected_impute = {
        "merge_subcluster_method": "mean",
        "subcluster_resolution": 3,
        "in_panel_subcluster_resolution": None,
        "prune": True,
        "n_neighbors": 1,
        "method": "mean",
        "graph_preprocess": True,
        "graph_n_pcs": 50,
    }
    assert requests["gene_panel"]["algorithm"] == {}
    assert requests["gene_dropout"]["algorithm"] == {}
    assert ENGINE_DEFAULTS["impute"] == expected_impute

    assert not any(key.startswith("sr_noise_") for key in ENGINE_DEFAULTS["sc"])
    adapter_source = Path("revise/backend/adapters.py").read_text()
    for removed in (
        "_inject_sr_spatial_leakage_noise",
        "sr_spatial_noise",
        "st_input_noisy",
    ):
        assert removed not in adapter_source


def test_benchmark_request_records_source_and_effective_hashes():
    path = Path("configs/benchmark/segmentation.yaml")
    request, metadata = _read_benchmark_request_with_metadata(str(path))

    assert metadata["source_path"] == str(path.resolve())
    assert metadata["source_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert metadata["effective_request"] == request
    assert metadata["effective_request_hash"] == hash_jsonable(request)
