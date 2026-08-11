from __future__ import annotations

import argparse
import hashlib
from importlib import resources
import json
import os
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import yaml

from revise.framework import REVISEPipeline
from revise.config.authority import ENGINE_DEFAULTS
from revise.utils.provenance import hash_jsonable

BENCHMARK_ROUTES = frozenset(
    {"segmentation", "bin2cell", "batch_effect", "spot_size", "gene_panel", "gene_dropout"}
)
SR_ROUTES = {"batch_effect", "spot_size"}
REMOVED_ASSIGNMENT_GUIDANCE_FLAGS = frozenset(
    {
        "--local-refinement-guidance",
        "--local-refinement-compatibility-mode",
        "--posterior-mode",
        "--posterior-strict",
        "--posterior-key",
        "--posterior-beta",
        "--posterior-min-affinity",
        "--posterior-cost-strength",
    }
)
_CLI_MIGRATION_ERROR = (
    "Assignment guidance options were removed; "
    "use --local-refinement-strength"
)
SR_REFINEMENT_PRESETS = {
    "none": {
        "sr_graph_agg_enabled": False,
        "sr_graph_agg_low_conf_only": False,
        "sr_graph_agg_anchor_only": False,
        "sr_graph_agg_conf_weighted_alpha": False,
    },
    "confidence_anchor": {
        "sr_graph_agg_enabled": True,
        "sr_graph_agg_low_conf_only": True,
        "sr_graph_agg_low_conf_quantile": 0.1,
        "sr_graph_agg_anchor_only": True,
        "sr_graph_agg_anchor_high_conf_quantile": 0.9,
        "sr_graph_agg_confidence_mode": "auto",
        "sr_graph_agg_conf_weighted_alpha": True,
        "sr_graph_agg_conf_alpha_min": 0.0,
        "sr_graph_agg_conf_alpha_max": 0.25,
        "sr_graph_agg_conf_alpha_power": 1.0,
    },
}


class _BenchmarkArgumentParser(argparse.ArgumentParser):
    def parse_known_args(self, args=None, namespace=None):
        raw_args = sys.argv[1:] if args is None else list(args)
        if any(
            str(token).split("=", 1)[0] in REMOVED_ASSIGNMENT_GUIDANCE_FLAGS
            for token in raw_args
        ):
            self.error(_CLI_MIGRATION_ERROR)
        return super().parse_known_args(raw_args, namespace)


def get_parser() -> argparse.ArgumentParser:
    parser = _BenchmarkArgumentParser("REVISE benchmark (unified interface wrapper)")
    parser.add_argument(
        "--platform",
        default="sim2real",
        choices=["sim2real"],
        help="Benchmark platform route",
    )
    parser.add_argument(
        "--config",
        required=True,
        type=str,
        help="Benchmark route YAML path or packaged template name",
    )
    parser.add_argument("--data-root", required=True, type=str, help="Benchmark data root path")
    parser.add_argument(
        "--dataset-task",
        default=None,
        type=str,
        help="Optional dataset task subdirectory under data-root",
    )
    parser.add_argument("--sample-name", required=True, type=str, help="Sample name (e.g. P2CRC/cut_part1)")
    parser.add_argument("--st-file", type=str, help="ST file name for routes that need an override")
    parser.add_argument(
        "--gt-svc-file",
        type=str,
        help="Ground-truth SVC file name for routes that need an override",
    )
    parser.add_argument(
        "--sc-ref-file",
        type=str,
        help="sc reference file name for routes that need an override",
    )
    parser.add_argument("--output-root", default="results_unified", type=str, help="Output root")
    parser.add_argument("--sample-size", type=int, default=None, help="Optional subsample size")
    parser.add_argument("--seed", type=int, default=42, help="Base random seed for reproducibility")
    parser.add_argument(
        "--seed-scope",
        type=str,
        choices=["process", "run"],
        default=None,
        help=(
            "Derive one reproducible seed per case (`process`) or reuse the "
            "base seed for every case (`run`)"
        ),
    )
    parser.add_argument(
        "--local-refinement-strength",
        type=float,
        default=None,
        help="Override the route default local OT conditioning strength",
    )
    parser.add_argument(
        "--sr-refinement-preset",
        choices=sorted(SR_REFINEMENT_PRESETS),
        default=None,
        help="SR-only graph refinement preset",
    )
    return parser


def get_args() -> argparse.Namespace:
    return get_parser().parse_args()


def _read_benchmark_request_with_metadata(
    config: str,
) -> tuple[dict[str, object], dict[str, object]]:
    path = Path(config)
    if path.is_file():
        raw = path.read_bytes()
        source_path = str(path.resolve())
    elif path.name == config and config in {f"{route}.yaml" for route in BENCHMARK_ROUTES}:
        raw = resources.files("revise.benchmark").joinpath("templates", config).read_bytes()
        source_path = f"package:revise.benchmark/templates/{config}"
    else:
        raise FileNotFoundError(f"Benchmark config not found: {config}")
    request = yaml.safe_load(raw)
    if not isinstance(request, dict) or request.get("schema_version") != 1:
        raise ValueError("Benchmark config must be a schema_version 1 mapping")
    if set(request) != {"schema_version", "route", "cases", "io", "algorithm"}:
        raise ValueError("Benchmark config must contain only schema_version/route/cases/io/algorithm")
    route = request["route"]
    if route not in BENCHMARK_ROUTES:
        raise ValueError(f"Unsupported benchmark route: {route}")
    if not all(isinstance(request[key], dict) for key in ("cases", "io", "algorithm")):
        raise ValueError("Benchmark cases/io/algorithm must be mappings")
    metadata = {
        "source_path": source_path,
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "effective_request": request,
        "effective_request_hash": hash_jsonable(request),
    }
    return request, metadata


def _read_benchmark_request(config: str) -> dict[str, object]:
    request, _ = _read_benchmark_request_with_metadata(config)
    return request


def _build_local_refinement_overrides(args: argparse.Namespace) -> Dict[str, object]:
    strength = getattr(args, "local_refinement_strength", None)
    if strength is None:
        return {}
    return {"local_refinement": {"strength": strength}}


def _build_sr_refinement_overrides(args: argparse.Namespace) -> Dict[str, object]:
    if args.sr_refinement_preset is None:
        return {}
    if args.route not in SR_ROUTES and args.sr_refinement_preset != "none":
        raise ValueError(
            "--sr-refinement-preset is only valid for batch_effect and spot_size benchmark routes"
        )
    return {"sc": dict(SR_REFINEMENT_PRESETS[args.sr_refinement_preset])}


def _build_algorithm_overrides(args: argparse.Namespace) -> Dict[str, object]:
    overrides = _build_local_refinement_overrides(args)
    overrides.update(_build_sr_refinement_overrides(args))
    return overrides


def _cli_overrides(args: argparse.Namespace) -> Dict[str, object]:
    keys = (
        "st_file",
        "gt_svc_file",
        "sc_ref_file",
        "sample_size",
        "local_refinement_strength",
        "sr_refinement_preset",
        "seed_scope",
    )
    return {
        key: value
        for key in keys
        if (value := getattr(args, key, None)) is not None
    }


def _aggregate_local_refinement(
    results: List[Dict[str, object]],
) -> List[Dict[str, object]]:
    return [
        dict(item["local_refinement"])
        for item in results
        if isinstance(item.get("local_refinement"), dict)
    ]


def _read_manifest_snapshot(
    path: Path,
) -> tuple[dict[str, int | str], Dict[str, object]] | None:
    try:
        before = path.stat()
        content = path.read_bytes()
        after = path.stat()
    except OSError:
        return None
    before_identity = (before.st_mtime_ns, before.st_size)
    after_identity = (after.st_mtime_ns, after.st_size)
    if before_identity != after_identity or len(content) != after.st_size:
        return None
    try:
        manifest = json.loads(content)
    except json.JSONDecodeError:
        return None
    if not isinstance(manifest, dict):
        return None
    return (
        {
            "mtime_ns": after.st_mtime_ns,
            "size": after.st_size,
            "sha256": hashlib.sha256(content).hexdigest(),
        },
        manifest,
    )


def _run_case(
    pipeline: REVISEPipeline,
    *,
    route: str,
    io_overrides: Dict[str, object],
    runtime_seed: int | None,
    algorithm_overrides: Dict[str, object] | None = None,
    benchmark_config_metadata: Dict[str, object] | None = None,
) -> Dict[str, object]:
    if runtime_seed is None:
        runtime_seed = int(ENGINE_DEFAULTS["runtime"]["seed"])
    if benchmark_config_metadata is not None:
        benchmark_config_metadata = dict(benchmark_config_metadata)
        effective_request = {
            "route": route,
            "io": dict(io_overrides),
            "algorithm": dict(algorithm_overrides or {}),
            "runtime": {"seed": runtime_seed},
        }
        benchmark_config_metadata["effective_request"] = effective_request
        benchmark_config_metadata["effective_request_hash"] = hash_jsonable(
            effective_request
        )
    try:
        svc = pipeline.run(
            svc_type=None,
            cf=route,
            runtime_overrides={"seed": runtime_seed},
            io_overrides=io_overrides,
            algorithm_overrides=algorithm_overrides,
            benchmark_config_metadata=benchmark_config_metadata,
            dry_run=False,
        )
        summary = svc.summary()
        run_dir = svc.provenance.get("run_dir")
        return {
            "ok": True,
            "profile": svc.provenance.get("profile"),
            "seed": runtime_seed,
            "run_dir": run_dir,
            "manifest_path": (
                str(Path(run_dir) / "provenance.json")
                if run_dir is not None
                else None
            ),
            "local_refinement": svc.provenance.get("local_refinement"),
            "summary": summary,
            "error": None,
        }
    except Exception as exc:  # pragma: no cover - wrapper behavior
        local_refinement = None
        failure_context = getattr(exc, "_revise_failure_context", None)
        run_dir = None
        manifest_path = None
        current_manifest = False
        if isinstance(failure_context, dict):
            run_dir = failure_context.get("run_dir")
            manifest_path = failure_context.get("manifest_path")
        if manifest_path is not None:
            snapshot = _read_manifest_snapshot(Path(manifest_path))
            if snapshot is not None:
                current_identity, manifest = snapshot
                current_manifest = (
                    failure_context.get("manifest_identity") == current_identity
                )
        if current_manifest:
            local_refinement = manifest.get("local_refinement")
        return {
            "ok": False,
            "profile": manifest.get("profile") if current_manifest else None,
            "seed": runtime_seed,
            "run_dir": str(run_dir) if current_manifest else None,
            "manifest_path": str(manifest_path) if current_manifest else None,
            "local_refinement": local_refinement,
            "summary": None,
            "error": f"{type(exc).__name__}: {exc}",
        }


def _base_io(args: argparse.Namespace) -> Dict[str, object]:
    data_root = os.path.join(args.data_root, args.dataset_task) if args.dataset_task else args.data_root
    output_task = args.dataset_task or args.route
    output_root = os.path.join(args.output_root, output_task)
    io_cfg: Dict[str, object] = {
        "data_root": data_root,
        "output_root": output_root,
        "sample_name": args.sample_name,
    }
    if args.sample_size is not None:
        io_cfg["sample_size"] = args.sample_size
    return io_cfg


def _append_result(results: List[Dict[str, object]], tag: str, run_result: Dict[str, object]) -> None:
    item = {"tag": tag}
    item.update(run_result)
    results.append(item)
    status = "OK" if run_result["ok"] else "FAIL"
    print(f"[{status}] {tag} -> {run_result.get('run_dir') or run_result.get('error')}")


def _resolve_seed_scope(args: argparse.Namespace) -> str:
    if args.seed_scope is not None:
        return args.seed_scope
    # Preserve distinct seeds across benchmark cases unless the caller asks to
    # reuse the same seed for every run.
    return "process"


def _runtime_seed_supplier(args: argparse.Namespace, seed_scope: str):
    if seed_scope == "run":
        return lambda: args.seed
    seed_stream = np.random.RandomState(args.seed)
    return lambda: int(seed_stream.randint(0, np.iinfo(np.int32).max))


def main(args: argparse.Namespace | None = None) -> None:
    cli_invocation = args is None
    if args is None:
        args = get_args()
        print(args)
    request, benchmark_config_metadata = _read_benchmark_request_with_metadata(
        args.config
    )
    benchmark_config_metadata["cli_overrides"] = _cli_overrides(args)
    args.route = request["route"]
    cases = request["cases"]
    request_io = request["io"]
    pipeline = REVISEPipeline()
    results: List[Dict[str, object]] = []
    base_io = _base_io(args)
    seed_scope = _resolve_seed_scope(args)
    next_runtime_seed = _runtime_seed_supplier(args, seed_scope)
    try:
        algorithm_overrides = dict(request["algorithm"])
        algorithm_overrides.update(_build_algorithm_overrides(args))
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    if args.route == "segmentation":
        st_file = args.st_file or request_io["st_file"]
        gt_svc_file = args.gt_svc_file or request_io["gt_svc_file"]
        sc_ref_file = args.sc_ref_file or request_io["sc_ref_file"]
        for seg_method in cases["segmentation_methods"]:
            io_cfg = dict(base_io)
            io_cfg.update(
                {
                    "st_file": st_file,
                    "gt_svc_file": gt_svc_file,
                    "sc_ref_file": sc_ref_file,
                    "seg_method": seg_method,
                }
            )
            run_result = _run_case(
                pipeline,
                route="segmentation",
                io_overrides=io_cfg,
                runtime_seed=next_runtime_seed(),
                algorithm_overrides=algorithm_overrides,
                benchmark_config_metadata=benchmark_config_metadata,
            )
            _append_result(results, f"segmentation:{seg_method}", run_result)

    elif args.route == "bin2cell":
        st_file = args.st_file or request_io["st_file"]
        gt_svc_file = args.gt_svc_file or request_io["gt_svc_file"]
        sc_ref_file = args.sc_ref_file or request_io["sc_ref_file"]
        for seg_method in cases["segmentation_methods"]:
            io_cfg = dict(base_io)
            io_cfg.update(
                {
                    "st_file": st_file,
                    "gt_svc_file": gt_svc_file,
                    "sc_ref_file": sc_ref_file,
                    "seg_method": seg_method,
                }
            )
            run_result = _run_case(
                pipeline,
                route="bin2cell",
                io_overrides=io_cfg,
                runtime_seed=next_runtime_seed(),
                algorithm_overrides=algorithm_overrides,
                benchmark_config_metadata=benchmark_config_metadata,
            )
            _append_result(results, f"bin2cell:{seg_method}", run_result)

    elif args.route == "batch_effect":
        st_file = args.st_file or request_io["st_file"]
        for spot_size in cases["spot_sizes"]:
            for batch in cases["batches"]:
                batch_num = batch["number"]
                gt_svc_file = batch["gt_svc_file"]
                sc_ref_file = batch["sc_ref_file"]
                io_cfg = dict(base_io)
                io_cfg.update(
                    {
                        "st_file": st_file,
                        "gt_svc_file": gt_svc_file,
                        "sc_ref_file": sc_ref_file,
                        "spot_size": spot_size,
                    }
                )
                run_result = _run_case(
                    pipeline,
                    route="batch_effect",
                    io_overrides=io_cfg,
                    runtime_seed=next_runtime_seed(),
                    algorithm_overrides=algorithm_overrides,
                    benchmark_config_metadata=benchmark_config_metadata,
                )
                _append_result(results, f"batch_effect:{spot_size}_{batch_num}", run_result)

    elif args.route == "spot_size":
        batch_num = cases["batch_number"]
        gt_svc_file = request_io["gt_svc_file"]
        sc_ref_file = request_io["sc_ref_file"]
        for spot_size in cases["spot_sizes"]:
            io_cfg = dict(base_io)
            io_cfg.update(
                {
                    "st_file": request_io["st_file"],
                    "gt_svc_file": gt_svc_file,
                    "sc_ref_file": sc_ref_file,
                    "spot_size": spot_size,
                }
            )
            run_result = _run_case(
                pipeline,
                route="spot_size",
                io_overrides=io_cfg,
                runtime_seed=next_runtime_seed(),
                algorithm_overrides=algorithm_overrides,
                benchmark_config_metadata=benchmark_config_metadata,
            )
            _append_result(results, f"spot_size:{spot_size}_{batch_num}", run_result)

    elif args.route == "gene_panel":
        io_cfg = dict(base_io)
        io_cfg.update(
            {
                "st_file": request_io["st_file"],
                "gt_svc_file": request_io["gt_svc_file"],
                "sc_ref_file": request_io["sc_ref_file"],
            }
        )
        run_result = _run_case(
            pipeline,
            route="gene_panel",
            io_overrides=io_cfg,
            runtime_seed=next_runtime_seed(),
            algorithm_overrides=algorithm_overrides,
            benchmark_config_metadata=benchmark_config_metadata,
        )
        _append_result(results, "gene_panel", run_result)

    elif args.route == "gene_dropout":
        io_cfg = dict(base_io)
        io_cfg.update(
            {
                "st_file": request_io["st_file"],
                "gt_svc_file": request_io["gt_svc_file"],
                "sc_ref_file": request_io["sc_ref_file"],
            }
        )
        run_result = _run_case(
            pipeline,
            route="gene_dropout",
            io_overrides=io_cfg,
            runtime_seed=next_runtime_seed(),
            algorithm_overrides=algorithm_overrides,
            benchmark_config_metadata=benchmark_config_metadata,
        )
        _append_result(results, "gene_dropout", run_result)

    else:
        raise NotImplementedError(f"Unsupported benchmark route: {args.route}")

    report = {
        "platform": args.platform,
        "route": args.route,
        "sample_name": args.sample_name,
        "dataset_task": args.dataset_task,
        "output_task": args.dataset_task or args.route,
        "seed": args.seed,
        "seed_scope": seed_scope,
        "local_refinement": _aggregate_local_refinement(results),
        "sr_refinement_preset": args.sr_refinement_preset,
        "ok": all(item["ok"] for item in results),
        "total_runs": len(results),
        "passed_runs": sum(1 for item in results if item["ok"]),
        "results": results,
    }
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if not report["ok"]:
        raise SystemExit(1)
    if cli_invocation:
        print("Done:", args.route)


if __name__ == "__main__":
    main()
