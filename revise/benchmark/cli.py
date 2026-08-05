from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np

from revise.framework import REVISEPipeline

# Base settings aligned with Sim2Real-ST benchmark convention.
SEG_METHODS = ["seg_1", "seg_2", "seg_3", "seg_4"]
BIN2CELL_METHODS = ["bin2cell"]
BATCH_NUMS = [1, 2, 3, 4]
SPOT_SIZES = [20, 50, 100, 200]
SR_CONFOUNDINGS = {"batch_effect", "spot_size"}
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
        "--confounding",
        required=True,
        type=str,
        choices=[
            "segmentation",
            "bin2cell",
            "batch_effect",
            "spot_size",
            "gene_panel",
            "gene_dropout",
        ],
        help="Confounding factor: segmentation/bin2cell/batch_effect/spot_size/gene_panel/gene_dropout",
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
    parser.add_argument("--config", default="revise/revise.yaml", help="Path to unified config")
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


def _build_local_refinement_overrides(args: argparse.Namespace) -> Dict[str, object]:
    strength = getattr(args, "local_refinement_strength", None)
    if strength is None:
        return {}
    return {"local_refinement": {"strength": strength}}


def _build_sr_refinement_overrides(args: argparse.Namespace) -> Dict[str, object]:
    if args.sr_refinement_preset is None:
        return {}
    if args.confounding not in SR_CONFOUNDINGS and args.sr_refinement_preset != "none":
        raise ValueError(
            "--sr-refinement-preset is only valid for batch_effect and spot_size benchmark routes"
        )
    return {"sc": dict(SR_REFINEMENT_PRESETS[args.sr_refinement_preset])}


def _build_algorithm_overrides(args: argparse.Namespace) -> Dict[str, object]:
    overrides = _build_local_refinement_overrides(args)
    overrides.update(_build_sr_refinement_overrides(args))
    return overrides


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
    confounding: str,
    io_overrides: Dict[str, object],
    runtime_seed: int | None,
    algorithm_overrides: Dict[str, object] | None = None,
) -> Dict[str, object]:
    try:
        svc = pipeline._execute_run(
            svc_type=None,
            cf=confounding,
            runtime_overrides={"seed": runtime_seed},
            io_overrides=io_overrides,
            algorithm_overrides=algorithm_overrides,
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
    output_task = args.dataset_task or args.confounding
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


def _discover_spot_sizes(*, data_root: str, sample_name: str, st_file: str) -> List[int]:
    sample_root = Path(data_root) / sample_name
    sizes: List[int] = []
    for spot_dir in sorted(sample_root.glob("spot_*")):
        if not spot_dir.is_dir():
            continue
        suffix = spot_dir.name.removeprefix("spot_")
        if not suffix.isdigit():
            continue
        if (spot_dir / st_file).exists():
            sizes.append(int(suffix))
    return sorted(set(sizes)) or list(SPOT_SIZES)


def _batch_spec(batch_num: int) -> tuple[str, str]:
    if batch_num == 1:
        return "selected_xenium.h5ad", "selected_xenium.h5ad"
    if batch_num == 2:
        return "selected_xenium.h5ad", "real_sc_ref_part.h5ad"
    if batch_num == 3:
        return "selected_xenium.h5ad", "real_sc_ref_all.h5ad"
    if batch_num == 4:
        return "selected_xenium.h5ad", "real_sc_ref_others.h5ad"
    raise NotImplementedError(f"batch_num {batch_num} not implemented")


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
    pipeline = REVISEPipeline(config_path=args.config)
    results: List[Dict[str, object]] = []
    base_io = _base_io(args)
    seed_scope = _resolve_seed_scope(args)
    next_runtime_seed = _runtime_seed_supplier(args, seed_scope)
    try:
        algorithm_overrides = _build_algorithm_overrides(args)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    if args.confounding == "segmentation":
        st_file = args.st_file or "xenium_spot.h5ad"
        gt_svc_file = args.gt_svc_file or "selected_xenium.h5ad"
        sc_ref_file = args.sc_ref_file or "real_sc_ref.h5ad"
        for seg_method in SEG_METHODS:
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
                confounding="segmentation",
                io_overrides=io_cfg,
                runtime_seed=next_runtime_seed(),
                algorithm_overrides=algorithm_overrides,
            )
            _append_result(results, f"segmentation:{seg_method}", run_result)

    elif args.confounding == "bin2cell":
        st_file = args.st_file or "xenium_spot.h5ad"
        gt_svc_file = args.gt_svc_file or "selected_xenium.h5ad"
        sc_ref_file = args.sc_ref_file or "real_sc_ref.h5ad"
        for seg_method in BIN2CELL_METHODS:
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
                confounding="bin2cell",
                io_overrides=io_cfg,
                runtime_seed=next_runtime_seed(),
                algorithm_overrides=algorithm_overrides,
            )
            _append_result(results, f"bin2cell:{seg_method}", run_result)

    elif args.confounding == "batch_effect":
        st_file = args.st_file or "xenium_spot.h5ad"
        spot_sizes = _discover_spot_sizes(
            data_root=str(base_io["data_root"]),
            sample_name=args.sample_name,
            st_file=st_file,
        )
        for spot_size in spot_sizes:
            for batch_num in BATCH_NUMS:
                gt_svc_file, sc_ref_file = _batch_spec(batch_num)
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
                    confounding="batch_effect",
                    io_overrides=io_cfg,
                    runtime_seed=next_runtime_seed(),
                    algorithm_overrides=algorithm_overrides,
                )
                _append_result(results, f"batch_effect:{spot_size}_{batch_num}", run_result)

    elif args.confounding == "spot_size":
        batch_num = 3
        gt_svc_file = "selected_xenium.h5ad"
        sc_ref_file = "real_sc_ref_all.h5ad"
        for spot_size in SPOT_SIZES:
            io_cfg = dict(base_io)
            io_cfg.update(
                {
                    "st_file": "xenium_spot.h5ad",
                    "gt_svc_file": gt_svc_file,
                    "sc_ref_file": sc_ref_file,
                    "spot_size": spot_size,
                }
            )
            run_result = _run_case(
                pipeline,
                confounding="spot_size",
                io_overrides=io_cfg,
                runtime_seed=next_runtime_seed(),
                algorithm_overrides=algorithm_overrides,
            )
            _append_result(results, f"spot_size:{spot_size}_{batch_num}", run_result)

    elif args.confounding == "gene_panel":
        io_cfg = dict(base_io)
        io_cfg.update(
            {
                "st_file": "selected_xenium.h5ad",
                "gt_svc_file": "selected_xenium.h5ad",
                "sc_ref_file": "real_sc_ref.h5ad",
            }
        )
        run_result = _run_case(
            pipeline,
            confounding="gene_panel",
            io_overrides=io_cfg,
            runtime_seed=next_runtime_seed(),
            algorithm_overrides=algorithm_overrides,
        )
        _append_result(results, "gene_panel", run_result)

    elif args.confounding == "gene_dropout":
        io_cfg = dict(base_io)
        io_cfg.update(
            {
                "st_file": "selected_xenium.h5ad",
                "gt_svc_file": "selected_xenium.h5ad",
                "sc_ref_file": "real_sc_ref.h5ad",
            }
        )
        run_result = _run_case(
            pipeline,
            confounding="gene_dropout",
            io_overrides=io_cfg,
            runtime_seed=next_runtime_seed(),
            algorithm_overrides=algorithm_overrides,
        )
        _append_result(results, "gene_dropout", run_result)

    else:
        raise NotImplementedError(f"Unsupported confounding: {args.confounding}")

    report = {
        "platform": args.platform,
        "confounding": args.confounding,
        "sample_name": args.sample_name,
        "dataset_task": args.dataset_task,
        "output_task": args.dataset_task or args.confounding,
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
        print("Done:", args.confounding)


if __name__ == "__main__":
    main()
