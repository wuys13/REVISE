from __future__ import annotations

import argparse
import json
import os
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
POSTERIOR_MODES = ["off", "cost", "reference"]
DEFAULT_POSTERIOR_MODE = "cost"

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


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser("REVISE benchmark (unified interface wrapper)")
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
        "--posterior-mode",
        choices=POSTERIOR_MODES,
        default=DEFAULT_POSTERIOR_MODE,
        help="Unified posterior-conditioning mode for every confounding route (default: cost)",
    )
    parser.add_argument(
        "--posterior-key",
        default=None,
        help="AnnData obsm key for posterior probabilities",
    )
    parser.add_argument(
        "--posterior-beta",
        type=float,
        default=None,
        help="Posterior affinity sharpness",
    )
    parser.add_argument(
        "--posterior-min-affinity",
        type=float,
        default=None,
        help="Posterior affinity floor",
    )
    parser.add_argument(
        "--posterior-cost-strength",
        type=float,
        default=None,
        help="Cost penalty strength for mode=cost",
    )
    parser.add_argument(
        "--posterior-strict",
        action="store_true",
        help="Fail if an active posterior mode cannot be applied to the local OT problem",
    )
    parser.add_argument(
        "--sr-refinement-preset",
        choices=sorted(SR_REFINEMENT_PRESETS),
        default=None,
        help="SR-only graph refinement preset; use confidence_anchor for controlled SR posterior-OT ablations",
    )
    return parser.parse_args()


def _build_posterior_overrides(args: argparse.Namespace) -> Dict[str, object]:
    values: Dict[str, object] = {}
    if args.posterior_mode is not None:
        values["enabled"] = args.posterior_mode != "off"
        values["mode"] = args.posterior_mode
    if args.posterior_key is not None:
        values["posterior_key"] = args.posterior_key
    if args.posterior_beta is not None:
        values["beta"] = args.posterior_beta
    if args.posterior_min_affinity is not None:
        values["min_affinity"] = args.posterior_min_affinity
    if args.posterior_cost_strength is not None:
        values["cost_strength"] = args.posterior_cost_strength
    if args.posterior_strict:
        values["strict"] = True
    return {"posterior_conditioning": values} if values else {}


def _build_sr_refinement_overrides(args: argparse.Namespace) -> Dict[str, object]:
    if args.sr_refinement_preset is None:
        return {}
    if args.confounding not in SR_CONFOUNDINGS and args.sr_refinement_preset != "none":
        raise ValueError(
            "--sr-refinement-preset is only valid for batch_effect and spot_size benchmark routes"
        )
    return {"sc": dict(SR_REFINEMENT_PRESETS[args.sr_refinement_preset])}


def _build_algorithm_overrides(args: argparse.Namespace) -> Dict[str, object]:
    if (
        args.posterior_mode in {"cost", "reference"}
        and args.confounding in SR_CONFOUNDINGS
        and args.sr_refinement_preset == "none"
    ):
        raise ValueError(
            "posterior-mode=cost/reference on SR benchmark routes requires graph aggregation; "
            "use the profile default or --sr-refinement-preset confidence_anchor"
        )

    overrides = _build_posterior_overrides(args)
    overrides.update(_build_sr_refinement_overrides(args))
    return overrides


def _run_case(
    pipeline: REVISEPipeline,
    *,
    platform: str,
    profile: str,
    confounding: str,
    io_overrides: Dict[str, object],
    runtime_seed: int | None,
    algorithm_overrides: Dict[str, object] | None = None,
) -> Dict[str, object]:
    try:
        runtime_cfg = {
            "platform": platform,
            "confounding": confounding,
            "seed": runtime_seed,
        }
        svc = pipeline._run_with_algorithm_overrides(
            profile=profile,
            runtime_overrides=runtime_cfg,
            io_overrides=io_overrides,
            algorithm_overrides=algorithm_overrides,
            dry_run=False,
        )
        summary = svc.summary()
        return {
            "ok": True,
            "profile": profile,
            "seed": runtime_seed,
            "run_dir": svc.provenance.get("run_dir"),
            "summary": summary,
            "error": None,
        }
    except Exception as exc:  # pragma: no cover - wrapper behavior
        return {
            "ok": False,
            "profile": profile,
            "seed": runtime_seed,
            "run_dir": None,
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
                platform=args.platform,
                profile="benchmark_seg",
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
                platform=args.platform,
                profile="benchmark_bin2cell",
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
                    platform=args.platform,
                    profile="benchmark_sr_batch",
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
                platform=args.platform,
                profile="benchmark_sr_spot_size",
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
            platform=args.platform,
            profile="benchmark_impute_panel",
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
            platform=args.platform,
            profile="benchmark_impute_dropout",
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
        "posterior_mode": args.posterior_mode,
        "posterior_strict": bool(args.posterior_strict),
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
