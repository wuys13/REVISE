from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
from typing import Dict, List

import numpy as np

from revise.backend.ops.assignment_guidance import AssignmentGuidanceCollector
from revise.config import merge_unified_config
from revise.framework import REVISEPipeline
from revise.utils import build_run_dir

# Base settings aligned with Sim2Real-ST benchmark convention.
SEG_METHODS = ["seg_1", "seg_2", "seg_3", "seg_4"]
BIN2CELL_METHODS = ["bin2cell"]
BATCH_NUMS = [1, 2, 3, 4]
SPOT_SIZES = [20, 50, 100, 200]
SR_CONFOUNDINGS = {"batch_effect", "spot_size"}
POSTERIOR_MODES = ["off", "cost", "reference"]
BENCHMARK_PROFILES = {
    "segmentation": "benchmark_seg",
    "bin2cell": "benchmark_bin2cell",
    "batch_effect": "benchmark_sr_batch",
    "spot_size": "benchmark_sr_spot_size",
    "gene_panel": "benchmark_impute_panel",
    "gene_dropout": "benchmark_impute_dropout",
}

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


def get_parser() -> argparse.ArgumentParser:
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
        "--local-refinement-guidance",
        choices=["off", "prefer", "require"],
        default=None,
        help="Optional assignment-guidance policy; omitted uses the route default",
    )
    parser.add_argument(
        "--local-refinement-compatibility-mode",
        choices=["cost", "reference"],
        default=None,
        help="Compatibility injection mode when local-refinement guidance is enabled",
    )
    parser.add_argument(
        "--posterior-mode",
        choices=POSTERIOR_MODES,
        default=None,
        help="Deprecated alias for local-refinement guidance and compatibility mode",
    )
    parser.add_argument(
        "--posterior-key",
        default=None,
        help=(
            "Deprecated: non-empty values are rejected; each route now uses "
            "route-provided Assignment State"
        ),
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
    return parser


def get_args() -> argparse.Namespace:
    return get_parser().parse_args()


def _build_posterior_overrides(args: argparse.Namespace) -> Dict[str, object]:
    new_values: Dict[str, object] = {}
    guidance = getattr(args, "local_refinement_guidance", None)
    compatibility_mode = getattr(
        args,
        "local_refinement_compatibility_mode",
        None,
    )
    if guidance is not None:
        new_values["guidance"] = guidance
    if compatibility_mode is not None:
        new_values["compatibility"] = {"mode": compatibility_mode}

    legacy_values: Dict[str, object] = {}
    posterior_mode = getattr(args, "posterior_mode", None)
    posterior_strict = bool(getattr(args, "posterior_strict", False))
    posterior_key = getattr(args, "posterior_key", None)
    if posterior_key not in (None, ""):
        raise ValueError(
            "--posterior-key is no longer supported; Assignment State is "
            "provided explicitly by each local-refinement route"
        )
    if posterior_mode == "off" and posterior_strict:
        raise ValueError("--posterior-mode off conflicts with --posterior-strict")
    if posterior_mode is not None:
        legacy_values["mode"] = posterior_mode
    if posterior_key == "":
        legacy_values["posterior_key"] = posterior_key
    for arg_name, key in (
        ("posterior_beta", "beta"),
        ("posterior_min_affinity", "min_affinity"),
        ("posterior_cost_strength", "cost_strength"),
    ):
        value = getattr(args, arg_name, None)
        if value is not None:
            legacy_values[key] = value
    if posterior_strict:
        legacy_values["strict"] = True

    legacy_guidance = None
    if posterior_mode == "off":
        legacy_guidance = "off"
    elif posterior_strict:
        legacy_guidance = "require"
    elif posterior_mode in {"cost", "reference"}:
        legacy_guidance = "prefer"
    if guidance is not None and legacy_guidance is not None and guidance != legacy_guidance:
        raise ValueError(
            "--local-refinement-guidance conflicts with deprecated posterior flags"
        )
    legacy_mode = posterior_mode if posterior_mode in {"cost", "reference"} else None
    if (
        compatibility_mode is not None
        and legacy_mode is not None
        and compatibility_mode != legacy_mode
    ):
        raise ValueError(
            "--local-refinement-compatibility-mode conflicts with --posterior-mode"
        )

    overrides: Dict[str, object] = {}
    if new_values:
        overrides["local_refinement"] = new_values
    if legacy_values:
        overrides["posterior_conditioning"] = legacy_values
    return overrides


def _build_sr_refinement_overrides(args: argparse.Namespace) -> Dict[str, object]:
    if args.sr_refinement_preset is None:
        return {}
    if args.confounding not in SR_CONFOUNDINGS and args.sr_refinement_preset != "none":
        raise ValueError(
            "--sr-refinement-preset is only valid for batch_effect and spot_size benchmark routes"
        )
    return {"sc": dict(SR_REFINEMENT_PRESETS[args.sr_refinement_preset])}


def _build_algorithm_overrides(args: argparse.Namespace) -> Dict[str, object]:
    overrides = _build_posterior_overrides(args)
    overrides.update(_build_sr_refinement_overrides(args))
    return overrides


def _resolved_report_aliases(
    *,
    raw_config: Dict[str, object],
    profile: str,
    platform: str,
    confounding: str,
    algorithm_overrides: Dict[str, object],
) -> Dict[str, object]:
    resolved = merge_unified_config(
        raw_config=raw_config,
        profile=profile,
        runtime_overrides={
            "platform": platform,
            "confounding": confounding,
        },
        io_overrides={},
        algorithm_overrides=algorithm_overrides,
    )
    refinement = resolved["local_refinement"]
    guidance = refinement["guidance"]
    compatibility_mode = refinement["compatibility"]["mode"]
    guidance_evidence = AssignmentGuidanceCollector(
        request_evidence=getattr(resolved, "request_evidence", {}),
        resolved_request=refinement,
    ).manifest()
    return {
        "assignment_guidance": guidance_evidence,
        "posterior_mode": (
            "off" if guidance == "off" else compatibility_mode
        ),
        "posterior_strict": guidance == "require",
    }


def _aggregate_assignment_guidance(
    results: List[Dict[str, object]],
    *,
    request_manifest: Dict[str, object],
) -> Dict[str, object]:
    leaves = [
        (case_ordinal, item, item["assignment_guidance"])
        for case_ordinal, item in enumerate(results, start=1)
        if isinstance(item.get("assignment_guidance"), dict)
    ]
    if not leaves:
        return copy.deepcopy(request_manifest)

    configured = copy.deepcopy(leaves[0][2]["configured"])
    resolved = copy.deepcopy(leaves[0][2]["resolved"])
    if (
        configured != request_manifest["configured"]
        or resolved != request_manifest["resolved"]
    ):
        raise ValueError(
            "inconsistent assignment-guidance resolved request and leaf evidence"
        )

    events: List[Dict[str, object]] = []
    for case_ordinal, item, manifest in leaves:
        if manifest.get("schema_version") != 2:
            raise ValueError("unsupported assignment-guidance leaf schema")
        if manifest.get("configured") != configured:
            raise ValueError(
                "inconsistent assignment-guidance configured leaf evidence"
            )
        if manifest.get("resolved") != resolved:
            raise ValueError(
                "inconsistent assignment-guidance resolved leaf evidence"
            )
        case_tag = str(item["tag"])
        for event_index, leaf_event in enumerate(
            manifest.get("events", []),
            start=1,
        ):
            event = copy.deepcopy(leaf_event)
            leaf_ordinal = event.get("ordinal", event_index)
            leaf_problem_key = str(event["problem_key"])
            event.update(
                {
                    "ordinal": len(events) + 1,
                    "problem_key": f"{case_tag}::{leaf_problem_key}",
                    "case_ordinal": case_ordinal,
                    "case_tag": case_tag,
                    "leaf_ordinal": leaf_ordinal,
                }
            )
            events.append(event)

    collector = AssignmentGuidanceCollector()
    collector.events = copy.deepcopy(events)
    return {
        "schema_version": 2,
        "configured": configured,
        "resolved": resolved,
        "events": events,
        "summary": collector.summary(),
    }


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
    platform: str,
    profile: str,
    confounding: str,
    io_overrides: Dict[str, object],
    runtime_seed: int | None,
    algorithm_overrides: Dict[str, object] | None = None,
) -> Dict[str, object]:
    expected_run_dir = None
    expected_manifest_path = None
    try:
        resolved = merge_unified_config(
            raw_config=pipeline.raw_config,
            profile=profile,
            runtime_overrides={
                "platform": platform,
                "confounding": confounding,
                "seed": runtime_seed,
            },
            io_overrides=io_overrides,
            algorithm_overrides=algorithm_overrides or {},
        )
        route_key = (
            f"{resolved['runtime']['platform']}:"
            f"{resolved['runtime']['confounding']}"
        )
        expected_run_dir = build_run_dir(
            output_root=resolved["io"]["output_root"],
            sample_name=resolved["io"]["sample_name"],
            route_key=route_key,
            io_cfg=resolved["io"],
        )
        expected_manifest_path = expected_run_dir / "provenance.json"
    except (AttributeError, KeyError, TypeError, ValueError):
        # A pre-context configuration error has no durable run envelope.
        pass
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
        run_dir = svc.provenance.get("run_dir")
        return {
            "ok": True,
            "profile": profile,
            "seed": runtime_seed,
            "run_dir": run_dir,
            "manifest_path": (
                str(Path(run_dir) / "provenance.json")
                if run_dir is not None
                else None
            ),
            "assignment_guidance": svc.provenance.get(
                "assignment_guidance"
            ),
            "summary": summary,
            "error": None,
        }
    except Exception as exc:  # pragma: no cover - wrapper behavior
        guidance_evidence = None
        failure_context = getattr(exc, "_revise_failure_context", None)
        current_identity = None
        manifest = None
        if expected_manifest_path is not None:
            snapshot = _read_manifest_snapshot(expected_manifest_path)
            if snapshot is not None:
                current_identity, manifest = snapshot
        current_manifest = bool(
            isinstance(failure_context, dict)
            and expected_run_dir is not None
            and expected_manifest_path is not None
            and failure_context.get("run_dir") == str(expected_run_dir)
            and failure_context.get("manifest_path")
            == str(expected_manifest_path)
            and failure_context.get("manifest_identity") == current_identity
        )
        if current_manifest:
            guidance_evidence = manifest.get("assignment_guidance")
        return {
            "ok": False,
            "profile": profile,
            "seed": runtime_seed,
            "run_dir": (
                str(expected_run_dir)
                if current_manifest
                else None
            ),
            "manifest_path": (
                str(expected_manifest_path)
                if current_manifest
                else None
            ),
            "assignment_guidance": guidance_evidence,
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
    profile = BENCHMARK_PROFILES.get(args.confounding)
    if profile is None:
        raise NotImplementedError(f"Unsupported confounding: {args.confounding}")
    report_aliases = _resolved_report_aliases(
        raw_config=pipeline.raw_config,
        profile=profile,
        platform=args.platform,
        confounding=args.confounding,
        algorithm_overrides=algorithm_overrides,
    )

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
                profile=profile,
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
                profile=profile,
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
                    profile=profile,
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
                profile=profile,
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
            profile=profile,
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
            profile=profile,
            confounding="gene_dropout",
            io_overrides=io_cfg,
            runtime_seed=next_runtime_seed(),
            algorithm_overrides=algorithm_overrides,
        )
        _append_result(results, "gene_dropout", run_result)

    else:
        raise NotImplementedError(f"Unsupported confounding: {args.confounding}")

    report_aliases["assignment_guidance"] = _aggregate_assignment_guidance(
        results,
        request_manifest=report_aliases["assignment_guidance"],
    )
    report = {
        "platform": args.platform,
        "confounding": args.confounding,
        "sample_name": args.sample_name,
        "dataset_task": args.dataset_task,
        "output_task": args.dataset_task or args.confounding,
        "seed": args.seed,
        "seed_scope": seed_scope,
        **report_aliases,
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
