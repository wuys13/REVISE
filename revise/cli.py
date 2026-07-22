#!/usr/bin/env python3
"""REVISE 统一重建脚本。

本脚本通过 ``revise-svc`` 提供的 :class:`REVISEPipeline`，使用同一套命令行
接口处理三类空间转录组输入：

* hST（高分辨率 ST）：执行 sp-SVC 重建；
* iST（成像型 ST）：执行 sc-SVC 重建，并将空间结果与表达结果合并为单个
  ``SVC.h5ad``；
* sST（spot-based ST）：执行超分辨率 sc-SVC 重建。

显式传入 ``--ot-method`` 时，它会同时控制全局注释和局部重建所使用的 OT
实现，可选 ``pot`` 或 ``tacco``；省略时保留配置中两阶段各自的 solver。

iST 的内部流程会生成空间和表达两个 AnnData。最终文件严格遵循以下规则：

1. 行数、``obs``、``obsm`` 和 ``obsp`` 来自空间结果；
2. 基因、``var``、``varm`` 和 ``varp`` 来自表达结果；
3. 两侧 ``SVC_cluster`` 集合必须完全一致，否则立即报错；
4. 表达映射可使用同 cluster 均值（默认）或带随机种子的随机单细胞表达。

最终文件统一写入 ``<output-root>/<sample-name>/SVC.h5ad``，运行 manifest
记录结果是 sp-SVC 还是 sc-SVC。Pipeline 自身的中间 H5AD 持久化会被关闭，
但运行日志和 provenance 仍由 REVISE 保存。
"""

from __future__ import annotations

import argparse
import copy
import json
import os
import sys
import tempfile
from contextlib import redirect_stdout
from pathlib import Path
from typing import Iterable

from revise._version import __version__
from revise.framework import REVISEPipeline
from revise.utils import completed_artifact

import numpy as np
from anndata import AnnData
from scipy import sparse


# 每个平台对应：配置 profile、confounding route、Pipeline 主要输出 artifact。
# iST 的第三项为 None，因为它需要同时读取 spatial/expr 两个 artifact 并合并。
ROUTES = {
    "hST": ("application_sp", "bin2cell", "sp_svc"),
    "iST": ("application_sc", "segmentation", None),
    "sST": ("application_sc_sst", "spot_size", "sc_svc_dec"),
}
PUBLIC_SVC_TYPES = {
    "sp": "sp-SVC",
    "sc": "sc-SVC",
}
PLATFORM_SVC_KINDS = {
    "hST": "sp",
    "iST": "sc",
    "sST": "sc",
}
# 这些配置由高层 CLI 参数统一管理。用户不能再通过 --set 覆盖它们，避免
# 命令行显示使用 POT、实际配置却被改成 TACCO 一类的歧义。
MANAGED_SET_KEYS = {
    "io.save_outputs",
    "runtime.seed",
}
IST_MANAGED_SET_KEYS = {
    "sc.select_ct",
    "columns.cell_type_col",
    "columns.sub_cell_type_col",
}


def _copy_value(value):
    """复制 AnnData 轴元数据，避免输出与 Pipeline 内部对象共享可变状态。"""
    return value.copy() if hasattr(value, "copy") else copy.deepcopy(value)


def _cluster_labels(adata: AnnData, cluster_col: str, source: str) -> np.ndarray:
    """读取 cluster 标签，并在进入映射前拒绝缺失列或空标签。"""
    if cluster_col not in adata.obs:
        raise KeyError(f"{source} is missing required obs column {cluster_col!r}")
    labels = adata.obs[cluster_col]
    if labels.isna().any():
        raise ValueError(f"{source} contains null {cluster_col!r} labels")
    return labels.to_numpy(dtype=object)


def _cluster_key(label) -> tuple[type, object]:
    """Keep type identity so values such as 1, True, and '1' stay distinct."""
    return type(label), label


def _cluster_keys(labels: np.ndarray) -> list[tuple[type, object]]:
    """将标签转换为同时包含类型和值的比较键。"""
    return [_cluster_key(label) for label in labels]


def _matching_indices(keys: list[tuple[type, object]], cluster) -> np.ndarray:
    """返回指定 typed cluster 在表达矩阵中的全部行号。"""
    return np.fromiter(
        (index for index, key in enumerate(keys) if key == cluster),
        dtype=np.int64,
    )


def _cluster_set_mismatch_message(
    spatial_labels: np.ndarray,
    expr_labels: np.ndarray,
) -> str | None:
    """比较空间侧和表达侧的完整 cluster 集合并生成可诊断错误信息。"""
    spatial_set = set(_cluster_keys(spatial_labels))
    expr_set = set(_cluster_keys(expr_labels))
    if spatial_set == expr_set:
        return None
    sort_key = lambda item: (item[0].__name__, repr(item[1]))
    spatial_only = sorted(spatial_set - expr_set, key=sort_key)
    expr_only = sorted(expr_set - spatial_set, key=sort_key)
    return (
        "iST SVC cluster sets do not match exactly: "
        f"spatial_only={spatial_only} (n={len(spatial_only)}), "
        f"expr_only={expr_only} (n={len(expr_only)})"
    )


def _mean_rows_by_cluster(matrix, labels, clusters):
    """按 cluster 计算表达均值，同时保留稀疏矩阵表示。"""
    rows = []
    if sparse.issparse(matrix):
        # 对稀疏输入逐 cluster 求和，避免先把全部 single-cell 表达转成稠密矩阵。
        matrix = matrix.tocsr()
        for cluster in clusters:
            mask = _matching_indices(labels, cluster)
            rows.append(matrix[mask].sum(axis=0) / float(mask.size))
        return sparse.vstack([sparse.csr_matrix(row) for row in rows], format="csr")

    dense = np.asarray(matrix)
    return np.vstack(
        [dense[_matching_indices(labels, cluster)].mean(axis=0) for cluster in clusters]
    )


def merge_ist_svc(
    spatial_adata: AnnData,
    expr_adata: AnnData,
    *,
    mode: str = "mean",
    seed: int = 42,
    cluster_col: str = "SVC_cluster",
) -> AnnData:
    """将 iST 表达结果映射到空间行，生成唯一的最终 AnnData。

    Parameters
    ----------
    spatial_adata
        决定最终细胞行、空间坐标和细胞级元数据的 sc-SVC spatial 结果。
    expr_adata
        提供完整基因表达、基因注释及 cluster 候选细胞的 sc-SVC expr 结果。
    mode
        ``mean`` 为每个空间细胞赋予同 cluster 的平均表达；``random`` 为每个
        空间细胞独立抽取一个同 cluster 的表达行。
    seed
        ``random`` 模式使用的 NumPy 随机种子，保证结果可复现。
    cluster_col
        两个输入中用于建立映射的 ``obs`` 列，默认 ``SVC_cluster``。

    Returns
    -------
    AnnData
        行轴来自 spatial、列轴来自 expr 的合并结果。
    """
    # 先检查模式和基因轴。重复 gene name 会导致 AnnData 列语义不唯一，不能
    # 等到写文件时再静默接受。
    if mode not in {"mean", "random"}:
        raise ValueError("iST mapping mode must be one of ['mean', 'random']")
    if not expr_adata.var_names.is_unique:
        duplicates = expr_adata.var_names[expr_adata.var_names.duplicated()].unique()
        raise ValueError(
            "expression SVC contains duplicate gene names: "
            f"{duplicates[:10].tolist()}"
        )

    # cluster 不仅要求显示值一致，也要求类型一致。例如整数 1 不等于字符串
    # "1"，从而避免不同数据源编码错误被静默合并。
    spatial_labels = _cluster_labels(spatial_adata, cluster_col, "spatial SVC")
    expr_labels = _cluster_labels(expr_adata, cluster_col, "expression SVC")
    mismatch = _cluster_set_mismatch_message(spatial_labels, expr_labels)
    if mismatch:
        raise ValueError(mismatch)

    spatial_keys = _cluster_keys(spatial_labels)
    expr_keys = _cluster_keys(expr_labels)
    clusters = sorted(
        set(expr_keys), key=lambda item: (item[0].__name__, repr(item[1]))
    )
    if mode == "mean":
        # 先生成 cluster × gene 均值矩阵，再按照每个空间行的 cluster 展开。
        # 这样最终行数严格等于 spatial_adata.n_obs。
        cluster_means = _mean_rows_by_cluster(expr_adata.X, expr_keys, clusters)
        cluster_index = {cluster: index for index, cluster in enumerate(clusters)}
        X = cluster_means[[cluster_index[key] for key in spatial_keys]]
    else:
        # 每个空间行独立抽取一个同 cluster 表达行；default_rng 避免修改全局
        # NumPy RNG 状态，并使相同 seed 下的映射可复现。
        rng = np.random.default_rng(seed)
        candidates = {
            cluster: _matching_indices(expr_keys, cluster) for cluster in clusters
        }
        selected = np.fromiter(
            (rng.choice(candidates[key]) for key in spatial_keys),
            dtype=np.int64,
            count=spatial_adata.n_obs,
        )
        X = expr_adata.X[selected].copy()

    # H5AD 可以写入 NaN/Inf，但下游分析会产生隐蔽错误，因此在组装前明确拒绝。
    values = X.data if sparse.issparse(X) else np.asarray(X)
    if not np.all(np.isfinite(values)):
        raise ValueError("mapped iST expression contains non-finite values")

    # X 的行来自 cluster 映射，obs/空间图来自 spatial，var/基因图来自 expr。
    # 不复制 spatial.X 或 expr.obs，因为它们分别不属于最终结果对应的轴语义。
    merged = AnnData(
        X=X,
        obs=spatial_adata.obs.copy(),
        var=expr_adata.var.copy(),
        uns=copy.deepcopy(spatial_adata.uns),
    )
    merged.obs_names = spatial_adata.obs_names.copy()
    merged.var_names = expr_adata.var_names.copy()
    for key, value in spatial_adata.obsm.items():
        merged.obsm[key] = _copy_value(value)
    for key, value in spatial_adata.obsp.items():
        merged.obsp[key] = _copy_value(value)
    for key, value in expr_adata.varm.items():
        merged.varm[key] = _copy_value(value)
    for key, value in expr_adata.varp.items():
        merged.varp[key] = _copy_value(value)

    # 写入足够的映射 provenance，便于结果文件脱离运行目录后仍可追溯。
    merged.uns["revise_reconstruction"] = {
        "platform": "iST",
        "mapping_mode": mode,
        "mapping_seed": int(seed),
        "cluster_col": cluster_col,
        "spatial_n_obs": int(spatial_adata.n_obs),
        "expression_n_obs": int(expr_adata.n_obs),
        "n_vars": int(expr_adata.n_vars),
        "cluster_set_match": True,
    }
    return merged


def _override_keys(overrides: Iterable[str]) -> set[str]:
    """从重复的 ``--set KEY=VALUE`` 参数中提取 dotted config key。"""
    return {
        item.split("=", 1)[0].strip()
        for item in overrides
        if "=" in item and item.split("=", 1)[0].strip()
    }


def _build_set_overrides(args: argparse.Namespace) -> list[str]:
    """构造 Pipeline 配置覆盖，并阻止高层参数与 --set 相互冲突。"""
    managed = set(MANAGED_SET_KEYS)
    values = ["io.save_outputs=false"]
    if args.ot_method is not None:
        managed.update({"ot.ga.solver", "ot.lr.solver"})
        values.extend(
            [
                f"ot.ga.solver={args.ot_method}",
                f"ot.lr.solver={args.ot_method}",
            ]
        )
    if args.platform == "iST":
        # iST 默认重建全部 cell type，同时把用户指定的列名传给 package。
        managed.update(IST_MANAGED_SET_KEYS)
        values.extend(
            [
                f"sc.select_ct={args.select_ct}",
                f"columns.cell_type_col={args.cell_type_col}",
                f"columns.sub_cell_type_col={args.sub_cell_type_col}",
            ]
        )
    user_keys = _override_keys(args.set_overrides)
    # 同时检查完全相同、父级和子级 key。例如 ot.ga={...} 也会覆盖
    # ot.ga.solver，因此必须和直接写 solver 一样被拒绝。
    conflicts = sorted(
        user_key
        for user_key in user_keys
        if any(
            user_key == managed_key
            or user_key.startswith(f"{managed_key}.")
            or managed_key.startswith(f"{user_key}.")
            for managed_key in managed
        )
    )
    if conflicts:
        raise ValueError(
            "Conflicting high-level CLI option and --set override for: "
            + ", ".join(conflicts)
        )
    return values + list(args.set_overrides)


def _run_pipeline(
    args: argparse.Namespace,
    *,
    dry_run: bool = False,
    finalize_callback=None,
):
    """Use the exact same route and input overrides for preflight and execution."""
    # ROUTES 是本脚本唯一的平台分派表，避免在后续流程重复判断 profile/route。
    profile, confounding, output_key = ROUTES[args.platform]
    pipeline = REVISEPipeline(config_path=args.config)
    # 配置解析顺序仍由 package 负责：defaults -> profile -> runtime/io -> --set。
    # save_outputs=false 只关闭中间 H5AD；本函数会在末尾自行写出规范文件名。
    svc = pipeline.run(
        profile=profile,
        runtime_overrides={
            "platform": args.platform,
            "confounding": confounding,
            "seed": args.seed,
        },
        io_overrides={
            "data_root": args.data_root,
            "output_root": args.output_root,
            "sample_name": args.sample_name,
            "st_file": args.st_file,
            "sc_ref_file": args.sc_ref_file,
            "patient_key": args.patient_key,
        },
        set_overrides=_build_set_overrides(args),
        dry_run=dry_run,
        finalize_callback=finalize_callback,
    )
    return profile, output_key, svc


def _build_public_result(args, profile, output_key, ctx) -> tuple[AnnData, Path]:
    """Build and publish the canonical result inside the finalize stage."""
    svc = ctx.svc
    if args.platform not in PLATFORM_SVC_KINDS:
        raise ValueError(f"Unsupported platform: {args.platform!r}")
    expected_kind = PLATFORM_SVC_KINDS[args.platform]
    if svc.svc_kind != expected_kind:
        raise ValueError(
            f"Platform {args.platform!r} requires SVC type {expected_kind!r}; "
            f"strategy returned {svc.svc_kind!r}"
        )
    result_type = PUBLIC_SVC_TYPES[expected_kind]
    # 所有平台都从 SVC 标准载体的 artifacts 取结果，不依赖 run_dir 中的临时文件。
    outputs = dict(svc.artifacts.get("outputs", {}))

    if args.platform == "iST":
        # iST 必须同时有空间与表达结果；缺少任一项时不能生成语义完整的文件。
        required = {"sc_svc_spatial", "sc_svc_expr"}
        missing = sorted(required - outputs.keys())
        if missing:
            raise RuntimeError(
                f"iST pipeline did not return required outputs {missing}; "
                f"available={sorted(outputs)}"
            )
        result = merge_ist_svc(
            outputs["sc_svc_spatial"],
            outputs["sc_svc_expr"],
            mode=args.ist_mapping,
            seed=args.seed,
        )
    else:
        # hST/sST 已由对应 strategy 生成单一主要 AnnData，只需校验并复制。
        if output_key not in outputs:
            raise RuntimeError(
                f"{args.platform} pipeline did not return required output {output_key!r}; "
                f"available={sorted(outputs)}"
            )
        result = outputs[output_key].copy()
        result.uns["revise_reconstruction"] = {
            "platform": args.platform,
            "seed": int(args.seed),
        }

    output_dir = Path(args.output_root) / args.sample_name
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "SVC.h5ad"
    relative_run_dir = Path(os.path.relpath(ctx.run_dir, start=output_dir)).as_posix()

    # 三个平台统一补充运行模式和可随 output tree 移动的 manifest 链接。
    provenance = result.uns.setdefault("revise_reconstruction", {})
    provenance.update(
        {
            "profile": profile,
            "run_dir": relative_run_dir,
            "run_manifest": (Path(relative_run_dir) / "provenance.json").as_posix(),
            "ot_method_override": args.ot_method,
            "ot_config": copy.deepcopy(ctx.merged_config["ot"]),
            "ot_events": json.dumps(
                ctx.ot_events,
                ensure_ascii=False,
                sort_keys=True,
            ),
        }
    )

    # 先完整写入同目录临时文件；失败时保留此前成功发布的 canonical 结果。
    temporary_path = None
    backup_path = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=output_dir,
            prefix=f".{output_path.stem}.",
            suffix=".h5ad",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
        result.write_h5ad(temporary_path)
        artifact = completed_artifact("public_result", temporary_path)
        artifact["path"] = str(output_path)
        result_record = {
            "filename": output_path.name,
            "type": result_type,
        }
        had_previous_result = "result" in ctx.provenance
        previous_result = copy.deepcopy(ctx.provenance.get("result"))
        had_previous_svc_result = "result" in ctx.svc.provenance
        previous_svc_result = copy.deepcopy(ctx.svc.provenance.get("result"))

        if output_path.exists():
            with tempfile.NamedTemporaryFile(
                dir=output_dir,
                prefix=f".{output_path.stem}.previous.",
                suffix=output_path.suffix,
                delete=False,
            ) as handle:
                backup_path = Path(handle.name)
            backup_path.unlink()

        def commit():
            if backup_path is not None:
                backup_path.unlink(missing_ok=True)

        def rollback():
            if backup_path is None:
                output_path.unlink(missing_ok=True)
            elif backup_path.exists():
                os.replace(backup_path, output_path)

            if not had_previous_result:
                ctx.provenance.pop("result", None)
            else:
                ctx.provenance["result"] = previous_result
            if not had_previous_svc_result:
                ctx.svc.provenance.pop("result", None)
            else:
                ctx.svc.provenance["result"] = previous_svc_result

            for index in range(len(ctx.artifact_records) - 1, -1, -1):
                if ctx.artifact_records[index] == artifact:
                    del ctx.artifact_records[index]
                    break

        ctx.set_pending_publication(commit=commit, rollback=rollback)
        try:
            if backup_path is not None:
                os.replace(output_path, backup_path)
            os.replace(temporary_path, output_path)
            ctx.provenance["result"] = result_record
            ctx.record_artifact(artifact)
        except BaseException:
            ctx.rollback_pending_publication()
            raise
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
    return result, output_path


def reconstruct(args: argparse.Namespace) -> tuple[AnnData, Path, dict]:
    """运行选定平台的 REVISE Pipeline，整理并写出唯一的最终 H5AD。"""
    published = {}

    def publish(ctx):
        profile, _, output_key = ROUTES[args.platform]
        published["result"], published["path"] = _build_public_result(
            args,
            profile,
            output_key,
            ctx,
        )

    profile, _, svc = _run_pipeline(args, finalize_callback=publish)
    summary = svc.summary()
    summary.update(profile=profile, route=svc.provenance.get("route_key"))
    return published["result"], published["path"], summary


def get_args() -> argparse.Namespace:
    """定义统一 CLI；平台专用参数在不适用的平台上会被安全忽略。"""
    parser = argparse.ArgumentParser(
        description="Reconstruct one SVC from hST, iST, or sST data through revise-svc"
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"revise-reconstruct {__version__}",
    )
    parser.add_argument("--platform", required=True, choices=list(ROUTES))
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--st-file", required=True)
    parser.add_argument("--sc-ref-file", required=True)
    parser.add_argument("--data-root", required=True)
    parser.add_argument("--output-root", default="output/reconstruct")
    parser.add_argument("--config", default="revise/revise.yaml")
    parser.add_argument("--ot-method", choices=["pot", "tacco"], default=None)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--patient-key", default="Patient")
    parser.add_argument("--select-ct", default="all")
    parser.add_argument("--cell-type-col", default="Level1")
    parser.add_argument("--sub-cell-type-col", default="Level2")
    parser.add_argument("--ist-mapping", choices=["mean", "random"], default="mean")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate resolved inputs and dependencies without reconstruction",
    )
    parser.add_argument(
        "--set",
        dest="set_overrides",
        action="append",
        default=[],
        metavar="KEY=VALUE",
    )
    return parser.parse_args()


def main() -> None:
    """CLI 入口：执行重建并向终端打印机器可读的运行摘要。"""
    args = get_args()
    if args.dry_run:
        with redirect_stdout(sys.stderr):
            profile, _, svc = _run_pipeline(args, dry_run=True)
        run_dir = Path(svc.provenance["run_dir"])
        summary = svc.summary()
        summary.update(profile=profile, route=svc.provenance.get("route_key"))
        print(
            json.dumps(
                {
                    "status": "ready",
                    "platform": args.platform,
                    "preflight": str(run_dir / "preflight.json"),
                    "pipeline": summary,
                },
                indent=2,
                ensure_ascii=False,
            )
        )
        return
    with redirect_stdout(sys.stderr):
        result, output_path, pipeline_summary = reconstruct(args)
    print(
        json.dumps(
            {
                "platform": args.platform,
                "output": str(output_path),
                "shape": list(result.shape),
                "pipeline": pipeline_summary,
            },
            indent=2,
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
