"""Local implementations for the explicit Moran and EMT notebook stages.

These functions preserve the original calculations; they do not orchestrate
partition, spatial diversity, or Region analysis.
"""
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy import sparse
from revise.analysis.paired_moran import compare_moran
from revise.analysis.advanced.aucell import get_aucell_provider_metadata, score_gene_set_aucell


def prepare_moran(
    RUNS,
    CONFIG,
    *,
    global_partition=None,
    raw_global=None,
    recon_global=None,
    moran_min_units=51,
    moran_n_neighs=6,
):
    """Validate paired Moran inputs and return the visible preparation audit."""

    scope_order = []
    sources = {}
    coordinate_sources = {}
    if global_partition is not None and raw_global is not None and recon_global is not None:
        scope_order.append("All")
        sources["All"] = (raw_global, recon_global, None)

    for parent in ("Fibroblast", "Mono_Macro", "T"):
        if parent not in RUNS:
            raise KeyError(f"Missing required Moran parent scope: {parent}")
        run = RUNS[parent]
        if "raw_expression" not in run or "recon_expression" not in run:
            raise KeyError(f"{parent}: RUNS must expose raw_expression and recon_expression")
        scope_order.append(parent)
        sources[parent] = (
            run["raw_expression"],
            run["recon_expression"],
            run.get("spatial_carrier"),
        )

    for scope in scope_order:
        raw_expr, recon_expr, carrier = sources[scope]
        if not raw_expr.obs_names.is_unique or not recon_expr.obs_names.is_unique:
            raise ValueError(f"{scope}: Moran observation IDs must be unique")
        if not raw_expr.obs_names.equals(recon_expr.obs_names):
            raise ValueError(f"{scope}: Raw and reconstruction observation IDs must match exactly")
        if raw_expr.n_obs == 0:
            raise ValueError(f"{scope}: Moran requires a nonempty paired observation axis")
        if not raw_expr.var_names.is_unique or not recon_expr.var_names.is_unique:
            raise ValueError(f"{scope}: Moran gene axes must be unique")

        if "spatial" in raw_expr.obsm:
            raw_xy = np.asarray(raw_expr.obsm["spatial"], dtype=float)
            coordinate_source = "raw_expression.obsm[spatial]"
        elif carrier is not None and "spatial" in carrier.obsm:
            if not carrier.obs_names.equals(raw_expr.obs_names):
                raise ValueError(f"{scope}: spatial carrier IDs do not match expression IDs")
            raw_xy = np.asarray(carrier.obsm["spatial"], dtype=float)
            coordinate_source = "spatial_carrier.obsm[spatial]"
        else:
            raise KeyError(f"{scope}: no Raw spatial coordinates are available")

        if "spatial" in recon_expr.obsm:
            recon_xy = np.asarray(recon_expr.obsm["spatial"], dtype=float)
        elif carrier is not None and "spatial" in carrier.obsm:
            if not carrier.obs_names.equals(recon_expr.obs_names):
                raise ValueError(f"{scope}: reconstruction carrier IDs do not match expression IDs")
            recon_xy = np.asarray(carrier.obsm["spatial"], dtype=float)
        else:
            raise KeyError(f"{scope}: no reconstruction spatial coordinates are available")

        if raw_xy.shape != recon_xy.shape or not np.isfinite(raw_xy).all() or not np.isfinite(recon_xy).all():
            raise ValueError(f"{scope}: Moran coordinates must be finite and shape-compatible")
        if not np.array_equal(raw_xy, recon_xy):
            raise ValueError(f"{scope}: Raw and reconstruction coordinates must agree exactly")
        if raw_xy.shape != (raw_expr.n_obs, 2):
            raise ValueError(f"{scope}: Moran expects two-dimensional spatial coordinates")

        raw_expr.obsm["spatial"] = raw_xy.copy()
        recon_expr.obsm["spatial"] = recon_xy.copy()
        coordinate_sources[scope] = coordinate_source

    spatial_config = CONFIG.get("spatial_region", {})
    coordinate_unit = spatial_config.get(
        "coordinate_unit", CONFIG.get("context", {}).get("coordinate_unit", "unknown")
    )
    coordinate_to_microns = spatial_config.get("microns_per_coordinate", np.nan)
    audit = pd.DataFrame([
        {
            "scope": scope,
            "units": sources[scope][0].n_obs,
            "raw_genes": sources[scope][0].n_vars,
            "reconstruction_genes": sources[scope][1].n_vars,
            "coordinate_axis": "paired Raw spatial IDs",
            "coordinate_unit": coordinate_unit,
            "coordinate_to_microns": coordinate_to_microns,
            "min_units": moran_min_units,
            "requested_neighbors": moran_n_neighs,
        }
        for scope in scope_order
    ])
    return scope_order, sources, coordinate_sources, coordinate_unit, coordinate_to_microns, audit


def summarize_moran(MORAN, MORAN_SOURCES, MORAN_SCOPE_ORDER):
    """Build full-gene status, all-valid, and shared-valid Moran summaries."""

    availability_rows = []
    summary_rows = []
    for scope in MORAN_SCOPE_ORDER:
        raw_expr, recon_expr, _carrier = MORAN_SOURCES[scope]
        scope_table = MORAN.loc[MORAN["scope"].eq(scope)].copy()
        raw_genes = set(raw_expr.var_names.astype(str))
        recon_genes = set(recon_expr.var_names.astype(str))
        availability = scope_table.loc[:, [
            "scope", "comparison_id", "gene_id", "raw_status", "raw_reason",
            "reconstruction_status", "reconstruction_reason", "comparison_status",
        ]].copy()
        availability["gene_id"] = availability["gene_id"].astype(str)
        availability["raw_provided"] = availability["gene_id"].isin(raw_genes)
        availability["reconstruction_provided"] = availability["gene_id"].isin(recon_genes)
        availability_rows.append(availability)

        raw_valid = scope_table.loc[scope_table["raw_status"].eq("computed"), "raw_moran"].dropna()
        recon_valid = scope_table.loc[
            scope_table["reconstruction_status"].eq("computed"), "reconstruction_moran"
        ].dropna()
        matched = scope_table.loc[
            scope_table["raw_status"].eq("computed")
            & scope_table["reconstruction_status"].eq("computed")
        ].copy()
        matched = matched.loc[
            np.isfinite(matched["raw_moran"]) & np.isfinite(matched["reconstruction_moran"])
        ].copy()
        paired_delta = matched["delta"].dropna().astype(float)
        summary_rows.append({
            "scope": scope,
            "comparison_id": scope_table["comparison_id"].iloc[0],
            "union_gene_n": int(len(scope_table)),
            "raw_nvalid": int(raw_valid.size),
            "raw_median": float(raw_valid.median()) if len(raw_valid) else np.nan,
            "raw_mean": float(raw_valid.mean()) if len(raw_valid) else np.nan,
            "raw_q1": float(raw_valid.quantile(0.25)) if len(raw_valid) else np.nan,
            "raw_q75": float(raw_valid.quantile(0.75)) if len(raw_valid) else np.nan,
            "reconstruction_nvalid": int(recon_valid.size),
            "reconstruction_median": float(recon_valid.median()) if len(recon_valid) else np.nan,
            "reconstruction_mean": float(recon_valid.mean()) if len(recon_valid) else np.nan,
            "reconstruction_q1": float(recon_valid.quantile(0.25)) if len(recon_valid) else np.nan,
            "reconstruction_q75": float(recon_valid.quantile(0.75)) if len(recon_valid) else np.nan,
            "matched_gene_n": int(len(matched)),
            "shared_raw_nvalid": int(len(matched)),
            "shared_raw_median": float(matched["raw_moran"].median()) if len(matched) else np.nan,
            "shared_raw_q1": float(matched["raw_moran"].quantile(0.25)) if len(matched) else np.nan,
            "shared_raw_q75": float(matched["raw_moran"].quantile(0.75)) if len(matched) else np.nan,
            "shared_reconstruction_nvalid": int(len(matched)),
            "shared_reconstruction_median": float(matched["reconstruction_moran"].median()) if len(matched) else np.nan,
            "shared_reconstruction_q1": float(matched["reconstruction_moran"].quantile(0.25)) if len(matched) else np.nan,
            "shared_reconstruction_q75": float(matched["reconstruction_moran"].quantile(0.75)) if len(matched) else np.nan,
            "paired_delta_n": int(len(paired_delta)),
            "paired_delta_mean": float(paired_delta.mean()) if len(paired_delta) else np.nan,
            "paired_delta_median": float(paired_delta.median()) if len(paired_delta) else np.nan,
            "paired_delta_q1": float(paired_delta.quantile(0.25)) if len(paired_delta) else np.nan,
            "paired_delta_q3": float(paired_delta.quantile(0.75)) if len(paired_delta) else np.nan,
            "raw_unmeasured_n": int((scope_table["raw_status"] == "unmeasured").sum()),
            "raw_not_computable_n": int((scope_table["raw_status"] == "not_computable").sum()),
            "raw_insufficient_support_n": int((scope_table["raw_status"] == "insufficient_support").sum()),
            "reconstruction_unmeasured_n": int((scope_table["reconstruction_status"] == "unmeasured").sum()),
            "reconstruction_not_computable_n": int((scope_table["reconstruction_status"] == "not_computable").sum()),
            "reconstruction_insufficient_support_n": int((scope_table["reconstruction_status"] == "insufficient_support").sum()),
        })

    gene_availability = pd.concat(availability_rows, ignore_index=True)
    summary = pd.DataFrame(summary_rows)
    display_rows = []
    for row in summary.to_dict("records"):
        for gene_set, prefix in (("all_valid", ""), ("shared_valid", "shared_")):
            for side, side_prefix in (("Raw", "raw_"), ("Reconstruction", "reconstruction_")):
                field_prefix = prefix + side_prefix
                display_rows.append({
                    "scope": row["scope"],
                    "comparison_id": row["comparison_id"],
                    "gene_set": gene_set,
                    "side": side,
                    "n_valid": row[field_prefix + "nvalid"],
                    "median": row[field_prefix + "median"],
                    "q1": row[field_prefix + "q1"],
                    "q75": row[field_prefix + "q75"],
                })
    return gene_availability, summary, pd.DataFrame(display_rows)


def build_moran_distribution(MORAN, MORAN_SCOPE_ORDER):
    """Build all-valid and shared-valid long-form Moran plotting rows."""

    rows = []
    for scope in MORAN_SCOPE_ORDER:
        table = MORAN.loc[MORAN["scope"].eq(scope)]
        for side, value, status in (("Raw", "raw_moran", "raw_status"), ("Reconstruction", "reconstruction_moran", "reconstruction_status")):
            valid = table.loc[table[status].eq("computed"), ["gene_id", value]].dropna()
            rows.extend({"scope": scope, "gene_set": "all_valid", "side": side, "gene_id": gene, "moran_i": float(value)} for gene, value in valid.itertuples(index=False, name=None))
        shared = table.loc[
            table["raw_status"].eq("computed") & table["reconstruction_status"].eq("computed"),
            ["gene_id", "raw_moran", "reconstruction_moran"],
        ].dropna()
        for side, value in (("Raw", "raw_moran"), ("Reconstruction", "reconstruction_moran")):
            rows.extend({"scope": scope, "gene_set": "shared_valid", "side": side, "gene_id": gene, "moran_i": float(value)} for gene, value in shared.loc[:, ["gene_id", value]].itertuples(index=False, name=None))
    return pd.DataFrame(rows, columns=["scope", "gene_set", "side", "gene_id", "moran_i"])


def build_moran_paired_delta(MORAN, MORAN_SCOPE_ORDER):
    """Build paired Moran deltas without changing the full-gene table."""

    rows = [
        MORAN.loc[
            MORAN["scope"].eq(scope)
            & MORAN["raw_status"].eq("computed")
            & MORAN["reconstruction_status"].eq("computed"),
            ["scope", "gene_id", "delta"],
        ].dropna()
        for scope in MORAN_SCOPE_ORDER
    ]
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["scope", "gene_id", "delta"])


def compute_emt(RUNS, ROOT, OUTPUT_DIR, CONFIG, *, AUCELL_AUC_THRESHOLD=0.01, AUCELL_SEED=42):
    # The pathway stage is deliberately explicit in the notebook.
    import hashlib
    import json
    from revise.analysis.basic.gene_set_scoring import read_gmt

    EMT_NAME = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
    GMT_PATH = ROOT / "reproduce" / "case" / "pathway" / "h.all.v2025.1.Hs.symbols.gmt"
    PATHWAY_OUTPUT_DIR = OUTPUT_DIR / "analysis" / "pathway_activity"
    PATHWAY_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not GMT_PATH.is_file():
        raise FileNotFoundError(f"AUCell resource does not exist: {GMT_PATH}")
    resource_sha256 = hashlib.sha256(GMT_PATH.read_bytes()).hexdigest()
    gene_sets = read_gmt(GMT_PATH)
    if EMT_NAME not in gene_sets:
        raise KeyError(f"AUCell resource is missing {EMT_NAME!r}")
    # Preserve resource order while removing duplicate symbols before matching.
    EMT_GENES = list(dict.fromkeys(str(gene) for gene in gene_sets[EMT_NAME]))
    if not EMT_GENES:
        raise ValueError(f"AUCell resource {EMT_NAME!r} contains no genes")

    # Provider metadata is part of the audit contract.  Provider/scorer failures
    # propagate; only known scientific non-computability is represented as NA.
    AUCell_PROVIDER = get_aucell_provider_metadata()
    AUCell_CUTOFF_AUDIT = []
    AUCell_AVAIL_ROWS = []
    AUCell_RESOURCE_GENE_ROWS = []
    AUCell_SCORE_FRAMES = []
    AUCell_SIDE_RESULTS = {}

    for scope, run in RUNS.items():
        raw_view = run["raw_expression"]
        reconstruction_view = run["recon_expression"]
        raw_ids = pd.Index(raw_view.obs_names)
        reconstruction_ids = pd.Index(reconstruction_view.obs_names)
        if not raw_ids.is_unique or not reconstruction_ids.is_unique:
            raise ValueError(f"{scope}: AUCell observation IDs must be unique")
        if not raw_ids.equals(reconstruction_ids):
            raise ValueError(f"{scope}: Raw and reconstruction AUCell axes must match exactly")

        scope_side_results = {}
        for side, view in (("raw", raw_view), ("reconstruction", reconstruction_view)):
            expression = view.X
            n_observations = int(view.n_obs)
            n_genes = int(view.n_vars)
            if n_observations < 1 or n_genes < 1:
                raise ValueError(f"{scope}/{side}: AUCell expression view must be non-empty")
            if not pd.Index(view.var_names).is_unique:
                raise ValueError(f"{scope}/{side}: AUCell gene IDs must be unique")

            # This is equivalent to OmicVerse's derive_auc_threshold, including
            # explicit sparse zeros.  No normalization or overlap clipping occurs.
            if sparse.issparse(expression):
                detected_counts = np.asarray((expression != 0).sum(axis=1)).ravel().astype(np.int64, copy=False)
                detected_anywhere = np.asarray((expression != 0).sum(axis=0)).ravel().astype(bool, copy=False)
            else:
                expression_array = np.asarray(expression)
                detected_counts = np.count_nonzero(expression_array, axis=1).astype(np.int64, copy=False)
                detected_anywhere = np.any(expression_array != 0, axis=0)
            detected_count_quantile = float(pd.Series(detected_counts).quantile(AUCELL_AUC_THRESHOLD))
            provider_auc_threshold = float(detected_count_quantile / float(n_genes))
            effective_rank_length = int(round(provider_auc_threshold * n_genes))
            rank_cutoff_zero_based = effective_rank_length - 1

            var_names = pd.Index(str(gene) for gene in view.var_names)
            available_genes = [gene for gene in EMT_GENES if gene in var_names]
            available_indices = [int(var_names.get_loc(gene)) for gene in available_genes]
            detected_signature_gene_count = int(detected_anywhere[available_indices].sum()) if available_indices else 0
            coverage = float(len(available_genes) / len(EMT_GENES))

            status = "computed"
            reason = ""
            score_values = None
            scored = None
            score_key = f"{EMT_NAME}_aucell"
            if not available_genes:
                status = "unmeasured"
                reason = "resource_genes_unmeasured"
            elif int(detected_counts.max()) == 0:
                status = "not_computable"
                reason = "zero_expression_no_detectable_signature"
            elif effective_rank_length <= 0:
                status = "not_computable"
                reason = "zero_detected_gene_rank_cutoff"
            else:
                # Direct scorer call: it receives the original full-gene X.  The
                # wrapper copies/CSR-converts only for the provider and does not
                # normalize or crop the expression view.
                scored, score_key = score_gene_set_aucell(
                    view,
                    available_genes,
                    score_name=EMT_NAME,
                    AUC_threshold=AUCELL_AUC_THRESHOLD,
                    seed=AUCELL_SEED,
                )
                if score_key not in scored.obs:
                    raise RuntimeError(f"{scope}/{side}: AUCell provider did not create {score_key!r}")
                score_series = pd.Series(
                    np.asarray(scored.obs[score_key], dtype=float),
                    index=pd.Index(scored.obs_names),
                ).reindex(raw_ids)
                score_values = score_series.to_numpy(dtype=float)
                # score_gene_set_aucell already validates finiteness; keep an
                # explicit shape check at the notebook boundary as well.
                if score_values.shape != (n_observations,) or not np.isfinite(score_values).all():
                    raise RuntimeError(f"{scope}/{side}: AUCell provider returned invalid scores")

            cutoff_row = {
                "scope": str(scope),
                "side": side,
                "n_observations": n_observations,
                "n_genes": n_genes,
                "detected_count_quantile": detected_count_quantile,
                "provider_auc_threshold": provider_auc_threshold,
                "effective_rank_length": effective_rank_length,
                "rank_cutoff_zero_based": rank_cutoff_zero_based,
                "cutoff_formula": "quantile(count_nonzero(X), 0.01) / n_genes",
                "rank_length_formula": "int(round(provider_auc_threshold * n_genes))",
                "status": status,
                "reason": reason,
            }
            AUCell_CUTOFF_AUDIT.append(cutoff_row)
            AUCell_AVAIL_ROWS.append({
                "scope": str(scope),
                "side": side,
                "pathway": EMT_NAME,
                "resource_gene_count": len(EMT_GENES),
                "available_gene_count": len(available_genes),
                "coverage": coverage,
                "detected_signature_gene_count": detected_signature_gene_count,
                "n_observations": n_observations,
                "n_genes": n_genes,
                "detected_count_quantile": detected_count_quantile,
                "provider_auc_threshold": provider_auc_threshold,
                "effective_rank_length": effective_rank_length,
                "rank_cutoff_zero_based": rank_cutoff_zero_based,
                "status": status,
                "reason": reason,
            })
            scope_side_results[side] = {
                "status": status,
                "reason": reason,
                "scores": score_values,
                "available_genes": available_genes,
                "coverage": coverage,
                "detected_signature_gene_count": detected_signature_gene_count,
            }

        # Save the complete resource membership and the side-specific measured
        # state.  This is per scope because the two full gene axes can differ.
        raw_result = scope_side_results["raw"]
        reconstruction_result = scope_side_results["reconstruction"]
        raw_var_names = pd.Index(str(gene) for gene in raw_view.var_names)
        reconstruction_var_names = pd.Index(str(gene) for gene in reconstruction_view.var_names)
        for resource_order, gene in enumerate(EMT_GENES):
            raw_available = gene in raw_var_names
            reconstruction_available = gene in reconstruction_var_names
            raw_index = int(raw_var_names.get_loc(gene)) if raw_available else None
            reconstruction_index = int(reconstruction_var_names.get_loc(gene)) if reconstruction_available else None
            if raw_index is None:
                raw_detected = False
            elif sparse.issparse(raw_view.X):
                raw_detected = bool(np.asarray((raw_view.X[:, raw_index] != 0).sum()).ravel()[0])
            else:
                raw_detected = bool(np.any(np.asarray(raw_view.X)[:, raw_index] != 0))
            if reconstruction_index is None:
                reconstruction_detected = False
            elif sparse.issparse(reconstruction_view.X):
                reconstruction_detected = bool(np.asarray((reconstruction_view.X[:, reconstruction_index] != 0).sum()).ravel()[0])
            else:
                reconstruction_detected = bool(np.any(np.asarray(reconstruction_view.X)[:, reconstruction_index] != 0))
            AUCell_RESOURCE_GENE_ROWS.append({
                "scope": str(scope),
                "pathway": EMT_NAME,
                "resource_gene_order": resource_order,
                "resource_gene": gene,
                "raw_available": raw_available,
                "reconstruction_available": reconstruction_available,
                "raw_detected_anywhere": raw_detected,
                "reconstruction_detected_anywhere": reconstruction_detected,
            })

        raw_scores = raw_result["scores"]
        reconstruction_scores = reconstruction_result["scores"]
        status_parts = [raw_result["status"], reconstruction_result["status"]]
        if all(part == "computed" for part in status_parts):
            comparison_status = "computed"
        elif any(part == "computed" for part in status_parts):
            comparison_status = "partial"
        else:
            comparison_status = "unavailable"
        reason_parts = []
        if raw_result["reason"]:
            reason_parts.append(f"raw:{raw_result['reason']}")
        if reconstruction_result["reason"]:
            reason_parts.append(f"reconstruction:{reconstruction_result['reason']}")
        scope_scores = pd.DataFrame({
            "scope": str(scope),
            "unit_id": raw_ids.astype(str),
            "raw_score": raw_scores if raw_scores is not None else np.full(len(raw_ids), np.nan),
            "reconstruction_score": reconstruction_scores if reconstruction_scores is not None else np.full(len(raw_ids), np.nan),
            "raw_status": raw_result["status"],
            "reconstruction_status": reconstruction_result["status"],
            "raw_reason": raw_result["reason"],
            "reconstruction_reason": reconstruction_result["reason"],
            "status": comparison_status,
            "reason": ";".join(reason_parts),
        })
        scope_scores["delta"] = scope_scores["reconstruction_score"] - scope_scores["raw_score"]
        AUCell_SCORE_FRAMES.append(scope_scores)
        AUCell_SIDE_RESULTS[str(scope)] = scope_side_results

    AUCELL_COLUMNS = [
        "scope", "unit_id", "raw_score", "reconstruction_score", "delta",
        "raw_status", "reconstruction_status", "raw_reason", "reconstruction_reason",
        "status", "reason",
    ]
    AUCELL = pd.concat(AUCell_SCORE_FRAMES, ignore_index=True) if AUCell_SCORE_FRAMES else pd.DataFrame(columns=AUCELL_COLUMNS)
    AUCELL = AUCELL.loc[:, AUCELL_COLUMNS]
    AUCELL["pathway"] = EMT_NAME
    AUCELL["sample_id"] = CONFIG["sample"]["id"]
    AVAIL = pd.DataFrame(AUCell_AVAIL_ROWS)
    AUCELL_AVAIL = AVAIL
    RESOURCE_GENES = pd.DataFrame(AUCell_RESOURCE_GENE_ROWS)
    CUTOFF_AUDIT = pd.DataFrame(AUCell_CUTOFF_AUDIT)

    summary_rows = []
    for scope, group in AUCELL.groupby("scope", sort=False):
        raw_valid = group["raw_score"].dropna().astype(float)
        reconstruction_valid = group["reconstruction_score"].dropna().astype(float)
        paired_delta = group["delta"].dropna().astype(float)
        summary_rows.append({
            "scope": scope,
            "pathway": EMT_NAME,
            "raw_status": group["raw_status"].iloc[0] if len(group) else "unavailable",
            "reconstruction_status": group["reconstruction_status"].iloc[0] if len(group) else "unavailable",
            "comparison_status": group["status"].iloc[0] if len(group) else "unavailable",
            "raw_n_valid": int(raw_valid.size),
            "raw_mean": float(raw_valid.mean()) if len(raw_valid) else np.nan,
            "raw_median": float(raw_valid.median()) if len(raw_valid) else np.nan,
            "raw_q1": float(raw_valid.quantile(0.25)) if len(raw_valid) else np.nan,
            "raw_q3": float(raw_valid.quantile(0.75)) if len(raw_valid) else np.nan,
            "reconstruction_n_valid": int(reconstruction_valid.size),
            "reconstruction_mean": float(reconstruction_valid.mean()) if len(reconstruction_valid) else np.nan,
            "reconstruction_median": float(reconstruction_valid.median()) if len(reconstruction_valid) else np.nan,
            "reconstruction_q1": float(reconstruction_valid.quantile(0.25)) if len(reconstruction_valid) else np.nan,
            "reconstruction_q3": float(reconstruction_valid.quantile(0.75)) if len(reconstruction_valid) else np.nan,
            "paired_n": int(paired_delta.size),
            "paired_delta_mean": float(paired_delta.mean()) if len(paired_delta) else np.nan,
            "paired_delta_median": float(paired_delta.median()) if len(paired_delta) else np.nan,
            "paired_delta_q1": float(paired_delta.quantile(0.25)) if len(paired_delta) else np.nan,
            "paired_delta_q3": float(paired_delta.quantile(0.75)) if len(paired_delta) else np.nan,
        })
    aucell_summary = pd.DataFrame(summary_rows)

    # Publish all reusable evidence under one pathway aspect directory.  The
    # in-memory AUCELL table remains the source for the later spatial-field cell.
    AUCELL.to_csv(PATHWAY_OUTPUT_DIR / "scores.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    AVAIL.to_csv(PATHWAY_OUTPUT_DIR / "availability.csv", index=False)
    aucell_summary.to_csv(PATHWAY_OUTPUT_DIR / "summary.csv", index=False)
    RESOURCE_GENES.to_csv(PATHWAY_OUTPUT_DIR / "resource_genes.csv", index=False)
    CUTOFF_AUDIT.to_csv(PATHWAY_OUTPUT_DIR / "cutoff_audit.csv", index=False)
    resource_payload = {
        "pathway": EMT_NAME,
        "resource_name": "localHallmark2025.1",
        "resource_version": "2025.1",
        "resource_species": "Hs",
        "gene_id_type": "symbol",
        "resource_path": str(GMT_PATH),
        "resource_sha256": resource_sha256,
        "resource_gene_count": len(EMT_GENES),
        "provider": AUCell_PROVIDER,
        "scorer": "AUCell",
        "auc_threshold_quantile": AUCELL_AUC_THRESHOLD,
        "seed": AUCELL_SEED,
        "cutoff_policy": "provider_detected_gene_quantile",
        "cutoff_formula": "quantile(count_nonzero(X), 0.01) / n_genes",
        "rank_length_formula": "int(round(provider_auc_threshold * n_genes))",
        "side_cutoffs": CUTOFF_AUDIT.to_dict(orient="records"),
        "output_files": {
            "scores": "scores.csv.gz",
            "availability": "availability.csv",
            "summary": "summary.csv",
            "resource_genes": "resource_genes.csv",
            "cutoff_audit": "cutoff_audit.csv",
            "resource": "resource.json",
        },
    }
    (PATHWAY_OUTPUT_DIR / "resource.json").write_text(
        json.dumps(resource_payload, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )



    return AUCELL, AVAIL, aucell_summary, EMT_NAME, GMT_PATH

def compute_moran(MORAN_SOURCES, MORAN_SCOPE_ORDER, RUNS, CONFIG, OUTPUT_DIR,
                  MORAN_COORDINATE_SOURCES, _coordinate_unit, _coordinate_to_microns,
                  *, MORAN_MIN_UNITS=51, MORAN_N_NEIGHS=6):
    # Normalize both full-gene views and compare them on one shared Raw graph.
    # Row normalization follows sq.gr.spatial_autocorr(transformation=True).
    import json
    from scipy.sparse.csgraph import connected_components

    MORAN_TABLES = []
    MORAN_GRAPH_AUDIT = []
    MORAN_DIR = OUTPUT_DIR / "analysis" / "moran"
    MORAN_GRAPH_DIR = MORAN_DIR / "graphs"
    MORAN_GRAPH_DIR.mkdir(parents=True, exist_ok=True)

    for _scope in MORAN_SCOPE_ORDER:
        _raw_expr, _recon_expr, _carrier = MORAN_SOURCES[_scope]
        _raw_work = _raw_expr.copy()
        _recon_work = _recon_expr.copy()
        _n_units = int(_raw_work.n_obs)

        sc.pp.normalize_total(_raw_work, target_sum=1e4, inplace=True)
        sc.pp.log1p(_raw_work, copy=False)
        sc.pp.normalize_total(_recon_work, target_sum=1e4, inplace=True)
        sc.pp.log1p(_recon_work, copy=False)

        # Squidpy/scanpy can leave a numerically tiny denominator for a constant
        # column after sparse normalization.  Audit constant genes explicitly so
        # they remain not_computable instead of becoming a spurious finite Moran.
        _constant_genes = {"raw": set(), "reconstruction": set()}
        for _side, _work in (("raw", _raw_work), ("reconstruction", _recon_work)):
            for _start in range(0, _work.n_vars, 64):
                _block = _work.X[:, _start : _start + 64]
                _values = np.asarray(_block.toarray() if sparse.issparse(_block) else _block, dtype=float)
                _scale = np.maximum(np.max(np.abs(_values), axis=0), 1.0)
                _constant = np.ptp(_values, axis=0) <= (1e-10 * _scale)
                _constant_genes[_side].update(
                    str(_gene)
                    for _gene, _is_constant in zip(_work.var_names[_start : _start + 64], _constant)
                    if _is_constant
                )

        _scope_dir = MORAN_GRAPH_DIR / str(_scope).replace("/", "_")
        _scope_dir.mkdir(parents=True, exist_ok=True)
        _graph_id = f"moran_graph_{_scope}"
        _actual_params = {}
        _graph_status = "computed"
        _raw_uns_has_spatial = "spatial" in _raw_work.uns

        # Do not construct a six-neighbor graph for under-supported groups.
        if _n_units < MORAN_MIN_UNITS:
            _graph_status = "not_built_too_few_units"
            _weights_raw = sparse.csr_matrix((_n_units, _n_units), dtype=float)
        else:
            # coord_type is intentionally omitted: Raw uns["spatial"] selects grid.
            sq.gr.spatial_neighbors(_raw_work, n_neighs=MORAN_N_NEIGHS)
            _weights_raw = sparse.csr_matrix(
                _raw_work.obsp["spatial_connectivities"], dtype=float, copy=True
            )
            _weights_raw.sum_duplicates()
            _weights_raw.setdiag(0.0)
            _weights_raw.eliminate_zeros()
            _actual_params = dict(_raw_work.uns.get("spatial_neighbors", {}).get("params", {}))

        _n_edges_raw = int(_weights_raw.nnz)
        _raw_row_sums = np.asarray(_weights_raw.sum(axis=1)).ravel()
        _n_isolates = int(np.count_nonzero(_raw_row_sums == 0))
        _n_components = int(
            connected_components(_weights_raw, directed=False, return_labels=False)
            if _n_units else 0
        )
        _inverse_row_sums = np.zeros_like(_raw_row_sums, dtype=float)
        np.divide(1.0, _raw_row_sums, out=_inverse_row_sums, where=_raw_row_sums > 0)
        _weights_row = sparse.diags(_inverse_row_sums).dot(_weights_raw).tocsr()
        _weights_row.sum_duplicates()
        _weights_row.eliminate_zeros()
        _n_edges_row = int(_weights_row.nnz)

        sparse.save_npz(_scope_dir / "weights_raw.npz", _weights_raw)
        sparse.save_npz(_scope_dir / "weights_row_normalized.npz", _weights_row)
        _axis = pd.DataFrame(
            {
                "row_index": np.arange(_n_units, dtype=int),
                "unit_id": _raw_work.obs_names.astype(str).to_numpy(),
                "x": _raw_work.obsm["spatial"][:, 0].astype(float),
                "y": _raw_work.obsm["spatial"][:, 1].astype(float),
            }
        )
        _axis.to_csv(_scope_dir / "observations.csv.gz", index=False, compression="gzip")

        _comparison_id = json.dumps(
            {
                "analysis": "moran",
                "edge": "raw_vs_reconstruction",
                "route": str(CONFIG.get("route_kind", "unknown")),
                "sample_id": str(CONFIG.get("sample", {}).get("id", "unknown")),
                "scope": str(_scope),
            },
            sort_keys=True,
            separators=(",", ":"),
        )
        _result = compare_moran(
            _raw_work,
            _recon_work,
            _weights_row,
            min_units=MORAN_MIN_UNITS,
            min_edges=1,
        )
        for _side in ("raw", "reconstruction"):
            _side_mask = (
                _result[f"{_side}_status"].eq("computed")
                & _result["gene_id"].astype(str).isin(_constant_genes[_side])
            )
            _result.loc[_side_mask, f"{_side}_moran"] = np.nan
            _result.loc[_side_mask, f"{_side}_status"] = "not_computable"
            _result.loc[_side_mask, f"{_side}_reason"] = "constant_expression"
        _both_computed = _result["raw_status"].eq("computed") & _result["reconstruction_status"].eq("computed")
        _result["delta"] = np.where(
            _both_computed,
            _result["reconstruction_moran"] - _result["raw_moran"],
            np.nan,
        )
        _result["comparison_status"] = np.where(_both_computed, "computed", "unavailable")
        _result["scope"] = _scope
        _result["comparison_id"] = _comparison_id
        _result["comparison_edge"] = "raw_vs_reconstruction"
        _result["graph_id"] = _graph_id
        _result["raw_expression_view"] = "native_full_gene_space"
        _reconstruction_view = RUNS.get(_scope, {}).get("expression_view")
        if _reconstruction_view is None:
            _reconstruction_view = (
                "cluster_mean_projection"
                if "spatial_carrier" in RUNS.get(_scope, {})
                else "native_full_gene_space"
            )
        _result["reconstruction_expression_view"] = _reconstruction_view
        _result["n_units"] = _n_units
        _result["n_edges_raw"] = _n_edges_raw
        _result["n_edges_row_normalized"] = _n_edges_row
        _result["support_status"] = (
            "too_few_units" if _n_units < MORAN_MIN_UNITS
            else "no_edges" if _n_edges_raw == 0
            else "eligible"
        )
        MORAN_TABLES.append(_result)
        MORAN_GRAPH_AUDIT.append(
            {
                "scope": _scope,
                "comparison_id": _comparison_id,
                "graph_id": _graph_id,
                "graph_status": _graph_status,
                "n_units": _n_units,
                "n_edges_raw": _n_edges_raw,
                "n_edges_row_normalized": _n_edges_row,
                "n_isolates": _n_isolates,
                "n_components": _n_components,
                "support_status": (
                    "too_few_units" if _n_units < MORAN_MIN_UNITS
                    else "no_edges" if _n_edges_raw == 0
                    else "eligible"
                ),
                "raw_uns_has_spatial": _raw_uns_has_spatial,
                "coordinate_source": MORAN_COORDINATE_SOURCES[_scope],
                "coordinate_unit": _coordinate_unit,
                "coordinate_to_microns": _coordinate_to_microns,
                "coord_type_requested": "auto_from_raw_uns_spatial",
                "coord_type_actual": str(_actual_params.get("coord_type", "not_built")),
                "n_neighbors_requested": MORAN_N_NEIGHS,
                "n_neighbors_actual": _actual_params.get("n_neighs", _actual_params.get("n_neighbors", np.nan)),
                "actual_graph_parameters": json.dumps(_actual_params, default=str, sort_keys=True),
                "weight_transformation": "row_l1_same_as_sq_gr_spatial_autocorr_transformation_true",
                "min_units": MORAN_MIN_UNITS,
                "min_edges": 1,
                "scanpy_version": getattr(sc, "__version__", "unknown"),
                "squidpy_version": getattr(sq, "__version__", "unknown"),
            }
        )

    MORAN = pd.concat(MORAN_TABLES, ignore_index=True)
    MORAN_GRAPH_AUDIT = pd.DataFrame(MORAN_GRAPH_AUDIT)
    MORAN.to_csv(MORAN_DIR / "moran_by_cell_type.csv", index=False)
    MORAN_GRAPH_AUDIT.to_csv(MORAN_DIR / "moran_graph_audit.csv", index=False)


    return MORAN, MORAN_GRAPH_AUDIT


def summarize_local_results(RUNS, CONFIG):
    def _distribution(values):
        values = pd.Series(values).dropna().astype(float)
        return {
            "n_observations": int(values.size),
            "mean": float(values.mean()) if len(values) else np.nan,
            "median": float(values.median()) if len(values) else np.nan,
            "q1": float(values.quantile(.25)) if len(values) else np.nan,
            "q3": float(values.quantile(.75)) if len(values) else np.nan,
        }

    LOCAL_DIVERSITY_ROWS = []
    LOCAL_WIDE_ROWS = []
    REGION_EXTENT_ROWS = []
    for parent, run in RUNS.items():
        spatial = run["spatial"]
        metrics = spatial["metrics"].copy()
        valid = metrics.loc[metrics["valid_window"]].copy()
        scale = float(spatial["side_um"])
        threshold = dict(spatial["state_threshold"])
        common = {
            "scope": parent,
            "scale_um": scale,
            "total_windows": int(len(metrics)),
            "n_valid_windows": int(len(valid)),
            "valid_units": int(valid["n_units"].sum()),
            "min_parent_units": int(CONFIG["spatial_region"]["min_parent_units"]),
            "rarefaction_draws": int(CONFIG["spatial_region"]["rarefaction_draws"]),
        }
        wide = dict(common)
        for metric in ("k_obs", "neff", "evenness"):
            reconstruction = _distribution(valid[f"{metric}_recon"])
            for name, baseline_column, delta_column in (
                ("raw_leiden", f"{metric}_raw", f"delta_{metric}_vs_raw_leiden"),
                ("raw_level2", f"{metric}_level2", f"delta_{metric}_vs_raw_level2"),
            ):
                baseline = _distribution(valid[baseline_column])
                delta = _distribution(valid[delta_column])
                LOCAL_DIVERSITY_ROWS.append({
                    **common,
                    "metric": metric,
                    "baseline": name,
                    "baseline_value_column": baseline_column,
                    "delta_value_column": delta_column,
                    **{f"baseline_{key}": value for key, value in baseline.items()},
                    **{f"reconstruction_{key}": value for key, value in reconstruction.items()},
                    **{f"delta_{key}": value for key, value in delta.items()},
                })
                wide[f"local_{metric}_{name}_baseline_median"] = baseline["median"]
                wide[f"local_{metric}_{name}_delta_median"] = delta["median"]
                wide[f"local_{metric}_{name}_delta_n"] = delta["n_observations"]
            wide[f"local_{metric}_reconstruction_median"] = reconstruction["median"]
        extent = spatial["extent"].copy()
        extent["scope"] = parent
        extent["scale_um"] = scale
        for key, value in threshold.items():
            extent[f"state_threshold_{key}"] = value
        REGION_EXTENT_ROWS.append(extent)
        wide["state_threshold_status"] = threshold.get("status")
        wide["state_threshold"] = threshold.get("threshold")
        wide["state_n_windows"] = threshold.get("n_windows")
        wide["state_n_valid_bootstrap"] = threshold.get("n_valid_bootstrap")
        overall_extent = extent.loc[extent["level1_region"].eq("Overall")]
        if not overall_extent.empty:
            overall = overall_extent.iloc[0]
            for column in (
                "region_available", "valid_windows", "region_windows", "region_area_um2",
                "region_area_mm2", "area_fraction", "valid_units", "region_units",
                "unit_fraction", "threshold_status",
            ):
                if column in overall:
                    wide[f"state_region_overall_{column}"] = overall[column]
        LOCAL_WIDE_ROWS.append(wide)

    LOCAL_DIVERSITY_CROSS_PARENT = pd.DataFrame(LOCAL_DIVERSITY_ROWS)
    REGION_EXTENT_CROSS_PARENT = pd.concat(REGION_EXTENT_ROWS, ignore_index=True)
    LOCAL_WIDE = pd.DataFrame(LOCAL_WIDE_ROWS)

    return LOCAL_DIVERSITY_CROSS_PARENT, REGION_EXTENT_CROSS_PARENT, LOCAL_WIDE
