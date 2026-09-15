"""Local presentation helpers; scientific stages remain explicit in the notebooks."""
import json
from pathlib import Path
from types import SimpleNamespace
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm
from revise.analysis.reconstruction_impact import map_raw_level2_labels

ANATOMY_ORDER = ["Tumor", "Normal", "Interface", "Other"]
ANATOMY_COLORS = {"Tumor": "#d73027", "Normal": "#4575b4", "Interface": "#984ea3", "Other": "#d9d9d9"}

def extent_table(impact, column):
    source = impact.region_extent_by_anatomy if column == "in_state_region" else impact.gain_region_extent_by_anatomy
    return source.loc[source["level1_region"].isin(["Overall", "Tumor", "Normal", "Interface"])]

def coordinates(adata_obj):
    return pd.DataFrame(adata_obj.obsm["spatial"], index=adata_obj.obs_names, columns=["x", "y"])

def raw_level2_mapping(raw_parent, parent, *, config, root):
    mapping_config = config["raw_level2_mapping"]
    level1 = mapping_config["level1_column"]
    reference_context = ad.read_h5ad(root / mapping_config["reference_h5ad"], backed="r")
    reference_mask = reference_context.obs[level1].astype(str).str.replace("/", "_", regex=False).eq(parent)
    reference_filter = mapping_config.get("reference_filter")
    if reference_filter:
        reference_mask &= reference_context.obs[reference_filter["column"]].astype(str).eq(str(reference_filter["value"]))
    reference_parent = reference_context[reference_mask].to_memory()
    reference_context.file.close()
    pot = mapping_config.get("pot", {}); tacco = mapping_config.get("tacco", {})
    return map_raw_level2_labels(raw_parent, reference_parent, parent_value=parent, method=mapping_config["method"], level1_col=level1, level2_col=mapping_config["level2_column"], reference_filter_column=(reference_filter or {}).get("column"), reference_filter_value=(reference_filter or {}).get("value"), pot_reg=pot.get("reg", .1), pot_reg_m=pot.get("reg_m", 0.), pot_reg_type=pot.get("reg_type", "entropy"), tacco_multi_center=tacco.get("multi_center"), tacco_lamb=tacco.get("lamb"))


PAIRING_AUDIT_COLUMNS = [
    "scope", "coordinate_source", "coordinate_unit", "rtol", "atol",
    "n_checked", "id_unique", "id_subset", "coordinates_match",
    "coordinate_finite",
]


def _pairing_row(
    scope,
    candidate_ids,
    candidate_coordinates,
    full_coordinates,
    coordinate_unit,
    coordinate_source,
):
    """Validate one existing Raw/carrier coordinate check and record its range."""

    candidate_ids = pd.Index(candidate_ids)
    id_unique = bool(candidate_ids.is_unique)
    id_subset = bool(candidate_ids.isin(full_coordinates.index).all())
    if not id_unique or not id_subset:
        raise ValueError("Spatial IDs must be a unique subset of Raw")
    expected = full_coordinates.loc[candidate_ids].to_numpy(float)
    observed = np.asarray(candidate_coordinates, dtype=float)
    coordinate_finite = bool(np.isfinite(expected).all() and np.isfinite(observed).all())
    coordinates_match = bool(
        coordinate_finite
        and observed.shape == expected.shape
        and np.allclose(observed, expected, rtol=0, atol=1e-8)
    )
    # Keep the original assertion as the authoritative check and preserve its
    # existing tolerance.  The audit row is written only after it succeeds.
    np.testing.assert_allclose(observed, expected, rtol=0, atol=1e-8)
    return {
        "scope": scope,
        "coordinate_source": coordinate_source,
        "coordinate_unit": coordinate_unit,
        "rtol": 0,
        "atol": 1e-8,
        "n_checked": int(len(candidate_ids)),
        "checked_basis": "full_spatial_carrier_before_cohort_selection",
        "id_unique": id_unique,
        "id_subset": id_subset,
        "coordinates_match": coordinates_match,
        "coordinate_finite": coordinate_finite,
    }


def prepare_route_inputs(
    config,
    root,
    raw_context,
    full_coordinates,
    full_level1,
    figures,
    source_files,
    input_rows,
    *,
    seed,
    sample_limit,
    use_full_visiumhd_cohort,
    parent_source,
    sample_fn,
):
    """Prepare Raw-defined route cohorts and save the existing pairing audits."""

    source_files = dict(source_files)
    runs = {}
    global_inputs = None
    pairing_rows = []
    coordinate_unit = config.get("spatial_region", {}).get(
        "coordinate_unit", config.get("context", {}).get("coordinate_unit", "unknown")
    )

    if config["route_kind"] == "sp_svc":
        comparison_config = config["partition_change"]["comparisons"][0]
        source_files["reconstruction"] = root / comparison_config["reconstructed_spatial_h5ad"]
        recon_context = ad.read_h5ad(source_files["reconstruction"], backed="r")
        reconstructed_ids = recon_context.obs_names
        pair_source = (
            "reconstruction.obsm[spatial] vs Raw.obsm[spatial] "
            "(full reconstruction carrier before cohort selection)"
        )
        pairing = _pairing_row(
            "All", reconstructed_ids, coordinates(recon_context), full_coordinates,
            coordinate_unit, pair_source,
        )
        pairing_rows.append(pairing)
        # The same full-carrier check is the coordinate authority for each
        # parent; cohort sizes remain in observations/carriers tables.
        pairing_rows.extend(
            dict(pairing, scope=parent)
            for parent in parent_source
        )
        global_ids = sample_fn(reconstructed_ids, sample_limit)
        global_inputs = (
            raw_context[global_ids].to_memory(),
            recon_context[global_ids].to_memory(),
        )
        input_rows.append({
            "role": "Reconstruction full carrier", "n_units": recon_context.n_obs,
            "n_genes": recon_context.n_vars, "spatial_axis": True,
        })
        for parent, source_label in parent_source.items():
            eligible = reconstructed_ids[full_level1.loc[reconstructed_ids].eq(source_label)]
            ids = sample_fn(eligible, sample_limit)
            runs[parent] = {
                "raw_expression": raw_context[ids].to_memory(),
                "recon_expression": recon_context[ids].to_memory(),
                "available": len(eligible), "sampled": len(ids),
                "expression_view": "native_reconstruction",
            }
            membership = full_coordinates.assign(
                unit_id=full_coordinates.index,
                raw_level1=full_level1,
                reconstruction_covered=full_coordinates.index.isin(reconstructed_ids),
            )
            membership["included"] = membership.index.isin(ids)
            membership["reason"] = np.select(
                [
                    ~membership.reconstruction_covered,
                    ~membership.raw_level1.eq(source_label),
                    ~membership.included,
                ],
                ["not_reconstructed", "outside_raw_level1", "not_sampled"],
                default="included",
            )
            figures.save_table(
                membership.reset_index(drop=True).assign(scope=parent),
                "inputs/observations.csv.gz", scope=parent,
            )
        global_membership = full_coordinates.assign(
            unit_id=full_coordinates.index, raw_level1=full_level1,
        )
        global_membership["included"] = global_membership.index.isin(global_ids)
        global_membership["reason"] = np.select(
            [
                ~global_membership.index.isin(reconstructed_ids),
                ~global_membership.included,
            ],
            ["not_reconstructed", "not_sampled"],
            default="included",
        )
        figures.save_table(
            global_membership.reset_index(drop=True).assign(scope="All"),
            "inputs/observations.csv.gz", scope="All",
        )
        recon_context.file.close()
    else:
        for spec in config["partition_change"]["comparisons"]:
            parent = spec["name"]
            source_files[parent + "_spatial"] = root / spec["reconstructed_spatial_h5ad"]
            source_files[parent + "_expression"] = root / spec["expression_h5ad"]
            spatial = ad.read_h5ad(source_files[parent + "_spatial"])
            pair_source = (
                "spatial_carrier.obsm[spatial] vs Raw.obsm[spatial] "
                "(full spatial carrier before parent filter)"
            )
            pairing_rows.append(_pairing_row(
                parent, spatial.obs_names, coordinates(spatial), full_coordinates,
                coordinate_unit, pair_source,
            ))
            eligible = spatial.obs_names[full_level1.loc[spatial.obs_names].eq(parent_source[parent])]
            raw_parent = raw_context[eligible].to_memory()
            carrier = spatial[eligible].copy()
            runs[parent] = {
                "raw_expression": raw_parent, "spatial_carrier": carrier,
                "expression_h5ad": spec["expression_h5ad"],
                "available": len(eligible), "sampled": len(eligible),
                "expression_view": "cluster_mean_projection",
            }
            membership = full_coordinates.assign(
                unit_id=full_coordinates.index, raw_level1=full_level1,
            )
            membership["included"] = membership.index.isin(eligible)
            membership["reason"] = np.select(
                [
                    ~membership.index.isin(spatial.obs_names),
                    ~membership.raw_level1.eq(parent_source[parent]),
                ],
                ["not_reconstructed", "raw_level1_mismatch"],
                default="included",
            )
            figures.save_table(
                membership.reset_index(drop=True).assign(scope=parent),
                "inputs/observations.csv.gz", scope=parent,
            )
            input_rows.append({
                "role": parent + " spatial carrier", "n_units": spatial.n_obs,
                "n_genes": spatial.n_vars,
                "excluded_raw_level1_mismatch": spatial.n_obs - len(eligible),
                "spatial_axis": True,
            })

    pairing_audit = pd.DataFrame(pairing_rows, columns=PAIRING_AUDIT_COLUMNS)
    figures.save_table(pairing_audit, "inputs/pairing_audit.csv")
    return {
        "source_files": source_files,
        "runs": runs,
        "global_inputs": global_inputs,
        "input_rows": input_rows,
        "pairing_audit": pairing_audit,
    }


def prepare_expression_views(runs, config, root, figures, source_files, input_rows):
    """Build complete-gene expression views and save their source audits."""

    from revise.application.ist_assembly import cluster_means

    gene_rows = []
    projection_audits = []
    for parent, run in runs.items():
        raw_expr = run["raw_expression"]
        if config["route_kind"] == "sc_svc":
            expression = ad.read_h5ad(source_files[parent + "_expression"])
            means, clusters = cluster_means(expression)
            spatial = run["spatial_carrier"]
            cluster_keys = spatial.obs["SVC_cluster"].tolist()
            if not set(cluster_keys).issubset(set(clusters)):
                raise ValueError("Every retained spatial cluster needs an expression mean")
            lookup = {key: i for i, key in enumerate(clusters)}
            positions = [lookup[key] for key in cluster_keys]
            projected = ad.AnnData(
                means[positions].copy(), obs=spatial.obs.copy(), var=expression.var.copy(),
                obsm={"spatial": spatial.obsm["spatial"].copy()},
            )
            projected.uns["expression_view"] = "cluster_mean_projection"
            run["recon_expression"] = projected
            mapping = pd.DataFrame({"unit_id": spatial.obs_names, "cluster": cluster_keys})
            mapping["reference_cells"] = mapping.cluster.map(expression.obs["SVC_cluster"].value_counts())
            figures.save_table(mapping.assign(scope=parent), "inputs/mapping.csv.gz", scope=parent)
            figures.save_table(
                pd.DataFrame({
                    "reference_cell_id": expression.obs_names,
                    "cluster": expression.obs["SVC_cluster"].to_numpy(),
                }),
                "inputs/reference_clusters.csv.gz", scope=parent,
            )
            input_rows.append({
                "role": parent + " expression carrier", "n_units": expression.n_obs,
                "n_genes": expression.n_vars, "spatial_axis": False,
            })
            projection_audits.append(
                mapping.groupby("cluster", sort=False).agg(
                    spatial_units=("unit_id", "size"),
                    reference_cells=("reference_cells", "first"),
                ).reset_index()
            )
        recon_expr = run["recon_expression"]
        if not raw_expr.var_names.is_unique or not recon_expr.var_names.is_unique:
            raise ValueError("Expression gene identities must be unique")
        genes = raw_expr.var_names.union(recon_expr.var_names, sort=False)
        for side, view in (("raw", raw_expr), ("reconstruction", recon_expr)):
            available = genes.isin(view.var_names)
            gene_rows.append(pd.DataFrame({
                "scope": parent, "gene_id": genes, "side": side,
                "provided": available,
                "status": np.where(available, "provided", "unmeasured"),
                "expression_view": "measured" if side == "raw" else run["expression_view"],
            }))
    gene_availability = pd.concat(gene_rows, ignore_index=True)
    figures.save_table(gene_availability, "inputs/gene_availability.csv")
    figures.save_table(pd.DataFrame(input_rows), "inputs/carriers.csv")
    return gene_availability, input_rows, projection_audits


def save_partition_evidence(
    partitions, global_inputs, runs, sample_id, seed, figures,
    *, complexity_kind="same_resolution",
):
    """Persist partition comparison tables without changing their calculations."""

    comparisons = []
    partition_tables = {
        name: [] for name in [
            "summary", "assignments", "mapping", "contingency",
            "resolution_sweep", "complexity_sweep",
        ]
    }
    for scope, part in partitions.items():
        raw_expr = global_inputs[0] if scope == "All" else runs[scope]["raw_expression"]
        comparison_groups = [
            ("matched_k", part.comparisons),
            (complexity_kind, part.complexity_comparisons),
        ]
        for kind_name, comparisons_by_edge in comparison_groups:
            for edge, comparison in comparisons_by_edge.items():
                comparison_id = json.dumps(
                    [sample_id, scope, kind_name, str(edge)], separators=(",", ":")
                )
                comparisons.append({
                    "comparison_id": comparison_id, "sample_id": sample_id,
                    "scope": scope, "analysis": "partition",
                    "comparison_edge": str(edge), "comparison_kind": kind_name,
                    "observation_basis": "Raw-defined partition cohort",
                    "gene_rule": part.audit.get("feature_selection"), "seed": seed,
                })
                context = {
                    "comparison_id": comparison_id, "scope": scope,
                    "comparison_kind": kind_name,
                }
                for table_name in ["summary", "assignments", "mapping"]:
                    partition_tables[table_name].append(
                        getattr(comparison, table_name).assign(**context)
                    )
                counts = comparison.contingency.rename_axis(
                    index="raw_cluster", columns="recon_cluster"
                )
                for normalization, matrix in [
                    ("count", counts),
                    ("raw_row_fraction", counts.div(counts.sum(axis=1), axis=0)),
                ]:
                    long = matrix.stack().rename("value").reset_index().assign(
                        normalization=normalization, **context
                    )
                    partition_tables["contingency"].append(long)
        partition_tables["resolution_sweep"].append(part.sweep.assign(scope=scope))
        partition_tables["complexity_sweep"].append(part.complexity_sweep.assign(scope=scope))
        retained = next(iter(part.comparisons.values())).assignments.index
        figures.save_table(
            pd.DataFrame({
                "unit_id": raw_expr.obs_names,
                "partition_included": raw_expr.obs_names.isin(retained),
                "reason": np.where(
                    raw_expr.obs_names.isin(retained), "retained", "Raw_QC"
                ),
            }).assign(scope=scope),
            "inputs/partition_observations.csv.gz", scope=scope,
        )
        figures.save_table(
            pd.DataFrame({
                "gene_id": raw_expr.var_names,
                "partition_feature": raw_expr.var_names.isin(part.feature_names),
            }).assign(scope=scope),
            "inputs/partition_features.csv", scope=scope,
        )
        figures.save_json(
            {
                "audit": part.audit,
                "representation_audit": part.representation_audit,
                "matched_cluster_status": part.matched_cluster_status,
            },
            f"reconstruction_impact/partition/{scope}_audit.json",
        )
    for table_name, tables in partition_tables.items():
        suffix = ".csv.gz" if table_name == "assignments" else ".csv"
        figures.save_table(
            pd.concat(tables, ignore_index=True),
            "reconstruction_impact/partition/" + table_name + suffix,
        )
    figures.save_table(pd.DataFrame(comparisons), "comparisons.csv")
    return comparisons


def save_anatomy_evidence(anatomy, full_level1, config, analysis_root, figures):
    """Persist anatomy context and return the unit-level anatomy assignment."""

    anatomy_dir = Path(analysis_root) / "anatomy"
    anatomy_dir.mkdir(parents=True, exist_ok=True)
    anatomy.anatomy_windows.to_csv(anatomy_dir / "anatomy_windows.csv", index=False)
    anatomy.anatomy_context_summary.to_csv(
        anatomy_dir / "anatomy_context_summary.csv", index=False
    )
    anatomy.support_sensitivity.to_csv(
        anatomy_dir / "support_sensitivity.csv", index=False
    )
    figures.save_json(
        {
            "scale": anatomy.scale_audit,
            "selection": anatomy.support_selection,
            "source_coordinate_unit": config["spatial_region"].get("coordinate_unit"),
            "saved_coordinate_unit": "micron",
        },
        "anatomy/window_decision.json",
    )
    anatomy_units = anatomy.anatomy_unit_assignments.copy()
    anatomy_units["raw_level1"] = full_level1.reindex(anatomy_units.index)
    anatomy_units["level1_region"] = anatomy_units.window_id.map(
        anatomy.anatomy_windows.set_index("window_id").level1_region
    )
    anatomy_units.reset_index(names="unit_id").to_csv(
        anatomy_dir / "unit_window_assignments.csv.gz", index=False
    )
    return anatomy_units


def save_spatial_feature_values(runs, aucell, emt_name, config, analysis_root):
    """Assemble and persist the EMT spatial field from already-scored rows."""

    feature_rows = []
    for parent, run in runs.items():
        ids = run["raw_expression"].obs_names
        xy = coordinates(run["raw_expression"][ids]) * config["spatial_region"]["microns_per_coordinate"]
        score_rows = aucell.loc[aucell.scope.eq(parent)].set_index("unit_id").reindex(ids)
        if score_rows.index.has_duplicates:
            raise ValueError(f"{parent}: AUCell spatial score axis is not unique")
        feature_rows.append(pd.DataFrame({
            "scope": parent,
            "unit_id": ids,
            "feature": emt_name,
            "x": xy.x.to_numpy(),
            "y": xy.y.to_numpy(),
            "raw_value": score_rows["raw_score"].to_numpy(float),
            "reconstruction_value": score_rows["reconstruction_score"].to_numpy(float),
            "delta": score_rows["delta"].to_numpy(float),
            "raw_status": score_rows["raw_status"].to_numpy(),
            "reconstruction_status": score_rows["reconstruction_status"].to_numpy(),
        }))
    features = pd.concat(feature_rows, ignore_index=True)
    features["expression_layer"] = "AUCell_original_X"
    features["coordinate_unit"] = "micron"
    field_dir = Path(analysis_root) / "spatial_fields"
    field_dir.mkdir(parents=True, exist_ok=True)
    features.to_csv(
        field_dir / "values.csv.gz", index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    return features


def prepare_local_state(
    parent, run, config, root, anatomy, analysis_root, figures, cell_equivalent_um
):
    """Map the two Raw baselines and create the primary paired windows."""

    from revise.analysis.basic.spatial_region import (
        assign_square_windows,
        convert_coordinates_to_microns,
        select_window_scale,
    )

    assignments = run["comparison"].assignments.set_index("unit_id")
    raw = run["raw_expression"][assignments.index].copy()
    level2 = raw_level2_mapping(raw, parent, config=config, root=root)
    paired_um = convert_coordinates_to_microns(
        coordinates(raw),
        microns_per_coordinate=config["spatial_region"]["microns_per_coordinate"],
    )
    origin = (
        float(anatomy.scale_audit["origin_x_um"]),
        float(anatomy.scale_audit["origin_y_um"]),
    )
    selection, support = select_window_scale(
        paired_um,
        candidate_window_sides=config["spatial_region"]["candidate_window_sides_um"],
        min_parent_units=4,
        origin=origin,
    )
    side_um = float(
        selection["window_side_length"]
        or min(config["spatial_region"]["candidate_window_sides_um"])
    )
    windows = assign_square_windows(
        paired_um, window_side_length=side_um, origin=origin
    )
    run["spatial"] = {
        "raw": raw,
        "assignments": assignments,
        "level2": level2,
        "paired_um": paired_um,
        "origin": origin,
        "windows": windows,
        "selection": selection,
        "support": support,
        "side_um": side_um,
    }
    local_dir = Path(analysis_root) / "local_state" / parent
    local_dir.mkdir(parents=True, exist_ok=True)
    level2.labels.rename("raw_level2").rename_axis("unit_id").to_csv(
        local_dir / "raw_level2_labels.csv"
    )
    level2.posterior.rename_axis("unit_id").to_csv(
        local_dir / "raw_level2_posterior.csv.gz",
        compression={"method": "gzip", "mtime": 0},
    )
    support.to_csv(local_dir / "support_sensitivity.csv", index=False)
    run["raw_level2"] = level2
    unit_assignments = windows.copy()
    unit_assignments["raw_cluster"] = assignments.raw_cluster
    unit_assignments["reconstructed_cluster"] = assignments.recon_cluster
    run["impact"] = SimpleNamespace(
        unit_assignments=unit_assignments,
        window_metrics=pd.DataFrame(),
        support_sensitivity=support,
        scale_audit={
            "main_window_side_um": side_um,
            "main_window_cells_per_side": side_um / cell_equivalent_um,
            "min_parent_units": 4,
            "rarefaction_draws": 200,
        },
    )
    figures.save_json(
        {
            "selection": selection,
            "scale": run["impact"].scale_audit,
            "raw_level2": level2.audit,
        },
        f"local_state/{parent}/window_decision.json",
    )
    return run


def save_local_state_evidence(parent, run, anatomy, analysis_root):
    """Join local metrics to anatomy and persist the existing summaries."""

    from revise.analysis.reconstruction_impact import (
        _map_to_anatomy_regions,
        _summarize_parent_window_anatomy,
    )
    from revise.analysis.basic.spatial_region import (
        summarize_cluster_change_by_anatomy,
        summarize_diversity_by_anatomy,
    )

    spatial = run["spatial"]
    unit_assignments = spatial["windows"].loc[:, [
        "x", "y", "window_id", "window_x_index", "window_y_index",
    ]].copy()
    assignments = spatial["assignments"]
    unit_assignments["raw_cluster"] = assignments.raw_cluster.reindex(
        unit_assignments.index
    ).astype(str)
    unit_assignments["reconstructed_cluster"] = assignments.recon_cluster.reindex(
        unit_assignments.index
    ).astype(str)
    unit_assignments["unit_changed"] = assignments.unit_changed.reindex(
        unit_assignments.index
    ).astype(bool)
    unit_assignments["raw_level2"] = spatial["level2"].labels.reindex(
        unit_assignments.index
    )
    context = _map_to_anatomy_regions(
        unit_assignments.loc[:, ["x", "y"]], anatomy
    )
    unit_assignments["anatomy_window_id"] = context.anatomy_window_id.to_numpy()
    unit_assignments["level1_region"] = context.level1_region.to_numpy()
    change_window = unit_assignments.groupby("window_id").unit_changed.agg(
        n_changed_units="sum", unit_change_fraction="mean"
    ).reset_index()
    metrics = spatial["metrics"].merge(
        change_window, on="window_id", validate="one_to_one"
    ).merge(
        _summarize_parent_window_anatomy(unit_assignments),
        on="window_id", how="left", validate="one_to_one",
    )
    spatial.update(unit_assignments=unit_assignments, metrics=metrics)
    local_dir = Path(analysis_root) / "local_state" / parent
    unit_assignments.reset_index(names="unit_id").to_csv(
        local_dir / "unit_window_assignments.csv.gz", index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    metrics.to_csv(local_dir / "window_metrics.csv", index=False)
    anatomy_summary = summarize_diversity_by_anatomy(metrics)
    anatomy_summary = anatomy_summary.rename(columns={
        f"{stat}_delta_{metric}": f"{stat}_delta_{metric}_vs_raw_leiden"
        for stat in ("median", "q1", "q3")
        for metric in ("k_obs", "entropy", "neff", "evenness")
    })
    anatomy_summary = anatomy_summary.assign(
        scope=parent, scale_um=spatial["side_um"]
    )
    anatomy_summary.to_csv(local_dir / "diversity_by_anatomy.csv", index=False)
    change_summary = summarize_cluster_change_by_anatomy(
        unit_assignments, metrics
    ).assign(
        scope=parent,
        scale_um=spatial["side_um"],
        matched_cluster_status=run["partition"].matched_cluster_status,
        headline_eligible=run["partition"].matched_cluster_status
        != "unmatched_cluster_complexity",
    )
    change_summary.to_csv(local_dir / "change_by_anatomy.csv", index=False)
    return change_summary, anatomy_summary


def save_final_evidence(
    runs,
    config,
    output_dir,
    comparisons,
    change_table,
    moran_summary,
    aucell_summary,
    sample_id,
    emt_name,
    gmt_path,
    source_files,
    figures,
):
    """Persist cross-parent, provenance, and artifact tables after all stages."""

    from revise.analysis.reconstruction_impact import file_sha256
    from reproduce.case.reconstruction_impact.notebook_analysis import summarize_local_results

    local_diversity, region_extent, local_wide = summarize_local_results(runs, config)
    figures.save_table(local_diversity, "cross_parent_local_diversity.csv")
    figures.save_table(region_extent, "cross_parent_state_region_extent.csv")
    cross_parent = (
        change_table
        .merge(moran_summary, on="scope", how="outer")
        .merge(aucell_summary, on="scope", how="outer", suffixes=("_moran", "_aucell"))
        .merge(local_wide, on="scope", how="outer")
    )
    figures.save_table(cross_parent, "cross_parent_summary.csv")
    seed = config["partition_change"]["random_state"]
    for row in moran_summary.to_dict("records"):
        comparisons.append({
            "comparison_id": row["comparison_id"], "sample_id": sample_id,
            "scope": row["scope"], "analysis": "moran",
            "comparison_edge": "raw_vs_reconstruction",
            "gene_rule": "full_space_per_side", "seed": seed,
        })
    for scope, run in runs.items():
        comparisons.append({
            "comparison_id": json.dumps([sample_id, scope, "aucell", emt_name]),
            "sample_id": sample_id, "scope": scope, "analysis": "aucell",
            "comparison_edge": "raw_vs_reconstruction",
            "gene_rule": "full_space_per_side", "seed": seed,
        })
        comparisons.append({
            "comparison_id": json.dumps([sample_id, scope, "local_diversity"]),
            "sample_id": sample_id, "scope": scope, "analysis": "local_diversity",
            "comparison_edge": "reconstruction_vs_raw_leiden_and_raw_level2",
            "scale_um": float(run["spatial"]["side_um"]),
            "min_parent_units": int(config["spatial_region"]["min_parent_units"]),
            "rarefaction_draws": int(config["spatial_region"]["rarefaction_draws"]),
            "seed": seed,
        })
        comparisons.append({
            "comparison_id": json.dumps([sample_id, scope, "state_region"]),
            "sample_id": sample_id, "scope": scope, "analysis": "state_region",
            "comparison_edge": "reconstructed_neff_state",
            "scale_um": float(run["spatial"]["side_um"]),
            "threshold_status": run["spatial"]["state_threshold"].get("status"),
            "seed": seed,
        })
    figures.save_table(pd.DataFrame(comparisons), "comparisons.csv")

    source_files["Hallmark_EMT_resource"] = gmt_path
    repo_root = Path(__file__).resolve().parents[3]
    for helper_name in (
        "build_notebooks.py", "notebook_helpers.py", "notebook_analysis.py",
        "content_contract.py",
    ):
        source_files[helper_name] = (
            repo_root / "reproduce" / "case" / "reconstruction_impact" / helper_name
        )
    source_files["notebook-block.md"] = (
        repo_root / "docs" / "design" / "reconstruction-impact" / "templates" /
        "notebook-block.md"
    )
    source_audit = [
        {
            "role": role, "path": str(path), "size_bytes": path.stat().st_size,
            "sha256": file_sha256(path),
        }
        for role, path in source_files.items()
    ]
    figures.save_table(pd.DataFrame(source_audit), "inputs/source_files.csv")
    figures.save_json(
        {
            "sample_id": sample_id, "route": config["route_kind"], "seed": seed,
            "cohort": "independent Raw-defined scope sampling before analysis-specific QC",
            "expression": "complete gene space per side",
            "region": "high reconstructed diversity, not maximum change",
            "inputs": source_audit,
        },
        "audit.json",
    )
    artifact_rows = []
    output_dir = Path(output_dir)
    for folder in (output_dir / "analysis", output_dir / "figures"):
        for path in sorted(folder.rglob("*")):
            if path.is_file() and path.name != "artifacts.csv":
                artifact_rows.append({
                    "path": str(path.relative_to(output_dir)),
                    "size_bytes": path.stat().st_size,
                    "sha256": file_sha256(path),
                })
    artifacts = pd.DataFrame(artifact_rows)
    figures.save_table(artifacts, "artifacts.csv")
    return {
        "cross_parent": cross_parent,
        "local_diversity": local_diversity,
        "region_extent": region_extent,
        "local_wide": local_wide,
        "artifacts": artifacts,
    }

def compact_iqr(values):
    values = pd.Series(values).dropna()
    if values.empty:
        return "NA"
    return f"{values.median():.3g} [{values.quantile(.25):.3g}, {values.quantile(.75):.3g}]"


def prepare_aucell_overall(scores, summary, scope_labels, scope_order):
    """Prepare the published AUCell rows for the overall display cell."""

    scores = scores.copy()
    scores["scope"] = scores["scope"].astype(str)
    scores["scope_label"] = scores["scope"].map(scope_labels).fillna(scores["scope"])
    long = scores.melt(
        id_vars=["scope", "scope_label", "unit_id"],
        value_vars=["raw_score", "reconstruction_score"],
        var_name="side",
        value_name="AUCell score",
    ).dropna(subset=["AUCell score"])
    long["side"] = long["side"].map({"raw_score": "Raw", "reconstruction_score": "Reconstructed"})
    delta = scores.dropna(subset=["delta"]).copy()
    summary = summary.copy()
    if not summary.empty:
        summary["scope"] = summary["scope"].astype(str)
        summary["scope_order"] = summary["scope"].map({value: index for index, value in enumerate(scope_order)})
        summary = summary.sort_values("scope_order", kind="stable")
    computed = summary.loc[summary["comparison_status"].eq("computed")] if not summary.empty else summary
    return long, delta, summary, computed


def prepare_aucell_parent(scores, summary, scope_key, scope_label):
    """Prepare one parent AUCell display without recomputing its scores."""

    parent = scores.loc[scores["scope"].astype(str).eq(scope_key)].copy()
    long = parent.melt(
        id_vars=["scope", "unit_id"],
        value_vars=["raw_score", "reconstruction_score"],
        var_name="side",
        value_name="AUCell score",
    ).dropna(subset=["AUCell score"])
    long["side"] = long["side"].map({"raw_score": "Raw", "reconstruction_score": "Reconstructed"})
    delta = parent.dropna(subset=["delta"])
    parent_summary = summary.loc[summary["scope"].astype(str).eq(scope_key)].copy()
    return parent, long, delta, parent_summary, scope_label


class NotebookFigures:
    """Keep the shared figure axes and output destination explicit."""

    def __init__(self, output_dir, xlim, ylim, origin, cell_equivalent_um=8.0):
        self.output_dir = Path(output_dir)
        self.xlim = xlim
        self.ylim = ylim
        self.origin = origin
        self.cell_equivalent_um = cell_equivalent_um

    def save_figure(self, fig, name):
        figure_dir = self.output_dir / "figures"
        figure_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(figure_dir / f"{name}.png", dpi=180, bbox_inches="tight")

    def save_table(self, frame, name, *, scope=None):
        destination = self.output_dir / "analysis"
        if scope is not None:
            destination = destination / scope
        destination = destination / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(destination, index=False)

    def set_spatial_axes(self, ax):
        ax.set(xlim=self.xlim, ylim=self.ylim, xlabel="x (µm)", ylabel="y (µm)")
        ax.set_aspect("equal", adjustable="box")

    def window_field(self, ax, frame, value, *, side_um, cmap="viridis", norm=None, title="", categorical=False):
        valid = frame.dropna(subset=[value])
        if valid.empty:
            self.set_spatial_axes(ax); ax.set_title(title + " (no valid windows)"); return None
        x_indices = np.arange(valid["window_x_index"].min(), valid["window_x_index"].max() + 1)
        y_indices = np.arange(valid["window_y_index"].min(), valid["window_y_index"].max() + 1)
        grid = np.full((len(y_indices), len(x_indices)), np.nan)
        for row in valid.itertuples():
            grid[int(row.window_y_index - y_indices[0]), int(row.window_x_index - x_indices[0])] = getattr(row, value)
        x_edges = self.origin[0] + np.arange(x_indices[0], x_indices[-1] + 2) * side_um
        y_edges = self.origin[1] + np.arange(y_indices[0], y_indices[-1] + 2) * side_um
        image = ax.pcolormesh(x_edges, y_edges, grid, shading="flat", cmap=cmap, norm=norm)
        self.set_spatial_axes(ax); ax.set_title(title)
        return image

    def anatomy_map(self, ax, anatomy, focus=None):
        regions = anatomy.anatomy_windows.copy()
        codes = regions["level1_region"].map({name: i for i, name in enumerate(ANATOMY_ORDER)})
        if focus is not None:
            codes = np.where(regions["level1_region"].eq(focus), ANATOMY_ORDER.index(focus), ANATOMY_ORDER.index("Other"))
        frame = regions.assign(anatomy_code=codes)
        cmap = ListedColormap([ANATOMY_COLORS[name] for name in ANATOMY_ORDER])
        image = self.window_field(ax, frame, "anatomy_code", side_um=float(anatomy.scale_audit["main_window_side_um"]), cmap=cmap, norm=BoundaryNorm(np.arange(-.5, 4.5), 4), title=focus or "Tumor / Normal / Interface / Other")
        if image is not None and focus is None:
            colorbar = plt.colorbar(image, ax=ax, ticks=range(4), fraction=.046, pad=.04)
            colorbar.ax.set_yticklabels(ANATOMY_ORDER)

    def plot_support_curve(self, impact, title):
        table = impact.support_sensitivity
        fig, axes = plt.subplots(1, 2, figsize=(10, 3.4))
        cells = table["window_side_length"] / self.cell_equivalent_um
        axes[0].plot(cells, table["retained_parent_unit_fraction"], marker="o")
        axes[0].axvline(impact.scale_audit["main_window_cells_per_side"], color="black", ls="--")
        axes[0].set(xlabel="cells per side", ylabel="retained-unit fraction", title="retention")
        axes[1].plot(cells, table["valid_window_fraction"], marker="o")
        axes[1].axvline(impact.scale_audit["main_window_cells_per_side"], color="black", ls="--")
        axes[1].set(xlabel="cells per side", ylabel="valid-window fraction", title="valid windows")
        fig.suptitle(title); fig.tight_layout(); return fig

    def plot_contingency(self, comparisons, title):
        fig, axes = plt.subplots(1, len(comparisons), figsize=(4.5 * len(comparisons), 3.8), squeeze=False)
        for ax, (scope, comparison) in zip(axes[0], comparisons.items()):
            table = comparison.contingency.div(comparison.contingency.sum(axis=1), axis=0).fillna(0)
            sns.heatmap(table, ax=ax, cmap="mako", vmin=0, vmax=1, cbar=True)
            ax.set(title=scope, xlabel="Reconstructed / Final", ylabel="Raw")
        fig.suptitle(title); fig.tight_layout(); return fig

    def plot_changed_units(self, runs, title):
        fig, axes = plt.subplots(1, len(runs), figsize=(4.5 * len(runs), 4), squeeze=False)
        for ax, (scope, run) in zip(axes[0], runs.items()):
            frame = run["partition_coordinates"].join(run["comparison"].assignments[["unit_changed"]])
            colors = np.where(frame["unit_changed"], "#d73027", "#d9d9d9")
            ax.scatter(frame["x"], frame["y"], c=colors, s=1, linewidths=0, rasterized=True)
            self.set_spatial_axes(ax); ax.set_title(scope)
        fig.suptitle(title); fig.tight_layout(); return fig

    def plot_cluster_pair(self, run, title):
        if run["partition"].matched_cluster_status == "unmatched_cluster_complexity":
            title += " [audit: unmatched K]"
        frame = run["impact"].unit_assignments
        mapping = run["comparison"].mapping.set_index("recon_cluster")["raw_cluster"].to_dict()
        raw_values = sorted(frame["raw_cluster"].unique())
        palette = {value: color for value, color in zip(raw_values, sns.color_palette("tab20", len(raw_values)))}
        recon_colors = [palette.get(mapping.get(value), "#111111") for value in frame["reconstructed_cluster"]]
        fig, axes = plt.subplots(1, 2, figsize=(9, 4))
        axes[0].scatter(frame["x"], frame["y"], c=frame["raw_cluster"].map(palette), s=1, linewidths=0, rasterized=True)
        axes[1].scatter(frame["x"], frame["y"], c=recon_colors, s=1, linewidths=0, rasterized=True)
        for ax, label in zip(axes, ["Raw Leiden", "Reconstructed / Final"]):
            self.set_spatial_axes(ax); ax.set_title(label)
        fig.suptitle(title + " — matched clusters share colors"); fig.tight_layout(); return fig

    def plot_reconstructed_clusters(self, run, title):
        frame = run["impact"].unit_assignments
        values = sorted(frame["reconstructed_cluster"].unique())
        palette = {value: color for value, color in zip(values, sns.color_palette("tab20", len(values)))}
        fig, ax = plt.subplots(figsize=(5.2, 4.4))
        ax.scatter(frame["x"], frame["y"], c=frame["reconstructed_cluster"].map(palette), s=1, linewidths=0, rasterized=True)
        self.set_spatial_axes(ax); ax.set_title(title); fig.tight_layout(); return fig

    def plot_metric_comparison(self, run, metric, title):
        metrics = run["impact"].window_metrics.loc[lambda x: x["valid_window"]].copy()
        side = float(run["impact"].scale_audit["main_window_side_um"])
        state_columns = [f"{metric}_raw", f"{metric}_recon", f"{metric}_level2"]
        state_max = max(1.0, float(metrics[state_columns].max().max()))
        state_min = 0.0 if metric == "evenness" else 1.0
        state_max = 1.0 if metric == "evenness" else state_max
        delta_columns = [f"delta_{metric}_vs_raw_leiden", f"delta_{metric}_vs_raw_level2"]
        delta_max = max(0.01, float(metrics[delta_columns].abs().max().max()))
        state_norm = plt.Normalize(state_min, state_max)
        delta_norm = TwoSlopeNorm(vcenter=0, vmin=-delta_max, vmax=delta_max)
        panels = [
            (f"{metric}_raw", "Raw Leiden", "viridis", state_norm),
            (f"{metric}_recon", "Reconstructed", "viridis", state_norm),
            (f"delta_{metric}_vs_raw_leiden", "Recon − Raw Leiden", "coolwarm", delta_norm),
            (f"{metric}_level2", "Raw Level2", "viridis", state_norm),
            (f"{metric}_recon", "Reconstructed", "viridis", state_norm),
            (f"delta_{metric}_vs_raw_level2", "Recon − Raw Level2", "coolwarm", delta_norm),
        ]
        fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.2), sharex=True, sharey=True)
        for ax, (column, label, cmap, norm) in zip(axes.flat, panels):
            image = self.window_field(ax, metrics, column, side_um=side, cmap=cmap, norm=norm, title=label)
            if image is not None:
                plt.colorbar(image, ax=ax, fraction=.046, pad=.04)
        fig.suptitle(title); fig.tight_layout(); return fig

    def plot_high_diversity(self, run, title):
        metrics = run["impact"].window_metrics.copy()
        side = float(run["impact"].scale_audit["main_window_side_um"])
        valid = metrics.loc[metrics["valid_window"]].copy()
        maximum = max(1.0, float(valid["neff_recon"].max()))
        fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), sharex=True, sharey=True)
        image = self.window_field(axes[0], valid, "neff_recon", side_um=side, cmap="viridis", norm=plt.Normalize(1, maximum), title="Reconstructed Neff")
        if image is not None:
            plt.colorbar(image, ax=axes[0], fraction=.046, pad=.04)
        threshold_ok = run["impact"].state_threshold.get("status") == "ok"
        mask = valid.assign(region_mask=valid["in_state_region"].fillna(False).astype(float))
        if threshold_ok:
            image = self.window_field(axes[1], mask, "region_mask", side_um=side, cmap=ListedColormap(["#efefef", "#fdae61"]), norm=BoundaryNorm([-0.5, 0.5, 1.5], 2), title="High-diversity Region")
            if image is not None:
                colorbar = plt.colorbar(image, ax=axes[1], ticks=[0, 1], fraction=.046, pad=.04)
                colorbar.ax.set_yticklabels(["outside", "inside"])
        else:
            self.window_field(axes[1], mask, "region_mask", side_um=side, cmap=ListedColormap(["#efefef"]), norm=BoundaryNorm([-0.5, 0.5], 1), title="No stable Region threshold")
        fig.suptitle(title); fig.tight_layout(); return fig

    def plot_moran_summary(self, summary):
        fig, ax = plt.subplots(figsize=(8.5, 3.8))
        long = summary.melt(id_vars="scope", value_vars=["raw_q75", "reconstruction_q75"], var_name="side", value_name="Moran I Q75")
        sns.barplot(data=long, x="scope", y="Moran I Q75", hue="side", ax=ax, palette=["#4c78a8", "#f58518"])
        ax.set(title="Full-gene Moran I by cell type", xlabel="")
        fig.tight_layout(); return fig

    def plot_moran_q75(self, summary, scope_order):
        """Plot all-valid and shared-valid Moran Q75 summaries on one scale."""

        all_valid = summary.set_index("scope").reindex(scope_order)[["raw_q75", "reconstruction_q75"]].rename(columns={"raw_q75": "Raw", "reconstruction_q75": "Reconstruction"})
        shared_valid = summary.set_index("scope").reindex(scope_order)[["shared_raw_q75", "shared_reconstruction_q75"]].rename(columns={"shared_raw_q75": "Raw", "shared_reconstruction_q75": "Reconstruction"})
        values = np.concatenate([all_valid.to_numpy(float).ravel(), shared_valid.to_numpy(float).ravel()])
        values = values[np.isfinite(values)]
        fig, axes = plt.subplots(1, 2, figsize=(10.4, max(2.8, 0.55 * len(all_valid))), sharey=True, constrained_layout=True)
        if not len(values):
            for ax in axes:
                ax.axis("off")
                ax.text(.5, .5, "No computed Moran values", ha="center", va="center")
        else:
            limit = max(float(np.abs(values).max()), 1e-3)
            norm = TwoSlopeNorm(vcenter=0, vmin=-limit, vmax=limit)
            image = None
            for ax, table, title in zip(axes, (all_valid, shared_valid), ("All valid genes", "Shared valid genes")):
                if table.dropna(how="all").empty:
                    ax.axis("off")
                    ax.text(.5, .5, f"{title}: no values", ha="center", va="center")
                    continue
                image = sns.heatmap(table, annot=True, fmt=".3f", cmap="RdBu_r", norm=norm, cbar=False, ax=ax, linewidths=.5, linecolor="white")
                ax.set(xlabel="", ylabel="scope", title=title)
                ax.tick_params(axis="y", labelrotation=0)
            if image is not None:
                fig.colorbar(image.collections[0], ax=axes, label="Moran I Q75", fraction=.046, pad=.04)
        fig.suptitle("Moran I Q75: all-valid and shared-valid gene sets")
        return fig

    def plot_moran_distribution(self, distribution, scope_order):
        """Plot all-valid and shared-valid Moran distributions."""

        fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
        if distribution.empty:
            for ax in axes:
                ax.axis("off")
                ax.text(.5, .5, "No computed Moran values", ha="center", va="center")
        else:
            values = distribution["moran_i"].to_numpy(float)
            values = values[np.isfinite(values)]
            limit = max(float(np.abs(values).max()), 1e-3)
            for ax, gene_set, title in zip(axes, ("all_valid", "shared_valid"), ("All valid genes", "Shared valid genes")):
                sns.boxplot(data=distribution.loc[distribution["gene_set"].eq(gene_set)], x="scope", y="moran_i", hue="side", order=scope_order, hue_order=["Raw", "Reconstruction"], palette=["#4c78a8", "#f58518"], showfliers=False, ax=ax)
                ax.set(title=title, xlabel="", ylabel="Moran I", ylim=(-limit, limit))
                if ax.legend_ is not None:
                    ax.legend_.remove()
            axes[0].legend(title="side", loc="best")
        fig.suptitle("Moran I distributions by scope")
        fig.tight_layout()
        return fig

    def plot_moran_shared_scatter(self, moran, scope_order):
        """Plot paired Moran values for genes computed on both sides."""

        scopes = [scope for scope in scope_order if not moran.loc[
            moran["scope"].eq(scope)
            & moran["raw_status"].eq("computed")
            & moran["reconstruction_status"].eq("computed")
        ].empty]
        n_panels = max(1, len(scopes))
        fig, axes = plt.subplots(1, n_panels, figsize=(4.3 * n_panels, 4.0), squeeze=False)
        for ax, scope in zip(axes[0], scopes):
            shared = moran.loc[
                moran["scope"].eq(scope)
                & moran["raw_status"].eq("computed")
                & moran["reconstruction_status"].eq("computed")
            ]
            x = shared["raw_moran"].to_numpy(float)
            y = shared["reconstruction_moran"].to_numpy(float)
            finite = np.isfinite(x) & np.isfinite(y)
            x, y = x[finite], y[finite]
            if len(x):
                lo, hi = float(min(x.min(), y.min())), float(max(x.max(), y.max()))
                pad = max((hi - lo) * .05, 1e-3)
                ax.scatter(x, y, s=7, alpha=.45, linewidths=0, rasterized=True)
                ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="black", ls="--", lw=1)
                ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad)
            ax.set(title=f"{scope} (n={len(x)})", xlabel="Raw Moran I", ylabel="Reconstruction Moran I")
        for ax in axes[0][len(scopes):]:
            ax.axis("off")
        if not scopes:
            axes[0][0].axis("off")
            axes[0][0].text(.5, .5, "No shared computed genes", ha="center", va="center")
        fig.suptitle("Shared-gene Moran I comparison")
        fig.tight_layout()
        return fig

    def plot_moran_paired_delta(self, paired_delta, scope_order):
        """Plot paired reconstruction-minus-Raw Moran deltas."""

        fig, ax = plt.subplots(figsize=(8.5, 4.2))
        if paired_delta.empty:
            ax.axis("off")
            ax.text(.5, .5, "No shared computed genes", ha="center", va="center")
        else:
            sns.boxplot(data=paired_delta, x="scope", y="delta", order=scope_order, color="#7b3294", showfliers=False, ax=ax)
            ax.axhline(0, color="black", linestyle="--", linewidth=.8)
            ax.set(title="Shared-gene paired Moran I delta", xlabel="", ylabel="Reconstruction − Raw Moran I")
        fig.tight_layout()
        return fig

    def plot_score_distribution(self, scores):
        long = scores.melt(id_vars=["scope", "unit_id"], value_vars=["raw_score", "reconstruction_score"], var_name="side", value_name="AUCell score").dropna()
        fig, ax = plt.subplots(figsize=(8.5, 3.8))
        sns.boxplot(data=long, x="scope", y="AUCell score", hue="side", ax=ax, palette=["#4c78a8", "#f58518"], fliersize=0)
        ax.set(title="EMT AUCell distribution by cell type", xlabel="")
        fig.tight_layout(); return fig

    def plot_aucell_overall(self, long, delta, scope_order, scope_labels, pathway):
        """Plot overall AUCell distributions and paired deltas."""

        fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
        order = [scope_labels.get(scope, scope) for scope in scope_order]
        if long.empty:
            axes[0].axis("off")
            axes[0].text(.5, .5, "No finite Raw/Reconstructed AUCell scores", ha="center", va="center")
        else:
            sns.boxplot(data=long, x="scope_label", y="AUCell score", hue="side", order=order, hue_order=["Raw", "Reconstructed"], ax=axes[0], palette=["#4c78a8", "#f58518"], fliersize=0)
            axes[0].set(title=f"{pathway}: Raw / reconstruction", xlabel="", ylabel="AUCell score")
        if delta.empty:
            axes[1].axis("off")
            axes[1].text(.5, .5, "No paired AUCell scores", ha="center", va="center")
        else:
            delta = delta.copy()
            delta["scope_label"] = delta["scope"].map(scope_labels).fillna(delta["scope"])
            sns.boxplot(data=delta, x="scope_label", y="delta", order=order, ax=axes[1], color="#7b3294", fliersize=0)
            axes[1].axhline(0, color="black", linestyle="--", linewidth=.8)
            axes[1].set(title="Paired reconstruction − Raw", xlabel="", ylabel="AUCell delta")
        fig.suptitle("Overall EMT AUCell comparison")
        fig.tight_layout()
        return fig

    def plot_aucell_parent(self, long, delta, scope_label, pathway):
        """Plot one parent's published AUCell distribution and delta."""

        fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))
        if long.empty:
            axes[0].axis("off")
            axes[0].text(.5, .5, "No finite Raw/Reconstructed AUCell scores", ha="center", va="center")
        else:
            sns.boxplot(data=long, x="side", y="AUCell score", order=["Raw", "Reconstructed"], hue="side", hue_order=["Raw", "Reconstructed"], ax=axes[0], palette=["#4c78a8", "#f58518"], fliersize=0, legend=False)
            axes[0].set(title=f"{scope_label}: Raw / reconstruction", xlabel="", ylabel="AUCell score")
        if delta.empty:
            axes[1].axis("off")
            axes[1].text(.5, .5, "No paired AUCell scores", ha="center", va="center")
        else:
            sns.boxplot(data=delta, y="delta", color="#7b3294", ax=axes[1], fliersize=0)
            axes[1].axhline(0, color="black", linestyle="--", linewidth=.8)
            axes[1].set(title="Paired reconstruction − Raw", xlabel="", ylabel="AUCell delta")
        fig.suptitle(f"{scope_label}: {pathway}")
        fig.tight_layout()
        return fig

    def plot_feature_field(self, frame, feature, title):
        subset = frame.loc[frame["feature"].eq(feature)].copy()
        fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharex=True, sharey=True)
        values = subset[["raw_value", "reconstruction_value"]].to_numpy(float)
        finite = values[np.isfinite(values)]
        norm = plt.Normalize(finite.min(), finite.max()) if finite.size and finite.min() != finite.max() else None
        delta = subset["delta"].to_numpy(float); delta_max = max(float(np.nanmax(np.abs(delta))), 0.01) if np.isfinite(delta).any() else 1.0
        panels = [("raw_value", "Raw", "viridis", norm), ("reconstruction_value", "Reconstructed", "viridis", norm), ("delta", "Recon − Raw", "coolwarm", TwoSlopeNorm(vcenter=0, vmin=-delta_max, vmax=delta_max))]
        for ax, (column, label, cmap, panel_norm) in zip(axes, panels):
            valid = subset.dropna(subset=[column])
            if valid.empty:
                self.set_spatial_axes(ax)
                ax.set_title(label + " (unavailable)")
                continue
            image = ax.scatter(valid.x, valid.y, c=valid[column], s=1, cmap=cmap, norm=panel_norm, linewidths=0, rasterized=True)
            self.set_spatial_axes(ax); ax.set_title(label); fig.colorbar(image, ax=ax, fraction=.046, pad=.04)
        fig.suptitle(title); fig.tight_layout(); return fig

    def save_json(self, value, name):
        destination = self.output_dir / "analysis" / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        # Pandas serializes missing numerical values as JSON null.
        destination.write_text(pd.Series(value, dtype=object).to_json(indent=2), encoding="utf-8")


def plot_threshold_reliability(runs, n_bootstrap):
    """Display existing threshold diagnostics before any binary Region map."""
    rows = []
    for parent, run in runs.items():
        spatial = run['spatial']
        threshold = spatial['state_threshold']
        valid = spatial['metrics'].loc[lambda frame: frame.valid_window, 'neff_recon'].dropna()
        span = valid.max() - valid.min() if len(valid) else np.nan
        lower, upper = threshold.get('ci_lower'), threshold.get('ci_upper')
        relative_width = ((upper - lower) / span if lower is not None and upper is not None and span > 0 else np.nan)
        rows.append({'scope': parent, 'scale_um': spatial['side_um'],
                     'valid_windows': len(valid), 'neff_range': span,
                     'relative_ci_width': relative_width,
                     'valid_bootstrap_fraction': threshold.get('n_valid_bootstrap', 0) / n_bootstrap,
                     'threshold': threshold.get('threshold'), 'status': threshold['status']})
    table = pd.DataFrame(rows)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    labels = table.scope.str.replace('_', '/', regex=False)
    for ax, column, reference, label in zip(axes,
            ('relative_ci_width', 'valid_bootstrap_fraction'), (.25, .8),
            ('Threshold CI width / Neff range', 'Valid bootstrap fraction')):
        ax.bar(labels, table[column], color='#4c78a8')
        ax.axhline(reference, color='black', linestyle='--', label=f'Criterion: {reference:.0%}')
        ax.set(ylabel=label)
        ax.legend()
        for index, value in enumerate(table[column]):
            if pd.notna(value): ax.annotate(f'{value:.1%}', (index, value), xytext=(0, 4), textcoords='offset points', ha='center')
        ax.margins(y=.2)
    fig.suptitle('State Region identifiability at the selected scale')
    fig.tight_layout()
    return table, fig
