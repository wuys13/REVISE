"""Bilingual, factual reading candidates for saved reconstruction-impact records.

The functions in this module are presentation helpers.  They consume the
already-built records and do not read scientific artefacts or calculate a new
metric.  Every sentence is deliberately phrased as a candidate for review so
that the existing record and section-summary review process remains the
authority for publication.
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Any, Iterable, Mapping

from .content_contract import NODE_SPECS


# These are the reader-facing nodes whose records need a bilingual batch
# interpretation.  The existing 2.3 and layer-end summary nodes remain
# represented by their current record/summary mechanisms.
NARRATIVE_NODE_IDS = (
    "2.1",
    "2.2",
    "2.4",
    "2.5",
    "3.1",
    "3.2",
    "3.3",
    "3.4",
    "3.5",
    "3.6",
    "3.7",
    "3.8",
    "4.1",
    "4.2",
    "4.3",
    "4.4",
)

_NARRATIVE_SPECS = {
    str(spec["node_id"]): spec
    for spec in NODE_SPECS
    if str(spec["node_id"]) in NARRATIVE_NODE_IDS
}
_SCOPE_ORDER = ("All", "Fibroblast", "Mono_Macro", "T")


def _record_list(records: Mapping[str, Any] | Iterable[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    """Accept either a package or its record list without mutating it."""

    if isinstance(records, Mapping):
        records = records.get("records", [])
    if not isinstance(records, Iterable) or isinstance(records, (str, bytes)):
        return []
    return [record for record in records if isinstance(record, Mapping)]


def _identity(record: Mapping[str, Any]) -> Mapping[str, Any]:
    value = record.get("identity", {})
    return value if isinstance(value, Mapping) else {}


def _results(record: Mapping[str, Any] | None) -> Mapping[str, Any]:
    if not isinstance(record, Mapping):
        return {}
    value = record.get("results", {})
    return value if isinstance(value, Mapping) else {}


def _number(value: Any) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _fmt(value: Any) -> str:
    """Keep compact numeric values readable while preserving NA explicitly."""

    number = _number(value)
    if number is None:
        if value is None:
            return "NA"
        if isinstance(value, bool):
            return "true" if value else "false"
        return str(value)
    if number == 0:
        return "0"
    if number.is_integer():
        return str(int(number))
    return f"{number:.4g}"


def _fmt_fraction(value: Any) -> str:
    number = _number(value)
    if number is None:
        return "NA"
    return f"{_fmt(number)} ({number * 100:.4g}%)"


def _fmt_percent(value: Any) -> str:
    number = _number(value)
    return "NA" if number is None else f"{number * 100:.4g}%"


def _fmt_bool(value: Any) -> str:
    if value is True:
        return "true"
    if value is False:
        return "false"
    return "NA"


def _scope(record: Mapping[str, Any] | None) -> str:
    value = _identity(record or {}).get("scope")
    return str(value) if value not in (None, "") else "unknown"


def _scope_order(values: Iterable[str]) -> list[str]:
    rank = {value: index for index, value in enumerate(_SCOPE_ORDER)}
    return sorted(set(values), key=lambda value: (rank.get(value, len(rank)), value))


def _question(record: Mapping[str, Any]) -> str:
    value = _identity(record).get("question")
    return str(value) if value not in (None, "") else ""


def _records_for_questions(record_list: Iterable[Mapping[str, Any]], questions: Iterable[str]) -> list[Mapping[str, Any]]:
    wanted = set(questions)
    return [record for record in record_list if _question(record) in wanted]


def _main_records(record_list: Iterable[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    """Keep parent-internal readings separate from the HD global cohort."""

    return [record for record in record_list if _scope(record) != "All"]


def _by_scope(record_list: Iterable[Mapping[str, Any]], question: str) -> dict[str, Mapping[str, Any]]:
    result: dict[str, Mapping[str, Any]] = {}
    for record in record_list:
        if _question(record) == question:
            result.setdefault(_scope(record), record)
    return result


def _by_scope_questions(record_list: Iterable[Mapping[str, Any]], questions: Iterable[str]) -> dict[str, dict[str, Mapping[str, Any]]]:
    result: dict[str, dict[str, Mapping[str, Any]]] = defaultdict(dict)
    wanted = set(questions)
    for record in record_list:
        question = _question(record)
        if question in wanted:
            result[_scope(record)].setdefault(question, record)
    return dict(result)


def _record_ids(records: Iterable[Mapping[str, Any]]) -> list[str]:
    result: list[str] = []
    for record in records:
        record_id = _identity(record).get("record_id")
        if record_id is not None:
            result.append(str(record_id))
    return list(dict.fromkeys(result))


def _nested(mapping: Mapping[str, Any], *keys: str) -> Any:
    value: Any = mapping
    for key in keys:
        if not isinstance(value, Mapping):
            return None
        value = value.get(key)
    return value


def _display_list(values: Any, limit: int = 8) -> str:
    if not isinstance(values, (list, tuple)):
        return "NA"
    items = [str(value) for value in values if value not in (None, "")]
    if len(items) > limit:
        return ", ".join(items[:limit]) + f", +{len(items) - limit} more"
    return ", ".join(items) if items else "NA"


def _direction(value: Any, *, zh: bool = False) -> str:
    number = _number(value)
    if number is None:
        return "未确定" if zh else "unknown"
    if number > 0:
        return "升高" if zh else "increased"
    if number < 0:
        return "降低" if zh else "decreased"
    return "无方向变化" if zh else "unchanged"


def _support_row(results: Mapping[str, Any]) -> Mapping[str, Any] | None:
    rows = results.get("support_sensitivity")
    main = _number(results.get("main_window_side_um"))
    if not isinstance(rows, list):
        return None
    for row in rows:
        if not isinstance(row, Mapping):
            continue
        side = _number(row.get("window_side_length"))
        if main is not None and side is not None and abs(main - side) <= 1e-9:
            return row
    return None


def _metric_values(record: Mapping[str, Any] | None, metric: str) -> Mapping[str, Any]:
    values = _results(record).get("metrics", {})
    if not isinstance(values, Mapping):
        return {}
    # Existing records use Kobs, Neff, and evenness.  The lower-case fallback
    # keeps this reader compatible with hand-built records in review fixtures.
    value = values.get(metric)
    if value is None:
        value = values.get(metric.lower())
    return value if isinstance(value, Mapping) else {}


def _baseline_metric(record: Mapping[str, Any] | None, metric: str, field: str) -> Any:
    return _metric_values(record, metric).get(field)


def _missing_line(scope: str, question: str, *, zh: bool) -> str:
    if zh:
        return f"{scope}：没有保存的 {question} record，解读保持 NA。"
    return f"{scope}: no saved {question} record was found; the interpretation remains NA."


def _complexity_lines(record_list: list[Mapping[str, Any]]) -> tuple[list[str], list[str]]:
    by_scope = _by_scope(record_list, "complexity")
    zh_lines: list[str] = []
    en_lines: list[str] = []
    for scope in _scope_order(by_scope):
        record = by_scope[scope]
        results = _results(record)
        kind = results.get("comparison_kind")
        ari = _fmt(results.get("ari"))
        status = _fmt(results.get("status"))
        if kind == "raw_reference_vs_fixed_final_clusters":
            raw = _fmt(results.get("raw_k_reference_resolution", results.get("raw_k")))
            recon = _fmt(results.get("fixed_final_cluster_k"))
            zh_lines.append(f"{scope}：reference-resolution Raw K={raw}，固定最终 cluster K={recon}，ARI={ari}，status={status}。")
            en_lines.append(f"{scope}: reference-resolution Raw K={raw}, fixed final-cluster K={recon}, ARI={ari}, status={status}.")
        else:
            raw = _fmt(results.get("raw_k"))
            recon = _fmt(results.get("reconstruction_k_at_raw_resolution"))
            zh_lines.append(f"{scope}：Raw K={raw}，重建 Raw-resolution K={recon}，ARI={ari}，status={status}。")
            en_lines.append(f"{scope}: Raw K={raw}, reconstruction-at-raw-resolution K={recon}, ARI={ari}, status={status}.")
    if not zh_lines:
        zh_lines.append(_missing_line("所有 scope", "complexity", zh=True))
        en_lines.append(_missing_line("all scopes", "complexity", zh=False))
    return zh_lines, en_lines


def _matched_lines(record_list: list[Mapping[str, Any]]) -> tuple[list[str], list[str]]:
    by_scope = _by_scope(record_list, "matched_k")
    zh_lines: list[str] = []
    en_lines: list[str] = []
    for scope in _scope_order(by_scope):
        results = _results(by_scope[scope])
        matched = results.get("matched")
        status = _fmt(results.get("status"))
        raw = _fmt(results.get("raw_k"))
        recon = _fmt(results.get("reconstruction_k"))
        fraction = _fmt_fraction(results.get("unit_change_fraction"))
        balanced = _fmt_fraction(results.get("balanced_change"))
        ari = _fmt(results.get("ari"))
        eligible = _fmt_bool(results.get("headline_eligible"))
        if matched is True:
            zh_lines.append(f"{scope}：matched-K 可用（matched=true，status={status}，headline_eligible={eligible}）；Raw K={raw}、Recon K={recon}，assignment change 比例={fraction}，balanced change={balanced}，ARI={ari}。")
            en_lines.append(f"{scope}: matched-K is available (matched=true, status={status}, headline_eligible={eligible}); Raw K={raw}, Recon K={recon}, assignment-change fraction={fraction}, balanced change={balanced}, ARI={ari}.")
        else:
            zh_lines.append(f"{scope}：matched-K 不可用（matched={_fmt_bool(matched)}，status={status}，headline_eligible={eligible}）；Raw K={raw}、Recon K={recon}，assignment change 比例={fraction}，该结果不作 headline 或下游 assignment 结论。")
            en_lines.append(f"{scope}: matched-K is unavailable (matched={_fmt_bool(matched)}, status={status}, headline_eligible={eligible}); Raw K={raw}, Recon K={recon}, assignment-change fraction={fraction}. It is not used for a headline or downstream assignment interpretation.")
    if not zh_lines:
        zh_lines.append(_missing_line("所有 scope", "matched_k", zh=True))
        en_lines.append(_missing_line("all scopes", "matched_k", zh=False))
    return zh_lines, en_lines


def _moran_lines(record_list: list[Mapping[str, Any]]) -> tuple[list[str], list[str]]:
    grouped = _by_scope_questions(record_list, ("moran_all_valid", "moran_shared_valid"))
    zh_lines: list[str] = []
    en_lines: list[str] = []
    for scope in _scope_order(grouped):
        all_results = _results(grouped[scope].get("moran_all_valid"))
        shared_results = _results(grouped[scope].get("moran_shared_valid"))
        raw_all = all_results.get("raw", {}) if isinstance(all_results.get("raw", {}), Mapping) else {}
        recon_all = all_results.get("reconstruction", {}) if isinstance(all_results.get("reconstruction", {}), Mapping) else {}
        raw_shared = shared_results.get("raw", {}) if isinstance(shared_results.get("raw", {}), Mapping) else {}
        recon_shared = shared_results.get("reconstruction", {}) if isinstance(shared_results.get("reconstruction", {}), Mapping) else {}
        all_delta = all_results.get("median_difference")
        shared_delta = _nested(shared_results, "paired_delta", "median")
        shared_n = _nested(shared_results, "paired_delta", "n")
        missing = all_results.get("missing_status_counts") or shared_results.get("missing_status_counts") or {}
        missing_items = []
        for key in ("raw_unmeasured", "raw_not_computable", "raw_insufficient_support", "reconstruction_unmeasured", "reconstruction_not_computable", "reconstruction_insufficient_support"):
            if isinstance(missing, Mapping) and missing.get(key) is not None:
                missing_items.append(f"{key}={_fmt(missing.get(key))}")
        missing_text_zh = f"；missing-status: {', '.join(missing_items)}" if missing_items else ""
        missing_text_en = f"; missing-status counts: {', '.join(missing_items)}" if missing_items else ""
        zh_lines.append(
            f"{scope}：all_valid（各自有效基因）Raw n={_fmt(raw_all.get('n_valid'))}、median={_fmt(raw_all.get('median'))}，Recon n={_fmt(recon_all.get('n_valid'))}、median={_fmt(recon_all.get('median'))}，marginal median difference={_fmt(all_delta)}；shared_valid（共同有效基因）paired n={_fmt(shared_n)}，Raw median={_fmt(raw_shared.get('median'))}，Recon median={_fmt(recon_shared.get('median'))}，paired median Δ={_fmt(shared_delta)}。两种视角分开解释{missing_text_zh}。"
        )
        en_lines.append(
            f"{scope}: all_valid (side-specific valid genes) has Raw n={_fmt(raw_all.get('n_valid'))}, median={_fmt(raw_all.get('median'))}, and Recon n={_fmt(recon_all.get('n_valid'))}, median={_fmt(recon_all.get('median'))}; marginal median difference={_fmt(all_delta)}. shared_valid (shared valid genes) has paired n={_fmt(shared_n)}, Raw median={_fmt(raw_shared.get('median'))}, Recon median={_fmt(recon_shared.get('median'))}, and paired median delta={_fmt(shared_delta)}. The two views remain separate{missing_text_en}."
        )
    if not zh_lines:
        zh_lines.append(_missing_line("所有 scope", "Moran", zh=True))
        en_lines.append(_missing_line("all scopes", "Moran", zh=False))
    return zh_lines, en_lines


def _emt_lines(record_list: list[Mapping[str, Any]]) -> tuple[list[str], list[str]]:
    grouped = _by_scope_questions(record_list, ("emt_coverage", "emt_score"))
    zh_lines: list[str] = []
    en_lines: list[str] = []
    for scope in _scope_order(grouped):
        coverage = _results(grouped[scope].get("emt_coverage"))
        score = _results(grouped[scope].get("emt_score"))
        raw_cov = coverage.get("raw", {}) if isinstance(coverage.get("raw", {}), Mapping) else {}
        recon_cov = coverage.get("reconstruction", {}) if isinstance(coverage.get("reconstruction", {}), Mapping) else {}
        raw_score = score.get("raw", {}) if isinstance(score.get("raw", {}), Mapping) else {}
        recon_score = score.get("reconstruction", {}) if isinstance(score.get("reconstruction", {}), Mapping) else {}
        paired = score.get("paired", {}) if isinstance(score.get("paired", {}), Mapping) else {}
        raw_cutoff = _nested(score, "cutoff_audit", "raw") or {}
        recon_cutoff = _nested(score, "cutoff_audit", "reconstruction") or {}
        if not isinstance(raw_cutoff, Mapping):
            raw_cutoff = {}
        if not isinstance(recon_cutoff, Mapping):
            recon_cutoff = {}
        carrier = score.get("carrier_condition", coverage.get("carrier_condition"))
        zh_lines.append(
            f"{scope}：EMT coverage 为 Raw={_fmt(raw_cov.get('coverage'))}（{_fmt_percent(raw_cov.get('coverage'))}，{_fmt(raw_cov.get('available_gene_count'))}/{_fmt(raw_cov.get('resource_gene_count'))} genes）到 Recon={_fmt(recon_cov.get('coverage'))}（{_fmt_percent(recon_cov.get('coverage'))}，{_fmt(recon_cov.get('available_gene_count'))}/{_fmt(recon_cov.get('resource_gene_count'))} genes）；score median 为 Raw={_fmt(raw_score.get('median'))}、Recon={_fmt(recon_score.get('median'))}，paired median Δ={_fmt(paired.get('median_delta'))}（{_direction(paired.get('median_delta'), zh=True)}）。cutoff effective rank length 为 Raw={_fmt(raw_cutoff.get('effective_rank_length'))}、Recon={_fmt(recon_cutoff.get('effective_rank_length'))}，provider AUC threshold 为 Raw={_fmt(raw_cutoff.get('provider_auc_threshold'))}、Recon={_fmt(recon_cutoff.get('provider_auc_threshold'))}；carrier={_fmt(carrier)}。coverage、cutoff 和 score 分开读取。"
        )
        en_lines.append(
            f"{scope}: EMT coverage is Raw={_fmt(raw_cov.get('coverage'))} ({_fmt_percent(raw_cov.get('coverage'))}, {_fmt(raw_cov.get('available_gene_count'))}/{_fmt(raw_cov.get('resource_gene_count'))} genes) to Recon={_fmt(recon_cov.get('coverage'))} ({_fmt_percent(recon_cov.get('coverage'))}, {_fmt(recon_cov.get('available_gene_count'))}/{_fmt(recon_cov.get('resource_gene_count'))} genes); score median is Raw={_fmt(raw_score.get('median'))} versus Recon={_fmt(recon_score.get('median'))}, with paired median delta={_fmt(paired.get('median_delta'))} ({_direction(paired.get('median_delta'))}). Cutoff effective rank length is Raw={_fmt(raw_cutoff.get('effective_rank_length'))} and Recon={_fmt(recon_cutoff.get('effective_rank_length'))}; provider AUC threshold is Raw={_fmt(raw_cutoff.get('provider_auc_threshold'))} and Recon={_fmt(recon_cutoff.get('provider_auc_threshold'))}; carrier={_fmt(carrier)}. Coverage, cutoffs, and scores remain separate."
        )
    if not zh_lines:
        zh_lines.append(_missing_line("所有 scope", "EMT coverage/score", zh=True))
        en_lines.append(_missing_line("all scopes", "EMT coverage/score", zh=False))
    return zh_lines, en_lines


def _anatomy_lines(record_list: list[Mapping[str, Any]]) -> tuple[list[str], list[str]]:
    by_scope = _by_scope(record_list, "anatomy")
    zh_lines: list[str] = []
    en_lines: list[str] = []
    for scope in _scope_order(by_scope):
        results = _results(by_scope[scope])
        regions = results.get("regions", [])
        support = results.get("support_scales", [])
        scales = [_number(row.get("window_side_length")) for row in support if isinstance(row, Mapping)] if isinstance(support, list) else []
        scales = [value for value in scales if value is not None]
        scale_text = f"{_fmt(min(scales))}–{_fmt(max(scales))} μm" if scales else "NA"
        names = [row.get("level1_region", row.get("level1")) for row in regions if isinstance(row, Mapping)] if isinstance(regions, list) else []
        names_text = _display_list(names)
        zh_lines.append(f"{scope}：保存的 anatomy context 包含 n_regions={_fmt(results.get('n_regions', len(regions) if isinstance(regions, list) else None))}（{names_text}），支持窗口尺度覆盖 {scale_text}；这是空间背景与支持条件。")
        en_lines.append(f"{scope}: the saved anatomy context contains n_regions={_fmt(results.get('n_regions', len(regions) if isinstance(regions, list) else None))} ({names_text}), with support scales spanning {scale_text}; it supplies spatial context and support conditions.")
    if not zh_lines:
        zh_lines.append(_missing_line("所有 scope", "anatomy", zh=True))
        en_lines.append(_missing_line("all scopes", "anatomy", zh=False))
    return zh_lines, en_lines


def _changed_lines(record_list: list[Mapping[str, Any]]) -> tuple[list[str], list[str]]:
    by_scope = _by_scope(record_list, "changed_units")
    zh_lines: list[str] = []
    en_lines: list[str] = []
    for scope in _scope_order(by_scope):
        results = _results(by_scope[scope])
        overall = results.get("overall", {}) if isinstance(results.get("overall", {}), Mapping) else {}
        changed = overall.get("changed_units")
        total = overall.get("paired_units", overall.get("total_units"))
        fraction = overall.get("change_fraction")
        scale = overall.get("scale_um")
        valid_windows = overall.get("n_valid_windows")
        median_window = overall.get("median_window_change")
        matched = overall.get("matched_cluster_status")
        eligible = overall.get("headline_eligible")
        strata = results.get("stratification")
        support = f"scale={_fmt(scale)} μm, valid windows={_fmt(valid_windows)}, median window change={_fmt_fraction(median_window)}" if any(value is not None for value in (scale, valid_windows, median_window)) else "window support=NA"
        zh_lines.append(f"{scope}：assignment change 数量={_fmt(changed)}/{_fmt(total)}，比例={_fmt_fraction(fraction)}；{support}；stratification={_fmt(strata)}，matched_cluster_status={_fmt(matched)}，headline_eligible={_fmt_bool(eligible)}。数量、比例和窗口支持分别保留。")
        en_lines.append(f"{scope}: assignment-change count={_fmt(changed)}/{_fmt(total)}, fraction={_fmt_fraction(fraction)}; {support}; stratification={_fmt(strata)}, matched_cluster_status={_fmt(matched)}, headline_eligible={_fmt_bool(eligible)}. Count, fraction, and window support are kept separate.")
    if not zh_lines:
        zh_lines.append(_missing_line("所有 scope", "changed_units", zh=True))
        en_lines.append(_missing_line("all scopes", "changed_units", zh=False))
    return zh_lines, en_lines


def _emt_spatial_lines(record_list: list[Mapping[str, Any]]) -> tuple[list[str], list[str]]:
    by_scope = _by_scope(record_list, "emt_spatial_fields")
    zh_lines: list[str] = []
    en_lines: list[str] = []
    for scope in _scope_order(by_scope):
        results = _results(by_scope[scope])
        metadata = results.get("carrier_metadata", {}) if isinstance(results.get("carrier_metadata", {}), Mapping) else {}
        metadata_flags = ", ".join(f"{key}={_fmt_bool(value)}" for key, value in metadata.items()) or "NA"
        zh_lines.append(f"{scope}：EMT spatial field metadata={_fmt_bool(results.get('available'))}，n_units={_fmt(results.get('n_units'))}，features={_display_list(results.get('features'))}，expression_layers={_display_list(results.get('expression_layers'))}，coordinate_units={_display_list(results.get('coordinate_units'))}，carrier={_fmt(results.get('carrier_condition'))}，carrier_metadata={metadata_flags}。当前只形成字段 metadata 记录，未形成区域定量判读，也不编富集结论。")
        en_lines.append(f"{scope}: EMT spatial-field metadata available={_fmt_bool(results.get('available'))}, n_units={_fmt(results.get('n_units'))}, features={_display_list(results.get('features'))}, expression_layers={_display_list(results.get('expression_layers'))}, coordinate_units={_display_list(results.get('coordinate_units'))}, carrier={_fmt(results.get('carrier_condition'))}, carrier_metadata={metadata_flags}. Only field metadata are recorded here; no region-level quantitative interpretation or enrichment claim is formed.")
    if not zh_lines:
        zh_lines.append(_missing_line("所有 scope", "emt_spatial_fields", zh=True))
        en_lines.append(_missing_line("all scopes", "emt_spatial_fields", zh=False))
    return zh_lines, en_lines


def _window_lines(record_list: list[Mapping[str, Any]]) -> tuple[list[str], list[str]]:
    by_scope = _by_scope(record_list, "window_support")
    zh_lines: list[str] = []
    en_lines: list[str] = []
    for scope in _scope_order(by_scope):
        results = _results(by_scope[scope])
        selected = _support_row(results) or {}
        valid = selected.get("n_valid_windows", results.get("n_valid_windows"))
        total = selected.get("n_tissue_windows", results.get("total_windows"))
        valid_fraction = selected.get("valid_window_fraction")
        retained = selected.get("retained_parent_units", results.get("valid_units"))
        scale = results.get("main_window_side_um")
        zh_lines.append(f"{scope}：选定窗口={_fmt(scale)} μm，有效窗口={_fmt(valid)}/{_fmt(total)}（{_fmt_percent(valid_fraction)}），保留 parent units={_fmt(retained)}；min parent units={_fmt(results.get('min_parent_units'))}，rarefaction draws={_fmt(results.get('rarefaction_draws'))}。")
        en_lines.append(f"{scope}: selected window={_fmt(scale)} μm, valid windows={_fmt(valid)}/{_fmt(total)} ({_fmt_percent(valid_fraction)}), retained parent units={_fmt(retained)}; minimum parent units={_fmt(results.get('min_parent_units'))}, rarefaction draws={_fmt(results.get('rarefaction_draws'))}.")
    if not zh_lines:
        zh_lines.append(_missing_line("所有 scope", "window_support", zh=True))
        en_lines.append(_missing_line("all scopes", "window_support", zh=False))
    return zh_lines, en_lines


def _local_metric_lines(record_list: list[Mapping[str, Any]], metric: str) -> tuple[list[str], list[str]]:
    grouped = _by_scope_questions(record_list, ("local_vs_raw_leiden", "local_vs_raw_level2"))
    label = metric
    zh_lines: list[str] = []
    en_lines: list[str] = []
    for scope in _scope_order(grouped):
        leiden = grouped[scope].get("local_vs_raw_leiden")
        level2 = grouped[scope].get("local_vs_raw_level2")
        if leiden is None and level2 is None:
            continue
        leiden_values = _metric_values(leiden, metric)
        level2_values = _metric_values(level2, metric)
        scale = _results(leiden or level2).get("scale_um")
        windows = _results(leiden or level2).get("n_valid_windows")
        zh_lines.append(
            f"{scope}：{label} 在 scale={_fmt(scale)} μm、valid windows={_fmt(windows)} 下，Raw Leiden baseline 为 {_fmt(leiden_values.get('baseline_median'))}→Recon {_fmt(leiden_values.get('reconstruction_median'))}（delta={_fmt(leiden_values.get('delta_median'))}，n={_fmt(leiden_values.get('delta_n_observations'))}）；Raw Level2 baseline 为 {_fmt(level2_values.get('baseline_median'))}→Recon {_fmt(level2_values.get('reconstruction_median'))}（delta={_fmt(level2_values.get('delta_median'))}，n={_fmt(level2_values.get('delta_n_observations'))}）。"
        )
        en_lines.append(
            f"{scope}: {label} at scale={_fmt(scale)} μm with valid windows={_fmt(windows)} is Raw Leiden baseline {_fmt(leiden_values.get('baseline_median'))} to Recon {_fmt(leiden_values.get('reconstruction_median'))} (delta={_fmt(leiden_values.get('delta_median'))}, n={_fmt(leiden_values.get('delta_n_observations'))}); Raw Level2 baseline {_fmt(level2_values.get('baseline_median'))} to Recon {_fmt(level2_values.get('reconstruction_median'))} (delta={_fmt(level2_values.get('delta_median'))}, n={_fmt(level2_values.get('delta_n_observations'))})."
        )
    if not zh_lines:
        zh_lines.append(_missing_line("所有 scope", label, zh=True))
        en_lines.append(_missing_line("all scopes", label, zh=False))
    return zh_lines, en_lines


def _anatomy_local_lines(record_list: list[Mapping[str, Any]]) -> tuple[list[str], list[str]]:
    grouped = _by_scope_questions(record_list, ("local_vs_raw_leiden", "local_vs_raw_level2"))
    zh_lines: list[str] = []
    en_lines: list[str] = []
    for scope in _scope_order(grouped):
        chunks_zh: list[str] = []
        chunks_en: list[str] = []
        for baseline, label in (("local_vs_raw_leiden", "Raw Leiden"), ("local_vs_raw_level2", "Raw Level2")):
            result = _results(grouped[scope].get(baseline))
            rows = result.get("anatomy_summary", [])
            if not isinstance(rows, list) or not rows:
                continue
            items_zh: list[str] = []
            items_en: list[str] = []
            for row in rows:
                if not isinstance(row, Mapping):
                    continue
                region = row.get("level1_region", "unknown")
                baseline_name = "raw_leiden" if baseline == "local_vs_raw_leiden" else "raw_level2"
                kobs_key = f"median_delta_k_obs_vs_{baseline_name}"
                neff_key = f"median_delta_neff_vs_{baseline_name}"
                evenness_key = f"median_delta_evenness_vs_{baseline_name}"
                items_zh.append(f"{region}: ΔKobs={_fmt(row.get(kobs_key))}, ΔNeff={_fmt(row.get(neff_key))}, Δevenness={_fmt(row.get(evenness_key))}, valid windows={_fmt(row.get('n_valid_windows'))}")
                items_en.append(f"{region}: ΔKobs={_fmt(row.get(kobs_key))}, ΔNeff={_fmt(row.get(neff_key))}, Δevenness={_fmt(row.get(evenness_key))}, valid windows={_fmt(row.get('n_valid_windows'))}")
            if items_zh:
                chunks_zh.append(f"{label}：" + "; ".join(items_zh))
                chunks_en.append(f"{label}: " + "; ".join(items_en))
        if chunks_zh:
            zh_lines.append(f"{scope}：" + "；".join(chunks_zh) + "。")
            en_lines.append(f"{scope}: " + "; ".join(chunks_en) + ".")
    if not zh_lines:
        zh_lines.append(_missing_line("所有 scope", "anatomy-stratified local results", zh=True))
        en_lines.append(_missing_line("all scopes", "anatomy-stratified local results", zh=False))
    return zh_lines, en_lines


def _threshold_lines(record_list: list[Mapping[str, Any]]) -> tuple[list[str], list[str]]:
    by_scope = _by_scope(record_list, "threshold_reliability")
    zh_lines: list[str] = []
    en_lines: list[str] = []
    for scope in _scope_order(by_scope):
        state = _results(by_scope[scope]).get("state", {})
        state = state if isinstance(state, Mapping) else {}
        zh_lines.append(f"{scope}：threshold status={_fmt(state.get('status'))}，analysis_status={_fmt(state.get('analysis_status'))}，threshold={_fmt(state.get('threshold'))}，CI=[{_fmt(state.get('ci_lower'))}, {_fmt(state.get('ci_upper'))}]，relative CI width={_fmt_percent(state.get('relative_ci_width'))}，valid bootstrap={_fmt(state.get('n_valid_bootstrap'))}/{_fmt(state.get('n_total_bootstrap'))}，n_windows={_fmt(state.get('n_windows'))}；可靠性先于区域范围。")
        en_lines.append(f"{scope}: threshold status={_fmt(state.get('status'))}, analysis_status={_fmt(state.get('analysis_status'))}, threshold={_fmt(state.get('threshold'))}, CI=[{_fmt(state.get('ci_lower'))}, {_fmt(state.get('ci_upper'))}], relative CI width={_fmt_percent(state.get('relative_ci_width'))}, valid bootstrap={_fmt(state.get('n_valid_bootstrap'))}/{_fmt(state.get('n_total_bootstrap'))}, n_windows={_fmt(state.get('n_windows'))}; reliability is assessed before extent.")
    if not zh_lines:
        zh_lines.append(_missing_line("所有 scope", "threshold_reliability", zh=True))
        en_lines.append(_missing_line("all scopes", "threshold_reliability", zh=False))
    return zh_lines, en_lines


def _region_extent_lines(record_list: list[Mapping[str, Any]], *, coverage: bool) -> tuple[list[str], list[str]]:
    by_scope = _by_scope(record_list, "region_extent")
    zh_lines: list[str] = []
    en_lines: list[str] = []
    for scope in _scope_order(by_scope):
        results = _results(by_scope[scope])
        overall = results.get("overall", {}) if isinstance(results.get("overall", {}), Mapping) else {}
        status = results.get("threshold_status", overall.get("threshold_status"))
        available = overall.get("region_available")
        if coverage:
            rows = results.get("by_region", [])
            if not isinstance(rows, list) or not rows:
                rows = [overall]
            region_text_zh: list[str] = []
            region_text_en: list[str] = []
            for row in rows:
                if not isinstance(row, Mapping):
                    continue
                region = row.get("level1_region", "Overall")
                region_text_zh.append(f"{region}: windows={_fmt(row.get('region_windows'))}/{_fmt(row.get('valid_windows'))}, area={_fmt(row.get('region_area_mm2'))} mm² ({_fmt_percent(row.get('region_area_fraction', row.get('area_fraction')))}), units={_fmt(row.get('region_units'))}/{_fmt(row.get('valid_units'))} ({_fmt_percent(row.get('region_unit_fraction', row.get('unit_fraction')))}), threshold={_fmt(row.get('threshold_status', status))}")
                region_text_en.append(f"{region}: windows={_fmt(row.get('region_windows'))}/{_fmt(row.get('valid_windows'))}, area={_fmt(row.get('region_area_mm2'))} mm² ({_fmt_percent(row.get('region_area_fraction', row.get('area_fraction')))}), units={_fmt(row.get('region_units'))}/{_fmt(row.get('valid_units'))} ({_fmt_percent(row.get('region_unit_fraction', row.get('unit_fraction')))}), threshold={_fmt(row.get('threshold_status', status))}")
            zh_lines.append(f"{scope}：threshold_status={_fmt(status)}，region_available={_fmt_bool(available)}；" + "；".join(region_text_zh) + "。State Region coverage 保留 area、window 和 unit 分母；它不是 EMT hotspot。")
            en_lines.append(f"{scope}: threshold_status={_fmt(status)}, region_available={_fmt_bool(available)}; " + "; ".join(region_text_en) + ". State Region coverage retains area, window, and unit denominators; it is not an EMT hotspot.")
        else:
            zh_lines.append(f"{scope}：State Region 由 Recon Neff 连续场的 threshold_status={_fmt(status)} 和可识别性决定，region_available={_fmt_bool(available)}，valid windows={_fmt(overall.get('valid_windows'))}，region windows={_fmt(overall.get('region_windows'))}；区域不可识别时范围保持 NA。State Region 不是 EMT hotspot，也不代表 EMT 富集。")
            en_lines.append(f"{scope}: the State Region is defined from the reconstructed Neff field with threshold_status={_fmt(status)} and identifiability region_available={_fmt_bool(available)}; valid windows={_fmt(overall.get('valid_windows'))}, region windows={_fmt(overall.get('region_windows'))}. Extent remains NA when the region is not identifiable. The State Region is not an EMT hotspot and does not represent EMT enrichment.")
    if not zh_lines:
        zh_lines.append(_missing_line("所有 scope", "region_extent", zh=True))
        en_lines.append(_missing_line("all scopes", "region_extent", zh=False))
    return zh_lines, en_lines


def _scale_lines(record_list: list[Mapping[str, Any]]) -> tuple[list[str], list[str]]:
    by_scope = _by_scope(record_list, "scale_sensitivity")
    zh_lines: list[str] = []
    en_lines: list[str] = []
    for scope in _scope_order(by_scope):
        results = _results(by_scope[scope])
        rows = results.get("rows", [])
        rows = rows if isinstance(rows, list) else []
        usable = [row for row in rows if isinstance(row, Mapping) and _number(row.get("window_side_length")) is not None]
        if usable:
            first = usable[0]
            last = usable[-1]
            first_scale = _fmt(first.get("window_side_length"))
            last_scale = _fmt(last.get("window_side_length"))
            first_leiden = first.get("median_delta_neff_vs_raw_leiden")
            last_leiden = last.get("median_delta_neff_vs_raw_leiden")
            first_level2 = first.get("median_delta_neff_vs_raw_level2")
            last_level2 = last.get("median_delta_neff_vs_raw_level2")
            leiden_values = [_number(row.get("median_delta_neff_vs_raw_leiden")) for row in usable]
            level2_values = [_number(row.get("median_delta_neff_vs_raw_level2")) for row in usable]
            leiden_values = [value for value in leiden_values if value is not None]
            level2_values = [value for value in level2_values if value is not None]
            leiden_direction = "一致为升高" if leiden_values and all(value > 0 for value in leiden_values) else "一致为降低" if leiden_values and all(value < 0 for value in leiden_values) else "方向混合或为零"
            level2_direction = "一致为升高" if level2_values and all(value > 0 for value in level2_values) else "一致为降低" if level2_values and all(value < 0 for value in level2_values) else "方向混合或为零"
            leiden_direction_en = "consistently increased" if leiden_values and all(value > 0 for value in leiden_values) else "consistently decreased" if leiden_values and all(value < 0 for value in leiden_values) else "mixed or zero"
            level2_direction_en = "consistently increased" if level2_values and all(value > 0 for value in level2_values) else "consistently decreased" if level2_values and all(value < 0 for value in level2_values) else "mixed or zero"
            zh_lines.append(f"{scope}：selected scale={_fmt(results.get('main_scale_um'))} μm；{first_scale}→{last_scale} μm 中，Neff vs Raw Leiden 的首末 delta={_fmt(first_leiden)}→{_fmt(last_leiden)}（{leiden_direction}），vs Raw Level2={_fmt(first_level2)}→{_fmt(last_level2)}（{level2_direction}）。")
            en_lines.append(f"{scope}: selected scale={_fmt(results.get('main_scale_um'))} μm; across {first_scale}→{last_scale} μm, Neff delta vs Raw Leiden changed from {_fmt(first_leiden)} to {_fmt(last_leiden)} ({leiden_direction_en}), while vs Raw Level2 changed from {_fmt(first_level2)} to {_fmt(last_level2)} ({level2_direction_en}).")
        else:
            zh_lines.append(f"{scope}：selected scale={_fmt(results.get('main_scale_um'))} μm；没有可用的尺度行，Neff 双 baseline 方向保持 NA。")
            en_lines.append(f"{scope}: selected scale={_fmt(results.get('main_scale_um'))} μm; no usable scale rows were saved, so Neff directions under both baselines remain NA.")
    if not zh_lines:
        zh_lines.append(_missing_line("所有 scope", "scale_sensitivity", zh=True))
        en_lines.append(_missing_line("all scopes", "scale_sensitivity", zh=False))
    return zh_lines, en_lines


def _build_node_text(node_id: str, record_list: list[Mapping[str, Any]]) -> tuple[str, list[str], str, list[str]]:
    """Return lead, bilingual interpretation lines, and next-step text."""

    if node_id == "2.1":
        zh, en = _complexity_lines(record_list)
        return "保存的 complexity 记录按 scope 给出实际 K、比较类型和 ARI，作为结构对照候选。", zh, "The saved complexity records retain the actual K, comparison kind, and ARI by scope as a structural reading candidate.", en
    if node_id == "2.2":
        zh, en = _matched_lines(record_list)
        return "先看 matched-K 是否成立，再把 assignment change 放在对应资格条件下阅读。", zh, "Read matched-K eligibility first, then place assignment change under its qualification condition.", en
    if node_id == "2.4":
        zh, en = _moran_lines(record_list)
        return "Moran 保留 all_valid 的两侧分布视角和 shared_valid 的逐基因配对视角。", zh, "Moran retains the side-specific all-valid view and the gene-wise paired shared-valid view.", en
    if node_id == "2.5":
        zh, en = _emt_lines(record_list)
        return "EMT 的 coverage、score 和 rank cutoff 按保存结果分开呈现，避免把覆盖差异写成生物学变化。", zh, "EMT coverage, scores, and rank cutoffs are shown separately from the saved results so coverage differences are not written as biological change.", en
    if node_id == "3.1":
        zh, en = _anatomy_lines(record_list)
        return "Anatomy 解读只整理已保存的区域名称和窗口支持背景。", zh, "The anatomy reading only organizes saved region labels and window-support context.", en
    if node_id == "3.2":
        zh, en = _changed_lines(record_list)
        return "Spatial assignment change 同时保留改变数量、改变比例和窗口支持分母。", zh, "The spatial assignment-change reading keeps changed counts, changed fractions, and window-support denominators together.", en
    if node_id == "3.3":
        zh, en = _emt_spatial_lines(record_list)
        return "EMT spatial 节点只从保存的字段 metadata 生成候选，不把它扩展为区域定量或富集结论。", zh, "The EMT spatial node generates a candidate from saved field metadata only and does not expand it into regional quantification or enrichment.", en
    if node_id == "3.4":
        zh, en = _window_lines(record_list)
        return "窗口支持候选把主尺度、有效窗口、保留 units、最小支持和 rarefaction 分母放在一起。", zh, "The window-support candidate keeps the selected scale, valid windows, retained units, minimum support, and rarefaction denominator together.", en
    if node_id == "3.5":
        zh, en = _local_metric_lines(record_list, "Kobs")
        return "Kobs 的候选解读分别对照 Raw Leiden 与 Raw Level2，并保留各自窗口分母。", zh, "The Kobs candidate is read separately against Raw Leiden and Raw Level2 with each window denominator retained.", en
    if node_id == "3.6":
        zh, en = _local_metric_lines(record_list, "Neff")
        return "Neff 的候选解读分别对照两个 Raw baseline；方向只作为有条件的保存结果描述。", zh, "The Neff candidate is read separately against both Raw baselines; direction is a conditioned description of saved results.", en
    if node_id == "3.7":
        zh, en = _local_metric_lines(record_list, "evenness")
        return "Evenness 只作组成辅助描述，不预设越高越好。", zh, "Evenness is used as a supporting composition description without assuming that higher is better.", en
    if node_id == "3.8":
        zh, en = _anatomy_local_lines(record_list)
        return "Anatomy 分层保留两条 Raw baseline 下的局部指标 delta 与有效窗口分母。", zh, "Anatomy stratification retains local metric deltas and valid-window denominators under both Raw baselines.", en
    if node_id == "4.1":
        zh, en = _threshold_lines(record_list)
        return "先从保存的 CI、bootstrap 和窗口数判断阈值可靠性，再决定区域范围是否可读。", zh, "Assess threshold reliability from saved CI, bootstrap, and window counts before reading region extent.", en
    if node_id == "4.2":
        zh, en = _region_extent_lines(record_list, coverage=False)
        return "State Region 候选来自 Recon Neff 连续场及其可识别性状态，与 EMT hotspot 语义分开。", zh, "The State Region candidate comes from the reconstructed Neff field and its identifiability state, separate from EMT-hotspot semantics.", en
    if node_id == "4.3":
        zh, en = _region_extent_lines(record_list, coverage=True)
        return "Region coverage 候选保留面积、窗口和 unit 三类分母，并沿用 threshold 状态。", zh, "The Region-coverage candidate retains area, window, and unit denominators together with threshold status.", en
    if node_id == "4.4":
        zh, en = _scale_lines(record_list)
        return "尺度敏感性候选逐个保留候选窗口及其相对两条 Raw baseline 的 Neff delta。", zh, "The scale-sensitivity candidate retains each candidate window and its Neff delta against both Raw baselines.", en
    raise ValueError(f"unsupported narrative node: {node_id}")


_NEXT = {
    "2.1": ("下一步检查 2.2 的 matched-K 资格，再读取 assignment change。", "Next, check matched-K eligibility in 2.2 before reading assignment change."),
    "2.2": ("随后在 2.4 分开阅读 Moran 的 all_valid 与 shared_valid 两种视角。", "Then read the two Moran views, all_valid and shared_valid, in 2.4."),
    "2.4": ("随后在 2.5 分开阅读 EMT coverage、cutoff 与 score。", "Then read EMT coverage, cutoffs, and scores separately in 2.5."),
    "2.5": ("随后进入 3.1–3.4 的 anatomy、assignment、EMT spatial 和窗口支持。", "Then move to anatomy, assignment, EMT spatial fields, and window support in 3.1–3.4."),
    "3.1": ("下一步在 3.2 定位 assignment change 的空间分布。", "Next, locate assignment change spatially in 3.2."),
    "3.2": ("下一步在 3.3 单独查看 EMT spatial field metadata。", "Next, inspect EMT spatial-field metadata separately in 3.3."),
    "3.3": ("下一步在 3.4 核对局部状态比较的窗口支持。", "Next, check window support for local-state comparison in 3.4."),
    "3.4": ("下一步在 3.5–3.7 分别阅读 Kobs、Neff 和 evenness。", "Next, read Kobs, Neff, and evenness separately in 3.5–3.7."),
    "3.5": ("下一步在 3.6 查看 Neff 的两个 baseline。", "Next, inspect Neff under the two baselines in 3.6."),
    "3.6": ("下一步在 3.7 用 evenness 补充组成描述。", "Next, use evenness for supporting composition description in 3.7."),
    "3.7": ("下一步在 3.8 查看 anatomy 分层结果。", "Next, inspect anatomy-stratified results in 3.8."),
    "3.8": ("下一步进入 4.1，先判断 State Region 阈值可靠性。", "Next, enter 4.1 and assess State Region threshold reliability first."),
    "4.1": ("下一步在 4.2 区分连续 State、可识别 mask 和 NA。", "Next, separate continuous state, identifiable mask, and NA in 4.2."),
    "4.2": ("下一步在 4.3 报告 Region 的面积、窗口和 unit 覆盖。", "Next, report Region area, window, and unit coverage in 4.3."),
    "4.3": ("下一步在 4.4 检查候选尺度敏感性。", "Next, check candidate-scale sensitivity in 4.4."),
    "4.4": ("最后回到现有 4.5 与层末 summary，汇总已审阅的条件化结果。", "Finally, return to the existing 4.5 and layer-end summaries for reviewed, conditioned results."),
}


def build_node_narratives(records: Mapping[str, Any] | Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Build bilingual candidate narratives from an existing record package.

    The returned list is intentionally detached from the input records.  It
    contains no scientific values beyond those already stored in ``results``
    and starts with a pending review marker; ``save_records`` attaches the
    package digest when persisting it.
    """

    record_list = _main_records(_record_list(records))
    narratives: list[dict[str, Any]] = []
    for node_id in NARRATIVE_NODE_IDS:
        spec = _NARRATIVE_SPECS[node_id]
        selected = _records_for_questions(record_list, spec.get("questions", ()))
        lead_zh, interpretations_zh, lead_en, interpretations_en = _build_node_text(node_id, record_list)
        next_zh, next_en = _NEXT[node_id]
        narratives.append({
            "node_id": node_id,
            "title": {"zh": str(spec["title_zh"]), "en": str(spec["title_en"])},
            "lead": {"zh": lead_zh, "en": lead_en},
            "interpretations": {"zh": list(interpretations_zh), "en": list(interpretations_en)},
            "next": {"zh": next_zh, "en": next_en},
            "record_ids": _record_ids(selected),
            "review": {"status": "pending"},
        })
    return narratives


__all__ = ["NARRATIVE_NODE_IDS", "build_node_narratives"]
