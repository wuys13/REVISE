"""Render saved reconstruction-impact conclusions as a static HTML report.

The renderer is deliberately a reader. It consumes the structured records,
CSV references, and figures already published by the route notebooks; it does
not load an AnnData object or calculate a scientific metric.
"""
from __future__ import annotations

import csv
import hashlib
import html
import json
import math
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping
from urllib.parse import quote

from .content_contract import LAYERS, NODE_SPECS, QUESTION_LABELS, QUESTION_LAYERS, SCOPES

TITLES = {
    "foundation": ("比较基础", "Comparison foundation"),
    "overall": ("总体变化", "Overall changes"),
    "localization": ("空间定位与局部状态", "Spatial localization and local states"),
    "region": ("Region 与综合解释", "Region and interpretation"),
}
LAYER_INTROS = {
    "foundation": "先固定比较载体、cohort、坐标和基因空间，后续结论只在这些条件下解释。",
    "overall": "在比较基础成立后，分别查看结构、空间自相关和 EMT 条件；各问题保留自己的分母。",
    "localization": "再把变化放回空间位置与局部状态，区分 assignment 定位和状态方向。",
    "region": "最后先判断 State Region 阈值是否可靠，再解释可识别范围；不能识别时保持 NA。",
}
GROUPED_QUESTIONS = {
    "foundation",
    "complexity",
    "matched_k",
    "moran_all_valid",
    "moran_shared_valid",
    "emt_coverage",
    "emt_score",
    "threshold_reliability",
}
FIGURE_ROLE_BY_NAME = {
    "matched_k_contingency.png": "matched_k",
    "moran_all_gene_distribution.png": "moran_all_valid",
    "moran_q75_heatmap.png": "moran_all_valid",
    "moran_shared_gene_delta_distribution.png": "moran_shared_valid",
    "moran_shared_gene_scatter.png": "moran_shared_valid",
    "threshold_reliability.png": "threshold_reliability",
    "anatomy_regions.png": "anatomy",
    "anatomy_window_decision.png": "anatomy",
    "changed_units.png": "changed_units",
    "level1_change_fraction.png": "changed_units",
}
FIGURE_TITLE_BY_NAME = {
    "matched_k_contingency.png": "Matched-K assignment contingency",
    "moran_all_gene_distribution.png": "Moran all-valid/shared-valid distribution comparison",
    "moran_q75_heatmap.png": "Moran all-valid/shared-valid Q75 comparison",
    "moran_shared_gene_delta_distribution.png": "Moran shared-valid paired-delta distribution",
    "moran_shared_gene_scatter.png": "Moran shared-valid paired-delta scatter",
    "threshold_reliability.png": "Threshold reliability and bootstrap support",
    "anatomy_regions.png": "Anatomy context and spatial localization",
    "anatomy_window_decision.png": "Anatomy context and window decision",
    "changed_units.png": "Changed-unit spatial localization",
    "level1_change_fraction.png": "Changed-unit fraction by region",
}
FIGURE_TITLE_BY_ROLE = {
    "complexity": "Reconstructed clustering complexity",
    "matched_k": "Matched-K assignment contingency",
    "moran_all_valid": "Moran all-valid/shared-valid comparison",
    "moran_shared_valid": "Moran shared-valid paired-delta comparison",
    "threshold_reliability": "Threshold reliability and bootstrap support",
    "region_extent": "Reconstructed Neff field and State Region",
    "window_support": "Window support and selected scale",
    "anatomy": "Anatomy context and spatial localization",
    "changed_units": "Changed-unit spatial localization",
    "emt_spatial_fields": "EMT spatial field and localization",
    "emt_score": "EMT AUCell distribution and paired delta",
    "local_vs_raw_leiden": "Local diversity evidence matrix",
    "local_vs_raw_level2": "Local diversity evidence matrix",
}


JUDGMENTS = {
    "foundation available": "比较基础已登记",
    "foundation incomplete": "比较基础不完整",
    "complexity diagnostic available": "复杂度诊断",
    "complexity diagnostic unavailable": "复杂度诊断不可用",
    "matched-K eligible": "可作 matched-K 比较",
    "matched-K unavailable": "matched-K 未成立",
    "side-specific distributions": "各自基因集合的分布对照",
    "paired directional shift": "共同基因的配对变化",
    "unstable / extent NA": "无稳定阈值；覆盖 NA",
    "not computable / NA": "阈值不可计算",
    "side-specific distributions; marginal comparison only": "各自基因集合的分布对照",
    "descriptive directional shift": "描述性方向变化",
    "descriptive score shift": "描述性评分变化",
    "coverage differs; comparability limited": "覆盖不同；可比性受限",
    "anatomy context available": "组织背景可用",
    "descriptive changed-unit fraction": "归属变化定位",
    "EMT spatial field available": "EMT 空间场可用",
    "window support available": "窗口支持已登记",
    "supports concentrated state direction": "支持局部状态更集中方向",
    "supports finer state direction": "支持更细内部状态方向",
    "does not support requested direction": "未支持该 baseline 的方向",
    "no directional support": "中位差为零",
    "direction unknown": "方向不可判定",
    "reliable": "可稳定识别",
    "unreliable / NA": "无稳定阈值",
    "extent available": "覆盖可计算",
    "extent NA": "覆盖为 NA",
    "conditional sensitivity": "按尺度检查方向",
    "sensitivity unavailable": "尺度敏感性不可用",
}
CONDITIONS = {
    "valid": "证据可用",
    "ok": "条件满足",
    "available": "结果可用",
    "review": "条件需审阅",
    "missing": "证据缺失",
    "invalid": "证据无效",
    "unmatched_cluster_complexity": "受 cluster complexity 限制",
}


def _load(package: Any) -> dict[str, Any]:
    if isinstance(package, (str, Path)):
        package = json.loads(Path(package).read_text(encoding="utf-8"))
    if isinstance(package, Mapping):
        result = dict(package)
        result["records"] = list(result.get("records", []))
        return result
    if isinstance(package, list):
        records = list(package)
        sample_id = str(records[0].get("identity", {}).get("sample_id", "unknown")) if records else "unknown"
        route_kind = str(records[0].get("identity", {}).get("route_kind", "unknown")) if records else "unknown"
        return {"sample_id": sample_id, "route_kind": route_kind, "records": records, "run_identity": {}}
    raise TypeError("package must be a saved report package, path, or record list")


def _labels() -> Mapping[str, tuple[str, str]]:
    return QUESTION_LABELS


def _number(value: Any, percent: bool = False) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, (float, int)):
        if not math.isfinite(value):
            return "NA"
        if percent:
            return f"{float(value):.1%}"
        return f"{float(value):,.0f}" if float(value).is_integer() else f"{float(value):.4g}"
    return str(value)


def _key_values(record: Mapping[str, Any]) -> str:
    """Select compact display facts while retaining the structured record."""

    results = record.get("results", {}) or {}
    question = record.get("identity", {}).get("question")
    pairs: list[tuple[str, Any]] = []
    if question == "foundation":
        pairs = [("Input n", results.get("input_units")), ("Partition n", results.get("paired_units")), ("Seed", results.get("seed"))]
    elif question == "complexity":
        fixed = results.get("comparison_kind") == "raw_reference_vs_fixed_final_clusters" or results.get("fixed_final_cluster_k") is not None
        pairs = [("Raw K", results.get("raw_k")), ("Fixed final K" if fixed else "Recon K", results.get("fixed_final_cluster_k") if fixed else results.get("reconstruction_k_at_raw_resolution"))]
    elif question == "matched_k":
        pairs = [("Raw / Recon K", f"{_number(results.get('raw_k'))} / {_number(results.get('reconstruction_k'))}"), ("Change", _number(results.get("unit_change_fraction"), True)), ("ARI", results.get("ari"))]
    elif str(question).startswith("moran_"):
        for side, label in (("raw", "Raw"), ("reconstruction", "Recon")):
            side_results = results.get(side, {}) or {}
            pairs.extend(((label + " Q75", side_results.get("q75")), (label + " n", side_results.get("n_valid"))))
        if question == "moran_shared_valid":
            pairs.append(("Paired median ΔI", (results.get("paired_delta", {}) or {}).get("median")))
    elif question == "emt_coverage":
        for side, label in (("raw", "Raw"), ("reconstruction", "Recon")):
            side_results = results.get(side, {}) or {}
            pairs.append((label, f"{_number(side_results.get('available_gene_count'))}/{_number(side_results.get('resource_gene_count'))} ({_number(side_results.get('coverage'), True)})"))
    elif question == "emt_score":
        pairs = [("Raw median", (results.get("raw", {}) or {}).get("median")), ("Recon median", (results.get("reconstruction", {}) or {}).get("median")), ("Paired median Δ", (results.get("paired", {}) or {}).get("median_delta")), ("Paired n", (results.get("paired", {}) or {}).get("n"))]
        for side, label in (("raw", "Raw"), ("reconstruction", "Recon")):
            cutoff = (results.get("cutoff_audit", {}) or {}).get(side, {}) or {}
            pairs.append((label + " rank length", cutoff.get("effective_rank_length")))
    elif str(question).startswith("local_vs_"):
        for metric, value in (results.get("metrics", {}) or {}).items():
            if not isinstance(value, Mapping):
                continue
            delta = value.get("delta_median")
            if delta is None:
                delta = (value.get("raw_leiden", {}) or {}).get("delta_median")
            pairs.append((f"{metric} median Δ", delta))
        pairs += [("Scale (μm)", results.get("scale_um")), ("Valid windows", results.get("n_valid_windows"))]
    elif question == "threshold_reliability":
        state = results.get("state", {}) or {}
        pairs = [("CI width / Neff range", _number(state.get("relative_ci_width"), True)), ("Valid bootstrap", _number(state.get("valid_bootstrap_fraction"), True)), ("Threshold", state.get("threshold")), ("Windows", state.get("n_windows"))]
    elif question == "region_extent":
        overall = results.get("overall", {}) or {}
        pairs = [("Area (mm²)", overall.get("region_area_mm2")), ("Area fraction", _number(overall.get("area_fraction"), True)), ("Unit fraction", _number(overall.get("unit_fraction"), True)), ("Valid units", overall.get("valid_units")), ("Valid windows", overall.get("valid_windows"))]
    elif question == "window_support":
        pairs = [("Scale (μm)", results.get("main_window_side_um")), ("Valid windows", results.get("n_valid_windows")), ("Minimum units", results.get("min_parent_units")), ("Paired draws", results.get("rarefaction_draws"))]
    elif question == "changed_units":
        overall = results.get("overall") or next((item for item in results.get("by_region", []) if item.get("level1_region") == "Overall"), {})
        denominator = overall.get("paired_units")
        if denominator is None:
            denominator = overall.get("total_units")
        pairs = [("Change", _number(overall.get("change_fraction"), True)), ("Changed / paired", f"{_number(overall.get('changed_units'))}/{_number(denominator)}")]
    elif question == "emt_spatial_fields":
        pairs = [("Units", results.get("n_units")), ("Coordinates", ", ".join(results.get("coordinate_units", [])))]
    elif question == "anatomy":
        pairs = [(item.get("level1_region", ""), _number(item.get("area_fraction"), True)) for item in results.get("regions", [])]
    elif question == "scale_sensitivity":
        pairs = [("Main scale (μm)", results.get("main_scale_um")), ("Scales checked", len(results.get("rows", [])))]
    return "<br>".join(f'<span class="metric">{html.escape(str(label))}: <strong>{html.escape(_number(value))}</strong></span>' for label, value in pairs)


def _headline_metric(record: Mapping[str, Any]) -> tuple[str, str]:
    """Return one saved fact for the overview, never a newly computed value."""

    identity = record.get("identity", {}) or {}
    results = record.get("results", {}) or {}
    question = identity.get("question")
    if question == "foundation":
        return "Input units", _number(results.get("input_units"))
    if question == "complexity":
        raw_k = _number(results.get("raw_k"))
        fixed_k = results.get("fixed_final_cluster_k")
        if fixed_k is not None:
            return "Raw → fixed K", f"{raw_k} → {_number(fixed_k)}"
        return "Raw → Recon K", f"{raw_k} → {_number(results.get('reconstruction_k_at_raw_resolution'))}"
    if question == "matched_k":
        if results.get("headline_eligible") is False or results.get("matched") is False:
            return "Headline", "matched-K unavailable"
        return "Change", _number(results.get("unit_change_fraction"), True)
    if question == "moran_all_valid":
        raw_q75 = (results.get("raw", {}) or {}).get("q75")
        recon_q75 = (results.get("reconstruction", {}) or {}).get("q75")
        return "Raw → Recon Q75", f"{_number(raw_q75)} → {_number(recon_q75)}"
    if question == "moran_shared_valid":
        return "Paired median ΔI", _number((results.get("paired_delta", {}) or {}).get("median"))
    if question == "emt_coverage":
        raw = results.get("raw", {}) or {}
        recon = results.get("reconstruction", {}) or {}
        raw_hits = f"{_number(raw.get('available_gene_count'))}/{_number(raw.get('resource_gene_count'))}"
        recon_hits = f"{_number(recon.get('available_gene_count'))}/{_number(recon.get('resource_gene_count'))}"
        return "Raw → Recon hits / total", f"{raw_hits} → {recon_hits}"
    if question == "emt_score":
        return "Paired median Δ", _number((results.get("paired", {}) or {}).get("median_delta"))
    if question == "anatomy":
        return "Regions", _number(len(results.get("regions", [])))
    if question == "changed_units":
        overall = results.get("overall") or next((item for item in results.get("by_region", []) if item.get("level1_region") == "Overall"), {})
        return "Change", _number(overall.get("change_fraction"), True)
    if question == "emt_spatial_fields":
        return "Spatial units", _number(results.get("n_units"))
    if question == "window_support":
        return "Chosen scale (μm)", _number(results.get("main_window_side_um"))
    if str(question).startswith("local_vs_"):
        metrics = results.get("metrics", {}) or {}
        neff = metrics.get("Neff", {}) or {}
        delta = neff.get("delta_median")
        if delta is None:
            baseline = results.get("baseline", "raw_leiden")
            delta = (neff.get(baseline, {}) or {}).get("delta_median")
        return "Neff median Δ", _number(delta)
    if question == "threshold_reliability":
        return "CI width / range", _number((results.get("state", {}) or {}).get("relative_ci_width"), True)
    if question == "region_extent":
        return "Area fraction", _number((results.get("overall", {}) or {}).get("area_fraction"), True)
    if question == "scale_sensitivity":
        rows = [row for row in results.get("rows", []) if isinstance(row, Mapping)]
        leiden = [row.get("median_delta_neff_vs_raw_leiden") for row in rows if row.get("median_delta_neff_vs_raw_leiden") is not None]
        level2 = [row.get("median_delta_neff_vs_raw_level2") for row in rows if row.get("median_delta_neff_vs_raw_level2") is not None]
        def span(values: list[Any]) -> str:
            return "NA" if not values else f"{_number(min(values))}…{_number(max(values))}"
        return "Neff Δ ranges", f"Leiden {span(leiden)}; Level2 {span(level2)}"
    return "Saved metric", _number(results.get("headline"))


def _display_scope(scope: Any) -> str:
    return str(scope).replace("Mono_Macro", "Mono/Macro")


def _judgment(record: Mapping[str, Any], language: str = "zh") -> str:
    judgment = str(record.get("performance_judgment", "NA"))
    condition = record.get("evidence_condition", {}) or {}
    status = str(condition.get("status", "pending"))
    if language == "zh":
        judgment = JUDGMENTS.get(judgment, judgment)
        status = CONDITIONS.get(status, status)
        if condition.get("matched_k_status") == "unmatched_cluster_complexity":
            status += "；matched-K 未成立，依赖其控制的解释受限"
    elif condition.get("matched_k_status") == "unmatched_cluster_complexity":
        status += "; matched-K is unavailable; dependent interpretation is limited"
    return f"{judgment} · {status}"


def _overview_conclusion(record: Mapping[str, Any], language: str = "zh") -> str:
    """Keep the overview concise by using saved judgment and condition only."""

    return _judgment(record, language)


def _brief_conclusion(record: Mapping[str, Any], language: str = "en", limit: int = 220) -> str:
    text = str((record.get("conclusion", {}) or {}).get(language, "") or "")
    if not text:
        return _overview_conclusion(record, language)
    text = re.split(r"(?<=[.!?。！？])\s+", text, maxsplit=1)[0]
    return text if len(text) <= limit else text[: limit - 1].rstrip() + "…"


def _figure_id(path: str) -> str:
    return "figure-" + hashlib.sha256(path.encode("utf-8")).hexdigest()[:12]


def _slug(value: Any) -> str:
    return re.sub(r"[^a-z0-9_-]+", "-", str(value).lower()).strip("-") or "item"


def _module_id(layer: str, question: str, all_scope: bool = False) -> str:
    return f"module-{_slug(question)}{'-all' if all_scope else ''}"


def _row_id(module_id: str, scope: str) -> str:
    return f"row-{module_id}-{_slug(scope)}"


def _metric_for_figure(path: str) -> str | None:
    name = Path(path).name.lower()
    for metric in ("kobs", "neff", "evenness"):
        if metric in name:
            return metric
    return None


def _record_layer(record: Mapping[str, Any]) -> str:
    identity = record.get("identity", {}) or {}
    question = str(identity.get("question", "unknown"))
    return str(QUESTION_LAYERS.get(question, "unknown"))


def _question_order(records: Iterable[Mapping[str, Any]]) -> list[str]:
    present = {str((record.get("identity", {}) or {}).get("question", "")) for record in records}
    labels = _labels()
    known = [question for question in labels if question in present]
    return known + sorted(present.difference(known))


def _record_map(records: Iterable[Mapping[str, Any]]) -> dict[tuple[str, str], Mapping[str, Any]]:
    return {(str((record.get("identity", {}) or {}).get("question", "")), str((record.get("identity", {}) or {}).get("scope", ""))): record for record in records}


def _output_root(package: Mapping[str, Any], destination: Path) -> Path:
    run_identity = package.get("run_identity", {}) or {}
    root = run_identity.get("output_root") if isinstance(run_identity, Mapping) else None
    return Path(str(root)).expanduser() if root else destination.parent


def _href(path: str, output_root: Path, destination_parent: Path) -> str:
    source = Path(path)
    if not source.is_absolute():
        source = output_root / source
    try:
        relative = os.path.relpath(source, destination_parent)
    except ValueError:
        relative = str(source)
    return quote(relative.replace(os.sep, "/"), safe="/:@?=&%")


def _repository_root(destination_parent: Path) -> Path:
    for candidate in (destination_parent, *destination_parent.parents):
        if (candidate / "docs" / "design" / "reconstruction-impact").is_dir():
            return candidate
    return Path(__file__).resolve().parents[3]


def _node_for_question(question: str) -> Mapping[str, Any] | None:
    return next((node for node in NODE_SPECS if question in node["questions"]), None)


FORMAL_NARRATIVE_REVIEW_STATUSES = {"reviewed", "approved"}


def _node_narrative(package: Mapping[str, Any], node_id: str) -> Mapping[str, Any] | None:
    """Return the package-owned narrative for a node, when one is present.

    Narrative text is deliberately kept separate from ``records``.  The
    package builder owns the text and its review metadata; the renderer only
    selects and escapes it.
    """

    narratives = package.get("node_narratives", [])
    if not isinstance(narratives, list):
        return None
    return next(
        (
            item
            for item in narratives
            if isinstance(item, Mapping) and str(item.get("node_id", "")) == str(node_id)
        ),
        None,
    )


def _narrative_reviewed(package: Mapping[str, Any], narrative: Mapping[str, Any] | None) -> bool:
    """Whether a node narrative is allowed into the formal reader output."""

    if narrative is None:
        return False
    review = narrative.get("review", {})
    if not isinstance(review, Mapping):
        return False
    status = str(review.get("status", "pending")).strip().lower()
    if status not in FORMAL_NARRATIVE_REVIEW_STATUSES:
        return False
    package_digest = package.get("record_digest")
    narrative_digest = review.get("record_digest")
    if package_digest:
        package_review = package.get("review", {})
        if not isinstance(package_review, Mapping):
            return False
        if str(package_review.get("status", "pending")).strip().lower() not in FORMAL_NARRATIVE_REVIEW_STATUSES:
            return False
        if not narrative_digest or str(narrative_digest) != str(package_digest):
            return False
    return True


def _narrative_text(narrative: Mapping[str, Any], key: str, language: str) -> str:
    value = narrative.get(key, {})
    if isinstance(value, Mapping):
        value = value.get(language, "")
    return str(value or "").strip()


def _narrative_interpretations(narrative: Mapping[str, Any], language: str) -> list[str]:
    value = narrative.get("interpretations", {})
    if isinstance(value, Mapping):
        value = value.get(language, [])
    if isinstance(value, str):
        value = [value]
    if not isinstance(value, (list, tuple)):
        return []
    return [str(item).strip() for item in value if str(item or "").strip()]


def _node_title(package: Mapping[str, Any], node: Mapping[str, Any], language: str = "zh", *, use_narrative: bool = True) -> str:
    """Use a reviewed package narrative title, otherwise the static contract title."""

    title_key = "title_zh" if language == "zh" else "title_en"
    if use_narrative:
        narrative = _node_narrative(package, str(node["node_id"]))
        if _narrative_reviewed(package, narrative):
            title = _narrative_text(narrative or {}, "title", language)
            if title:
                return title
    return str(node.get(title_key, node.get("title_zh", node["node_id"])))


def _narrative_record_attrs(narrative: Mapping[str, Any]) -> str:
    record_ids = narrative.get("record_ids", [])
    if isinstance(record_ids, (list, tuple)):
        value = ",".join(str(item) for item in record_ids if item)
    else:
        value = str(record_ids or "")
    return f' data-record-ids="{html.escape(value, quote=True)}"' if value else ""


def _narrative_lead_panel(package: Mapping[str, Any], node: Mapping[str, Any], language: str = "zh") -> str:
    """Render the reviewed lead immediately after a node heading."""

    narrative = _node_narrative(package, str(node["node_id"]))
    if narrative is None:
        return ""
    node_id = html.escape(str(node["node_id"]), quote=True)
    if not _narrative_reviewed(package, narrative):
        return f'<aside class="node-narrative pending" data-node-id="{node_id}"><span class="narrative-review">Review: pending</span></aside>'
    lead = _narrative_text(narrative, "lead", language)
    if not lead:
        return ""
    attrs = _narrative_record_attrs(narrative)
    return f'<aside class="node-narrative reviewed narrative-lead-panel" data-node-id="{node_id}"{attrs}><p class="narrative-lead">{html.escape(lead)}</p></aside>'


def _narrative_tail_panel(package: Mapping[str, Any], node: Mapping[str, Any], language: str = "zh") -> str:
    """Render reviewed interpretations and next step after the evidence."""

    narrative = _node_narrative(package, str(node["node_id"]))
    if narrative is None or not _narrative_reviewed(package, narrative):
        return ""
    interpretations = _narrative_interpretations(narrative, language)
    next_text = _narrative_text(narrative, "next", language)
    if not interpretations and not next_text:
        return ""
    node_id = html.escape(str(node["node_id"]), quote=True)
    attrs = _narrative_record_attrs(narrative)
    parts = [f'<aside class="node-narrative reviewed narrative-tail-panel" data-node-id="{node_id}"{attrs}>']
    if interpretations:
        parts.append('<ul class="narrative-interpretations">')
        parts.extend(f'<li>{html.escape(item)}</li>' for item in interpretations)
        parts.append("</ul>")
    if next_text:
        parts.append(f'<p class="narrative-next"><strong>Next</strong> · {html.escape(next_text)}</p>')
    parts.append("</aside>")
    return "".join(parts)


def _narrative_panel(package: Mapping[str, Any], node: Mapping[str, Any], language: str = "zh", phase: str = "full") -> str:
    """Render narrative phases without repeating the node title.

    The article owns the dynamic heading.  ``phase`` lets the renderer place
    the lead before evidence and the interpretation/next block after it;
    ``full`` remains a small compatibility wrapper for callers outside the
    article renderer.
    """

    if phase == "lead":
        return _narrative_lead_panel(package, node, language)
    if phase == "tail":
        return _narrative_tail_panel(package, node, language)
    return _narrative_lead_panel(package, node, language) + _narrative_tail_panel(package, node, language)


def _notebook_narrative(package: Mapping[str, Any], node: Mapping[str, Any], language: str = "en") -> str:
    """Return a compact reviewed narrative block for executed notebooks."""

    narrative = _node_narrative(package, str(node["node_id"]))
    if narrative is None:
        return ""
    if not _narrative_reviewed(package, narrative):
        return '<p class="node-narrative pending">Review: pending</p>'
    title = _narrative_text(narrative, "title", language) or _node_title(package, node, language)
    lead = _narrative_text(narrative, "lead", language)
    interpretations = _narrative_interpretations(narrative, language)
    next_text = _narrative_text(narrative, "next", language)
    parts = [f'<section class="notebook-narrative"><h3>{html.escape(title)}</h3>']
    if lead:
        parts.append(f'<p>{html.escape(lead)}</p>')
    if interpretations:
        parts.append("<ul>" + "".join(f"<li>{html.escape(item)}</li>" for item in interpretations) + "</ul>")
    if next_text:
        parts.append(f'<p><strong>Next</strong> · {html.escape(next_text)}</p>')
    return "".join(parts) + "</section>"


def _method_link(question: str, destination_parent: Path) -> str:
    node = _node_for_question(question)
    doc, anchor = node["method"] if node is not None else ("outputs-and-test-plan.md", "")
    target = _repository_root(destination_parent) / "docs" / "design" / "reconstruction-impact" / doc
    href = _href(str(target), _repository_root(destination_parent), destination_parent) + "#" + anchor
    return f'<a class="method-link" href="{html.escape(href, quote=True)}">方法与判读规范</a>'


def _tables(record: Mapping[str, Any]) -> list[tuple[str, bool]]:
    evidence = record.get("evidence", {}) or {}
    raw = evidence.get("tables", []) if isinstance(evidence, Mapping) else []
    result = []
    for item in raw:
        if isinstance(item, Mapping):
            path = item.get("path")
            exists = item.get("exists", True)
        else:
            path = item
            exists = True
        if path:
            result.append((str(path), bool(exists)))
    return result


def _figure_paths(record: Mapping[str, Any]) -> list[str]:
    evidence = record.get("evidence", {}) or {}
    raw = evidence.get("figures", []) if isinstance(evidence, Mapping) else []
    return [str(path) for path in raw if path]


def _preferred_figure_question(path: str, questions: set[str]) -> str | None:
    name = Path(path).name.lower()
    preferred = FIGURE_ROLE_BY_NAME.get(name)
    scope_slugs = {"fibroblast", "mono_macro", "t"}
    for scope in scope_slugs:
        if name == f"{scope}_reconstructed_clusters.png":
            preferred = "complexity"
        elif name == f"{scope}_matched_clusters.png":
            preferred = "matched_k"
        elif name == f"{scope}_hallmark_epithelial_mesenchymal_transition_spatial.png":
            preferred = "emt_spatial_fields"
        elif name == f"{scope}_window_decision.png":
            preferred = "window_support"
        elif name == f"{scope}_state_region.png":
            preferred = "region_extent"
        elif name == f"{scope}_kobs_matrix.png" or name == f"{scope}_neff_matrix.png" or name == f"{scope}_evenness_matrix.png":
            preferred = "local_vs_raw_leiden"
        elif name == f"aucell_{scope}_distribution_and_delta.png":
            preferred = "emt_score"
    if name == "aucell_overall_distribution_and_delta.png":
        preferred = "emt_score"
    return preferred if preferred in questions else None


def _is_common_figure(path: str, scopes: set[str]) -> bool:
    parent_scopes = scopes.difference({"All"})
    if len(parent_scopes) < 2:
        return False
    name = Path(path).name.lower()
    scope_slugs = {"fibroblast", "mono_macro", "t"}
    scoped_suffixes = {
        "reconstructed_clusters.png",
        "matched_clusters.png",
        "hallmark_epithelial_mesenchymal_transition_spatial.png",
        "window_decision.png",
        "state_region.png",
        "kobs_matrix.png",
        "neff_matrix.png",
        "evenness_matrix.png",
    }
    if any(name == f"{scope}_{suffix}" for scope in scope_slugs for suffix in scoped_suffixes):
        return False
    if any(name == f"aucell_{scope}_distribution_and_delta.png" for scope in scope_slugs):
        return False
    return True


def _figure_title(path: str, role: str) -> str:
    """Return the finite reader-facing title used by captions and shared links."""

    name = Path(path).name.lower()
    if name.endswith("_reconstructed_clusters.png"):
        return "Reconstructed partition used by downstream local analysis; not a same-resolution K comparison"
    if name in FIGURE_TITLE_BY_NAME:
        return FIGURE_TITLE_BY_NAME[name]
    if role.startswith("local_vs_"):
        metric = _metric_for_figure(path)
        if metric:
            return f"Local {metric.capitalize()} 2×3 evidence matrix"
    return FIGURE_TITLE_BY_ROLE.get(role, "Saved evidence figure")


def _build_figure_registry(records: Iterable[Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
    occurrences: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for record in records:
        for path in _figure_paths(record):
            occurrences[path].append(record)
    registry: dict[str, dict[str, Any]] = {}
    for path, items in occurrences.items():
        questions = {str((item.get("identity", {}) or {}).get("question", "")) for item in items}
        preferred = _preferred_figure_question(path, questions)
        owner_items = [item for item in items if (item.get("identity", {}) or {}).get("question") == preferred] if preferred else []
        if not owner_items:
            owner_items = items
        # When the same figure is registered for HD All and parent scopes,
        # the primary parent batch row owns the visual. HD All keeps the
        # numeric record and links back to that shared visual.
        parent_items = [item for item in owner_items if (item.get("identity", {}) or {}).get("scope") != "All"]
        owner = parent_items[0] if parent_items else owner_items[0]
        owner_identity = owner.get("identity", {}) or {}
        owner_question = str(owner_identity.get("question", "unknown"))
        owner_layer = _record_layer(owner)
        scopes = {str((item.get("identity", {}) or {}).get("scope", "")) for item in items}
        common = _is_common_figure(path, scopes)
        all_scope = str(owner_identity.get("scope")) == "All"
        module = _module_id(owner_layer, owner_question, all_scope=all_scope)
        if owner_question in GROUPED_QUESTIONS:
            owner_row = _row_id(module, "batch")
        elif common:
            owner_candidates = [item for item in owner_items if (item.get("identity", {}) or {}).get("scope") != "All"]
            row_scope = str((owner_candidates[0] if owner_candidates else owner).get("identity", {}).get("scope", "All"))
            owner_row = _row_id(module, row_scope)
        else:
            owner_row = _row_id(module, str(owner_identity.get("scope", "All")))
            metric = _metric_for_figure(path) if owner_question.startswith("local_vs_") else None
            if metric:
                module = f"module-local-diversity-{_slug(metric)}"
                owner_row = f"row-{module}-{_slug(str(owner_identity.get('scope', 'All')))}"
        figure_role = preferred or owner_question
        registry[path] = {
            "path": path,
            "figure_id": _figure_id(path),
            "owner_module": module,
            "owner_row": owner_row,
            "common": common,
            "scopes": tuple(sorted(scopes)),
            "title": _figure_title(path, figure_role),
        }
    return registry


def _caption(path: str, info: Mapping[str, Any]) -> str:
    name = str(info.get("title") or "Saved evidence figure")
    scopes = ", ".join(_display_scope(scope) for scope in info.get("scopes", ()) if scope)
    if info.get("common") and scopes:
        suffix = f" · shared scopes: {scopes}"
        if "All" in info.get("scopes", ()):
            suffix += " · HD All uses a separate denominator"
    elif scopes:
        suffix = f" · scope: {scopes}"
    else:
        suffix = ""
    return f"{name}{suffix}"


def _render_figure(path: str, info: Mapping[str, Any], output_root: Path, destination_parent: Path) -> str:
    src = _href(path, output_root, destination_parent)
    figure_id = info["figure_id"]
    figure_class = "evidence-figure shared" if info.get("common") else "evidence-figure"
    return f'<figure class="{figure_class}" id="{figure_id}"><a href="{html.escape(src, quote=True)}"><img loading="lazy" src="{html.escape(src, quote=True)}" alt="{html.escape(Path(path).stem)}"></a><figcaption>{html.escape(_caption(path, info))} · <a href="{html.escape(src, quote=True)}">Original figure</a></figcaption></figure>'


def _figure_cell_records(records: Iterable[Mapping[str, Any]], row_id: str, registry: Mapping[str, Mapping[str, Any]], output_root: Path, destination_parent: Path) -> str:
    rendered: list[str] = []
    seen: set[str] = set()
    for record in records:
        for path in _figure_paths(record):
            if row_id.startswith("row-module-local-diversity-"):
                metric = _metric_for_figure(path)
                if not metric or not row_id.startswith(f"row-module-local-diversity-{_slug(metric)}-"):
                    continue
            if path in seen:
                continue
            seen.add(path)
            info = registry.get(path)
            if info is None:
                continue
            if info.get("owner_row") == row_id:
                rendered.append(_render_figure(path, info, output_root, destination_parent))
            else:
                owner = info.get("owner_row", "")
                title = str(info.get("title") or "Saved evidence figure")
                rendered.append(f'<a class="figure-reference" href="#{html.escape(owner, quote=True)}">Shared evidence: {html.escape(title)}</a>')
    if not rendered:
        return '<td class="evidence-figure empty">No figure recorded</td>'
    return '<td class="evidence-figure"><div class="figure-grid">' + "".join(rendered) + '</div></td>'


def _figure_cell(record: Mapping[str, Any] | None, row_id: str, registry: Mapping[str, Mapping[str, Any]], output_root: Path, destination_parent: Path) -> str:
    if record is None:
        return '<td class="evidence-figure empty">—</td>'
    return _figure_cell_records([record], row_id, registry, output_root, destination_parent)


def _proof_content(record: Mapping[str, Any] | None, output_root: Path, destination_parent: Path) -> str:
    if record is None:
        return 'No saved record'
    condition = record.get("evidence_condition", {}) or {}
    status = CONDITIONS.get(str(condition.get("status", "pending")), str(condition.get("status", "pending")))
    issues = condition.get("issues", [])
    if isinstance(issues, str):
        issues = [issues]
    condition_html = f'<div class="proof-condition">{html.escape(status)}'
    if issues:
        condition_html += " · " + html.escape("; ".join(str(issue) for issue in issues))
    condition_html += "</div>"
    parts = [f'<div class="proof-metrics">{_key_values(record) or "Saved result fields: NA"}</div>', condition_html]
    tables = _tables(record)
    if tables:
        links = []
        for path, exists in tables:
            href = _href(path, output_root, destination_parent)
            label = Path(path).name + (" [missing]" if not exists else "")
            links.append(f'<a class="csv-link" href="{html.escape(href, quote=True)}" title="{html.escape(path, quote=True)}">{html.escape(label)}</a>')
        parts.append('<div class="proof-links">' + " · ".join(links) + "</div>")
    else:
        parts.append('<div class="proof-links empty">No saved CSV link</div>')
    return "".join(parts)


def _local_metric_proof(record: Mapping[str, Any], metric: str, output_root: Path, destination_parent: Path) -> str:
    results = record.get("results", {}) or {}
    baseline = str((record.get("identity", {}) or {}).get("baseline") or results.get("baseline") or "raw_leiden")
    metric_values = (results.get("metrics", {}) or {}).get(metric, {}) or {}
    values = metric_values.get(baseline, metric_values) if isinstance(metric_values, Mapping) else {}
    compact = "<br>".join(
        f'<span class="metric">{html.escape(label)}: <strong>{html.escape(_number(value))}</strong></span>'
        for label, value in (
            ("Scale (μm)", results.get("scale_um")),
            ("Valid windows", results.get("n_valid_windows")),
            ("Baseline windows", values.get("baseline_n_observations")),
            ("Baseline median", values.get("baseline_median")),
            ("Recon median", values.get("reconstruction_median")),
            ("Delta median", values.get("delta_median")),
            ("Delta Q1", values.get("delta_q1")),
            ("Delta Q3", values.get("delta_q3")),
        )
    )
    proof = _proof_content(record, output_root, destination_parent)
    proof = re.sub(r'<div class="proof-metrics">.*?</div>', f'<div class="proof-metrics">{compact}</div>', proof, count=1, flags=re.DOTALL)
    return proof


LOCAL_DIVERSITY_METRICS = (("kobs", "Kobs"), ("neff", "Neff"), ("evenness", "evenness"))


def _local_metric_values(record: Mapping[str, Any] | None, metric: str) -> Mapping[str, Any]:
    if record is None:
        return {}
    values = ((record.get("results", {}) or {}).get("metrics", {}) or {}).get(metric, {}) or {}
    baseline = str((record.get("identity", {}) or {}).get("baseline") or (record.get("results", {}) or {}).get("baseline") or "")
    return values.get(baseline, values) if isinstance(values, Mapping) else {}


def _local_diversity_table(*, metric: str, label: str, records: Iterable[Mapping[str, Any]], registry: Mapping[str, Mapping[str, Any]], output_root: Path, destination_parent: Path, all_scope: bool = False) -> str:
    records = list(records)
    by_scope = _record_map(records)
    scopes = ("All",) if all_scope else SCOPES
    parts = []
    for question, baseline in (("local_vs_raw_leiden", "Raw Leiden"), ("local_vs_raw_level2", "Raw Level2")):
        rows = []
        for scope in scopes:
            record = by_scope.get((question, scope))
            if record is None:
                continue
            values = _local_metric_values(record, label)
            results = record.get("results", {}) or {}
            rows.append((_display_scope(scope), values.get("baseline_median"), values.get("reconstruction_median"), values.get("delta_median"), results.get("n_valid_windows"), results.get("scale_um")))
        parts.append(f"<h4>{html.escape(label)} · {baseline}</h4>" + _foundation_table(("Cell type", "Raw median", "Recon median", "Paired Δ median", "Valid windows", "Scale (μm)"), rows))
    parts.append('<div class="table-scroll"><table class="figure-evidence-table"><tbody>')
    for scope in scopes:
        row_id = f"row-module-local-diversity-{_slug(metric)}-{_slug(scope)}{'-all' if all_scope else ''}"
        selected = [by_scope[(q, scope)] for q in ("local_vs_raw_leiden", "local_vs_raw_level2") if (q, scope) in by_scope]
        anchors = "".join(f'<span id="row-module-{_slug(q)}-{_slug(scope)}-{_slug(metric)}"></span>' for q in ("local_vs_raw_leiden", "local_vs_raw_level2"))
        parts.append(f'<tr id="{row_id}"><th>{anchors}{html.escape(_display_scope(scope))} · {html.escape(label)}</th>{_figure_cell_records(selected, row_id, registry, output_root, destination_parent)}</tr>')
    return "".join(parts) + "</tbody></table></div>"


def _local_diversity_module(*, metric: str, label: str, records: Iterable[Mapping[str, Any]], registry: Mapping[str, Mapping[str, Any]], output_root: Path, destination_parent: Path, all_scope: bool = False) -> str:
    module_id = f"module-local-diversity-{_slug(metric)}{'-all' if all_scope else ''}"
    intro = {
        "kobs": "Kobs records observed label richness. It is descriptive and does not set a higher-is-better direction.",
        "neff": "Neff keeps the saved baseline-specific directional interpretation and support condition.",
        "evenness": "Evenness describes balance conditional on richness and does not define a higher-is-better direction.",
    }[metric]
    return f'<article class="report-module" id="{module_id}"><p class="module-intro">{html.escape(intro)}</p>{_local_diversity_table(metric=metric, label=label, records=records, registry=registry, output_root=output_root, destination_parent=destination_parent, all_scope=all_scope)}</article>'


def _scale_proof_content(record: Mapping[str, Any], output_root: Path, destination_parent: Path) -> str:
    rows = [row for row in (record.get("results", {}) or {}).get("rows", []) if isinstance(row, Mapping)]
    minimum_units = []
    for row in rows:
        value = row.get("min_parent_units")
        if value is not None and _number(value) not in minimum_units:
            minimum_units.append(_number(value))
    minimum_text = ", ".join(minimum_units) if minimum_units else "NA"
    note = '<div class="scale-proof-note">RL = Raw Leiden · L2 = Raw Level2 · Minimum parent units: ' + html.escape(minimum_text) + '</div>'
    table = ['<table class="scale-proof-table"><thead><tr><th>Scale (μm)</th><th>ΔNeff RL</th><th>ΔNeff L2</th><th>Valid windows</th></tr></thead><tbody>']
    for row in rows:
        table.append(
            "<tr>"
            f"<td>{html.escape(_number(row.get('window_side_length')))}</td>"
            f"<td>{html.escape(_number(row.get('median_delta_neff_vs_raw_leiden')))}</td>"
            f"<td>{html.escape(_number(row.get('median_delta_neff_vs_raw_level2')))}</td>"
            f"<td>{html.escape(_number(row.get('n_valid_windows')))}</td>"
            "</tr>"
        )
    if not rows:
        table.append('<tr><td colspan="4">No saved scale rows</td></tr>')
    table.append("</tbody></table>")
    proof = _proof_content(record, output_root, destination_parent)
    proof = re.sub(r'<div class="proof-metrics">.*?</div>', '<div class="proof-metrics">' + note + "".join(table) + "</div>", proof, count=1, flags=re.DOTALL)
    return proof


def _complexity_comparison_rows(records: Iterable[Mapping[str, Any]]) -> list[tuple[Any, ...]]:
    rows: list[tuple[Any, ...]] = []
    for record in records:
        identity = record.get("identity", {}) or {}
        results = record.get("results", {}) or {}
        _kind, target = _complexity_kind(results)
        rows.append(
            (
                _display_scope(identity.get("scope", "")),
                results.get("raw_k"),
                target,
            )
        )
    return rows


def _complexity_table(records: Iterable[Mapping[str, Any]], registry: Mapping[str, Mapping[str, Any]], output_root: Path, destination_parent: Path, *, all_scope: bool = False) -> str:
    """Keep 2.1 as a small Raw-K versus route-specific comparison table."""

    records = list(records)
    row_id = "row-module-complexity-all" if all_scope else "row-module-complexity-batch"
    comparison_rows = _complexity_comparison_rows(records)
    fixed_final = any(
        _complexity_kind(record.get("results", {}) or {})[0] == "Raw → X fixed final K"
        for record in records
    )
    target_label = "Recon fixed final-cluster K" if fixed_final else "Recon same-resolution K"
    route_note = (
        "Recon K is the fixed final-cluster K for this Xenium comparison."
        if fixed_final
        else "Recon K is measured at the Raw resolution for this same-resolution comparison."
    )
    body = [f'<p class="compact-note">{html.escape(route_note)}</p><table class="evidence-table judgment-table complexity-table"><thead><tr><th>Cell type</th><th>Raw K</th><th>{html.escape(target_label)}</th></tr></thead><tbody>']
    if records:
        for scope, raw_k, target in comparison_rows:
            body.append(
                "<tr>"
                f"<td>{html.escape(_number(scope))}</td>"
                f"<td>{html.escape(_number(raw_k))}</td>"
                f"<td>{html.escape(_number(target))}</td>"
                "</tr>"
            )
    else:
        body.append('<tr><td colspan="3" class="empty">No saved record</td></tr>')
    body.append("</tbody></table>")
    figures = (
        f'<details class="figure-audit"><summary>Supporting partition figures (folded; not the same-resolution K evidence)</summary><table class="figure-evidence-table"><thead><tr><th>Parent / shared figure</th><th>图形证据</th></tr></thead><tbody><tr id="{row_id}"><th>Saved complexity figures</th>'
        + _figure_cell_records(records, row_id, registry, output_root, destination_parent)
        + "</tr></tbody></table></details>"
        if records
        else ""
    )
    return '<div class="table-scroll">' + "".join(body) + figures + '</div>'


def _matched_k_table(records: Iterable[Mapping[str, Any]], registry: Mapping[str, Mapping[str, Any]], output_root: Path, destination_parent: Path, *, all_scope: bool = False) -> str:
    """Render matched-K diagnostics with their eligibility restriction visible."""

    records = list(records)
    row_id = "row-module-matched_k-all" if all_scope else "row-module-matched_k-batch"
    body = ['<table class="evidence-table judgment-table matched-k-table"><thead><tr><th>Parent</th><th>Raw / Recon K</th><th>Matched-K status</th><th>Diagnostic change</th><th>Diagnostic ARI</th></tr></thead><tbody>']
    if records:
        for record in records:
            identity = record.get("identity", {}) or {}
            results = record.get("results", {}) or {}
            condition = record.get("evidence_condition", {}) or {}
            unavailable = results.get("matched") is False or results.get("headline_eligible") is False or condition.get("matched_k_status") == "unmatched_cluster_complexity"
            status = "matched-K unavailable · headline restricted" if unavailable else _judgment(record, "zh")
            change = results.get("unit_change_fraction")
            ari = results.get("ari")
            if unavailable:
                change_cell = '<details><summary>NA · diagnostic only</summary>' + html.escape(_number(change, percent=True)) + "</details>"
                ari_cell = '<details><summary>NA · diagnostic only</summary>' + html.escape(_number(ari)) + "</details>"
            else:
                change_cell = html.escape(_number(change, percent=True))
                ari_cell = html.escape(_number(ari))
            body.append(
                "<tr>"
                f"<td>{html.escape(_display_scope(identity.get('scope', '')))}</td>"
                f"<td>{html.escape(_number(results.get('raw_k')))} / {html.escape(_number(results.get('reconstruction_k')))}</td>"
                f"<td>{html.escape(status)}</td>"
                f"<td>{change_cell}</td>"
                f"<td>{ari_cell}</td>"
                "</tr>"
            )
    else:
        body.append('<tr><td colspan="5" class="empty">No saved record</td></tr>')
    body.append("</tbody></table>")
    figures = (
        f'<table class="figure-evidence-table"><thead><tr><th>Parent / shared figure</th><th>图形证据</th></tr></thead><tbody><tr id="{row_id}"><th>Saved matched-K figures</th>'
        + _figure_cell_records(records, row_id, registry, output_root, destination_parent)
        + "</tr></tbody></table>"
        if records
        else ""
    )
    return '<div class="table-scroll">' + "".join(body) + figures + '<p class="compact-note">Assignment difference is a matched-K diagnostic; an unavailable match limits headline and dependent interpretation.</p></div>'


def _notebook_details(record: Mapping[str, Any]) -> str:
    """Keep finite saved issues and limitations visible in short notebook rows."""

    condition = record.get("evidence_condition", {}) or {}
    details: list[str] = []
    issues = condition.get("issues", []) if isinstance(condition, Mapping) else []
    if isinstance(issues, str):
        issues = [issues]
    details.extend(str(item) for item in issues if item)
    limitations = record.get("limitations", {}) or {}
    if isinstance(limitations, Mapping):
        values = limitations.get("en") or limitations.get("zh") or []
    else:
        values = limitations
    if isinstance(values, str):
        values = [values]
    details.extend(str(item) for item in values if item)
    return " · ".join(details)


def _evidence_table(*, layer: str, question: str, records: Iterable[Mapping[str, Any]], all_scope: bool, registry: Mapping[str, Mapping[str, Any]], output_root: Path, destination_parent: Path) -> str:
    module = _module_id(layer, question, all_scope=all_scope)
    by_scope = _record_map(records)
    scopes = ("All",) if all_scope else SCOPES
    if question in {"moran_all_valid", "moran_shared_valid", "emt_coverage", "emt_score", "changed_units"}:
        selected = [by_scope[(question, scope)] for scope in scopes if (question, scope) in by_scope]
        parts, numeric_rows = [], []
        for record in selected:
            scope = _display_scope(record["identity"]["scope"])
            results = record.get("results", {})
            raw, recon = results.get("raw", {}), results.get("reconstruction", {})
            if question.startswith("moran_"):
                headers = ("Cell type", "Raw valid N", "Recon valid N", "Raw Q75", "Recon Q75", "Paired ΔI median")
                numeric_rows.append((scope, raw.get("n_valid"), recon.get("n_valid"), raw.get("q75"), recon.get("q75"), results.get("paired_delta", {}).get("median")))
                if question == "moran_all_valid":
                    headers = headers[:-1]
                    numeric_rows[-1] = numeric_rows[-1][:-1]
            elif question == "emt_coverage":
                headers = ("Cell type", "Raw resource hits / total", "Recon resource hits / total")
                numeric_rows.append((scope, f'{raw.get("available_gene_count")} / {raw.get("resource_gene_count")}', f'{recon.get("available_gene_count")} / {recon.get("resource_gene_count")}'))
            elif question == "emt_score":
                headers = ("Cell type", "Raw median", "Recon median", "Paired Δ median", "Paired units", "Raw → Recon rank length")
                cutoff = results.get("cutoff_audit", {})
                ranks = " → ".join(_number(cutoff.get(side, {}).get("effective_rank_length")) for side in ("raw", "reconstruction"))
                numeric_rows.append((scope, raw.get("median"), recon.get("median"), results.get("paired", {}).get("median_delta"), results.get("paired", {}).get("n"), ranks))
            else:
                rows = [(r.get("level1_region"), f'{r.get("changed_units")} / {r.get("paired_units")}', _number(r.get("change_fraction"), percent=True), r.get("n_valid_windows")) for r in results.get("by_region", [])]
                condition = "" if results.get("overall", {}).get("headline_eligible") else " · matched-K 未成立：仅作诊断"
                parts.append(f'<h4>{html.escape(scope + condition)}</h4>' + _foundation_table(("Anatomy", "Changed / paired units", "Change fraction", "Valid windows"), rows))
        if numeric_rows:
            parts.append(_foundation_table(headers, numeric_rows))
        groups = [("batch", selected)] if question in GROUPED_QUESTIONS else [(r["identity"]["scope"], [r]) for r in selected]
        parts.append('<div class="table-scroll"><table class="figure-evidence-table"><tbody>')
        for scope, group in groups:
            row_id = _row_id(module, scope)
            parts.append(f'<tr id="{row_id}"><th>{html.escape(_labels().get(question, (question, question))[0])}</th>{_figure_cell_records(group, row_id, registry, output_root, destination_parent)}</tr>')
        return "".join(parts) + "</tbody></table></div>"
    rows = ['<table class="evidence-table"><thead><tr><th class="evidence-question">科学问题／观察方向</th><th class="evidence-conclusions">本批结论</th><th class="evidence-proof">数值证明与条件</th><th class="evidence-figure">图形证据</th></tr></thead><tbody>']
    if question.startswith("local_vs_") and not all_scope:
        label = _labels().get(question, (question, question))[0]
        for scope in scopes:
            record = by_scope.get((question, scope))
            if record is None:
                row_id = _row_id(module, scope)
                rows.append(f'<tr id="{row_id}"><th class="evidence-question"><span class="scope">{html.escape(_display_scope(scope))}</span></th><td class="evidence-conclusions empty">No saved record</td><td class="evidence-proof empty">No saved record</td><td class="evidence-figure empty">—</td></tr>')
                continue
            for metric in ("kobs", "neff", "evenness"):
                row_id = _row_id(module, scope) + "-" + metric
                metric_name = metric.upper() if metric == "kobs" else metric.capitalize()
                direction = _judgment(record, "zh") if metric == "neff" else "描述性指标；方向判读保留在 Neff 行"
                conclusion = str((record.get("conclusion", {}) or {}).get("zh", "")) if metric == "neff" else "同一 baseline 的保存指标图；方向判读保留在 Neff 行。"
                evidence = dict(record.get("evidence", {}) or {})
                evidence["figures"] = [path for path in _figure_paths(record) if _metric_for_figure(path) == metric]
                metric_record = dict(record)
                metric_record["evidence"] = evidence
                metric_key = {"kobs": "Kobs", "neff": "Neff", "evenness": "evenness"}[metric]
                rows.append(f'<tr class="metric-subrow" id="{row_id}"><th class="evidence-question"><span class="question-label">{html.escape(label)} · {html.escape(metric_name)}</span><span class="scope">{html.escape(_display_scope(scope))}</span><span class="direction">{html.escape(direction)}</span></th><td class="evidence-conclusions">{html.escape(conclusion)}</td><td class="evidence-proof">{_local_metric_proof(record, metric_key, output_root, destination_parent)}</td>{_figure_cell(metric_record, row_id, registry, output_root, destination_parent)}</tr>')
        return '<div class="table-scroll">' + "".join(rows) + "</tbody></table></div>"
    if question in GROUPED_QUESTIONS:
        grouped_records = [by_scope.get((question, scope)) for scope in scopes]
        grouped_records = [record for record in grouped_records if record is not None]
        row_id = _row_id(module, "batch")
        label = _labels().get(question, (question, question))[0]
        if not grouped_records:
            rows.append(f'<tr id="{row_id}"><th class="evidence-question"><span class="question-label">{html.escape(label)}</span><span class="scope">Batch</span></th><td class="evidence-conclusions empty">No saved record</td><td class="evidence-proof empty">No saved record</td><td class="evidence-figure empty">—</td></tr>')
        else:
            directions = "".join(
                f'<span class="scope-block"><strong>{html.escape(_display_scope((record.get("identity", {}) or {}).get("scope", "")))}</strong> · {html.escape(_judgment(record, "zh"))}</span>'
                for record in grouped_records
            )
            conclusions = "".join(
                f'<div class="scope-block"><strong>{html.escape(_display_scope((record.get("identity", {}) or {}).get("scope", "")))}</strong><br>{html.escape(str((record.get("conclusion", {}) or {}).get("zh", "")))}</div>'
                for record in grouped_records
            )
            proof = "".join(
                f'<div class="scope-block"><strong>{html.escape(_display_scope((record.get("identity", {}) or {}).get("scope", "")))}</strong>{_proof_content(record, output_root, destination_parent)}</div>'
                for record in grouped_records
            )
            rows.append(f'<tr id="{row_id}"><th class="evidence-question"><span class="question-label">{html.escape(label)}</span><span class="scope">Batch · {html.escape(" / ".join(_display_scope((record.get("identity", {}) or {}).get("scope", "")) for record in grouped_records))}</span><span class="direction">{directions}</span></th><td class="evidence-conclusions">{conclusions}</td><td class="evidence-proof">{proof}</td>{_figure_cell_records(grouped_records, row_id, registry, output_root, destination_parent)}</tr>')
        return '<div class="table-scroll">' + "".join(rows) + "</tbody></table></div>"
    for scope in scopes:
        record = by_scope.get((question, scope))
        row_id = _row_id(module, scope)
        if record is None:
            rows.append(f'<tr id="{row_id}"><th class="evidence-question"><span class="scope">{html.escape(_display_scope(scope))}</span></th><td class="evidence-conclusions empty">No saved record</td><td class="evidence-proof empty">No saved record</td><td class="evidence-figure empty">—</td></tr>')
            continue
        label = _labels().get(question, (question, question))[0]
        identity = record.get("identity", {}) or {}
        baseline = identity.get("baseline")
        baseline_text = f" · baseline={baseline}" if baseline else ""
        direction = _judgment(record, "zh")
        proof = _scale_proof_content(record, output_root, destination_parent) if question == "scale_sensitivity" else _proof_content(record, output_root, destination_parent)
        rows.append(f'<tr id="{row_id}"><th class="evidence-question"><span class="question-label">{html.escape(label)}</span><span class="scope">{html.escape(_display_scope(scope))}</span><span class="direction">{html.escape(direction + baseline_text)}</span></th><td class="evidence-conclusions">{html.escape(str((record.get("conclusion", {}) or {}).get("zh", "")))}</td><td class="evidence-proof">{proof}</td>{_figure_cell(record, row_id, registry, output_root, destination_parent)}</tr>')
    return '<div class="table-scroll">' + "".join(rows) + "</tbody></table></div>"


def _module(*, layer: str, question: str, records: Iterable[Mapping[str, Any]], all_scope: bool, registry: Mapping[str, Mapping[str, Any]], output_root: Path, destination_parent: Path, next_module: str | None = None, next_label: str | None = None) -> str:
    module_id = _module_id(layer, question, all_scope=all_scope)
    label = _labels().get(question, (question, question))[0]
    node = _node_for_question(question)
    intro = str(node["purpose_zh"]) if node is not None else "这一项读取已保存结果和证据条件。"
    table = _evidence_table(layer=layer, question=question, records=records, all_scope=all_scope, registry=registry, output_root=output_root, destination_parent=destination_parent)
    bridge = ""
    if next_module and next_label:
        bridge = f'<p class="module-next">下一问题：<a href="#{html.escape(next_module, quote=True)}">{html.escape(next_label)}</a></p>'
    return f'<article class="report-module" id="{module_id}"><h3>{html.escape(label)}</h3><p class="module-intro">{html.escape(intro)} · {_method_link(question, destination_parent)}</p>{table}{bridge}</article>'


def _overview_cell(record: Mapping[str, Any] | None) -> str:
    if record is None:
        return '<td class="overview-cell" data-scope="unknown">—</td>'
    label, value = _headline_metric(record)
    scope = _display_scope((record.get("identity", {}) or {}).get("scope", ""))
    return f'<td class="overview-cell" data-scope="{html.escape(scope, quote=True)}"><span class="overview-judgment">{html.escape(_overview_conclusion(record))}</span><span class="overview-headline">{html.escape(label)}: <strong>{html.escape(value)}</strong></span></td>'


def _overview_layer(layer: str, records: Iterable[Mapping[str, Any]]) -> str:
    subset = [record for record in records if _record_layer(record) == layer and (record.get("identity", {}) or {}).get("scope") != "All"]
    questions = _question_order(subset)
    by_scope = _record_map(subset)
    rows = [f'<section class="overview-layer" id="overview-{_slug(layer)}"><h3>{html.escape(TITLES[layer][0])}</h3><p class="layer-intro">{html.escape(LAYER_INTROS[layer])}</p><div class="table-scroll"><table class="overview-table"><thead><tr><th>科学问题</th>']
    rows.append("".join(f'<th data-scope="{html.escape(_display_scope(scope), quote=True)}">{html.escape(_display_scope(scope))}</th>' for scope in SCOPES))
    rows.append("</tr></thead><tbody>")
    for question in questions:
        label = _labels().get(question, (question, question))[0]
        module = _module_id(layer, question)
        cells = "".join(_overview_cell(by_scope.get((question, scope))) for scope in SCOPES)
        rows.append(f'<tr class="overview-question"><th><a class="module-link" href="#{module}">{html.escape(label)}</a></th>{cells}</tr>')
    if not questions:
        rows.append('<tr class="overview-question"><th colspan="4">No saved question records</th></tr>')
    return "".join(rows) + "</tbody></table></div></section>"


def _overview_hd(records: Iterable[Mapping[str, Any]]) -> str:
    subset = [record for record in records if (record.get("identity", {}) or {}).get("scope") == "All"]
    if not subset:
        return ""
    rows = ['<section class="overview-hd" id="overview-hd-all"><h3>HD All：独立队列与分母</h3><p class="layer-intro">HD All 单独保留其 cohort 与分母，不混入三个 parent 的内部比较。</p>']
    for layer in LAYERS:
        layer_records = [record for record in subset if _record_layer(record) == layer]
        questions = _question_order(layer_records)
        if not questions:
            continue
        by_scope = _record_map(layer_records)
        rows.append(f'<h4>{html.escape(TITLES[layer][0])}</h4><div class="table-scroll"><table class="overview-table overview-hd-table"><thead><tr><th>科学问题</th><th>All</th></tr></thead><tbody>')
        for question in questions:
            label = _labels().get(question, (question, question))[0]
            module = _module_id(layer, question, all_scope=True)
            rows.append(f'<tr class="overview-question"><th><a class="module-link" href="#{module}">{html.escape(label)}</a></th>{_overview_cell(by_scope.get((question, "All")))}</tr>')
        rows.append("</tbody></table></div>")
    return "".join(rows) + "</section>"


OVERVIEW_GROUPS = (
    ("partition", "分群", "Partition", ("complexity", "matched_k"), "2.1"),
    ("moran", "Moran", "Moran", ("moran_all_valid", "moran_shared_valid"), "2.4"),
    ("emt", "EMT", "EMT", ("emt_coverage", "emt_score"), "2.5"),
    ("local", "局部双 baseline", "Local state · two Raw baselines", ("local_vs_raw_leiden", "local_vs_raw_level2"), "3.5"),
    ("region", "Region", "Region", ("threshold_reliability", "region_extent", "scale_sensitivity"), "4.1"),
)


def _overview_group_records(records: Iterable[Mapping[str, Any]], questions: Iterable[str], scope: str) -> dict[str, Mapping[str, Any]]:
    wanted = set(questions)
    return {
        str((record.get("identity", {}) or {}).get("question", "")): record
        for record in records
        if str((record.get("identity", {}) or {}).get("scope", "")) == scope
        and str((record.get("identity", {}) or {}).get("question", "")) in wanted
    }


def _complexity_kind(results: Mapping[str, Any]) -> tuple[str, Any]:
    comparison_kind = str(results.get("comparison_kind", "")).lower()
    fixed = results.get("fixed_final_cluster_k")
    if fixed is not None or "fixed_final" in comparison_kind:
        return "Raw → X fixed final K", fixed
    reconstruction = results.get("reconstruction_k_at_raw_resolution")
    if reconstruction is not None or "same_resolution" in comparison_kind or "same-resolution" in comparison_kind:
        return "Raw → Recon same-resolution K", reconstruction
    return "Raw → Recon K", reconstruction


def _overview_group_metric(group: str, grouped: Mapping[str, Mapping[str, Any]]) -> tuple[str, str]:
    """Select compact saved facts for the five overview questions."""

    if group == "partition":
        complexity = grouped.get("complexity")
        matched = grouped.get("matched_k")
        if complexity is None and matched is None:
            return "Saved fact", "NA"
        pairs: list[str] = []
        if complexity is not None:
            results = complexity.get("results", {}) or {}
            label, target = _complexity_kind(results)
            pairs.append(f"K {_number(results.get('raw_k'))} → {_number(target)}")
            comparison_label = "fixed final" if label == "Raw → X fixed final K" else "same-resolution"
        else:
            comparison_label = "saved"
        if matched is not None:
            results = matched.get("results", {}) or {}
            condition = matched.get("evidence_condition", {}) or {}
            unavailable = (
                results.get("matched") is False
                or results.get("headline_eligible") is False
                or condition.get("matched_k_status") == "unmatched_cluster_complexity"
            )
            status = "matched-K unavailable" if unavailable else "matched-K eligible"
            change = _number(results.get("unit_change_fraction"), percent=True)
            change_label = "change diagnostic" if unavailable else "change"
            pairs.append(f"{status}; {change_label} {change}")
        return f"K / matched-K ({comparison_label})", "; ".join(pairs)
    if group == "moran":
        all_valid = grouped.get("moran_all_valid")
        shared = grouped.get("moran_shared_valid")
        if all_valid is None and shared is None:
            return "Saved fact", "NA"
        pairs: list[str] = []
        if all_valid is not None:
            results = all_valid.get("results", {}) or {}
            raw_q75 = (results.get("raw", {}) or {}).get("q75")
            recon_q75 = (results.get("reconstruction", {}) or {}).get("q75")
            pairs.append(f"Q75 {_number(raw_q75)} → {_number(recon_q75)}")
        if shared is not None:
            delta = ((shared.get("results", {}) or {}).get("paired_delta", {}) or {}).get("median")
            pairs.append(f"paired ΔI {_number(delta)}")
        return "Moran", "; ".join(pairs)
    if group == "emt":
        coverage = grouped.get("emt_coverage")
        score = grouped.get("emt_score")
        if coverage is None and score is None:
            return "Saved fact", "NA"
        pairs = []
        if coverage is not None:
            results = coverage.get("results", {}) or {}
            raw = results.get("raw", {}) or {}
            recon = results.get("reconstruction", {}) or {}
            pairs.append(
                f"hits {_number(raw.get('available_gene_count'))}/{_number(raw.get('resource_gene_count'))} → "
                f"{_number(recon.get('available_gene_count'))}/{_number(recon.get('resource_gene_count'))}"
            )
        if score is not None:
            delta = (((score.get("results", {}) or {}).get("paired", {}) or {}).get("median_delta"))
            pairs.append(f"paired Δ {_number(delta)}")
        return "EMT", "; ".join(pairs)
    if group == "local":
        pairs = []
        for question, label in (("local_vs_raw_leiden", "Leiden"), ("local_vs_raw_level2", "Level2")):
            record = grouped.get(question)
            if record is None:
                continue
            values = _local_metric_values(record, "Neff")
            pairs.append(f"{label} ΔNeff {_number(values.get('delta_median'))}")
        return "Neff median Δ", "; ".join(pairs) if pairs else "NA"
    if group == "region":
        extent = grouped.get("region_extent")
        threshold = grouped.get("threshold_reliability")
        if extent is not None:
            return "Area fraction", _number(((extent.get("results", {}) or {}).get("overall", {}) or {}).get("area_fraction"), percent=True)
        if threshold is not None:
            return _headline_metric(threshold)
        return "Area fraction", "NA"
    return "Saved fact", "NA"


def _overview_group_cell(group: str, grouped: Mapping[str, Mapping[str, Any]]) -> str:
    # Partition status is governed by matched-K eligibility when that record is
    # available; the saved complexity record still supplies the K transition.
    record = grouped.get("matched_k") if group == "partition" and grouped.get("matched_k") is not None else next(iter(grouped.values()), None)
    if record is None:
        return '<td class="overview-cell" data-scope="unknown">—</td>'
    label, value = _overview_group_metric(group, grouped)
    # Keep the cell judgment tied to the saved representative record.  The
    # detailed module carries the complete per-question conditions.
    judgment = _overview_conclusion(record)
    scope = _display_scope((record.get("identity", {}) or {}).get("scope", ""))
    return (
        f'<td class="overview-cell" data-scope="{html.escape(scope, quote=True)}">'
        f'<span class="overview-judgment">{html.escape(judgment)}</span>'
        f'<span class="overview-headline">{html.escape(label)}: <strong>{html.escape(value)}</strong></span></td>'
    )


def _overview_reorganized(records: Iterable[Mapping[str, Any]]) -> str:
    """Render the approved five-row reader overview for narrative packages."""

    records = list(records)
    parent_records = [record for record in records if (record.get("identity", {}) or {}).get("scope") != "All"]
    rows = ['<section class="overview-main" id="overview-main"><h3>五项判断</h3><p class="layer-intro">按分群、Moran、EMT、局部双 baseline 和 Region 并列查看三个 parent；每格只保留保存的关键事实。</p><div class="table-scroll"><table class="overview-table"><thead><tr><th>判断</th>']
    rows.append("".join(f'<th data-scope="{html.escape(_display_scope(scope), quote=True)}">{html.escape(_display_scope(scope))}</th>' for scope in SCOPES))
    rows.append("</tr></thead><tbody>")
    for group, title_zh, _title_en, questions, node_id in OVERVIEW_GROUPS:
        cells = "".join(_overview_group_cell(group, _overview_group_records(parent_records, questions, scope)) for scope in SCOPES)
        rows.append(
            f'<tr class="overview-question overview-group" id="overview-group-{_slug(group)}">'
            f'<th><a class="module-link" href="#node-{_slug(node_id)}">{html.escape(title_zh)}</a></th>{cells}</tr>'
        )
    rows.append("</tbody></table></div></section>")

    all_records = [record for record in records if (record.get("identity", {}) or {}).get("scope") == "All"]
    if all_records:
        rows.append('<section class="overview-hd" id="overview-hd-all"><h3>HD global：独立 cohort 与分母</h3><p class="layer-intro">HD global 单独保留 cohort、分母和保存条件，不混入三个 parent 的内部比较。</p><div class="table-scroll"><table class="overview-table overview-hd-table"><thead><tr><th>判断</th><th>HD global</th></tr></thead><tbody>')
        for group, title_zh, _title_en, questions, node_id in OVERVIEW_GROUPS:
            grouped = _overview_group_records(all_records, questions, "All")
            rows.append(
                f'<tr class="overview-question overview-group" id="overview-hd-{_slug(group)}">'
                f'<th><a class="module-link" href="#node-{_slug(node_id)}-all">{html.escape(title_zh)}</a></th>{_overview_group_cell(group, grouped)}</tr>'
            )
        rows.append("</tbody></table></div></section>")
    return "".join(rows)


SUMMARY_NODE_IDS = {"1.6", "2.6", "3.9", "4.5", "summary"}
FOUNDATION_NODE_IDS = {"1.1", "1.2", "1.3", "1.4", "1.5"}
COMPACT_NODE_IDS = {"2.3", "3.8", "4.3"}
LOCAL_NODE_METRICS = {"3.5": ("kobs", "Kobs"), "3.6": ("neff", "Neff"), "3.7": ("evenness", "evenness")}
LEGACY_ANCHOR_NODE = {
    "foundation": "1.1", "complexity": "2.1", "matched_k": "2.2",
    "moran_all_valid": "2.4", "moran_shared_valid": "2.4", "emt_coverage": "2.5",
    "emt_score": "2.5", "anatomy": "3.1", "changed_units": "3.2",
    "emt_spatial_fields": "3.3", "window_support": "3.4",
    "local_vs_raw_leiden": "3.5", "local_vs_raw_level2": "3.5",
    "threshold_reliability": "4.1", "region_extent": "4.2", "scale_sensitivity": "4.4",
}


def _node_method_link(node: Mapping[str, Any], destination_parent: Path) -> str:
    doc, anchor = node["method"]
    target = _repository_root(destination_parent) / "docs" / "design" / "reconstruction-impact" / str(doc)
    href = _href(str(target), _repository_root(destination_parent), destination_parent)
    if anchor:
        href += "#" + str(anchor)
    return f'<a class="method-link" href="{html.escape(href, quote=True)}">方法与判读规范</a>'


def _node_records(node: Mapping[str, Any], records: Iterable[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    questions = set(node["questions"])
    return [record for record in records if str((record.get("identity", {}) or {}).get("question", "")) in questions]


def _section_summary(package: Mapping[str, Any], node_id: str) -> Mapping[str, Any] | None:
    summaries = package.get("section_summaries", [])
    if not isinstance(summaries, list):
        return None
    return next((item for item in summaries if isinstance(item, Mapping) and item.get("node_id") == node_id), None)


def _summary_review_status(package: Mapping[str, Any], summary: Mapping[str, Any] | None) -> str:
    if summary is None:
        return "pending"
    review = summary.get("review", {}) if isinstance(summary.get("review", {}), Mapping) else {}
    digest = package.get("record_digest")
    if review.get("status") == "reviewed" and digest and review.get("record_digest") == digest:
        return "reviewed"
    return "pending"


def _summary_block(package: Mapping[str, Any], node: Mapping[str, Any]) -> str:
    summary = _section_summary(package, str(node["node_id"]))
    status = _summary_review_status(package, summary)
    if summary is None or status != "reviewed":
        return '<p class="section-summary pending">Review: pending · No reviewed section summary.</p>'
    conclusion = summary.get("conclusion", {}) if isinstance(summary.get("conclusion", {}), Mapping) else {}
    text = str(conclusion.get("zh", ""))
    return f'<p class="section-summary {html.escape(status, quote=True)}">Review: {html.escape(status)} · {html.escape(text or "No saved section summary.")}</p>'


def _foundation_table(headers: Iterable[str], rows: Iterable[Iterable[Any]]) -> str:
    body = ["<div class=\"table-scroll\"><table class=\"foundation-table\"><thead><tr>"]
    body.extend(f"<th>{html.escape(header)}</th>" for header in headers)
    body.append("</tr></thead><tbody>")
    has_rows = False
    for row in rows:
        has_rows = True
        body.append("<tr>")
        body.extend(f"<td>{html.escape(_number(value))}</td>" for value in row)
        body.append("</tr>")
    if not has_rows:
        body.append('<tr><td colspan="8">No saved foundation fact</td></tr>')
    return "".join(body) + "</tbody></table></div>"


def _source_path(results: Mapping[str, Any], side: str) -> Any:
    source_files = results.get("source_files", {}) or {}
    rows = source_files.get("rows", []) if isinstance(source_files, Mapping) else []
    roles = {"raw"} if side == "raw" else {"reconstruction", f"{results.get('scope')}_expression"}
    for row in rows:
        if isinstance(row, Mapping) and str(row.get("role", "")) in roles:
            return row.get("path")
    return None


def _preprocessing_detail(key: str, fact: Mapping[str, Any]) -> str:
    fields = fact.get("fields", {}) if isinstance(fact.get("fields", {}), Mapping) else {}
    if key == "partition":
        pairs = (("mode", fields.get("mode")), ("resolution", fields.get("within_level1_resolution")), ("top genes", fields.get("n_top_genes")), ("units", fields.get("n_units")), ("features", fields.get("n_features")))
    elif key == "moran":
        pairs = (("coordinates", fields.get("coordinate_source")), ("unit", fields.get("coordinate_unit")), ("graph units", fields.get("n_units")), ("neighbors", fields.get("n_neighbors_actual")), ("min units", fields.get("min_units")))
    else:
        cutoffs = fields.get("cutoff_by_side", {}) if isinstance(fields.get("cutoff_by_side", {}), Mapping) else {}
        raw = cutoffs.get("raw", {}) if isinstance(cutoffs.get("raw", {}), Mapping) else {}
        reconstruction = cutoffs.get("reconstruction", {}) if isinstance(cutoffs.get("reconstruction", {}), Mapping) else {}
        pairs = (("scorer", fields.get("scorer")), ("Raw rank", raw.get("effective_rank_length")), ("Recon rank", reconstruction.get("effective_rank_length")))
    values = [f"{label}={_number(value)}" for label, value in pairs if value is not None]
    return "; ".join(values) or "NA"


def _carrier_presentation(item: Mapping[str, Any]) -> tuple[Any, Any, str, str]:
    role = str(item.get("role", ""))
    spatial = item.get("spatial_axis")
    if "raw full context" in role.lower():
        return "Spatial unit", "Raw expression, Level1 labels and coordinates", "Cohort selection, Raw baselines and anatomy", "Input provenance does not establish biological truth"
    if "reconstruction full carrier" in role.lower():
        return "Reconstructed spatial unit", "Native reconstructed expression and coordinates", "Recon partition, Moran, EMT and local state", "Use Raw-defined scopes and separately audited pairing"
    if spatial is True:
        return "Spatial unit", "Spatial coordinates and final cluster assignments", "Spatial comparison axis and projection destination", "Expression donor IDs are not spatial IDs"
    return "Reference cell", "Expression grouped by SVC cluster", "Cluster-mean projection before Moran and EMT", "Projected field; not a per-spatial-cell measurement"


def _foundation_block(node: Mapping[str, Any], records: Iterable[Mapping[str, Any]]) -> str:
    node_id = str(node["node_id"])
    record_list = list(records)
    if node_id == "1.1":
        rows = []
        seen = set()
        for record in record_list:
            results = record.get("results", {}) or {}
            carrier = results.get("carrier", {}) or {}
            for item in carrier.get("rows", []) if isinstance(carrier, Mapping) else []:
                if isinstance(item, Mapping):
                    role = item.get("role")
                    if role in seen:
                        continue
                    seen.add(role)
                    units, content, use, boundary = _carrier_presentation(item)
                    rows.append((role, units, content, use, boundary))
        sources = (record_list[0].get("results", {}) or {}).get("source_files", {}) if record_list else {}
        if any(item.get("role") == "raw_level2_reference" for item in sources.get("rows", [])):
            rows.append(("Raw Level2 reference", "Reference cell", "Expression and Level2 labels", "Map the second Raw baseline", "Separate from the reconstructed expression carrier"))
        flow = '<p class="source-flow">Raw spatial units → Raw-defined comparison scopes ← Recon spatial carrier.'
        if any("expression carrier" in str(row[0]) for row in rows):
            flow += ' Xenium expression carrier → cluster mean → spatial expression field.'
        flow += '</p>'
        return flow + _foundation_table(("Carrier", "Observation units", "Content", "Use", "Boundary"), rows)
    if node_id == "1.2":
        selection_rows = []
        denominator_rows = []
        for record in record_list:
            results = record.get("results", {}) or {}
            denominator = results.get("denominator", {}) or {}
            selection = results.get("selection", {}) or {}
            excluded = denominator.get("excluded_units", {}) or {}
            scope = _display_scope((record.get("identity", {}) or {}).get("scope", ""))
            exclusions = selection.get("stage_exclusions", {})
            selection_rows.append((scope, selection.get("raw_parent_label"), selection.get("raw_parent_n"), selection.get("reconstruction_covered_parent_n"), selection.get("sampled_n"), f"Not reconstructed: {_number(exclusions.get('not_reconstructed'))}; not sampled: {_number(exclusions.get('not_sampled'))}"))
            preprocessing = results.get("preprocessing", {}) or {}
            moran = preprocessing.get("moran", {}) if isinstance(preprocessing, Mapping) else {}
            emt = preprocessing.get("emt", {}) if isinstance(preprocessing, Mapping) else {}
            moran_fields = moran.get("fields", {}) if isinstance(moran, Mapping) else {}
            emt_fields = emt.get("fields", {}) if isinstance(emt, Mapping) else {}
            cutoffs = emt_fields.get("cutoff_by_side", {}) if isinstance(emt_fields, Mapping) else {}
            raw_cutoff = cutoffs.get("raw", {}) if isinstance(cutoffs, Mapping) else {}
            recon_cutoff = cutoffs.get("reconstruction", {}) if isinstance(cutoffs, Mapping) else {}
            denominator_rows.append((scope, denominator.get("paired_units"), moran_fields.get("n_units") if isinstance(moran_fields, Mapping) else None, raw_cutoff.get("n_observations") if isinstance(raw_cutoff, Mapping) else None, recon_cutoff.get("n_observations") if isinstance(recon_cutoff, Mapping) else None, f"{scope}: separate scope; Raw-QC excluded: {_number(excluded.get('raw_qc'))}."))
        return '<h4>Parent selection from full-Raw observations</h4>' + _foundation_table(("Scope", "Raw parent", "Full-Raw parent n", "Recon covered", "Sampled", "Selection exclusions"), selection_rows) + '<h4>Saved analysis denominators</h4><p>Partition QC does not restrict Moran or EMT. Partition n counts retained assignments; Moran n counts the shared spatial graph; EMT n counts scoring inputs. Valid genes and scores are reported in their result nodes.</p>' + _foundation_table(("Scope", "Partition paired", "Moran graph n", "EMT Raw n", "EMT Recon n", "Conclusion applicability"), denominator_rows)
    if node_id == "1.3":
        rows = []
        for record in record_list:
            pairing = (record.get("results", {}) or {}).get("pairing_audit", {}) or {}
            checked = pairing.get("n_checked")
            if checked is not None:
                checked = f"{_number(checked)} (full-spatial carrier check)"
            rows.append((_display_scope((record.get("identity", {}) or {}).get("scope", "")), pairing.get("coordinate_source"), pairing.get("coordinate_unit"), checked, pairing.get("id_unique"), pairing.get("id_subset"), pairing.get("coordinates_match"), pairing.get("coordinate_finite")))
        return _foundation_table(("Scope", "Coordinate source", "Unit", "Checked", "ID unique", "ID subset", "Coordinates match", "Finite"), rows)
    if node_id == "1.4":
        rows = []
        for record in record_list:
            results = record.get("results", {}) or {}
            gene = results.get("gene", {}) or {}
            by_side = gene.get("by_side", {}) if isinstance(gene, Mapping) else {}
            for side, label in (("raw", "Raw"), ("reconstruction", "Recon")):
                facts = by_side.get(side, {}) if isinstance(by_side, Mapping) else {}
                total = facts.get("provided_n") if isinstance(facts, Mapping) else None
                source = facts.get("source") if isinstance(facts, Mapping) else None
                rows.append((_display_scope((record.get("identity", {}) or {}).get("scope", "")), f"{label} total genes", total, source or _source_path(results, side) or gene.get("path"), gene.get("status")))
        return _foundation_table(("Scope", "Side", "Total genes", "Source", "Availability status"), rows)
    rows = []
    for record in record_list:
        preprocessing = (record.get("results", {}) or {}).get("preprocessing", {}) or {}
        for key, label in (("partition", "Partition"), ("moran", "Moran"), ("emt", "EMT")):
            fact = preprocessing.get(key, {}) if isinstance(preprocessing, Mapping) else {}
            rows.append((_display_scope((record.get("identity", {}) or {}).get("scope", "")), label, _preprocessing_detail(key, fact) if isinstance(fact, Mapping) else "NA", (fact.get("path") or fact.get("source")) if isinstance(fact, Mapping) else None, fact.get("status") if isinstance(fact, Mapping) else "missing"))
    return _foundation_table(("Scope", "Analysis", "Saved method / parameters", "Source", "Status"), rows)


def _compact_block(records: Iterable[Mapping[str, Any]], node_id: str = "2.3", language: str = "zh") -> str:
    records = list(records)
    if node_id == "3.8":
        grouped: dict[tuple[str, str], list[tuple[Any, ...]]] = {}
        for record in records:
            identity = record.get("identity", {})
            scope = _display_scope(identity.get("scope", ""))
            baseline = str(identity.get("baseline", "") or "NA")
            group = grouped.setdefault((scope, baseline), [])
            for row in record.get("results", {}).get("anatomy_summary", []):
                group.append((row.get("level1_region"), row.get("n_valid_windows"), row.get("scale_um"), *[row.get(f"median_delta_{metric}_vs_{baseline}") for metric in ("k_obs", "neff", "evenness")]))
        blocks = []
        for (scope, baseline), rows in grouped.items():
            blocks.append(f'<section class="compact-subtable"><h4>{html.escape(scope)} · {html.escape(baseline)}</h4>{_foundation_table(("Anatomy", "Valid windows", "Scale (μm)", "Median ΔKobs", "Median ΔNeff", "Median Δevenness"), rows)}</section>')
        return "".join(blocks) or _foundation_table(("Anatomy", "Valid windows", "Scale (μm)", "Median ΔKobs", "Median ΔNeff", "Median Δevenness"), [])
    if node_id == "4.3":
        grouped: dict[str, list[tuple[Any, ...]]] = {}
        for record in records:
            scope = _display_scope(record.get("identity", {}).get("scope", ""))
            group = grouped.setdefault(scope, [])
            for row in record.get("results", {}).get("by_region", []):
                group.append((row.get("level1_region"), row.get("valid_windows"), row.get("region_windows"), row.get("region_area_mm2"), _number(row.get("area_fraction"), percent=True), row.get("valid_units"), _number(row.get("unit_fraction"), percent=True), row.get("threshold_status")))
        blocks = []
        for scope, rows in grouped.items():
            blocks.append(f'<section class="compact-subtable"><h4>{html.escape(scope)}</h4>{_foundation_table(("Anatomy", "Valid windows", "Region windows", "Area (mm²)", "Area fraction", "Valid units", "Unit fraction", "Threshold status"), rows)}</section>')
        return "".join(blocks) or _foundation_table(("Anatomy", "Valid windows", "Region windows", "Area (mm²)", "Area fraction", "Valid units", "Unit fraction", "Threshold status"), [])
    lines = []
    for record in records:
        identity = record.get("identity", {}) or {}
        question = identity.get("question", "")
        label = _labels().get(question, (question, question))[0 if language == "zh" else 1]
        lines.append(f'<tr><td>{html.escape(_display_scope(identity.get("scope", "")))}</td><td>{html.escape(label)}</td><td>{html.escape(_brief_conclusion(record, language))}</td><td>{_key_values(record)}<br>{html.escape(_judgment(record, language))}</td></tr>')
    return '<div class="table-scroll"><table class="foundation-table"><thead><tr><th>Scope</th><th>Question</th><th>Saved conclusion</th><th>Evidence and condition</th></tr></thead><tbody>' + "".join(lines) + "</tbody></table></div>"


def _node_article(package: Mapping[str, Any], node: Mapping[str, Any], records: list[Mapping[str, Any]], registry: Mapping[str, Mapping[str, Any]], output_root: Path, destination_parent: Path, *, all_scope: bool = False) -> str:
    node_id = str(node["node_id"])
    selected = _node_records(node, records)
    generic = node_id not in SUMMARY_NODE_IDS | FOUNDATION_NODE_IDS | COMPACT_NODE_IDS | set(LOCAL_NODE_METRICS)
    suffix = "-all" if all_scope else ""
    anchors = "" if generic else "".join(f'<span id="module-{_slug(question)}{suffix}"></span>' for question in node["questions"] if LEGACY_ANCHOR_NODE.get(question) == node_id)
    # All-scope evidence is its own denominator.  Package narratives describe
    # the parent batch and must not be copied into the HD global article.
    use_narrative = not all_scope
    heading = f'{node_id}. {_node_title(package, node, use_narrative=use_narrative)}'
    opening = f'<article class="report-module node-module" id="node-{_slug(node_id)}{suffix}">{anchors}<h3>{html.escape(heading)}</h3>'
    lead = _narrative_panel(package, node, phase="lead") if use_narrative else ""
    method = f'<p class="module-intro">{_node_method_link(node, destination_parent)}</p>'
    tail = _narrative_panel(package, node, phase="tail") if use_narrative else ""
    if node_id in SUMMARY_NODE_IDS:
        return opening + lead + method + _summary_block(package, node) + tail + "</article>"
    if node_id in FOUNDATION_NODE_IDS:
        facts_by_node = {"1.1": ("carrier", "source_files"), "1.2": ("cohort_audit",), "1.3": ("pairing_audit",), "1.4": ("gene",), "1.5": ("preprocessing",)}
        paths = set()
        for record in selected:
            for key in facts_by_node[node_id]:
                fact = record.get("results", {}).get(key, {})
                if key == "cohort_audit":
                    fact = fact.get("observation", {})
                items = fact.values() if key == "preprocessing" else (fact,)
                for item in items:
                    if isinstance(item, Mapping) and item.get("path"):
                        paths.add(str(item["path"]))
        links = []
        for path in sorted(paths):
            if (output_root / path).is_file():
                label = path.removeprefix("analysis/")
                links.append(f'<a class="csv-link" href="{html.escape(_href(path, output_root, destination_parent), quote=True)}">{html.escape(label)}</a>')
        return opening + lead + method + _foundation_block(node, selected) + '<p class="proof-links">' + " · ".join(links) + "</p>" + tail + "</article>"
    if node_id in LOCAL_NODE_METRICS:
        metric, label = LOCAL_NODE_METRICS[node_id]
        return opening + lead + method + _local_diversity_module(metric=metric, label=label, records=selected, registry=registry, output_root=output_root, destination_parent=destination_parent, all_scope=all_scope) + tail + "</article>"
    if node_id in COMPACT_NODE_IDS:
        paths = sorted({path for record in selected for path, exists in _tables(record) if exists and ((node_id == "3.8" and "diversity_by_anatomy" in path) or (node_id == "4.3" and "state_region_extent" in path))})
        links = [f'<a href="{html.escape(_href(path, output_root, destination_parent), quote=True)}">{html.escape(path.removeprefix("analysis/"))}</a>' for path in paths]
        return opening + lead + method + _compact_block(selected, node_id=node_id) + '<p class="proof-links">' + " · ".join(links) + "</p>" + tail + "</article>"
    blocks = []
    questions = list(node["questions"])
    for index, question in enumerate(questions):
        question_records = [record for record in selected if (record.get("identity", {}) or {}).get("question") == question]
        if question_records or (node_id in {"2.1", "2.2"} and ("node_narratives" in package or node_id in {"2.1", "2.2"})):
            next_question = questions[index + 1] if index + 1 < len(questions) else None
            bridge = f'<p class="module-next">下一问题：<a href="#module-{_slug(next_question)}{suffix}">{html.escape(_labels().get(next_question, (next_question, next_question))[0])}</a></p>' if next_question else ""
            if "node_narratives" in package and node_id == "2.1" and question == "complexity":
                evidence = _complexity_table(question_records, registry, output_root, destination_parent, all_scope=all_scope)
            elif "node_narratives" in package and node_id == "2.2" and question == "matched_k":
                evidence = _matched_k_table(question_records, registry, output_root, destination_parent, all_scope=all_scope)
            else:
                evidence = _evidence_table(layer=str(node["layer"]), question=question, records=question_records, all_scope=all_scope, registry=registry, output_root=output_root, destination_parent=destination_parent)
            blocks.append(f'<section class="legacy-question" id="module-{_slug(question)}{suffix}"><h4>{html.escape(_labels().get(question, (question, question))[0])}</h4>{evidence}{bridge}</section>')
    return opening + lead + method + ("".join(blocks) or '<p class="empty">No saved record</p>') + tail + "</article>"


def _node_placement(node: Mapping[str, Any]) -> str:
    """Read the placement owned by the content contract."""

    return str(node.get("reading_placement", "main"))


def _reference_node_article(package: Mapping[str, Any], node: Mapping[str, Any], records: list[Mapping[str, Any]], registry: Mapping[str, Mapping[str, Any]], output_root: Path, destination_parent: Path, *, all_scope: bool = False) -> str:
    """Keep reference nodes available without repeating a main evidence table."""

    node_id = str(node["node_id"])
    suffix = "-all" if all_scope else ""
    selected = _node_records(node, records)
    use_narrative = not all_scope
    title = f'{node_id}. {_node_title(package, node, use_narrative=use_narrative)}'
    lead = _narrative_panel(package, node, phase="lead") if use_narrative else ""
    tail = _narrative_panel(package, node, phase="tail") if use_narrative else ""
    content = _compact_block(selected, node_id=node_id, language="zh") if selected else '<p class="empty">No saved reference record</p>'
    return (
        f'<details class="reference-details" id="node-{_slug(node_id)}{suffix}">'
        + f'<summary>{html.escape(title)} · reference view</summary>'
        + f'<p class="module-intro">{_node_method_link(node, destination_parent)}</p>'
        + lead
        + content
        + tail
        + "</details>"
    )


def _body_reorganized(package: Mapping[str, Any], records: list[Mapping[str, Any]], registry: Mapping[str, Mapping[str, Any]], output_root: Path, destination_parent: Path) -> str:
    """Keep the web reading order separate from the scientific node tree."""
    return _web_body(package, records, registry, output_root, destination_parent)


def _body(package: Mapping[str, Any], records: list[Mapping[str, Any]], registry: Mapping[str, Mapping[str, Any]], output_root: Path, destination_parent: Path) -> str:
    if "node_narratives" in package:
        return _body_reorganized(package, records, registry, output_root, destination_parent)
    chunks: list[str] = []
    main_records = [record for record in records if (record.get("identity", {}) or {}).get("scope") != "All"]
    all_records = [record for record in records if (record.get("identity", {}) or {}).get("scope") == "All"]
    for layer in LAYERS:
        chunks.append(f'<section id="layer-{_slug(layer)}" class="body-layer"><h2>{html.escape(TITLES[layer][0])}</h2><p class="transition">{html.escape(LAYER_INTROS[layer])}</p>')
        for node in (item for item in NODE_SPECS if item["layer"] == layer):
            chunks.append(_node_article(package, node, main_records, registry, output_root, destination_parent))
        chunks.append("</section>")
    if all_records:
        chunks.append('<section id="hd-all" class="hd-all"><h2>HD All：独立 cohort 与分母</h2><p class="transition">HD All 的总体分母单独展示；它不进入 parent-internal 三列。</p>')
        for question in _question_order(all_records):
            subset = [record for record in all_records if (record.get("identity", {}) or {}).get("question") == question]
            chunks.append(_module(layer=_record_layer(subset[0]), question=question, records=subset, all_scope=True, registry=registry, output_root=output_root, destination_parent=destination_parent))
        chunks.append("</section>")
    summary = next(node for node in NODE_SPECS if node["node_id"] == "summary")
    chunks.append(f'<section id="summary" class="body-layer">{_node_article(package, summary, main_records, registry, output_root, destination_parent)}</section>')
    return "".join(chunks)


# These are reading positions, not new scientific nodes or analyses.
WEB_READING_PLAN = (
    ("layer-overall", "总体差异", (
        ("2.1", "分群数量与整体一致性"), ("2.2", "控制群数后的成员归属"),
        ("2.3", "跨细胞类型对照"), ("2.4", "基因提供范围与空间结构"), ("2.5", "EMT 覆盖与评分"),
    )),
    ("layer-local", "局部邻域差异", (
        ("3.4", "比较范围与窗口支持"), ("local-baselines", "两个 Raw 对照下的邻域状态"),
        ("3.5", "状态数量 · Kobs 证据"), ("3.6", "有效状态数 · Neff 证据"),
        ("3.7", "均衡程度 · evenness 证据"), ("4.4", "局部结果的尺度敏感性"),
    )),
    ("layer-localization", "空间定位", (
        ("3.1", "组织背景"), ("partition-maps", "分群与重建状态地图"),
        ("3.2", "成员归属变化的位置"), ("3.3", "EMT 变化的空间分布"),
        ("3.8", "不同组织区域的局部差异"), ("4.1", "State Region 可识别性"),
        ("4.2", "连续状态场与 State Region"), ("4.3", "State Region 的区域覆盖"),
        ("gain-audit", "Gain Region · 独立诊断"),
    )),
    ("summary", "综合解释", (
        ("2.6", "总体差异归纳"), ("3.9", "局部与空间结果归纳"),
        ("4.5", "Region 解释范围"), ("summary", "跨类型综合解释"),
    )),
)
WEB_NODE_TITLES = {key: title for _, _, items in WEB_READING_PLAN for key, title in items}
WEB_BASELINES = (("local_vs_raw_leiden", "Raw Leiden"), ("local_vs_raw_level2", "Raw Level2"))


def _web_anchor(key: str) -> str:
    return "node-" + _slug(key) if key in {str(n["node_id"]) for n in NODE_SPECS} else key


def _web_sorted(records: Iterable[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    order = {scope: index for index, scope in enumerate((*SCOPES, "All"))}
    return sorted(records, key=lambda r: order.get(r.get("identity", {}).get("scope"), len(order)))


def _web_unmatched(record: Mapping[str, Any]) -> bool:
    result = record.get("results", {})
    condition = record.get("evidence_condition", {})
    return result.get("matched") is False or result.get("headline_eligible") is False or condition.get("matched_k_status") == "unmatched_cluster_complexity"


def _web_sentence(record: Mapping[str, Any], metric: str | None = None) -> str:
    """Describe saved values without computing a new effect or judgment."""
    q = record.get("identity", {}).get("question", "")
    r = record.get("results", {})
    raw, recon = r.get("raw", {}) or {}, r.get("reconstruction", {}) or {}
    if q == "complexity":
        kind, k = _complexity_kind(r)
        context = "Raw 参考分群 → Recon 最终分群" if "fixed final" in kind else "Raw 分辨率下的两侧分群"
        return f"{context}：群数 {_number(r.get('raw_k'))} → {_number(k)}，该对照的 ARI 为 {_number(r.get('ari'))}。"
    if q == "matched_k":
        status = "控制 K 条件未成立，以下为诊断值" if _web_unmatched(record) else _judgment(record)
        return f"Raw / Recon K = {_number(r.get('raw_k'))} / {_number(r.get('reconstruction_k'))}；{status}。成员改变率 {_number(r.get('unit_change_fraction'), True)}，ARI {_number(r.get('ari'))}。"
    if q == "moran_all_valid":
        return f"两侧各自有效基因 {_number(raw.get('n_valid'))} → {_number(recon.get('n_valid'))}；Moran Q75 {_number(raw.get('q75'))} → {_number(recon.get('q75'))}。"
    if q == "moran_shared_valid":
        return f"在 {_number(raw.get('n_valid'))} 个共同有效基因上，Moran 的配对 ΔI 中位数为 {_number(r.get('paired_delta', {}).get('median'))}。"
    if q == "emt_coverage":
        return f"EMT 可用基因 {_number(raw.get('available_gene_count'))}/{_number(raw.get('resource_gene_count'))} → {_number(recon.get('available_gene_count'))}/{_number(recon.get('resource_gene_count'))}。"
    if q == "emt_score":
        return f"评分中位数 {_number(raw.get('median'))} → {_number(recon.get('median'))}；配对差值中位数 {_number(r.get('paired', {}).get('median_delta'))}。覆盖与排名截断条件见上表。"
    if q.startswith("local_vs_"):
        values = _local_metric_values(record, metric or "Neff")
        baseline = "Raw Leiden" if q == "local_vs_raw_leiden" else "Raw Level2"
        return f"{baseline} → Recon：{metric or 'Neff'} {_number(values.get('baseline_median'))} → {_number(values.get('reconstruction_median'))}，配对 Δ 中位数 {_number(values.get('delta_median'))}。"
    if q == "window_support":
        return f"主窗口尺度 {_number(r.get('main_window_side_um'))} μm，有效窗口 {_number(r.get('n_valid_windows'))} 个；每窗最低支持 {_number(r.get('min_parent_units'))} 个单元。"
    if q == "region_extent":
        e = r.get("overall", {})
        return f"State Region 覆盖有效窗口面积的 {_number(e.get('area_fraction'), True)}，面积 {_number(e.get('region_area_mm2'))} mm²。该范围描述重建后的状态分布，不表示改善面积。"
    if q == "changed_units":
        e = r.get("overall") or next((x for x in r.get("by_region", []) if x.get("level1_region", x.get("level1")) == "Overall"), {})
        total = e.get("paired_units") if e.get("paired_units") is not None else e.get("total_units")
        condition = "控制 K 条件未成立，此处仅作诊断。" if _web_unmatched(record) else ""
        stratification = "Raw Level1 类型" if record.get("identity", {}).get("scope") == "All" else "组织区域"
        return f"{condition}配对范围内 {_number(e.get('changed_units'))}/{_number(total)} 个单元的归属发生变化（{_number(e.get('change_fraction'), True)}）；按{stratification}列出的数量与比例见表。"
    return _brief_conclusion(record, "zh", limit=360)


def _web_result_list(records: Iterable[Mapping[str, Any]], metric: str | None = None) -> str:
    items = []
    for record in _web_sorted(records):
        scope = _display_scope(record.get("identity", {}).get("scope", ""))
        items.append(f'<li><strong>{html.escape(scope)}</strong> · {html.escape(_web_sentence(record, metric))}</li>')
    return '<ul class="scope-results">' + "".join(items) + "</ul>" if items else ""


def _web_links(records: Iterable[Mapping[str, Any]], output_root: Path, destination_parent: Path) -> str:
    paths = {path: exists for record in records for path, exists in _tables(record)}
    if not paths:
        return ""
    links = []
    for path, exists in sorted(paths.items()):
        label = path.removeprefix("analysis/") + (" · 文件缺失" if not exists else "")
        links.append(f'<a class="csv-link" href="{html.escape(_href(path, output_root, destination_parent), quote=True)}">{html.escape(label)}</a>')
    return '<details class="source-details"><summary>查看本问题的数据与条件</summary><div class="proof-links">' + " · ".join(links) + "</div></details>"


def _web_values_table(question: str, records: Iterable[Mapping[str, Any]]) -> str:
    records = _web_sorted(records)
    if not records:
        return '<p class="empty">未保存该问题的结果。</p>'
    rows: list[tuple[Any, ...]] = []
    headers: tuple[str, ...] = ()
    for record in records:
        scope = _display_scope(record["identity"]["scope"])
        r = record.get("results", {})
        raw, recon = r.get("raw", {}) or {}, r.get("reconstruction", {}) or {}
        if question == "complexity":
            kind, target = _complexity_kind(r)
            fixed = "fixed final" in kind
            headers = ("细胞类型", "Raw 参考 K" if fixed else "Raw K", "Recon 最终 K" if fixed else "同分辨率 Recon K", "此组对照 ARI")
            rows.append((scope, r.get("raw_k"), target, r.get("ari")))
        elif question == "matched_k":
            headers = ("细胞类型", "Raw / Recon K", "控制条件", "成员改变率", "1 − macro-F1", "此组对照 ARI")
            rows.append((scope, f"{_number(r.get('raw_k'))} / {_number(r.get('reconstruction_k'))}", "未成立 · 诊断值" if _web_unmatched(record) else _judgment(record), _number(r.get("unit_change_fraction"), True), r.get("balanced_change"), r.get("ari")))
        elif question.startswith("moran_"):
            headers = ("细胞类型", "Raw 有效基因", "Recon 有效基因", "Raw Q75", "Recon Q75")
            row = (scope, raw.get("n_valid"), recon.get("n_valid"), raw.get("q75"), recon.get("q75"))
            if question == "moran_shared_valid":
                headers = ("细胞类型", "Raw 共同有效数", "Recon 共同有效数", "Raw Q75", "Recon Q75", "配对 ΔI 中位数")
                row += (r.get("paired_delta", {}).get("median"),)
            rows.append(row)
        elif question == "emt_coverage":
            headers = ("细胞类型", "Raw 可用 / 资源基因", "Recon 可用 / 资源基因")
            rows.append((scope, f"{_number(raw.get('available_gene_count'))} / {_number(raw.get('resource_gene_count'))}", f"{_number(recon.get('available_gene_count'))} / {_number(recon.get('resource_gene_count'))}"))
        elif question == "emt_score":
            cutoff = r.get("cutoff_audit", {})
            headers = ("细胞类型", "Raw 中位数", "Recon 中位数", "配对 Δ 中位数", "配对单元数", "Raw → Recon 排名长度")
            rows.append((scope, raw.get("median"), recon.get("median"), r.get("paired", {}).get("median_delta"), r.get("paired", {}).get("n"), " → ".join(_number(cutoff.get(side, {}).get("effective_rank_length")) for side in ("raw", "reconstruction"))))
        elif question == "window_support":
            headers = ("细胞类型", "窗口尺度 μm", "有效窗口", "最低单元数 / 窗", "配对抽样次数")
            rows.append((scope, r.get("main_window_side_um"), r.get("n_valid_windows"), r.get("min_parent_units"), r.get("rarefaction_draws")))
        elif question == "threshold_reliability":
            state = r.get("state", {})
            headers = ("细胞类型", "可识别状态", "阈值", "有效 bootstrap", "CI 宽度 / Neff 范围", "窗口数")
            rows.append((scope, _judgment(record), state.get("threshold"), _number(state.get("valid_bootstrap_fraction"), True), _number(state.get("relative_ci_width"), True), state.get("n_windows")))
        elif question == "region_extent":
            e = r.get("overall", {})
            headers = ("细胞类型", "Region / 有效窗口", "面积 mm²", "有效面积覆盖", "Region / 有效单元", "单元覆盖")
            rows.append((scope, f"{_number(e.get('region_windows'))} / {_number(e.get('valid_windows'))}", e.get("region_area_mm2"), _number(e.get("area_fraction"), True), f"{_number(e.get('region_units'))} / {_number(e.get('valid_units'))}", _number(e.get("unit_fraction"), True)))
        elif question == "changed_units":
            if scope == "All":
                headers = ("Raw Level1 类型", "改变 / 配对单元", "改变比例", "95% Wilson 区间", "低样本量标记")
                for e in r.get("by_region", []):
                    denominator = e.get("paired_units") if e.get("paired_units") is not None else e.get("total_units")
                    flag = {True: "是", False: "否"}.get(e.get("low_sample_size"), "NA")
                    rows.append((e.get("level1", e.get("level1_region")), f"{_number(e.get('changed_units'))} / {_number(denominator)}", _number(e.get("change_fraction"), True), f"{_number(e.get('wilson_ci_lower'), True)}–{_number(e.get('wilson_ci_upper'), True)}", flag))
                continue
            headers = ("细胞类型", "组织区域", "改变 / 配对单元", "改变比例", "有效窗口", "解释条件")
            for e in r.get("by_region", []):
                denominator = e.get("paired_units") if e.get("paired_units") is not None else e.get("total_units")
                status = "仅作诊断" if record.get("evidence_condition", {}).get("matched_k_status") == "unmatched_cluster_complexity" else "归属差异定位"
                rows.append((scope, e.get("level1_region"), f"{_number(e.get('changed_units'))} / {_number(denominator)}", _number(e.get("change_fraction"), True), e.get("n_valid_windows"), status))
        elif question == "anatomy":
            headers = ("细胞类型上下文", "组织区域", "Raw 背景面积占比")
            rows.extend((scope, e.get("level1_region"), _number(e.get("area_fraction"), True)) for e in r.get("regions", []))
        elif question == "emt_spatial_fields":
            headers = ("细胞类型", "空间单元数", "坐标单位")
            rows.append((scope, r.get("n_units"), ", ".join(r.get("coordinate_units", []))))
        elif question == "scale_sensitivity":
            headers = ("细胞类型", "尺度 μm", "主尺度", "ΔNeff / Raw Leiden", "ΔNeff / Raw Level2", "有效窗口", "最低单元数 / 窗")
            for e in r.get("rows", []):
                scale = e.get("window_side_length")
                rows.append((scope, scale, "主尺度" if scale == r.get("main_scale_um") else "", e.get("median_delta_neff_vs_raw_leiden"), e.get("median_delta_neff_vs_raw_level2"), e.get("n_valid_windows"), e.get("min_parent_units")))
    return _foundation_table(headers, rows) if headers else ""


def _web_figure_registry(records: list[Mapping[str, Any]], registry: Mapping[str, Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
    result = {}
    old_rows: set[str] = set()
    for path, original in registry.items():
        info = dict(original)
        all_only = set(info.get("scopes", ())) == {"All"}
        suffix = "-all" if all_only else ""
        name = Path(path).name.lower()
        occurrences = [r for r in records if path in _figure_paths(r)]
        questions = {r["identity"]["question"] for r in occurrences}
        question = _preferred_figure_question(path, questions) or occurrences[0]["identity"]["question"]
        metric = _metric_for_figure(path)
        if name.endswith(("_reconstructed_clusters.png", "_matched_clusters.png")):
            owner = "partition-maps" + suffix
        elif metric and question.startswith("local_vs_"):
            owner = "node-" + {"kobs": "3-5", "neff": "3-6", "evenness": "3-7"}[metric] + suffix
        else:
            owner = "module-" + _slug(question) + suffix
        info["web_owner"] = owner
        info["records"] = occurrences
        old_row = str(info.get("owner_row", ""))
        info["legacy_row"] = old_row if old_row not in old_rows else ""
        old_rows.add(old_row)
        result[path] = info
    return result


def _web_figures(owner: str, registry: Mapping[str, Mapping[str, Any]], output_root: Path, destination_parent: Path) -> str:
    parts, aliases = [], set()
    for path, info in registry.items():
        if info.get("web_owner") != owner:
            continue
        old_row = str(info.get("legacy_row", ""))
        if old_row and old_row not in aliases:
            parts.append(f'<span id="{html.escape(old_row, quote=True)}"></span>')
            aliases.add(old_row)
        records = info.get("records", [])
        name = Path(path).name.lower()
        metric = _metric_for_figure(path)
        if name.endswith("_reconstructed_clusters.png"):
            explanation = '<p>最终 Recon 分群的空间分布，用于理解后续邻域分析中的状态。该图只显示 Recon；群数变化请看总体层的两侧数值对照。</p>'
        else:
            preferred = _preferred_figure_question(path, {r["identity"]["question"] for r in records})
            selected = records if metric else [r for r in records if r["identity"]["question"] == preferred]
            # Scope-specific plots stay with their scope; shared plots retain All.
            chosen = selected or records
            metric_label = {"kobs": "Kobs", "neff": "Neff", "evenness": "evenness"}.get(metric)
            parents = [r for r in chosen if r["identity"].get("scope") != "All"]
            globals_ = [r for r in chosen if r["identity"].get("scope") == "All"]
            explanation = _web_result_list(parents, metric_label)
            if globals_:
                explanation += '<div class="global-figure-note"><p><a href="#hd-all">HD All · 独立范围与分母</a></p>' + _web_result_list(globals_, metric_label) + '</div>'
        scopes = " / ".join(_display_scope(s) for s in (*SCOPES, "All") if s in info.get("scopes", ()))
        parts.append('<div class="figure-entry"><div class="figure-explanation"><h4>' + html.escape(scopes or "共同证据") + '</h4>' + explanation + '</div>' + _render_figure(path, info, output_root, destination_parent) + '</div>')
    return "".join(parts)


def _web_tail(package: Mapping[str, Any], node: Mapping[str, Any]) -> str:
    narrative = _node_narrative(package, str(node["node_id"]))
    if not _narrative_reviewed(package, narrative):
        return ""
    items = _narrative_interpretations(narrative or {}, "zh")
    return '<ul class="narrative-interpretations">' + "".join(f'<li>{html.escape(s)}</li>' for s in items) + '</ul>' if items else ""


def _web_gene_range(records: list[Mapping[str, Any]]) -> str:
    rows = []
    notes = []
    by_scope = _record_map(records)
    for scope in (*SCOPES, "All"):
        foundation = by_scope.get(("foundation", scope))
        if foundation is None:
            continue
        gene = foundation.get("results", {}).get("gene", {}).get("by_side", {})
        moran = by_scope.get(("moran_all_valid", scope), {}).get("results", {})
        rows.append((_display_scope(scope), gene.get("raw", {}).get("provided_n"), gene.get("reconstruction", {}).get("provided_n"), moran.get("raw", {}).get("n_valid"), moran.get("reconstruction", {}).get("n_valid")))
        if any(gene.get(side, {}).get("provided_n") is None for side in ("raw", "reconstruction")):
            notes.append(f'{_display_scope(scope)} 的侧别基因提供审计未保存，提供数保留 NA。')
            carrier_rows = foundation.get("results", {}).get("carrier", {}).get("rows", [])
            carrier_counts = {r.get("role"): r.get("n_genes") for r in carrier_rows}
            if "Raw full context" in carrier_counts and "Reconstruction full carrier" in carrier_counts:
                notes.append(f'已有完整载体的基因范围为 Raw {_number(carrier_counts["Raw full context"])} → Recon {_number(carrier_counts["Reconstruction full carrier"])}；这是载体范围。')
    note_html = '<p class="compact-note">' + html.escape(" ".join(notes)) + '</p>' if notes else ""
    return '<section class="web-question" id="gene-availability' + ('-all' if rows and all(row[0] == "All" for row in rows) else '') + '"><h4>输入提供范围与 Moran 可计算范围</h4><p>提供数量反映两侧输入的基因范围；有效数量反映当前分析中可计算 Moran 的范围。Recon 值保留其重建或投影视图含义。</p>' + note_html + _foundation_table(("细胞类型", "Raw 提供基因", "Recon 提供基因", "Raw Moran 有效基因", "Recon Moran 有效基因"), rows) + '</section>'


def _web_dual_baselines(records: list[Mapping[str, Any]]) -> str:
    parts = ['<article class="report-module" id="local-baselines"><h3>两个 Raw 对照下的邻域状态</h3><p>分别比较同一 Recon 与两个 Raw baseline。表中的 Δ 是已保存的配对差值中位数。</p><div class="baseline-grid">']
    for q, label in WEB_BASELINES:
        selected = _web_sorted(r for r in records if r["identity"]["question"] == q)
        parts.append(f'<section class="baseline-panel" id="module-{q}"><h4>相对 {label}</h4>')
        for r in selected:
            result = r.get("results", {})
            scope = _display_scope(r["identity"]["scope"])
            parts.append(f'<div class="baseline-scope"><h5>{html.escape(scope)}</h5><p>Neff 判读：{html.escape(_judgment(r))}</p>')
            rows = []
            for _, metric in LOCAL_DIVERSITY_METRICS:
                value = _local_metric_values(r, metric)
                rows.append((metric, value.get("baseline_median"), value.get("reconstruction_median"), value.get("delta_median")))
            parts.append(_foundation_table(("指标", "Raw 中位数", "Recon 中位数", "配对 Δ 中位数"), rows))
            parts.append(f'<p class="compact-note">有效窗口 {_number(result.get("n_valid_windows"))} · 尺度 {_number(result.get("scale_um"))} μm</p></div>')
        if not selected:
            parts.append('<p class="empty">未保存该 baseline 的比较结果。</p>')
        parts.append('</section>')
    return ''.join(parts) + '</div></article>'


def _web_question(question: str, records: list[Mapping[str, Any]], registry: Mapping[str, Mapping[str, Any]], output_root: Path, destination_parent: Path, all_scope: bool = False) -> str:
    selected = _web_sorted(r for r in records if r["identity"]["question"] == question)
    if not selected:
        return ""
    owner = "module-" + _slug(question) + ("-all" if all_scope else "")
    label = {
        "complexity": "原分群数量与一致性", "matched_k": "控制 K 后的归属对照",
        "moran_all_valid": "两侧各自全部有效基因", "moran_shared_valid": "共同有效同名基因",
        "emt_coverage": "可用资源基因", "emt_score": "评分与配对变化",
    }.get(question, QUESTION_LABELS.get(question, (question, question))[0])
    if question == "changed_units" and all_scope:
        label = "按 Raw Level1 类型分层的成员改变"
    note = ""
    if question == "complexity":
        fixed = any("fixed final" in _complexity_kind(r.get("results", {}))[0] for r in selected)
        note = "本路线比较 Raw 参考分群与 Recon 固定最终分群。下一个问题使用另一组控制 K 的对照。" if fixed else "两侧使用 Raw 分辨率作群数与 ARI 对照；下一个问题另行检查控制 K 后的成员变化。"
    elif question == "matched_k":
        note = "成员改变率按配对单元统计；balanced change 为 1 − macro-F1。控制条件未成立的行保留诊断值，其 ARI 不与上一组分群对照混用。"
    elif question == "emt_score":
        note = "结合上方基因覆盖及实际排名长度解释评分变化；配对 Δ 中位数来自已保存的单元配对结果。"
    elif question == "moran_shared_valid":
        note = "同名且两侧均可计算的基因逐一配对；此处的基因范围与上一表两侧各自的有效范围不同。"
    content = f'<section class="web-question" id="{owner}"><h4>{html.escape(label)}</h4>'
    if note:
        content += f'<p class="compact-note">{html.escape(note)}</p>'
    if all_scope:
        content += _web_result_list(selected)
    if question in {"emt_score", "emt_spatial_fields"}:
        views = {
            "cluster_mean_projection": "Recon 使用 cluster-mean 投影视图",
            "reconstructed_expression_spatial": "Recon 使用重建表达的空间视图",
        }
        conditions = list(dict.fromkeys(str(r.get("results", {}).get("carrier_condition")) for r in selected if r.get("results", {}).get("carrier_condition")))
        if conditions:
            content += '<p class="compact-note">表达条件：' + html.escape("；".join(views.get(value, value) for value in conditions)) + '。评分变化同时受资源覆盖和排名截断影响。</p>'
    content += _web_values_table(question, selected)
    if question in {"complexity", "matched_k"} and not all_scope:
        content += '<p class="compact-note"><a href="#partition-maps">查看 Raw / Recon 分群与重建状态地图</a></p>'
    figures = _web_figures(owner, registry, output_root, destination_parent)
    if figures:
        content += figures
    else:
        refs = {path for r in selected for path in _figure_paths(r)}
        links = []
        for path in sorted(refs):
            if path not in registry:
                continue
            info = registry[path]
            scopes = " / ".join(_display_scope(s) for s in (*SCOPES, "All") if s in info.get("scopes", ()))
            label = str(info["title"]) + (f" · {scopes}" if scopes else "")
            links.append(f'<a class="figure-reference" href="#{html.escape(str(info["figure_id"]), quote=True)}">{html.escape(label)}</a>')
        if links:
            content += '<p class="compact-note">相关图形证据</p>' + ''.join(links)
    return content + _web_links(selected, output_root, destination_parent) + '</section>'


def _web_node(package: Mapping[str, Any], node_id: str, records: list[Mapping[str, Any]], registry: Mapping[str, Mapping[str, Any]], output_root: Path, destination_parent: Path, all_scope: bool = False) -> str:
    node = next(n for n in NODE_SPECS if n["node_id"] == node_id)
    selected = _web_sorted(_node_records(node, records))
    suffix = "-all" if all_scope else ""
    if all_scope and (not selected or node_id in SUMMARY_NODE_IDS):
        return ""
    if node_id in FOUNDATION_NODE_IDS:
        return _node_article(package, node, records, registry, output_root, destination_parent, all_scope=all_scope)
    title = WEB_NODE_TITLES.get(node_id, str(node["title_zh"]))
    if all_scope and node_id == "3.2":
        title = "不同 Raw Level1 类型的成员改变"
    anchor = _web_anchor(node_id) + suffix
    if node_id == "2.3":
        return f'<details class="reference-details" id="{anchor}"><summary>{html.escape(title)}</summary>' + _compact_block(selected, node_id, "zh") + _web_links(selected, output_root, destination_parent) + '</details>'
    content = f'<article class="report-module" id="{anchor}"><h3>{html.escape(title)}</h3>'
    if node_id == "4.1" and not all_scope:
        content += '<span id="layer-region"></span><p>先看当前阈值能否识别，再阅读连续状态与覆盖。State Region 表示重建后的高多样性状态范围，不表示改善区域。</p>'
    if not all_scope:
        content += _narrative_lead_panel(package, node)
    if node_id in SUMMARY_NODE_IDS:
        content += _summary_block(package, node)
    elif node_id in LOCAL_NODE_METRICS:
        metric, _ = LOCAL_NODE_METRICS[node_id]
        content += f'<span id="module-local-diversity-{metric}{suffix}"></span>'
        content += '<p class="compact-note">两种 Raw baseline 的完整证据 panel 在此共同展示；上方并列表保留各自的配对差异。</p>'
        content += _web_figures(anchor, registry, output_root, destination_parent)
        content += _web_links(selected, output_root, destination_parent)
    elif node_id in {"3.8", "4.3"}:
        if node_id == "3.8":
            selected = [r for r in selected if r["identity"]["question"] in {q for q, _ in WEB_BASELINES}]
        content += _compact_block(selected, node_id, "zh") + _web_links(selected, output_root, destination_parent)
    else:
        if node_id == "2.4":
            content += _web_gene_range(records)
        for question in node["questions"]:
            content += _web_question(question, selected, registry, output_root, destination_parent, all_scope)
    if not all_scope:
        content += _web_tail(package, node)
    content += f'<p class="compact-note">{_node_method_link(node, destination_parent)}</p>'
    return content + '</article>'


def _web_csv_number(value: str | None) -> Any:
    if value is None or value.strip().lower() in {"", "nan", "na", "null", "none"}:
        return None
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except ValueError:
        return value


def _web_gain(output_root: Path, destination_parent: Path) -> str:
    content = ['<article class="report-module gain-audit" id="gain-audit"><h3>Gain Region · 独立诊断</h3><p>这里保留相对 Raw Leiden 的局部正向 ΔNeff 分布及其阈值诊断，与上方 Recon 状态范围分别阅读；不作为综合收益判断。</p>']
    for scope in SCOPES:
        folder = output_root / "analysis" / "local_state" / scope
        threshold_path = folder / "gain_threshold_audit.json"
        extent_path = folder / "gain_region_extent_audit.csv"
        bootstrap_path = folder / "gain_threshold_bootstrap_audit.csv"
        content.append(f'<section><h4>{html.escape(_display_scope(scope))}</h4>')
        if not threshold_path.is_file() or not extent_path.is_file():
            content.append('<p class="empty">Gain 审计文件未保存或缺失。</p></section>')
            continue
        try:
            saved = json.loads(threshold_path.read_text(encoding="utf-8"))
            threshold = saved[0] if isinstance(saved, list) and saved else saved
            with extent_path.open(encoding="utf-8", newline="") as handle:
                extent = list(csv.DictReader(handle))
        except (OSError, ValueError, csv.Error) as exc:
            content.append(f'<p class="empty">Gain 审计读取失败：{html.escape(str(exc))}</p></section>')
            continue
        if not isinstance(threshold, Mapping):
            content.append('<p class="empty">Gain 阈值记录格式无法识别。</p></section>')
            continue
        state = str(threshold.get("status", "未保存状态"))
        content.append(f'<p>阈值状态：{html.escape(state)}；阈值 {_number(threshold.get("threshold"))}；bootstrap 区间 {_number(threshold.get("ci_lower"))}–{_number(threshold.get("ci_upper"))}。</p>')
        content.append(f'<p class="compact-note">阈值审计窗口 {_number(threshold.get("n_windows"))} 个；有效 bootstrap {_number(threshold.get("n_valid_bootstrap"))} 次。覆盖的分母见下表。</p>')
        rows = []
        for row in extent:
            available = str(row.get("region_available", "")).lower() == "true"
            def measure(key: str) -> Any:
                return _web_csv_number(row.get(key)) if available else None
            rows.append((row.get("level1_region"), f"{_number(measure('region_windows'))} / {_number(_web_csv_number(row.get('valid_windows')))}", measure("region_area_mm2"), _number(measure("area_fraction"), True), f"{_number(measure('region_units'))} / {_number(_web_csv_number(row.get('valid_units')))}", _number(measure("unit_fraction"), True)))
        content.append(_foundation_table(("组织区域", "Gain / 有效窗口", "面积 mm²", "有效面积覆盖", "Gain / 有效单元", "单元覆盖"), rows))
        links = [f'<a href="{html.escape(_href(str(path), output_root, destination_parent), quote=True)}">{html.escape(path.name)}</a>' for path in (threshold_path, extent_path, bootstrap_path) if path.is_file()]
        content.append('<details class="source-details"><summary>Gain 原始诊断记录</summary>' + ' · '.join(links) + '</details></section>')
    return ''.join(content) + '</article>'


def _web_body(package: Mapping[str, Any], records: list[Mapping[str, Any]], registry: Mapping[str, Mapping[str, Any]], output_root: Path, destination_parent: Path) -> str:
    main_records = [r for r in records if r["identity"].get("scope") != "All"]
    all_records = [r for r in records if r["identity"].get("scope") == "All"]
    parts = []
    for section, title, items in WEB_READING_PLAN:
        parts.append(f'<section id="{section}" class="reading-layer"><h2>{html.escape(title)}</h2>')
        for key, label in items:
            if key == "local-baselines":
                parts.append(_web_dual_baselines(main_records))
            elif key == "partition-maps":
                parts.append('<article class="report-module" id="partition-maps"><h3>分群与重建状态地图</h3><p>先看 Raw / Recon 分群在组织中的分布，再结合下一节的改变比例与区域数量定位归属变化。</p>' + _web_figures(key, registry, output_root, destination_parent) + '</article>')
            elif key == "gain-audit":
                parts.append(_web_gain(output_root, destination_parent))
            else:
                parts.append(_web_node(package, key, main_records, registry, output_root, destination_parent))
        parts.append('</section>')
    if all_records:
        parts.append('<section id="hd-all" class="hd-all reading-layer"><h2>HD All · 独立全局分析</h2><p>All 使用独立 cohort 与分母。当前保存了分群、Moran 和成员改变结果；这里仅展开这些已存在的分析。</p>')
        foundation = next((r.get("results", {}) for r in all_records if r["identity"]["question"] == "foundation"), {})
        selection, denominator = foundation.get("selection", {}), foundation.get("denominator", {})
        if foundation:
            parts.append(f'<p>Raw 总范围 {_number(selection.get("raw_parent_n"))} 个单元，其中 Recon 覆盖 {_number(selection.get("reconstruction_covered_parent_n"))} 个；All 独立抽样 {_number(selection.get("sampled_n"))} 个。分群配对 {_number(denominator.get("paired_units"))} 个，Moran 使用 {_number(denominator.get("moran_units"))} 个。</p><p class="compact-note"><a href="#foundation-audit-hd">查看 All 的比较基础与预处理</a></p>')
        for key in ("2.1", "2.2", "2.3", "2.4", "3.2"):
            parts.append(_web_node(package, key, all_records, registry, output_root, destination_parent, True))
        if any(info.get("web_owner") == "partition-maps-all" for info in registry.values()):
            parts.append('<article id="partition-maps-all"><h3>All 分群地图</h3>' + _web_figures("partition-maps-all", registry, output_root, destination_parent) + '</article>')
        parts.append('</section>')
    parts.append('<details class="foundation-audit reading-layer" id="foundation-audit"><summary>比较基础 · 输入、配对、表达范围与预处理</summary><section id="layer-foundation" class="audit-layer">')
    for key in ("1.1", "1.2", "1.3", "1.4", "1.5", "1.6"):
        parts.append(_web_node(package, key, main_records, registry, output_root, destination_parent))
    parts.append('</section></details>')
    if all_records:
        parts.append('<details class="foundation-audit" id="foundation-audit-hd"><summary>HD All 比较基础 · 独立范围与分母</summary><section id="layer-foundation-hd">')
        for key in ("1.1", "1.2", "1.3", "1.4", "1.5"):
            parts.append(_web_node(package, key, all_records, registry, output_root, destination_parent, True))
        parts.append('</section></details>')
    return ''.join(parts)


def _web_nav(records: list[Mapping[str, Any]]) -> str:
    parts = ['<ul class="toc-list"><li><a href="#overview-main">样本总览</a></li>']
    for section, title, items in WEB_READING_PLAN:
        parts.append(f'<li class="toc-group"><a href="#{section}">{html.escape(title)}</a><ul>')
        for key, label in items:
            parts.append(f'<li><a href="#{_web_anchor(key)}">{html.escape(label)}</a>')
            if key == "local-baselines":
                parts.append('<ul>' + ''.join(f'<li><a href="#module-{q}">{baseline} 对照</a></li>' for q, baseline in WEB_BASELINES) + '</ul>')
            parts.append('</li>')
        parts.append('</ul></li>')
    if any(r["identity"].get("scope") == "All" for r in records):
        parts.append('<li><a href="#hd-all">HD All · 独立全局分析</a></li>')
    parts.append('<li><a href="#foundation-audit">比较基础</a></li><li><a href="#evidence-index">数据与原图入口</a></li></ul>')
    return ''.join(parts)


def _web_overview(records: list[Mapping[str, Any]]) -> str:
    labels = {"partition": "分群", "moran": "基因信息与 Moran", "emt": "EMT", "local": "局部双 baseline", "region": "State Region"}
    rows = ['<section class="overview-main" id="overview-main"><h3>五项固定对照</h3><div class="table-scroll"><table class="overview-table"><thead><tr><th>比较问题</th>']
    rows.extend(f'<th>{html.escape(_display_scope(scope))}</th>' for scope in SCOPES)
    rows.append('</tr></thead><tbody>')
    for group, _, _, questions, node in OVERVIEW_GROUPS:
        anchor = 'local-baselines' if group == 'local' else _web_anchor(node)
        rows.append(f'<tr class="overview-question" id="overview-group-{group}"><th><a href="#{anchor}">{labels[group]}</a></th>')
        for scope in SCOPES:
            grouped = _overview_group_records(records, questions, scope)
            lines = []
            if group == 'partition':
                if 'complexity' in grouped:
                    r = grouped['complexity']['results']; kind, target = _complexity_kind(r)
                    label = 'Raw 参考 → Recon 最终 K' if 'fixed final' in kind else '同分辨率 K'
                    lines.append(f'{label}：{_number(r.get("raw_k"))} → {_number(target)}')
                if 'matched_k' in grouped:
                    r = grouped['matched_k']['results']
                    status = '控制条件未成立 · 诊断' if _web_unmatched(grouped['matched_k']) else '控制 K 后'
                    lines.append(f'{status}：{_number(r.get("raw_k"))}/{_number(r.get("reconstruction_k"))}，成员改变 {_number(r.get("unit_change_fraction"), True)}')
            elif group == 'moran':
                foundation = _record_map(records).get(('foundation', scope), {}).get('results', {})
                gene = foundation.get('gene', {}).get('by_side', {})
                lines.append('提供基因：' + ' → '.join(_number(gene.get(side, {}).get('provided_n')) for side in ('raw', 'reconstruction')))
                _, value = _overview_group_metric(group, grouped); lines.append(value)
            elif group == 'local':
                for q, baseline in WEB_BASELINES:
                    if q in grouped:
                        value = _local_metric_values(grouped[q], 'Neff')
                        lines.append(f'{baseline}：ΔNeff {_number(value.get("delta_median"))}；{_judgment(grouped[q])}')
            else:
                label, value = _overview_group_metric(group, grouped)
                lines.append(f'{label}：{value}')
                if group == 'region' and 'threshold_reliability' in grouped:
                    lines.append(_judgment(grouped['threshold_reliability']))
                elif group == 'emt':
                    lines.append('覆盖与评分条件分别解释')
            rows.append(f'<td class="overview-cell" data-scope="{html.escape(_display_scope(scope), quote=True)}">' + ''.join(f'<span class="overview-headline">{html.escape(line)}</span>' for line in lines) + '</td>')
        rows.append('</tr>')
    rows.append('</tbody></table></div></section>')
    if any(r['identity'].get('scope') == 'All' for r in records):
        rows.append('<p id="overview-hd-all" class="compact-note"><a href="#hd-all">HD All 的已有全局分析单独展示</a>，使用独立范围与分母。</p>')
    return ''.join(rows)


def _template_path() -> Path:
    return Path(__file__).resolve().parents[3] / "docs" / "design" / "reconstruction-impact" / "templates" / "report.html"


def _render_shell(*, title: str, identity: str, nav: str, overview: str, body: str, evidence_links: str) -> str:
    template = _template_path().read_text(encoding="utf-8")
    replacements = {
        "{{TITLE}}": html.escape(title),
        "{{IDENTITY}}": identity,
        "{{NAV}}": nav,
        "{{OVERVIEW}}": overview,
        "{{BODY}}": body,
        "{{EVIDENCE_LINKS}}": evidence_links,
    }
    for token, value in replacements.items():
        template = template.replace(token, value)
    return template


def format_notebook_records(package: Any, layer: int | str | None = None, node_id: str | None = None) -> str:
    """Render short English question records without repeating the HTML matrix."""

    package = _load(package)
    if node_id is not None:
        node = next((item for item in NODE_SPECS if item["node_id"] == node_id), None)
        if node is None:
            raise ValueError(f"unknown content node: {node_id}")
        method_file, method_anchor = node["method"]
        narrative_block = _notebook_narrative(package, node, language="en")
        lead = f'**{html.escape(_node_title(package, node, "en"))}**\n\n[Method](../../../docs/design/reconstruction-impact/{method_file}#{method_anchor})\n\n'
        summary = _section_summary(package, node_id)
        if _summary_review_status(package, summary) == "reviewed":
            review = summary.get("review", {}) if isinstance(summary.get("review", {}), Mapping) else {}
            conclusion = summary.get("conclusion", {}) if isinstance(summary.get("conclusion", {}), Mapping) else {}
            return lead + narrative_block + f'<h3>Saved conclusions: {html.escape(_node_title(package, node, "en"))}</h3><p>Review: {html.escape(str(review.get("status", "pending")))}</p><p class="section-summary">{html.escape(str(conclusion.get("en", "No saved section summary.")))}</p>'
        if node_id in SUMMARY_NODE_IDS:
            return lead + narrative_block + f'<h3>Saved conclusions: {html.escape(_node_title(package, node, "en"))}</h3><p>Review: pending</p><p class="section-summary pending">No saved section summary.</p>'
        records = _node_records(node, package.get("records", []))
        if node_id in COMPACT_NODE_IDS:
            parent_records = [r for r in records if r.get("identity", {}).get("scope") != "All"]
            content = _compact_block(parent_records, node_id=node_id, language="en")
            global_records = [r for r in records if r.get("identity", {}).get("scope") == "All"]
            if global_records:
                content += '<h4>HD global — separate cohort and denominator</h4>' + _compact_block(global_records, node_id=node_id, language="en")
            return lead + narrative_block + content
        if node_id in FOUNDATION_NODE_IDS:
            parent_records = [r for r in records if r.get("identity", {}).get("scope") != "All"]
            content = _foundation_block(node, parent_records)
            global_records = [r for r in records if r.get("identity", {}).get("scope") == "All"]
            if node_id == "1.2" and global_records:
                content += '<h4>HD global — separate cohort and denominator</h4>' + _foundation_block(node, global_records)
            return lead + narrative_block + f'<h3>Saved evidence: {html.escape(_node_title(package, node, "en"))}</h3>' + content
        if node_id in LOCAL_NODE_METRICS and node_id != "3.6":
            metric, label = LOCAL_NODE_METRICS[node_id]
            rows = []
            for record in records:
                identity = record.get("identity", {})
                values = _local_metric_values(record, label)
                rows.append((identity.get("scope"), identity.get("baseline"), values.get("baseline_median"), values.get("reconstruction_median"), values.get("delta_median"), _notebook_details(record)))
            return lead + narrative_block + f'<h3>Saved evidence: {html.escape(label)}</h3><p>Descriptive support; no standalone improvement judgment.</p>' + _foundation_table(("Scope", "Baseline", "Raw median", "Recon median", "Delta median", "Conditions"), rows)
        return lead + narrative_block + format_notebook_records({**package, "records": records}, layer=node["layer"])
    selected = LAYERS[layer - 1] if isinstance(layer, int) else layer
    records = [record for record in package.get("records", []) if record.get("performance_judgment") != "not_applicable" and (selected is None or _record_layer(record) == selected)]
    title = "Saved conclusions: Overview" if selected is None else "Saved conclusions: " + TITLES.get(str(selected), (str(selected), str(selected)))[1]
    result = [f'<h3>{html.escape(title)}</h3><p>Review: {html.escape(str((package.get("review", {}) or {}).get("status", "pending")))}</p>']

    def append_layer_sections(section_records: list[Mapping[str, Any]], heading: str | None = None) -> None:
        grouped: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
        for record in section_records:
            grouped[_record_layer(record)].append(record)
        if heading:
            result.append(f'<section class="notebook-hd-all"><h4>{html.escape(heading)}</h4><p>HD All records retain their separate cohort and denominator.</p>')
        for current_layer in LAYERS:
            if current_layer not in grouped:
                continue
            result.append(f'<section class="notebook-layer"><h4>{html.escape(TITLES[current_layer][1])}</h4>')
            if "node_narratives" in package:
                for node in NODE_SPECS:
                    if node.get("layer") != current_layer or not _narrative_reviewed(package, _node_narrative(package, str(node["node_id"]))):
                        continue
                    result.append(_notebook_narrative(package, node, language="en"))
            result.append("<ul>")
            for record in grouped[current_layer]:
                identity = record.get("identity", {}) or {}
                question = str(identity.get("question", "unknown"))
                label = _labels().get(question, (question, question))[1]
                scope = _display_scope(identity.get("scope", ""))
                metric_label, metric_value = _headline_metric(record)
                judgment = _judgment(record, "en")
                details = _notebook_details(record)
                saved_details = f' · {html.escape(details)}' if details else ""
                result.append(
                    f'<li class="notebook-question"><strong>{html.escape(label)}</strong> · {html.escape(scope)} · '
                    f'{html.escape(_brief_conclusion(record, "en"))} '
                    f'<span class="notebook-judgment">[{html.escape(judgment)}{saved_details}]</span> '
                    f'<span class="notebook-metric">[{html.escape(metric_label)}: {html.escape(metric_value)}]</span></li>'
                )
            result.append("</ul></section>")
        if heading:
            result.append("</section>")

    append_layer_sections([record for record in records if (record.get("identity", {}) or {}).get("scope") != "All"])
    hd_records = [record for record in records if (record.get("identity", {}) or {}).get("scope") == "All"]
    if hd_records:
        append_layer_sections(hd_records, "HD global — separate cohort and denominator")
    return "".join(result)


def render_report(package: Any, destination: str | os.PathLike[str]) -> Path:
    package = _load(package)
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    records = [record for record in package.get("records", []) if "node_narratives" in package or record.get("performance_judgment") != "not_applicable"]
    output_root = _output_root(package, destination)
    registry = _build_figure_registry(records)
    title = str(package.get("sample_id", "sample")) + " · Reconstruction impact"
    review = str((package.get("review", {}) or {}).get("status", "pending"))
    identity = f'Route: {html.escape(str(package.get("route_kind", "NA")))} · Spec: {html.escape(str(package.get("schema_version", "NA")))} · Review: {html.escape(review)}'
    if "node_narratives" in package:
        title = str(package.get("sample_id", "sample"))
        registry = _web_figure_registry(records, registry)
        nav = _web_nav(records)
        overview = [_web_overview(records)]
        identity = 'Raw → Recon · Fibroblast / Mono/Macro / T 分别比较'
    else:
        nav = "".join(f'<a href="#layer-{layer}">{i + 1}. {html.escape(TITLES[layer][0])}</a>' for i, layer in enumerate(LAYERS))
        overview = []
        for layer in LAYERS:
            overview.append(_overview_layer(layer, records))
        overview.append(_overview_hd(records))
    body = _body(package, records, registry, output_root, destination.parent)
    evidence_links = '<p><a href="analysis/conclusions/report_records.json">完整结论记录 JSON</a> · <a href="analysis/conclusions/report_records.csv">结论 CSV</a> · <a href="analysis/conclusions/review.json">审阅记录</a> · <a href="analysis/artifacts.csv">科学产物索引</a> · <a href="notebook/execution.json">执行审计</a></p>'
    destination.write_text(_render_shell(title=title, identity=identity, nav=nav, overview="".join(overview), body=body, evidence_links=evidence_links), encoding="utf-8")
    return destination
