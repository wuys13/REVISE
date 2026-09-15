"""The reading report keeps fixed matrices and real evidence without computation."""
import importlib
from pathlib import Path
from html.parser import HTMLParser


def _record(*, question, layer, scope='Fibroblast', conclusion='A saved conclusion.', metric=None, figures=(), tables=()):
    results = {}
    if metric is not None:
        results = metric
    return {
        'identity': {'question': question, 'scope': scope, 'layer': layer, 'task_cell_type': scope},
        'performance_judgment': 'reliable',
        'results': results,
        'conclusion': {'zh': conclusion, 'en': 'A saved conclusion.'},
        'evidence_condition': {'status': 'valid'},
        'limitations': {'zh': [], 'en': []},
        'evidence': {
            'figures': list(figures),
            'tables': [{'path': path, 'exists': True} for path in tables],
        },
    }


def _package(records, tmp_path):
    return {
        'sample_id': 'fixture',
        'route_kind': 'sp_svc',
        'schema_version': '1',
        'run_identity': {'output_root': str(tmp_path)},
        'review': {'status': 'pending'},
        'records': records,
    }


def _web_package(records, tmp_path):
    package = _package(records, tmp_path)
    # The saved narrative key selects the HTML-only reader layout.  An empty
    # list is sufficient here because these tests exercise saved records and
    # do not need to rebuild reviewed prose.
    package['node_narratives'] = []
    return package


def _between(text, start, end):
    return text[text.index(start):text.index(end, text.index(start))]


def test_web_partition_keeps_two_k_ari_groups_and_unmatched_diagnostic(tmp_path):
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    complexity = _record(
        question='complexity', layer='overall', scope='Fibroblast',
        metric={'comparison_kind': 'raw_reference_vs_fixed_final_clusters', 'raw_k': 5, 'fixed_final_cluster_k': 10, 'ari': .42},
    )
    matched = _record(
        question='matched_k', layer='overall', scope='Fibroblast',
        metric={'raw_k': 10, 'reconstruction_k': 10, 'unit_change_fraction': .25, 'balanced_change': .80, 'ari': .91, 'matched': False},
    )
    matched['performance_judgment'] = 'matched-K unavailable'
    matched['evidence_condition'] = {'status': 'review', 'matched_k_status': 'unmatched_cluster_complexity'}
    text = renderer.render_report(_web_package([complexity, matched], tmp_path), tmp_path / 'report.html').read_text(encoding='utf-8')
    complexity_block = _between(text, 'id="module-complexity"', 'id="module-matched_k"')
    matched_block = _between(text, 'id="module-matched_k"', 'id="node-2-3"')

    assert '<h4>原分群数量与一致性</h4>' in complexity_block
    assert '<th>Raw 参考 K</th>' in complexity_block and '<th>Recon 最终 K</th>' in complexity_block
    assert '<td>5</td>' in complexity_block and '<td>10</td>' in complexity_block and '<td>0.42</td>' in complexity_block
    assert '<h4>控制 K 后的归属对照</h4>' in matched_block
    assert '未成立 · 诊断值' in matched_block
    assert '<td>10 / 10</td>' in matched_block and '<td>0.91</td>' in matched_block


def test_web_moran_keeps_provided_n_and_shared_paired_delta(tmp_path):
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    foundation = _record(
        question='foundation', layer='foundation', scope='Fibroblast',
        metric={'gene': {'by_side': {'raw': {'provided_n': 123}, 'reconstruction': {'provided_n': 234}}}},
    )
    all_valid = _record(
        question='moran_all_valid', layer='overall', scope='Fibroblast',
        metric={'raw': {'n_valid': 100, 'q75': .10}, 'reconstruction': {'n_valid': 200, 'q75': .20}},
    )
    shared_valid = _record(
        question='moran_shared_valid', layer='overall', scope='Fibroblast',
        metric={'raw': {'n_valid': 80, 'q75': .10}, 'reconstruction': {'n_valid': 70, 'q75': .20}, 'paired_delta': {'n': 65, 'median': -.03}},
    )
    text = renderer.render_report(_web_package([foundation, all_valid, shared_valid], tmp_path), tmp_path / 'report.html').read_text(encoding='utf-8')
    gene_block = _between(text, 'id="gene-availability"', 'id="module-moran_all_valid"')
    shared_block = _between(text, 'id="module-moran_shared_valid"', '</section>')

    assert '<th>Raw 提供基因</th>' in gene_block and '<th>Recon 提供基因</th>' in gene_block
    assert '<td>123</td>' in gene_block and '<td>234</td>' in gene_block
    assert '<th>Raw 共同有效数</th>' in shared_block
    assert '<th>配对 ΔI 中位数</th>' in shared_block and '<td>-0.03</td>' in shared_block
    assert '<td>100</td>' in text and '<td>200</td>' in text


def test_web_dual_baselines_and_state_preserve_na_and_zero(tmp_path):
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    metrics = {
        'Kobs': {'baseline_median': 2.0, 'reconstruction_median': 1.5, 'delta_median': -.5},
        'Neff': {'baseline_median': 1.8, 'reconstruction_median': 1.4, 'delta_median': -.4},
        'evenness': {'baseline_median': .9, 'reconstruction_median': .95, 'delta_median': .05},
    }
    local_records = []
    for question, baseline in (('local_vs_raw_leiden', 'raw_leiden'), ('local_vs_raw_level2', 'raw_level2')):
        local = _record(
            question=question, layer='localization', scope='Fibroblast',
            metric={'baseline': baseline, 'scale_um': 40, 'n_valid_windows': 12, 'metrics': metrics},
        )
        local['identity']['baseline'] = baseline
        local_records.append(local)
    threshold = _record(
        question='threshold_reliability', layer='region', scope='Fibroblast',
        metric={'state': {'threshold': None, 'valid_bootstrap_fraction': .50, 'relative_ci_width': 1.2, 'n_windows': 4}},
    )
    threshold['performance_judgment'] = 'unreliable / NA'
    extent = _record(
        question='region_extent', layer='region', scope='Fibroblast',
        metric={'overall': {'region_windows': 0, 'valid_windows': 4, 'region_area_mm2': 0.0, 'area_fraction': 0.0, 'region_units': 0, 'valid_units': 10, 'unit_fraction': 0.0}},
    )
    text = renderer.render_report(_web_package([*local_records, threshold, extent], tmp_path), tmp_path / 'report.html').read_text(encoding='utf-8')
    baselines = _between(text, 'id="local-baselines"', 'id="node-3-5"')
    leiden = _between(baselines, 'id="module-local_vs_raw_leiden"', 'id="module-local_vs_raw_level2"')
    level2 = _between(baselines, 'id="module-local_vs_raw_level2"', '</article>')
    state = _between(text, 'id="node-4-1"', 'id="node-4-2"')
    extent_block = _between(text, 'id="node-4-2"', 'id="node-4-3"')

    assert 'id="module-local_vs_raw_leiden"' in baselines and 'id="module-local_vs_raw_level2"' in baselines
    assert 'Raw Leiden' in baselines and 'Raw Level2' in baselines
    assert '<td>-0.4</td>' in leiden and '<td>-0.4</td>' in level2
    assert '无稳定阈值' in state and '<td>NA</td>' in state
    assert '<td>0.0%</td>' in extent_block


def test_web_anchors_shared_figure_once_and_hd_all_only_saved_questions(tmp_path):
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    figure = 'figures/moran_all_gene_distribution.png'
    (tmp_path / 'figures').mkdir()
    (tmp_path / figure).write_bytes(b'fixture')

    def foundation(scope):
        return _record(
            question='foundation', layer='foundation', scope=scope,
            metric={'gene': {'by_side': {'raw': {'provided_n': 10}, 'reconstruction': {'provided_n': 20}}}},
        )

    def moran(scope):
        return _record(
            question='moran_all_valid', layer='overall', scope=scope,
            metric={'raw': {'n_valid': 8, 'q75': .1}, 'reconstruction': {'n_valid': 9, 'q75': .2}},
            figures=(figure,),
        )

    records = [foundation(scope) for scope in ('Fibroblast', 'Mono_Macro', 'T', 'All')]
    records.extend(moran(scope) for scope in ('Fibroblast', 'Mono_Macro', 'T', 'All'))
    text = renderer.render_report(_web_package(records, tmp_path), tmp_path / 'report.html').read_text(encoding='utf-8')

    class Parse(HTMLParser):
        def __init__(self):
            super().__init__()
            self.ids = []
            self.anchors = []

        def handle_starttag(self, tag, attrs):
            attributes = dict(attrs)
            if 'id' in attributes:
                self.ids.append(attributes['id'])
            if tag == 'a' and attributes.get('href', '').startswith('#'):
                self.anchors.append(attributes['href'][1:])

    dom = Parse()
    dom.feed(text)
    all_section = _between(text, 'id="hd-all"', 'id="foundation-audit"')

    assert len(dom.ids) == len(set(dom.ids))
    assert set(dom.anchors) <= set(dom.ids)
    assert text.count(f'src="{figure}"') == 1
    assert 'id="node-2-4-all"' in all_section
    assert 'id="node-2-1-all"' not in all_section
    assert 'id="node-2-2-all"' not in all_section
    assert 'id="node-2-5-all"' not in all_section


def test_overview_is_grouped_by_layer_with_one_headline_and_module_links(tmp_path):
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    records = []
    for layer, question in (
        ('foundation', 'foundation'),
        ('overall', 'matched_k'),
        ('localization', 'changed_units'),
        ('region', 'region_extent'),
    ):
        for scope in ('Fibroblast', 'Mono_Macro', 'T'):
            metric = {'overall': {'area_fraction': .25, 'unit_fraction': .40}} if question == 'region_extent' else {'headline': .123}
            records.append(_record(question=question, layer=layer, scope=scope, metric=metric))
    output = renderer.render_report(_package(records, tmp_path), tmp_path / 'report.html')
    text = output.read_text(encoding='utf-8')
    overview = _between(text, '<section id="overview"', '<section id="body"')

    assert overview.count('class="overview-layer"') == 4
    assert overview.count('class="overview-question"') == 4
    assert '<img ' not in overview
    assert 'href="#module-foundation"' in overview
    assert 'href="#figure-' not in overview
    assert '25.0%' in overview
    assert '40.0%' not in overview
    for scope in ('Fibroblast', 'Mono/Macro', 'T'):
        assert overview.count(f'class="overview-cell" data-scope="{scope}"') == 4


def test_body_uses_four_column_evidence_rows_and_embeds_each_figure_once(tmp_path):
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    common = 'figures/common_overview.png'
    fibroblast = 'figures/fibroblast_spatial.png'
    (tmp_path / 'figures').mkdir()
    (tmp_path / 'figures/common_overview.png').write_bytes(b'common')
    (tmp_path / 'figures/fibroblast_spatial.png').write_bytes(b'parent')
    records = [
        _record(question='matched_k', layer='overall', scope=scope, figures=(common,), tables=('analysis/matched_k.csv',))
        for scope in ('Fibroblast', 'Mono_Macro', 'T')
    ]
    records.append(_record(question='emt_spatial_fields', layer='localization', scope='Fibroblast', figures=(fibroblast,), tables=('analysis/spatial.csv',)))
    output = renderer.render_report(_package(records, tmp_path), tmp_path / 'report.html')
    text = output.read_text(encoding='utf-8')

    assert text.count('class="evidence-table"') >= 2
    assert '科学问题／观察方向' in text
    assert '本批结论' in text
    assert '数值证明与条件' in text
    assert '图形证据' in text
    assert text.count('src="figures/common_overview.png"') == 1
    assert text.count('src="figures/fibroblast_spatial.png"') == 1
    assert text.count('id="row-module-matched_k-batch"') == 1
    assert 'id="row-module-emt_spatial_fields-fibroblast"' in text
    assert 'figures/common_overview.png' in text
    assert 'Original figure' in text
    assert 'analysis/matched_k.csv' in text
    assert 'common_overview.png' in text[text.index('id="module-matched_k"'):]


def test_common_figures_with_hd_all_are_owned_by_main_batch_row(tmp_path):
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    common = 'figures/moran_all_gene_distribution.png'
    (tmp_path / 'figures').mkdir()
    (tmp_path / common).write_bytes(b'common')
    records = [
        _record(question='moran_all_valid', layer='overall', scope=scope, figures=(common,))
        for scope in ('Fibroblast', 'Mono_Macro', 'T')
    ]
    records.append(_record(question='moran_all_valid', layer='overall', scope='All', figures=(common,)))
    text = renderer.render_report(_package(records, tmp_path), tmp_path / 'report.html').read_text(encoding='utf-8')
    main_row = text[text.index('id="row-module-moran_all_valid-batch"'):text.index('</tr>', text.index('id="row-module-moran_all_valid-batch"'))]
    hd_module = text[text.index('id="module-moran_all_valid-all"'):]
    assert '<figure ' in main_row
    assert 'Shared evidence: Moran all-valid/shared-valid distribution comparison' in hd_module
    assert 'moran_all_gene_distribution.png' not in hd_module
    assert text.count('src="figures/moran_all_gene_distribution.png"') == 1
    assert 'Moran all-valid/shared-valid distribution comparison' in text


def test_modules_link_to_the_next_question(tmp_path):
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    records = [
        _record(question='moran_all_valid', layer='overall', scope='Fibroblast'),
        _record(question='moran_shared_valid', layer='overall', scope='Fibroblast'),
    ]
    text = renderer.render_report(_package(records, tmp_path), tmp_path / 'report.html').read_text(encoding='utf-8')
    module = text[text.index('id="module-moran_all_valid"'):text.index('id="module-moran_shared_valid"')]
    assert '下一问题' in module
    assert 'href="#module-moran_shared_valid"' in module


def test_figure_roles_use_registered_names_only():
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    assert renderer._preferred_figure_question('figures/moran_all_gene_distribution.png', {'moran_all_valid'}) == 'moran_all_valid'
    assert renderer._preferred_figure_question('figures/mystery_moran_plot.png', {'moran_all_valid'}) is None


def test_scale_sensitivity_proof_keeps_saved_baseline_rows(tmp_path):
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    record = _record(
        question='scale_sensitivity', layer='region', scope='Fibroblast',
        metric={'main_scale_um': 40, 'rows': [
            {'window_side_length': 24, 'median_delta_neff_vs_raw_leiden': -0.4, 'median_delta_neff_vs_raw_level2': 0.2, 'n_valid_windows': 12, 'min_parent_units': 4},
        ]},
    )
    text = renderer.render_report(_package([record], tmp_path), tmp_path / 'report.html').read_text(encoding='utf-8')
    assert '<th>ΔNeff RL</th>' in text
    assert '<th>ΔNeff L2</th>' in text
    assert '<th>Minimum parent units</th>' not in text
    assert 'RL = Raw Leiden · L2 = Raw Level2 · Minimum parent units: 4' in text
    assert '-0.4' in text and '0.2' in text and '12' in text


def test_figure_captions_have_one_original_link_and_moran_panel_titles(tmp_path):
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    figure = 'figures/moran_all_gene_distribution.png'
    q75 = 'figures/moran_q75_heatmap.png'
    (tmp_path / 'figures').mkdir()
    (tmp_path / figure).write_bytes(b'distribution')
    (tmp_path / q75).write_bytes(b'q75')
    records = [
        _record(question='moran_all_valid', layer='overall', scope=scope, figures=(figure, q75))
        for scope in ('Fibroblast', 'Mono_Macro', 'T')
    ]
    text = renderer.render_report(_package(records, tmp_path), tmp_path / 'report.html').read_text(encoding='utf-8')
    assert 'Moran all-valid/shared-valid distribution comparison' in text
    assert 'Moran all-valid/shared-valid Q75 comparison' in text
    assert text.count('Original figure') == 2
    assert 'Original figure / 原图 · Original figure' not in text


def test_local_evenness_proof_reads_saved_lowercase_metric_key(tmp_path):
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    record = _record(
        question='local_vs_raw_leiden', layer='localization', scope='Fibroblast',
        metric={'metrics': {'evenness': {
            'baseline_median': 0.94,
            'baseline_n_observations': 2478,
            'reconstruction_median': 0.96,
            'delta_median': 0.01,
        }}, 'scale_um': 40, 'n_valid_windows': 12, 'baseline': 'raw_leiden'},
    )
    text = renderer.render_report(_package([record], tmp_path), tmp_path / 'report.html').read_text(encoding='utf-8')
    evenness_row = _between(text, 'id="node-3-7"', '</article>')
    assert '<h4>evenness · Raw Leiden</h4>' in evenness_row
    assert '<th>Raw median</th><th>Recon median</th><th>Paired Δ median</th>' in evenness_row
    assert '<td>0.94</td><td>0.96</td><td>0.01</td>' in evenness_row
    assert '<td>12</td><td>40</td>' in evenness_row


def test_changed_units_uses_total_units_when_paired_units_is_absent(tmp_path):
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    record = _record(
        question='changed_units', layer='localization', scope='All',
        metric={
            'overall': {'change_fraction': 0.44, 'changed_units': 11119, 'total_units': 25060},
            'by_region': [
                {'level1': 'Overall', 'change_fraction': 0.44, 'changed_units': 11119, 'paired_units': 25060, 'wilson_ci_lower': .40, 'wilson_ci_upper': .48, 'low_sample_size': False},
                {'level1': 'B', 'change_fraction': 0.50, 'changed_units': 5, 'paired_units': 10, 'wilson_ci_lower': .20, 'wilson_ci_upper': .80, 'low_sample_size': True},
            ],
        },
    )
    proof = renderer._key_values(record)
    assert '11,119/25,060' in proof
    assert '11,119/NA' not in proof
    text = renderer.render_report(_web_package([record], tmp_path), tmp_path / 'report.html').read_text(encoding='utf-8')
    all_section = _between(text, 'id="hd-all"', 'id="foundation-audit"')
    assert 'Raw Level1 类型' in all_section
    assert '<td>Overall</td>' in all_section and '<td>B</td>' in all_section
    assert '40.0%–48.0%' in all_section and '20.0%–80.0%' in all_section
    assert '<td>否</td>' in all_section and '<td>是</td>' in all_section


def test_region_reliability_precedes_extent_and_overview_uses_area_fraction(tmp_path):
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    records = [
        _record(
            question='threshold_reliability', layer='region', scope=scope,
            conclusion='Threshold is reliable before extent.',
            metric={'state': {'relative_ci_width': .10, 'threshold': 1.8}},
        )
        for scope in ('Fibroblast', 'Mono_Macro', 'T')
    ]
    records.extend(
        _record(
            question='region_extent', layer='region', scope=scope,
            metric={'overall': {'area_fraction': .25, 'unit_fraction': .40}},
        )
        for scope in ('Fibroblast', 'Mono_Macro', 'T')
    )
    output = renderer.render_report(_package(records, tmp_path), tmp_path / 'report.html')
    text = output.read_text(encoding='utf-8')
    assert text.index('module-threshold_reliability') < text.index('module-region_extent')
    overview = _between(text, '<section id="overview"', '<section id="body"')
    assert 'Area fraction' in overview
    assert '25.0%' in overview
    assert '40.0%' not in overview


def test_notebook_formatter_emits_short_question_records(tmp_path):
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    records = [
        _record(question='matched_k', layer='overall', scope=scope, conclusion='Short saved conclusion.')
        for scope in ('Fibroblast', 'Mono_Macro', 'T')
    ]
    formatted = renderer.format_notebook_records(_package(records, tmp_path))
    assert 'Saved conclusions: Overview' in formatted
    assert 'Partition complexity' not in formatted
    assert 'A saved conclusion.' in formatted
    assert 'Matched K' in formatted
    assert formatted.count('class="notebook-question"') == 3
    assert '<table' not in formatted


def test_notebook_short_records_keep_judgment_condition_and_hd_records(tmp_path):
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    records = [
        _record(
            question='matched_k', layer='overall', scope='Fibroblast',
            conclusion='A short saved conclusion. A second detail must stay out of the summary.',
        ),
        _record(
            question='matched_k', layer='overall', scope='All',
            conclusion='HD saved conclusion.',
        ),
    ]
    records[0]['performance_judgment'] = 'matched-K unavailable'
    records[0]['evidence_condition'] = {'status': 'review', 'matched_k_status': 'unmatched_cluster_complexity'}
    records[1]['conclusion']['en'] = 'HD saved conclusion.'
    formatted = renderer.format_notebook_records(_package(records, tmp_path))
    assert 'matched-K is unavailable' in formatted
    assert 'review' in formatted.lower() or 'condition' in formatted.lower()
    assert 'HD global — separate cohort and denominator' in formatted
    assert 'HD saved conclusion.' in formatted
    assert formatted.count('class="notebook-question"') == 2


def test_notebook_formatter_reads_saved_pending_section_summary_by_node(tmp_path):
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    package = _package([
        _record(question='matched_k', layer='overall', scope='Fibroblast'),
    ], tmp_path)
    package['section_summaries'] = [{
        'node_id': '2.6',
        'conclusion': {
            'zh': '已保存的总体观察。',
            'en': 'Saved overall observation.',
        },
        'record_ids': ['fixture:overall:Fibroblast'],
        'review': {'status': 'pending'},
    }]

    formatted = renderer.format_notebook_records(package, node_id='2.6')

    assert 'Saved overall observation.' not in formatted
    assert 'Review: pending' in formatted
    assert 'Matched K' not in formatted


def test_foundation_nodes_render_saved_roles_pairing_side_totals_and_preprocessing(tmp_path):
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    record = _record(question='foundation', layer='foundation', scope='Fibroblast')
    record['results'] = {
        'carrier': {
            'status': 'valid',
            'rows': [
                {'role': 'Raw full context', 'n_units': 20, 'n_genes': 100, 'spatial_axis': True},
                {'role': 'Reconstruction full carrier', 'n_units': 20, 'n_genes': 150, 'spatial_axis': True},
            ],
        },
        'source_files': {'status': 'valid', 'rows': [{'role': 'raw', 'path': 'raw.h5ad'}]},
        'cohort': {'status': 'valid', 'label': 'Raw-defined paired scope'},
        'selection': {
            'status': 'valid', 'raw_parent_label': 'Fibroblast',
            'raw_parent_n': 24, 'reconstruction_covered_parent_n': 20,
            'sampled_n': 18,
        },
        'denominator': {'input_units': 20, 'paired_units': 18, 'excluded_units': {'raw_qc': 2}},
        'pairing_audit': {
            'status': 'valid', 'coordinate_source': 'saved_observations', 'coordinate_unit': 'um',
            'n_checked': 18, 'id_unique': True, 'id_subset': True,
            'coordinates_match': True, 'coordinate_finite': True,
        },
        'gene': {
            'status': 'valid', 'gene_space': 'full',
            'by_side': {
                'raw': {'n_genes': 999, 'provided_n': 100, 'source': 'raw.h5ad'},
                'reconstruction': {'n_genes': 999, 'provided_n': 150, 'source': 'recon.h5ad'},
            },
        },
        'preprocessing': {
            'partition': {'status': 'valid', 'path': 'partition_audit.json', 'fields': {'mode': 'level1_ari', 'n_features': 2000}},
            'moran': {'status': 'valid', 'path': 'moran_graph_audit.csv', 'fields': {'n_units': 20, 'n_neighbors_actual': 6, 'min_units': 5}},
            'emt': {'status': 'valid', 'path': 'cutoff_audit.csv', 'fields': {'scorer': 'AUCell', 'cutoff_by_side': {'raw': {'effective_rank_length': 20}, 'reconstruction': {'effective_rank_length': 150}}}},
        },
    }

    text = renderer.render_report(_package([record], tmp_path), tmp_path / 'report.html').read_text(encoding='utf-8')

    assert 'Saved fact' not in text
    assert 'Raw full context' in text and 'Reconstruction full carrier' in text
    assert 'Observation units' in text and 'Boundary' in text
    assert 'Parent selection from full-Raw observations' in text
    assert 'Full-Raw parent n' in text and 'Fibroblast' in text and 'Saved analysis denominators' in text
    assert 'Moran graph n' in text and '20' in text
    assert 'saved_observations' in text and 'Coordinates match' in text
    assert '18 (full-spatial carrier check)' in text
    assert 'Raw total genes' in text and 'Recon total genes' in text
    assert '<td>100</td>' in text and '<td>150</td>' in text and '<td>999</td>' not in text
    assert 'Partition' in text and 'Moran' in text and 'EMT' in text
    assert 'features=2,000' in text and 'min units=5' in text and 'Raw rank=20' in text and 'Recon rank=150' in text


def test_local_diversity_groups_each_metric_with_two_saved_baselines(tmp_path):
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    metrics = {
        'Kobs': {'baseline_median': 2.0, 'baseline_n_observations': 12, 'reconstruction_median': 1.5, 'delta_median': -0.5},
        'Neff': {'baseline_median': 1.8, 'baseline_n_observations': 12, 'reconstruction_median': 1.4, 'delta_median': -0.4},
        'evenness': {'baseline_median': 0.9, 'baseline_n_observations': 12, 'reconstruction_median': 0.95, 'delta_median': 0.05},
    }
    records = []
    for question, baseline in (
        ('local_vs_raw_leiden', 'raw_leiden'),
        ('local_vs_raw_level2', 'raw_level2'),
    ):
        record = _record(
            question=question,
            layer='localization',
            scope='Fibroblast',
            metric={'baseline': baseline, 'scale_um': 40, 'n_valid_windows': 12, 'metrics': metrics},
        )
        record['identity']['baseline'] = baseline
        record['performance_judgment'] = 'supports concentrated state direction'
        record['evidence_condition'] = {'status': 'valid', 'baseline': baseline}
        records.append(record)

    text = renderer.render_report(_package(records, tmp_path), tmp_path / 'report.html').read_text(encoding='utf-8')

    for metric in ('kobs', 'neff', 'evenness'):
        module = _between(text, f'id="module-local-diversity-{metric}"', '</article>')
        assert 'Raw Leiden' in module
        assert 'Raw Level2' in module
        assert '<th>Raw median</th><th>Recon median</th><th>Paired Δ median</th>' in module
        expected_row = {
            'kobs': '<td>2</td><td>1.5</td><td>-0.5</td>',
            'neff': '<td>1.8</td><td>1.4</td><td>-0.4</td>',
            'evenness': '<td>0.9</td><td>0.95</td><td>0.05</td>',
        }[metric]
        assert module.count(expected_row) == 2
        assert module.count('<td>40</td>') == 2


def test_html_renderer_keeps_layers_and_unique_shared_figures(tmp_path):
    module_path = Path(__file__).resolve().parents[2] / 'reproduce/case/reconstruction_impact/report_html.py'
    assert module_path.exists(), 'Local report renderer is not implemented'
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    assert renderer._template_path().is_file(), 'Runtime report shell template is missing'
    figures = tmp_path/'figures'; figures.mkdir()
    (figures/'threshold_reliability.png').write_bytes(b'fixture')
    records=[]
    for scope in ('Fibroblast','Mono_Macro','T'):
        records.append({'identity':{'question':'threshold_reliability','scope':scope,'layer':'region'},
                        'performance_judgment':'reliable','results':{'state':{'threshold':1.78,'relative_ci_width':.233}},
                        'conclusion':{'zh':'当前条件下可识别。','en':'Identifiable under the current conditions.'},
                        'evidence_condition':{'status':'valid'},'limitations':{'zh':[],'en':[]},
                        'evidence':{'figures':['figures/threshold_reliability.png'],'tables':[]}})
    package={'sample_id':'fixture','route_kind':'sp_svc','schema_version':'1',
             'run_identity':{'output_root':str(tmp_path)},'review':{'status':'pending'},'records':records}
    output=renderer.render_report(package,tmp_path/'report.html')
    text=output.read_text()
    class Parse(HTMLParser):
        def __init__(self): super().__init__();self.ids=[];self.images=[];self.anchors=[]
        def handle_starttag(self,tag,attrs):
            a=dict(attrs)
            if 'id' in a:self.ids.append(a['id'])
            if tag=='img':self.images.append(a['src'])
            if tag=='a' and a.get('href','').startswith('#'):self.anchors.append(a['href'][1:])
    dom=Parse();dom.feed(text)
    assert len(dom.images)==1
    assert len(dom.ids)==len(set(dom.ids))
    assert all(a in dom.ids for a in dom.anchors)
    assert [dom.ids.index('layer-'+x) for x in ('foundation','overall','localization','region')]==sorted(dom.ids.index('layer-'+x) for x in ('foundation','overall','localization','region'))
    assert '当前条件下可识别' in text and '23.3%' in text
    english=renderer.format_notebook_records(package,layer=4)
    assert 'Identifiable under' in english and '当前条件' not in english


def test_display_preserves_route_specific_k_cutoff_and_missing_extent():
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    def record(question, results):
        return {'identity': {'question': question}, 'results': results}
    complexity = renderer._key_values(record('complexity', {'raw_k': 5, 'fixed_final_cluster_k': 10}))
    assert 'Fixed final K' in complexity and '>10<' in complexity
    hd = renderer._key_values(record('complexity', {'raw_k': 5, 'fixed_final_cluster_k': None, 'reconstruction_k_at_raw_resolution': 11}))
    assert 'Recon K' in hd and '>11<' in hd and 'Fixed final K' not in hd
    emt = renderer._key_values(record('emt_score', {'cutoff_audit': {'raw': {'effective_rank_length': 26}, 'reconstruction': {'effective_rank_length': 9008}}}))
    assert '>26<' in emt and '>9,008<' in emt
    extent = renderer._key_values(record('region_extent', {'overall': {'region_area_mm2': None, 'area_fraction': None, 'valid_units': 100}}))
    assert 'Area (mm²): <strong>NA' in extent


def test_foundation_roles_are_types_and_not_repeated_counts():
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    from reproduce.case.reconstruction_impact.content_contract import NODE_SPECS
    record = _record(question='foundation', layer='foundation', scope='Fibroblast')
    record['results'] = {'carrier': {'rows': [{'role': 'Raw full context', 'n_units': 123456, 'n_genes': 18085, 'spatial_axis': True}]}}
    block = renderer._foundation_block(NODE_SPECS[0], [record])
    assert '123456' not in block and '123,456' not in block and '18,085' not in block
    assert 'spatial unit' in block.lower()
    assert 'source-flow' in block


def test_local_matrix_image_belongs_to_metric_row(tmp_path):
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    record = _record(question='local_vs_raw_leiden', layer='localization', scope='Fibroblast')
    record['identity']['baseline'] = 'raw_leiden'
    record['evidence']['figures'] = ['figures/Fibroblast_neff_matrix.png']
    path = tmp_path / 'figures/Fibroblast_neff_matrix.png'
    path.parent.mkdir(); path.write_bytes(b'fixture')
    text = renderer.render_report(_package([record], tmp_path), tmp_path/'report.html').read_text()
    module = _between(text, 'id="node-3-6"', '</article>')
    assert '<img ' in module
    assert text.count('<img ') == 1


def test_notebook_foundation_nodes_select_their_own_saved_facts(tmp_path):
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    record = _record(question='foundation', layer='foundation', scope='Fibroblast')
    record['results'] = {
        'carrier': {'rows': [{'role': 'Raw full context', 'n_units': 200, 'spatial_axis': True}]},
        'gene': {'status': 'valid', 'by_side': {'raw': {'provided_n': 100}, 'reconstruction': {'provided_n': 150}}},
    }
    package = _package([record], tmp_path)
    carrier = renderer.format_notebook_records(package, node_id='1.1')
    genes = renderer.format_notebook_records(package, node_id='1.4')
    assert 'Raw full context' in carrier and 'Raw total genes' not in carrier
    assert 'Raw total genes' in genes and '<td>100</td>' in genes
    assert 'Raw full context' not in genes


def test_anatomy_and_extent_nodes_show_saved_stratum_denominators():
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    local = _record(question='local_vs_raw_leiden', layer='localization', scope='Fibroblast')
    local['identity']['baseline'] = 'raw_leiden'
    local['results']['anatomy_summary'] = [{'level1_region': 'Tumor', 'n_valid_windows': 12, 'scale_um': 40, 'median_delta_neff_vs_raw_leiden': -0.25}]
    anatomy = renderer._compact_block([local], node_id='3.8')
    assert 'Tumor' in anatomy and '<td>12</td>' in anatomy and '-0.25' in anatomy
    region = _record(question='region_extent', layer='region', scope='Fibroblast')
    region['results']['by_region'] = [{'level1_region': 'Tumor', 'valid_windows': 12, 'region_windows': None, 'area_fraction': None, 'threshold_status': 'ci_too_wide'}]
    extent = renderer._compact_block([region], node_id='4.3')
    assert 'Tumor' in extent and '<td>12</td>' in extent and '<td>NA</td>' in extent
    assert 'ci_too_wide' in extent


def test_notebook_cell_type_summary_keeps_hd_global_separate(tmp_path):
    renderer = importlib.import_module('reproduce.case.reconstruction_impact.report_html')
    records = [_record(question='matched_k', layer='overall', scope=scope) for scope in ('Fibroblast', 'All')]
    text = renderer.format_notebook_records(_package(records,tmp_path),node_id='2.3')
    parent, global_section = text.split('HD global — separate cohort and denominator')
    assert '<td>Fibroblast</td>' in parent and '<td>All</td>' not in parent
    assert '<td>All</td>' in global_section
