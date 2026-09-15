"""Overview compression must preserve the approved comparison, not select a side."""
from reproduce.case.reconstruction_impact import report_html as report


def record(question, results):
    return {'identity': {'question': question, 'scope': 'Fibroblast'}, 'results': results}


def test_overview_complexity_retains_both_cluster_counts():
    label, value = report._headline_metric(record('complexity', {
        'raw_k': 3, 'reconstruction_k_at_raw_resolution': 11, 'fixed_final_cluster_k': None,
    }))
    assert '3' in value and '11' in value and '→' in value


def test_overview_moran_full_space_uses_q75_pair_not_median_difference():
    label, value = report._headline_metric(record('moran_all_valid', {
        'raw': {'q75': .125}, 'reconstruction': {'q75': .625}, 'median_difference': .999,
    }))
    assert 'Q75' in label and '0.125' in value and '0.625' in value
    assert '.999' not in value


def test_overview_emt_coverage_retains_both_resource_hit_denominators():
    label, value = report._headline_metric(record('emt_coverage', {
        'raw': {'available_gene_count': 31, 'resource_gene_count': 200, 'coverage': .155},
        'reconstruction': {'available_gene_count': 187, 'resource_gene_count': 200, 'coverage': .935},
    }))
    assert all(token in value for token in ('31', '187', '200', '→'))


def test_overview_window_support_names_selected_scale():
    label, value = report._headline_metric(record('window_support', {'main_window_side_um': 40, 'n_valid_windows': 2409}))
    assert value == '40'
    assert 'μm' in label or 'µm' in label


def test_unmatched_k_fraction_stays_out_of_headline():
    label, value = report._headline_metric(record('matched_k', {
        'matched': False, 'headline_eligible': False, 'status': 'unmatched_cluster_complexity',
        'unit_change_fraction': .57069,
    }))
    assert '57' not in value


def test_notebook_summary_does_not_split_decimal_numbers_as_sentences():
    item = record('moran_shared_valid', {})
    item['conclusion'] = {'en': 'Raw median=0.125; paired median delta=0.625. Compare identical genes.'}
    text = report._brief_conclusion(item)
    assert '0.125' in text and '0.625' in text
