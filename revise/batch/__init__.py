"""Standard ST directories, hierarchical configuration and sequential tasks."""

__all__ = ['run_reconstruction_task', 'run_analysis_task']


def run_reconstruction_task(config_path, sample_id, *, cell_type=None):
    """Reconstruct exactly one configured sample (and one type for iST)."""
    from .runner import run_reconstruction_task as execute

    return execute(config_path, sample_id, cell_type=cell_type)


def run_analysis_task(config_path, sample_id, aspect, *, cell_type=None):
    """Run exactly one configured analysis aspect for a reconstruction task."""
    from .analysis import run_analysis_task as execute

    return execute(config_path, sample_id, aspect, cell_type=cell_type)
