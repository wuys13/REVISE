"""Command-line entrypoints for batch reconstruction and analysis."""
import argparse
import json


def main(argv=None):
    parser = argparse.ArgumentParser(description='Reconstruct standard ST samples with hierarchical batch.yaml configuration')
    parser.add_argument('--config', required=True, help='Project batch.yaml (schema_version: 2, input_root and output_root)')
    args = parser.parse_args(argv)
    from .runner import run_batch
    try:
        result = run_batch(args.config)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(2, f'{exc}\n')
    print(json.dumps(result['summary']))
    if result['summary']['failed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()


def analysis_main(argv=None):
    parser = argparse.ArgumentParser(description='Run configured analysis aspects for reconstructed samples')
    parser.add_argument('--config', required=True, help='Batch YAML with input_root, output_root and analysis')
    args = parser.parse_args(argv)
    from .analysis import run_analysis_batch
    try:
        result = run_analysis_batch(args.config)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(2, f'{exc}\n')
    print(json.dumps(result['summary']))
    if result['status'] != 'completed':
        raise SystemExit(1)
