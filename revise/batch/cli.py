"""Command-line entrypoints for sample preparation and batch reconstruction."""
import argparse
import json


def main(argv=None):
    parser = argparse.ArgumentParser(description='Reconstruct all sample packages under a data root')
    parser.add_argument('--config', required=True, help='Batch YAML (schema_version and data_root)')
    args = parser.parse_args(argv)
    from .runner import run_batch
    try:
        result = run_batch(args.config)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(2, f'{exc}\n')
    print(json.dumps(result['summary']))
    if result['summary']['failed']:
        raise SystemExit(1)


def prepare_main(argv=None):
    parser = argparse.ArgumentParser(description='Prepare one standard sample package')
    parser.add_argument('--config', required=True, help='sample.yaml')
    args = parser.parse_args(argv)
    from .sample import prepare_sample
    try:
        sample = prepare_sample(args.config)
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(2, f'{exc}\n')
    print(f'Prepared {sample.sample_id}: {sample.st_path}')


if __name__ == '__main__':
    main()
