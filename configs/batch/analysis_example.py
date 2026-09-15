"""Copyable adapter: inventory actual carriers, without a scientific comparison."""
import csv


def run(context):
    """Write carrier roles and axes using the verified reconstruction handoff."""
    with (context.output_dir / 'carriers.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['role', 'observation_role', 'observations', 'genes', 'source_path'])
        for role, carrier in context.reconstruction['outputs'].items():
            writer.writerow([role, carrier['observation_role'], *carrier['shape'], carrier['path']])
    return {
        'artifacts': {
            'carrier_inventory': {
                'path': 'carriers.csv',
                'description': 'Published carrier axes and source paths; not a biological result',
            },
        },
        'calculation': {
            'input_view': 'native_carriers',
            'parameters': {},
            'comparison_basis': 'inventory_only_no_scientific_comparison',
        },
    }
