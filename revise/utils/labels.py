"""Small, dependency-free label helpers."""


def normalize_cell_type_label(value: str) -> str:
    """Trim a cell-type label and replace slash separators with underscores."""
    return value.strip().replace("/", "_")
