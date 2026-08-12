from dataclasses import dataclass

@dataclass
class BaseConf:
    # runtime parameters
    sample_name: str
    raw_data_path: str
    result_root_path: str

    # annotate column keys
    cell_type_col: str
    confidence_col: str
    unknown_key: str
