#!/usr/bin/env bash
set -euo pipefail

raw_data_path="${RAW_DATA_PATH:-./raw_data/Real_application}"
sample_name="${SAMPLE_NAME:-P2CRC}"
data_type="${DATA_TYPE:-Xenium}"
sc_ref_file="${SC_REF_FILE:-adata_sc_all_reanno.h5ad}"
save_path="${SAVE_PATH:-output/sc_SVC_case}"
config_path="${CONFIG_PATH:-revise/revise.yaml}"

# Selected resolution for the target cell type.
select_ct="${SELECT_CT:-Fibroblast}"

echo "Start sample_name: ${sample_name} ....."
python -u application_sc_SVC_recon.py \
    --config "${config_path}" \
    --sample-name "${sample_name}" \
    --st-file "${data_type}.h5ad" \
    --data-root "${raw_data_path}" \
    --output-root "${save_path}" \
    --sc-ref-file "${sc_ref_file}" \
    --select-ct "${select_ct}" \
    --compatibility-mode

echo "Finished. Check outputs under ${save_path}."
