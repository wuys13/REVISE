#!/usr/bin/env bash
set -euo pipefail

raw_data_path="${RAW_DATA_PATH:-./raw_data/Real_application}"
save_path="${SAVE_PATH:-output/sp_SVC_case}"
config_path="${CONFIG_PATH:-revise/revise.yaml}"
sample_patients=(${SAMPLE_PATIENTS:-P1CRC})
data_type="${DATA_TYPE:-HD}"

for sample_patient in "${sample_patients[@]}"; do
    mkdir -p "0_records/${sample_patient}_application"
    echo "Start sample_name: ${sample_patient} ....."
    python -u application_sp_SVC_recon.py \
        --config "${config_path}" \
        --data-root "${raw_data_path}" \
        --sample-name "${sample_patient}" \
        --st-file "${data_type}.h5ad" \
        --sc-ref-file "adata_sc_all_reanno.h5ad" \
        --output-root "${save_path}" \
        > "0_records/${sample_patient}_application/sp_SVC.log" 2>&1 &
done

echo "All tasks are running. Check logs under 0_records/ and outputs under ${save_path}."
