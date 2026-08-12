#!/usr/bin/env bash
raw_data_path="./raw_data/Real_application"
sample_name="P2CRC"
data_type="Xenium"

sc_ref_file="adata_sc_all_reanno.h5ad"
save_path="output/sc_SVC_case"

 # selected resolution for Fibroblast, you can set it by your preference or ref to our automated pipeline
select_ct="Fibroblast"
select_resolution=0.8

echo "Start sample_name: ${sample_name} ....."
python -u application_sc_SVC_recon.py --sample_name ${sample_name} \
                         --data_type ${data_type} \
                         --raw_data_path ${raw_data_path} \
                         --save_path  ${save_path} \
                         --sc_ref_file ${sc_ref_file}\
                         --select_ct ${select_ct} \
                         --select_resolution ${select_resolution}

echo "All tasks are running! Check results in the 'logs' folder until all tasks are done."
