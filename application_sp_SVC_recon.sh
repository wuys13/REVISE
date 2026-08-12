raw_data_path="./raw_data/Real_application"
sample_patients=(P1CRC) # add others if you want
data_type="HD"

for sample_patient in "${sample_patients[@]}"; do
    mkdir -p 0_records/${sample_patient}_application
    echo "Start patient_id: ${sample_patient} ....."
    python -u application_sp_SVC_recon.py --raw_data_path ${raw_data_path} \
                         --sample_name ${sample_patient} \
                         --st_file ${data_type}.h5ad \
                         --sc_ref_file adata_sc_all_reanno.h5ad \
                         > 0_records/${sample_patient}_application/sp_SVC.log 2>&1 &
done
echo "All tasks are running! Check results in the 'output' folder until all tasks are done."
