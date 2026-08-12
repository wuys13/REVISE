raw_data_path="./raw_data/Sim2Real-ST"
# segmentation/bin2cell/batch_effect/spot_size/gene_panel/gene_dropout"
sample_patient=P2CRC
sample_parts=(part1 part2 part3)


cf="segmentation"
for sample_part in "${sample_parts[@]}"; do
    mkdir -p 0_records/${sample_patient}_${sample_part}
    echo "Start patient_id: ${sample_patient}; part ${sample_part}; task ${cf} ....."
    nohup python -u benchmark_main.py --cf ${cf} \
                             --raw_data_path ${raw_data_path} \
                             --task ${cf} \
                             --st_file xenium_spot.h5ad \
                             --gt_svc_file selected_xenium.h5ad \
                             --sc_ref_file real_sc_ref.h5ad \
                             --sample_name ${sample_patient}/cut_${sample_part} \
                             > 0_records/${sample_patient}_${sample_part}/${cf}.log 2>&1 &
done
wait