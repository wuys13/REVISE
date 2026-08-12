# segmentation/bin2cell/batch_effect/spot_size/gene_panel/gene_dropout"
cfs=(segmentation bin2cell batch_effect spot_size gene_panel gene_dropout)

for cf in "${cfs[@]}"; do
    echo "Start task ${cf} ....."
    bash ./reproduce/benchmark/benchmark_${cf}.sh

done

echo "All tasks are running! Check results in the 'logs' folder until all tasks are done."