CLUSTER_CMD=("bsub -n {threads} -R {cluster.resources} -M {cluster.memory} -o {cluster.output} -J {cluster.name} -W {cluster.time}") 
snakemake --configfile config.yaml \
    --use-singularity \
    --singularity-args "--bind /juno" \
    --jobs 4 \
    --skip-script-cleanup \
    --cluster-config ./cluster.yaml \
    --cluster "${CLUSTER_CMD}" \
    --rerun-incomplete \
    --allowed-rules pvacseq_expression_filter
