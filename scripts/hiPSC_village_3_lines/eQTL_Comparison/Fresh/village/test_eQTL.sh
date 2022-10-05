

SNAKEFILE="iPSC_Village_Publication/scripts/hiPSC_village_3_lines/eQTL_Comparison/Fresh/village/test_eQTL.smk" # available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication
LOG="/path/to/output/eQTL_Comparison/Fresh/uni_village/village/logs/"

mkdir -p $LOG



nohup \
  snakemake \
    --snakefile $SNAKEFILE \
    --jobs 200 \
    --use-singularity \
    --restart-times 1 \
    --keep-going \
    --cluster \
        "qsub -S /bin/bash \
        -q short.q \
        -r yes \
        -pe smp {threads} \
        -l tmp_requested={resources.disk_per_thread_gb}G \
        -l mem_requested={resources.mem_per_thread_gb}G \
        -e $LOG \
        -o $LOG \
        -j y \
        -V" \
  > $LOG/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &

