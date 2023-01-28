#!/bin/bash

SNAKEFILE="/path/to/scripts/hiPSC_village_18_lines/iPSC_maintenance/eQTL_Comparison/test_eQTL.smk"
LOG="/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/eQTL_Comparison/deboever/gene_separate/logs/"

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


