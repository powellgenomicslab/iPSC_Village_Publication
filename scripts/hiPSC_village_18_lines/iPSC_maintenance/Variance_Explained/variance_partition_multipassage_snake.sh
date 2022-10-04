
cd /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/multi-passage/Variance/

SNAKEFILE="iPSC_Village_Publication/scripts/hiPSC_village_18_lines/iPSC_maintenance/Variance_Explained/variance_partition_multipassage.smk" # available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication
LOG="/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/Variance//logs/"

mkdir -p $LOG


nohup \
  snakemake \
    --snakefile $SNAKEFILE \
    --jobs 100 \
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

