#!/bin/bash

SNAKEFILE=/path/to/scripts/hiPSC_village_18_lines/Growth_Rates/Snakefile

OUTDIR=/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/Growth_Rates/
LOGS=$OUTDIR/logs


mkdir $LOGS



cd /path/to/scripts/hiPSC_village_18_lines/Growth_Rates/




snakemake \
    --snakefile $SNAKEFILE \
    --rerun-incomplete \
    --jobs 20 \
    --use-conda \
    --keep-going \
    --cluster \
        "qsub -S /bin/bash \
        -q short.q \
        -r yes \
        -pe smp 4 \
        -l tmp_requested=8G \
        -l mem_requested=8G \
        -e $LOGS \
        -o $LOGS \
        -j y \
        -V"



        