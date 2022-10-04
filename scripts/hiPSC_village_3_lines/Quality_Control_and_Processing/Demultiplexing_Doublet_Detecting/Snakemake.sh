#!/bin/bash
##### Information about the use of this snakemake submission file file
## Author: Drew Neavin
## Date: 2022-10-03
## Description: This script is to subit the snakefile to run demultiplexing and doublet detecting on the 3-line hiPSC villages (iPSC_Village_Publication/scripts/hiPSC_village_3_lines/Quality_Control/Demultiplexing_Doublet_Detecting/Snakefile from github: https://github.com/powellgenomicslab/iPSC_Village_Publication)

export SINGULARITY_BINDPATH="/path/"
SNAKEFILE="iPSC_Village_Publication/scripts/hiPSC_village_3_lines/Quality_Control/Demultiplexing/Snakemake.sh" # available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication
LOGS="/path/to/output/demultiplexing_snakemake/logs"
mkdir -p $LOGS

nohup snakemake --snakefile $SNAKEFILE \
    --rerun-incomplete --jobs 100 --use-singularity --restart-times 4 --keep-going \
    --cluster "qsub -S /bin/bash -q short.q -r yes -pe smp {threads} -l tmp_requested={resources.disk_per_thread_gb}G -l mem_requested={resources.mem_per_thread_gb}G -e $LOGS -o $LOGS -j y -V" \
    > $LOGS/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &
