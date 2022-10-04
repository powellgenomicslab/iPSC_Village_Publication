#!/bin/bash
##############################################################################
# Author: Drew Neavin
# Date: 2022-10-03
# Description: This script is to run souporcell on the 18-line hiPSC village with Demuxafy (https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/). 
##############################################################################


SAMPLE_INFO="iPSC_Village_Publication/scripts/hiPSC_village_18_lines/Pool_meta_data.tsv" # (accessible from the iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication)
SAMPLE=$(tail -n +2 $SAMPLE_INFO | tr '\n' ':' | cut -d ':' -f $SGE_TASK_ID)

COUNTS="$DATA/$SAMPLE/outs/filtered_feature_bc_matrix" ### Available by processing fastq data through cellranger or on request from d.neavin @ garvan.org.au


mkdir -p mkdir -p $OUT/$SAMPLE

singularity exec Demuxafy.sif scds.R -o $OUT/$SAMPLE -t $COUNTS
