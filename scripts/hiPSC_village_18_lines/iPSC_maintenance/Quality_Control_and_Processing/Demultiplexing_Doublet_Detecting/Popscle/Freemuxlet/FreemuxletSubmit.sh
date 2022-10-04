#!/bin/bash
##############################################################################
# Author: Drew Neavin
# Date: 2022-10-03
# Description: This script is to run Freemuxlet on the 18-line hiPSC village with "Freemuxlet.sh" script
##############################################################################


DATA="/path/to/cellranger/gene_expression/"
OUT="/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/Popscle/Freemuxlet"
LOG="$OUT/logs"
PIPELINE="iPSC_Village_Publication/scripts/hiPSC_village_18_lines/iPSC_maintenance/Demultiplexing_Doublet_Detecting/Popscle/Freemuxlet/Freemuxlet.sh" # available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication

mkdir -p $LOG



qsub -S /bin/bash \
    -q short.q \
    -r yes \
    -t 1-3 \
    -tc 3 \
    -l mem_requested=10G \
    -l tmp_requested=10G \
    -pe smp 16 \
    -N popscle_freemuxlet \
    -cwd \
    -j y \
    -e $LOG \
    -o $LOG \
    -v T=20,DATA=$DATA,OUT=$OUT \
    -C '' $PIPELINE

