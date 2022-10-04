#!/bin/bash
##############################################################################
# Author: Drew Neavin
# Date: 2022-10-04
# Description: This script is to run Demuxafy combination of demultiplexing methods for 18-line hiPSC village
##############################################################################


INDIR="/path/to/output/hiPSC_village_18_lines/iPSC_maintenance"
OUT="/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/Combined"


SAMPLE_INFO="iPSC_Village_Publication/scripts/hiPSC_village_18_lines/iPSC_maintenance/Sample_meta.tsv" # (accessible from the iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication)



for $SAMPLE in `tail -n +2 $SAMPLE_INFO | tr '\n' ' '`
do
    singularity exec Demuxafy.sif Combine_Results.R \
        -o $OUT/$SAMPLE \
        -d $INDIR/Popscle/Demuxlet/$SAMPLE \
        -f $INDIR/Popscle/Freemuxlet/$SAMPLE \
        -u $INDIR/Souporcell/$SAMPLE \
        -v $INDIR/Vireo/$SAMPLE \
        --method "AtLeastHalfSinglet"
done
