#!/bin/bash
## Date: 4 October, 2022
## Author: Drew Neavin
## Reason: RNA velocity on iPSC village samples for pseudotime following the Allevin tutorial https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/


SCRIPT="iPSC_Village_Publication/scripts/hiPSC_village_3_lines/RNA_Velocity_Pseudotime/salmon_alevin.sh" # available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication
OUT="/path/to/output/RNA_Velocity_Pseudotime/preprocess"
LOG=$OUT/logs

mkdir -p $LOG

T=10

for pool in DRENEA_{1..6} Village_B_1_week Village_A_Baseline
do
	qsub -S /bin/bash \
        -q short.q \
        -r yes \
        -l mem_requested=4G \
        -l tmp_requested=4G \
        -N $pool\_alevin \
		-pe smp $T \
        -cwd \
        -j y \
        -e $LOG \
        -o $LOG \
        -V \
        -v OUT=$OUT,pool=$pool,T=$T \
        -C '' $SCRIPT
done