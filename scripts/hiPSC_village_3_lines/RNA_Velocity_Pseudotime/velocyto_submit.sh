#!/bin/bash
## Date: 4 October, 2022
## Author: Drew Neavin
## Reason: RNA velocity on iPSC village samples for pseudotime following the Allevin tutorial https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/


SCRIPT="iPSC_Village_Publication/scripts/hiPSC_village_3_lines/RNA_Velocity_Pseudotime/velocyto.sh"
TENxDIR="/path/to/cellranger/gene_expression/" ### Available by processing fastq data through cellranger or on request from d.neavin @ garvan.org.au
OUT="/path/to/output/RNA_Velocity_Pseudotime/preprocess"
GTF="/path/to/gencode.v38.annotation.gtf" # available from https://www.gencodegenes.org/human/release_38.html
LOG=$OUT/logs

T=32


for pool in DRENEA_{1..6} Village_B_1_week Village_A_Baseline
do
	qsub -S /bin/bash \
        -q short.q \
        -r yes \
        -l mem_requested=4G \
        -l tmp_requested=4G \
		-pe smp $T \
        -N $pool\_velocyto \
        -cwd \
        -j y \
        -e $LOG \
        -o $LOG \
        -V \
        -v TENxDIR=$TENxDIR,OUT=$OUT,pool=$pool,GTF=$GTF,T=$T \
        -C '' $SCRIPT
done