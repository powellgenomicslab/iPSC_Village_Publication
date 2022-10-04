## Date: 4 October, 2022
## Author: Drew Neavin
## Reason: RNA velocity on iPSC village samples for pseudotime following the Allevin tutorial https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/



SIDX="/path/to/output/RNA_Velocity_Pseudotime/preprocess/gencode.v38.annotation.expanded.sidx" ## Produced by iPSC_Village_Publication/scripts/hiPSC_village_3_lines/RNA_Velocity_Pseudotime/scVelo_preprocess.sh script available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication
tx2gene="/path/to/output/RNA_Velocity_Pseudotime/preprocess/gencode.v38.annotation.expanded.tx2gene.tsv" ## Produced by iPSC_Village_Publication/scripts/hiPSC_village_3_lines/RNA_Velocity_Pseudotime/scVelo_preprocess.sh script available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication



salmon alevin -l ISR -i $SIDX \
	-1 $OUT/$pool/$pool\_L001_R1_001.fastq.gz \
	-2 $OUT/$pool/$pool\_L001_R2_001.fastq.gz \
	-o $OUT/$pool/alevin_out --tgMap $tx2gene \
	--chromium --dumpFeatures --expectCells 20000 -p $T

