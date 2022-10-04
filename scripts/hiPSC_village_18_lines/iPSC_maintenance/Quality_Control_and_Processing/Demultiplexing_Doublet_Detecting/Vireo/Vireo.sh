#!/bin/bash
##############################################################################
# Author: Drew Neavin
# Date: 2022-10-03
# Description: This script is to run vireo on the 18-line hiPSC village with Demuxafy (https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/). 
##############################################################################


SAMPLE_INFO="iPSC_Village_Publication/scripts/hiPSC_village_18_lines/iPSC_maintenance/Sample_meta.tsv" # (accessible from the iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication)
SAMPLE=$(tail -n +2 $SAMPLE_INFO | tr '\n' ':' | cut -d ':' -f $SGE_TASK_ID)

VCF="/path/to/output/vcf/snps_18_hiPSC_r2_maf0.3_exons_hg38.vcf" # Exons only
FASTA="/path/to/GRCh38/genome.fa"

BAM=$DATA/$SAMPLE/outs/possorted_genome_bam.bam ### Available by processing fastq data through cellranger or on request from d.neavin @ garvan.org.au
BARCODES="$DATA/$SAMPLE/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" ### Available by processing fastq data through cellranger or on request from d.neavin @ garvan.org.au


mkdir -p mkdir -p $OUT/$SAMPLE


singularity exec Demuxafy.sif cellsnp-lite -s $BAM -b $BARCODES -O $OUT/$SAMPLE -R $VCF -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

singularity exec Demuxafy.sif bcftools view $VCF -R $OUT/$SAMPLE/cellSNP.base.vcf.gz -Ov -o $OUT/$SAMPLE/donor_subset.vcf
singularity exec Demuxafy.sif vireo -c $OUT/$SAMPLE -d $OUT/$SAMPLE/donor_subset.vcf -o $OUT/$SAMPLE -t "GP"