#!/bin/bash
##############################################################################
# Author: Drew Neavin
# Date: 2022-10-03
# Description: This script is to run demuxlet on the 18-line hiPSC village with Demuxafy (https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/). 
##############################################################################


SAMPLE_INFO="iPSC_Village_Publication/scripts/hiPSC_village_18_lines/iPSC_maintenance/Sample_meta.tsv" # (accessible from the iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication)
SAMPLE=$(tail -n +2 $SAMPLE_INFO | tr '\n' ':' | cut -d ':' -f $SGE_TASK_ID)

VCF="/path/to/output/vcf/snps_18_hiPSC_r2_maf0.3_exons_hg38.vcf" # Exons only

BAM=$DATA/$SAMPLE/outs/possorted_genome_bam.bam ### Available by processing fastq data through cellranger or on request from d.neavin @ garvan.org.au
BARCODES="$DATA/$SAMPLE/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" ### Available by processing fastq data through cellranger or on request from d.neavin @ garvan.org.au


mkdir -p mkdir -p $OUT/$SAMPLE


singularity exec Demuxafy.sif popscle demuxlet --plp $OUT/$SAMPLE/pileup --vcf $VCF --field "GP" --group-list $BARCODES --geno-error-coeff 1.0 --geno-error-offset 0.05 --out $OUT/$SAMPLE/demuxlet