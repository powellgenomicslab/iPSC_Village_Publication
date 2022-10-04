#!/bin/bash
##############################################################################
# Author: Drew Neavin
# Date: 2022-10-03
# Description: This script is to filter reference SNP genotypes (Imputed) by R^2, MAF, overlappying exons with Demuxafy (https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/). 
#   and then to lift to hg38 with the SNP_imputation_1000g_hg38.sif image from the SNP imputation pipeline on hg38 (https://github.com/powellgenomicslab/SNP_imputation_1000g_hg38/wiki/SNP-Genotype-Imputation-Using-1000G-hg38-Reference)
##############################################################################


VCF="iPSC_Village_Publication/data/vcf/snps_3_hiPSC.vcf.gz" # available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication
OUTDIR="/path/to/output/vcf"

mkdir -p $OUTDIR


##### Filter the vcf for:
### R^2 >= 0.05
### MAF >= 0.05
### biallelic snps
### overlapping indesl
singularity exec Demuxafy.sif bcftools filter --include 'MAF>=0.05 & R2>=0.3' --regions-file hg19exonsUCSC.bed $VCF | \
    singularity exec Demuxafy.sif bcftools view --exclude-types indels --max-alleles 2 -Oz --output $OUTDIR/snps_3_hiPSC_r2_maf0.3_exons.vcf.gz


##### Lift to hg38 using CrossMap
singularity exec SNP_imputation_1000g_hg38.sif CrossMap.py /opt/GRCh37_to_GRCh38.chain $OUTDIR/snps_3_hiPSC_r2_maf0.3_exons.vcf.gz $OUTDIR/snps_3_hiPSC_r2_maf0.3_exons_hg38.vcf
