#!/bin/bash

OUTDIR="/path/to/output/hiPSC_village_18_lines/iPSC_maintenance"
mkdir -p $OUTDIR


deboever_genes="$OUTDIR/deboever_gene_snp_list.tsv"
variance_genes="$OUTDIR/Variance/gene_separated/residuals4qtl/"


ls $variance_genes | sed "s/_residuals4qtl.rds//g" | sort -u > $OUTDIR/residual_genes.tsv
awk '{print $2}' $deboever_genes | sort -u > $OUTDIR/deboever_genes.tsv


# comm -12 $OUTDIR/residual_genes.tsv head $OUTDIR/kilpinen_genes.tsv > $OUTDIR/genes_residual_and_kilpinen.tsv ## NONE
comm -12 $OUTDIR/residual_genes.tsv $OUTDIR/deboever_genes.tsv > $OUTDIR/genes_residual_and_deboever.tsv


### Get finalized list of genes that were eQTLs and demonstrate line effects
grep -F -f $OUTDIR/genes_residual_and_deboever.tsv $deboever_genes > $OUTDIR/finalized_deboever_gene_snp_list.tsv



### Filter vcf for just these snps ###
awk '{print $1}' $OUTDIR/finalized_deboever_gene_snp_list.tsv | sort -u | sed 's/:/\t/g' | sed 's/$/\t/g' > $OUTDIR/finalized_deboever_snps.tsv

grep "#" $OUTDIR/snps_18_hiPSC_noCNV_INV_filtered_diff_genotypes_hg19.vcf > $OUTDIR/deboever_finalized_snps.vcf
grep -f $OUTDIR/finalized_deboever_snps.tsv $OUTDIR/snps_18_hiPSC_noCNV_INV_filtered_diff_genotypes_hg19.vcf >> $OUTDIR/deboever_finalized_snps.vcf
