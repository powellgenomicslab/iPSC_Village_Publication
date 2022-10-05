#!/bin/bash

DIR="/path/to/output/eQTL_Comparison/uni_village/Fresh"
OUTDIR="$DIR/uniculture"

mkdir -p $OUTDIR

deboever_genes="$DIR/deboever_gene_snp_list.tsv"
variance_genes="/path/to/output/eQTL_Comparison/variance_partition_fresh/uniculture/gene_separated/residuals4qtl"


ls $variance_genes | sed "s/_residuals4qtl.rds//g" | sort -u > $OUTDIR/residual_genes.tsv
awk '{print $2}' $deboever_genes | sort -u > $OUTDIR/deboever_genes.tsv


comm -12 $OUTDIR/residual_genes.tsv $OUTDIR/deboever_genes.tsv > $OUTDIR/genes_residual_and_deboever.tsv


### Get finalized list of genes that were eQTLs and demonstrate line effects
grep -F -f $OUTDIR/genes_residual_and_deboever.tsv $deboever_genes > $OUTDIR/finalized_deboever_gene_snp_list.tsv



### Filter vcf for just these snps ###
awk '{print $1}' $OUTDIR/finalized_deboever_gene_snp_list.tsv | sort -u > $OUTDIR/finalized_deboever_snps.tsv


conda activate vcftools

    vcftools --vcf $DIR/snps_3_hiPSC_R2_0.3_filtered_diff_genotypes.vcf --snps $OUTDIR/finalized_deboever_snps.tsv --recode --recode-INFO-all --out $OUTDIR/deboever_finalized_snps

conda deactivate

