#!/bin/bash


vcf="iPSC_Village_Publication/data/vcf/snps_3_hiPSC.vcf.gz" # available from publication
deboever="iPSC_Village_Publication/scripts/hiPSC_village_3_lines/eQTL_Comparison/DeBoeveriPSCeQTLs.txt" ## Downloaded from DeBoever et al publication (2017, Cell Stem Cell) and provided in the iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication

OUT="/path/to/output/eQTL_Comparison"


INTERSECT_OUT="$OUT/uni_village"

mkdir -p $INTERSECT_OUT


### Filter on individuals and MAF for our lines ###
singularity exec Demuxafy.sif bcftools filter --include 'R2>=0.3' $vcf > /path/to/output/vcf/snps_3_hiPSC_R2_0.3.vcf

VCF_R2="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/SNPgenotypes/"

sed '/<CN/d' /path/to/output/vcf/snps_3_hiPSC_R2_0.3.vcf > $INTERSECT_OUT/snps_3_hiPSC_R2_0.3.vcf
sed -i '/<INV>/d' $INTERSECT_OUT/snps_3_hiPSC_R2_0.3.vcf


conda activate vcftools

    vcftools --vcf $INTERSECT_OUT/snps_3_hiPSC_R2_0.3.vcf --recode --recode-INFO-all --mac 1 --out $INTERSECT_OUT/snps_3_hiPSC_R2_0.3_filtered

    bcftools view -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES || COUNT(GT="RA")=N_SAMPLES' $INTERSECT_OUT/snps_3_hiPSC_R2_0.3_filtered.recode.vcf > $INTERSECT_OUT/snps_3_hiPSC_R2_0.3_filtered_diff_genotypes.vcf

conda deactivate



##### DEBOEVER DATA #####
### Intersect with these SNPs ##$#
conda activate bedtools

    bedtools intersect -a $OUT/deboever.bed -b $INTERSECT_OUT/snps_3_hiPSC_R2_0.3_filtered_diff_genotypes.vcf -wa -wb > $INTERSECT_OUT/deboever_imputed_overlapping_filtered.bed

conda deactivate


head -n 1 $OUT/deboever.bed > $INTERSECT_OUT/deboever_header_bed.tsv
grep "#CHROM" $INTERSECT_OUT/snps_3_hiPSC_R2_0.3_filtered_diff_genotypes.vcf > $INTERSECT_OUT/deboever_header_vcf.tsv
paste -d"\t" $INTERSECT_OUT/deboever_header_bed.tsv $INTERSECT_OUT/deboever_header_vcf.tsv > $INTERSECT_OUT/deboever_combined_header.tsv

cat $INTERSECT_OUT/deboever_combined_header.tsv > $INTERSECT_OUT/deboever_imputed_overlapping_filtered_header.bed
cat $INTERSECT_OUT/deboever_imputed_overlapping_filtered.bed >> $INTERSECT_OUT/deboever_imputed_overlapping_filtered_header.bed


### Use this file to test for trends in gene
awk 'BEGIN{FS=OFS="\t"}{print($31, $15, $25)}' $INTERSECT_OUT/deboever_imputed_overlapping_filtered_header.bed | sed 's/\..*\t/\t/g' > $INTERSECT_OUT/deboever_gene_snp_list.tsv



