library(vcfR)
library(data.table)
library(tidyverse)


##### Set up variables #####
dir <- "/path/to/output/eQTL_Comparison/Cryopreserved/uni_village/"


##### Read in data #####
vcf <- read.vcfR(paste0(dir, "uniculture/deboever_finalized_snps.recode.vcf"))
eqtls_overlap_dt <- fread(paste0(dir,"deboever_imputed_overlapping_filtered_header.bed"), sep = "\t")

geno <- data.table(extract.gt(element = "GT",vcf, IDtoRowNames = F))


if (!all(colSums(is.na(geno)) == nrow(geno))){
    message("Found GT genotype format in cluster vcf. Will use that metric for cluster correlation.")
    format_clust = "GT"

    if (any(grepl("\\|",geno[,1]))){
        separator = "|"
        message("Detected | separator for GT genotype format in cluster vcf")
    } else if (any(grepl("/",geno[,1]))) {
        separator = "/"
        message("Detected / separator for GT genotype format in cluster vcf")
    } else {
        format_clust = NA
        message("Can't identify a separator for the GT field in cluster vcf, moving on to using GP.")
    }
    if (!is.na(format_clust)){
        geno <- data.table(as_tibble(lapply(geno, function(x) {gsub(paste0("0\\",separator,"0"),0, x)}) %>%
                                lapply(., function(x) {gsub(paste0("0\\",separator,"1"),1, x)}) %>%
                                lapply(., function(x) {gsub(paste0("1\\",separator,"0"),1, x)}) %>%
                                lapply(., function(x) {gsub(paste0("1\\",separator,"1"),2, x)})))
    }
}

#### Add ID of the variants ####
geno$ID <- vcf@fix[,'ID']
geno$ID_ref_alt <- paste0(vcf@fix[,'CHROM'],":", vcf@fix[,'POS'],"_", vcf@fix[,'REF'],"_", vcf@fix[,'ALT'])

eqtls_overlap_dt <- eqtls_overlap_dt[ID %in% unique(geno$ID)]
eqtls_overlap_dt$ID_ref_alt <- paste0(eqtls_overlap_dt$ID, "_", eqtls_overlap_dt$REF, "_", eqtls_overlap_dt$ALT)


#### Only keep snps that are non-identical across the hiPSC lines ####
temp <- list()

for (gene in unique(eqtls_overlap_dt$gene_id)){
    print(gene)
    temp[[gene]] <- eqtls_overlap_dt[gene_id == gene]
    k = 1
    while (k < nrow(temp[[gene]])){
        test_snp <- temp[[gene]]$ID_ref_alt[k]
        for (snp in temp[[gene]][(k+1):nrow(temp[[gene]]),]$ID_ref_alt){
            if (abs(cor(as.numeric(geno[ID_ref_alt == test_snp,1:3]), as.numeric(geno[ID_ref_alt == snp,1:3]))) == 1){
                temp[[gene]] <- temp[[gene]][ID_ref_alt != snp]
            } 
        }
        k = k + 1
    }
}


eqtls_overlap_dt_subset <- do.call(rbind, temp)
eqtls_overlap_dt_subset$gene_id <- gsub("\\..+", "", eqtls_overlap_dt_subset$gene_id)


##### Keep only the genes whose variance was explained by the Line effect #####
resid_files <- list.files("/path/to/output/eQTL_Comparison/variance_partition_cryo/uniculture/gene_separated/residuals4qtl")
ensg_list <- unique(gsub("_residuals4qtl.rds", "", resid_files))

eqtls_overlap_dt_subset <- eqtls_overlap_dt_subset[gene_id %in% ensg_list]


fwrite(eqtls_overlap_dt_subset, paste0(dir, "uniculture/deboever_imputed_overlapping_filtered_header_pruned.bed"), sep = "\t")


eqtls_overlap_dt_subset_gene_snp <- eqtls_overlap_dt_subset[,c("gene_id","ID_ref_alt")]

fwrite(eqtls_overlap_dt_subset_gene_snp, paste0(dir, "uniculture/deboever_imputed_overlapping_filtered_header_pruned_snp_gene.tsv"), sep = "\t")


