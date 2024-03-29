library(Seurat)
library(tivdyverse)
library(data.table)


dir.create("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Kilpinen_eQTLs/")
dir.create("/path/to/output/eQTL_Comparison/uni_village/Cryopreserved/")



seurat <- readRDS( "/path/to/output/variance_partition_cryo/seurat_integrated_cryo_1pct_expressing.rds") ## Prepared by 'iPSC_Village_Publication/scripts/hiPSC_village_3_lines/Variance_Explained/prepare_data.R' script from the github: https://github.com/powellgenomicslab/iPSC_Village_Publication


### Make DF for modeling ###
df_hier_unscale <- data.frame("Village" = as.factor(ifelse(seurat@meta.data$Time == "Baseline", 0, 1)), "Line" = seurat@meta.data$Final_Assignment, "Replicate" = as.factor(gsub("[A-Z][a-z]+", "", seurat@meta.data$MULTI_ID)), "Cryopreserved" = seurat$Cryopreserved)
df_hier_unscale$Barcode <- rownames(df_hier_unscale)

fwrite(df_hier_unscale, "/path/to/output/eQTL_Comparison/uni_village/Cryopreserved/cell_meta.tsv", sep = "\t")

