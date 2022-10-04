library(Seurat)
liberary(tidyverse)
library(dplyr)
library(data.table)



##### Set up Directories #####
dir <- '/path/to/output/'



##### Read in data #####
seurat_integrated <- readRDS(paste0(dir, "Integrated/village_3_hiPSC_lines_high_quality_cells_integrated.rds"))


##### Create Variance Explained plot for fresh samples with at least 1% expressing #####
seurat_fresh <- subset(seurat_integrated, subset = Cryopreservation != "Cryopreserved")
seurat_fresh_1pct <- subset(seurat_fresh, features = rownames(seurat_fresh)[which(rowSums(seurat_fresh[["SCT"]]@counts > 0)/ncol(seurat_fresh[["SCT"]]@counts) >= 0.01)])
seurat_fresh_1pct <- SCTransform(seurat_fresh_1pct, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt", "percent.rb"), return.only.var.genes = FALSE)

saveRDS(seurat_fresh_1pct, paste0(dir,"variance_partition_fresh/seurat_integrated_noncryo_1pct_expressing.rds"))
fwrite(data.table(Gene = rownames(seurat_fresh_1pct)), paste0(dir,"variance_partition_fresh/seurat_integrated_noncryo_1pct_expressing_genes.tsv", sep = "\t"))


##### Create Variance Explained plot with at least 1% expressing #####
seurat_cryo <- subset(seurat_integrated, subset = Location == "Sydney")
seurat_cryo_1pct <- subset(seurat_cryo, features = rownames(seurat_cryo)[which(rowSums(seurat_cryo[["SCT"]]@counts > 0)/ncol(seurat_cryo[["SCT"]]@counts) >= 0.01)])
seurat_cryo_1pct <- SCTransform(seurat_cryo_1pct, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt", "percent.rb"), return.only.var.genes = FALSE)

saveRDS(seurat_cryo_1pct, paste0(dir,"variance_partition_cryo/seurat_integrated_cryopreserved_1pct_expressing.rds"))
fwrite(data.table(Gene = rownames(seurat_cryo_1pct)), paste0(dir,"variance_partition_cryo/seurat_integrated_cryopreserved_1pct_expressing_genes.tsv", sep = "\t"))
