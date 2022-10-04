library(tidyverse)
library(haven)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(Seurat)
library(data.table)




dir <- "/path/to/output/"
outdir <- paste0(dir,"RNA_Velocity_Pseudotime/Variance_Explained/")
dir.create(outdir, recursive = TRUE)



##### Read in Data #####
### seurat object ###
seurat <- readRDS(paste0(dir,"village_3_hiPSC_lines_high_quality_cells_integrated.rds"))

### velocity metadata ###
velo_meta <- fread("/path/to/output/RNA_Velocity_Pseudotime/scvelo/metadata.csv", sep = ",")
rownames(velo_meta) <- velo_meta$V1
velo_meta$V1 <- NULL

velo_meta_sub <- velo_meta[,c("n_unspliced_counts", "latent_time")]
rownames(velo_meta_sub) <- rownames(velo_meta)

seurat <- AddMetaData(seurat, velo_meta_sub)


seurat@meta.data$Location <- gsub("_Baseline", "", seurat@meta.data$Location) %>% gsub("_Village.+", "", .) %>% gsub("Thawed", "Cryopreserved",.)

seurat_noNA <- subset(seurat, subset = latent_time >= 0)


seurat_noNA@meta.data$Location <- gsub("_Baseline", "", seurat_noNA@meta.data$Location) %>% gsub("_Village.+", "", .) %>% gsub("Thawed", "Cryopreserved",.)
seurat_noNA@meta.data$Time <- gsub("Thawed Village Day 0", "Baseline", seurat_noNA@meta.data$Time) %>% gsub("Thawed Village Day 7", "Village", .) %>% gsub(" Day 4", "", .)
seurat_noNA@meta.data$Cryopreserved <-ifelse(seurat_noNA@meta.data$Location == "Sydney_Cryopreserved", "Cryopreserved", "Fresh")
seurat_noNA@meta.data$Location <- gsub("_Cryopreserved", "", seurat_noNA@meta.data$Location)


saveRDS(seurat_noNA, paste0(outdir, "village_3_hiPSC_lines_high_quality_cells_integrated_1pct_expressing.rds"))
seurat_noNA <- readRDS(paste0(outdir, "village_3_hiPSC_lines_high_quality_cells_integrated_1pct_expressing.rds"))


fwrite(data.table(Gene = rownames(seurat_noNA)), paste0(outdir,"village_3_hiPSC_lines_high_quality_cells_integrated_1pct_expressing.tsv"))