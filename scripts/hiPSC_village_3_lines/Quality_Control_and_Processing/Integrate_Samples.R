##############################################################################
# Author: Drew Neavin
# Date: 2022-10-04
# Description: Integrate the samples from different times and sites for visualization.
##############################################################################


##### Load libraries #####
library(Seurat)
library(tidyverse)
library(ggplot2)
library(data.table)



##### Set up directories #####
dir <- '/path/to/output/'
outdir <- paste0(dir,"Integrated/")
dir.create(outdir)



##### Read in Data #####
seurat_qc <- readRDS(paste0(dir, "QC/village_3_hiPSC_lines_high_quality_cells.rds")) ## Producerd by "hiPSC_village_3_lines_QC.R" from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication


##### Cluster to identify poor cells to additionally remove
##### Separate by cell line at each site #####
seurat_subset_list <- list()
for (site_time in unique(seurat_qc@meta.data$Location_Time)){
    seurat_sub <- subset(seurat_qc, subset = Location_Time == site_time)
	for (cell_line in unique(seurat_sub@meta.data$Final_Assignment)){
		seurat_subset_list[[paste0(site_time,"_",cell_line)]] <- subset(seurat_sub, subset = Final_Assignment == cell_line)
	}
}


seurat_subset_list <- lapply(seurat_subset_list, function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})


features <- SelectIntegrationFeatures(object.list = seurat_subset_list)
seurat_subset_list <- lapply(X = seurat_subset_list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt"))
    x <- RunPCA(x, features = features, verbose = FALSE)
})


anchors <- FindIntegrationAnchors(object.list = seurat_subset_list, reference = c(1, 2), reduction = "rpca", 
    dims = 1:30)
seurat_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)


seurat_integrated <- ScaleData(seurat_integrated, verbose = FALSE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt"))
seurat_integrated <- RunPCA(seurat_integrated, verbose = FALSE)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:30)
seurat_integrated <- FindNeighbors(seurat_integrated, resolution = 0.08)

### Cluster 4 and 5 have high and low MT% -> remove these clusters for downstream analysis ###
seurat_integrated <- subset(seurat_integrated, idents < 4) ### removes 584 cells



##### Rerun scaling and Normalization #####
seurat_integrated <- ScaleData(seurat_integrated, verbose = FALSE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt"))
seurat_integrated <- RunPCA(seurat_integrated, verbose = FALSE)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:30)
seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "pca", assay = "integrated", dims = 1:20)
seurat_integrated <- SCTransform(seurat_integrated, verbose = TRUE, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt"), return.only.var.genes = FALSE)


saveRDS(seurat_integrated, paste0(outdir, "village_3_hiPSC_lines_high_quality_cells_integrated.rds"))