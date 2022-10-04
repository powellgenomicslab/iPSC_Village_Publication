library(Seurat)
library(SeuratDisk)


##### Setting up Directories
datadir <- "/path/to/output/Integrated/"
outdir <- "/path/to/output/RNA_Velocity_Pseudotime/preprocess"

dir.create(outdir, recursive = TRUE)


##### Read in seurat object
seurat <- readRDS(paste0(datadir,"village_3_hiPSC_lines_high_quality_cells_integrated.rds")) ## Produced with "hiPSC_village_3_lines_QC.R" script from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication

write.csv(Cells(seurat), file = paste0(outdir, "cellID_obs.csv"), row.names = FALSE)
write.csv(Embeddings(seurat, reduction = "umap"), file = paste0(outdir,"cell_embeddings.csv"))
write.csv(seurat@meta.data, file = paste0(outdir,"metadata.csv"))

dir.create(paste0(outdir,"h5ad/"), recursive = TRUE)
SaveH5Seurat(seurat, filename = paste0(outdir,"village_3_hiPSC_lines_high_quality_cells_integrated.h5Seurat"))

Convert(paste0(outdir,"village_3_hiPSC_lines_high_quality_cells_integrated.h5Seurat"), dest = "h5ad")

