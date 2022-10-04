##############################################################################
# Author: Drew Neavin
# Date: 2022-10-04
# Description: This script is to QC the 18-line hiPSC village using the demultiplexing results, doublet detecting results and typical QC metrics (mt%)
##############################################################################

#### loading libraries ####
library(scran)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(Seurat)
library(jcolors)
library(cowplot)
library(RColorBrewer)
library(readr)
library(purrr)
library(clustree)
library(reticulate)
library(BiocParallel)
library(coop)
library(harmony)
library(awtools)
library(Nebulosa)

#### Setting up directories #####
print("Set up directories")
dir <- '/path/to/output/hiPSC_village_18_lines/iPSC_maintenance'
outdir <- paste0(outdir, "QC/")


### Set variables ###
sample_meta <- fread("iPSC_Village_Publication/scripts/hiPSC_village_18_lines/Pool_meta_data.tsv", sep = "\t") # available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication
pools <- sample_meta$Pool


## Read in expression data
counts_list <- lapply(pools, function(x){
    print(x)
    Read10X(paste0("/path/to/Expression/", pools, "/filtered_feature_bc_matrix/"), gene.column = 1)
})
names(counts_list) <- pools


## Add poolnames to cell names so can easily match back if there are any cells with same barcoe
counts_list <- lapply(names(counts_list), function(x){
    colnames(counts_list[[x]]) <- gsub("-1", "", paste0(x, "_", colnames(counts_list[[x]])))
    return(counts_list[[x]])
})
names(counts_list) <- pools


## Make seurat object ##
seurat_list <- lapply(names(counts_list), function(x){
    temp <- CreateSeuratObject(counts = counts_list[[x]])
    temp@meta.data$Passage <- x
    return(temp)
})
names(seurat_list) <- pools


seurat <- merge(seurat_list[[1]], y = seurat_list[2:length(seurat_list)],project = "Village_Phase1_multi-passage")



### Add gene ID data to gene level metadata ###
features <- fread(paste0(datadir,pools[1], "/GE/",pools[1],"/outs/per_sample_outs/",pools[1],"/count/sample_feature_bc_matrix/features.tsv.gz"), col.names = c("ENSG", "Gene_ID", "Assay"))
features$Assay <- NULL

features_df <- data.frame(features)
rownames(features_df) <- features$ENSG
features_df$ENSG <- NULL

seurat[["RNA"]] <- AddMetaData(seurat[["RNA"]], features_df)


## Add QC metrics
RbGeneList <- read.delim(file = "RibosomalGeneList_GeneID_ENSG.txt") # available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication
MtGeneList <- read.delim(file = "MtGeneList_GeneID_ENSG.txt") # available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication
print("Calculating Mt %")
seurat <- PercentageFeatureSet(seurat, features = MtGeneList$ENSG, col.name = "percent.mt")
print("Calculating Rb %")
RbGeneList <- RbGeneList[which(RbGeneList$ENSG %in% rownames(seurat)),]
seurat <- PercentageFeatureSet(seurat, features = RbGeneList$ENSG, col.name = "percent.rb")


## Add Demultiplexing Results ##
demultiplexing_list <- lapply(pools, function(pool){
    temp <- data.frame(fread(paste0(dir, "output/multi-passage/demultiplexed/with_demuxlet/updated_2022_06_26/atleasthalf_singlet/", pool, "/atleasthalf_singlet_w_combined_assignments.tsv")))
    rownames(temp) <- gsub("-1", "", paste0(pool, "_", temp$Barcode))
    return(temp)
})

demultiplexing_dt <- do.call(rbind,demultiplexing_list)


seurat <- AddMetaData(seurat, demultiplexing_dt)


print("Starting cell cycle step")
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))

countsENSG <- seurat[["RNA"]]@counts # Create a counts matrix that has the ENSG gene clasifiers so can determine cell cycle phase
print("The size of the counts object in genes x cells is:")
print(dim(countsENSG))

### Run the cell cycle identification ###
print("Starting cell cycle determination")
assigned <- cyclone(countsENSG, pairs=hs.pairs)  #Note, this takes hours
table(assigned$phases)
write.table(assigned, file = paste0(outdir,"CellCycleProportions.txt"), quote = F, sep = "\t") #Save so that can read in and don't have to wait to recompute again


assigned <- read.table(paste0(outdir,"CellCycleProportions.txt"), sep = "\t")
assigned <- as.data.frame(assigned)
rownames(assigned) <- colnames(seurat)
write.table(assigned, file = paste0(outdir,"CellCycleProportions.txt"), quote = F, sep = "\t") #Save so that can read in and don't have to wait to recompute again
assigned <- read.table(paste0(outdir,"CellCycleProportions.txt"), sep = "\t")



##### Subset seurat object by cells that have no info for cell assignments (removed by dropletQC) #####
seurat_sub <- subset(seurat, subset = AtLeastHalfSinglet_Individual_Assignment != "doublet")


##### Remove high mitochondrial % cells #####
seurat_sub_mt <- subset(seurat_sub, subset = percent.mt < 25)


### Integrate different time-points
seurat_sub_mt_list <- SplitObject(seurat_sub_mt, split.by = "Pool")
seurat_sub_mt_list <- lapply(X = seurat_sub_mt_list, function(x) SCTransform(x, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt", "percent.rb"), return.only.var.genes = FALSE))
features <- SelectIntegrationFeatures(object.list = seurat_sub_mt_list, nfeatures = 3000)
seurat_sub_mt_list <- PrepSCTIntegration(object.list = seurat_sub_mt_list, anchor.features = features)

seurat_sub_mt_anchors <- FindIntegrationAnchors(object.list = seurat_sub_mt_list, normalization.method = "SCT",
    anchor.features = features)
combined_sct <- IntegrateData(anchorset = seurat_sub_mt_anchors, normalization.method = "SCT")

combined_sct <- RunPCA(combined_sct, verbose = FALSE)
combined_sct <- RunUMAP(combined_sct, reduction = "pca", dims = 1:30)

saveRDS(combined_sct, paste0(outdir,"seurat_qc_integrated.rds"))



### Filter for 1% expressing ###
feats <- rownames(combined_sct[["SCT"]]@counts)[which((rowSums(combined_sct[["SCT"]]@counts > 0)/ncol(combined_sct[["SCT"]]@counts)) >= 0.01)]
combined_sct_filt <- subset(combined_sct, features = feats)

saveRDS(combined_sct_filt, paste0(outdir, "time-integrated_filtered_seurat_1pct_expressing.rds"))
fwrite(data.table(Gene = rownames(combined_sct_filt)), paste0(outdir, "time-integrated_filtered_seurat_1pct_expressing_genes.tsv"), sep = "\t")
