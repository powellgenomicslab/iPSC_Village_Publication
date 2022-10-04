##############################################################################
# Author: Drew Neavin
# Date: 2022-10-03
# Description: This script is to QC the 3-line hiPSC village using the hashing barcodes, demultiplexing results, doublet detecting results and typical QC metrics (mt%)
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
dir <- '/path/to/output/'
outdir <- paste0(outdir, "QC/")


### Set variables ###
sample_meta <- fread("iPSC_Village_Publication/scripts/hiPSC_village_3_lines/Quality_Control_and_Processing/Demultiplexing_Doublet_Detecting/Sample_meta.tsv", sep = "\t") # available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication
pools <- sample_meta$Pool
hiPSC_lines <- c("FSA0006","MBE1006","TOB0421")



### Set up Functions ###
mad_function <- function(seurat, column, number_mad){
    mad <- mad(seurat@meta.data[,column])
    low <- median(seurat@meta.data[,column]) - number_mad*mad
    high <- median(seurat@meta.data[,column]) + number_mad*mad
    print("The mad is:")
    print(mad)
    print("The lower bound is:")
    print(low)
    print("The upper bound is:")
    print(high)
    seurat@meta.data[,paste0(column,"_mad")] <- ifelse((seurat@meta.data[,column] > low & seurat@meta.data[,column] < high),"NotOutlier", "Outlier")
    return(seurat)
}


##### Read in and aggregate gene expression and hashing data #####
dirs10x <- paste0("/path/to/Expression/", pools, "/filtered_feature_bc_matrix/")
dirs_hash <- paste0("path/to/Hashing/", pools, "/umi_count/")



## Read in expression data
counts_list <- lapply(dirs10x, function(x){
    Read10X(x, gene.column = 1)
})
names(counts_list) <- pools

## Add poolnames to cell names so can easily match back if there are any cells with same barcoe
counts_list <- lapply(names(counts_list), function(x){
    colnames(counts_list[[x]]) <- paste0(x, "_", colnames(counts_list[[x]]))
    return(counts_list[[x]])
})
names(counts_list) <- pools


### Read in the hashing data ###
hash_list <- lapply(dirs_hash, function(x){
    Read10X(x, gene.column =1)
})
names(hash_list) <- pools

hash_list <- lapply(hash_list, function(x){
    x <- as.data.frame(t(as.data.frame(x)))
    x$Barcode <- rownames(x)
    return(x)
})
names(hash_list) <- pools

hash_list <- lapply(names(hash_list), function(x){
    tmp <- as.data.frame(t(hash_list[[x]]))
    colnames(tmp) <- hash_list[[x]]$Barcode
    tmp <- tmp[1:(nrow(tmp)-1),]
})
names(hash_list) <- pools

### Add pool names to barcodes so all are unique ###
hash_list <- lapply(names(hash_list), function(x){
    colnames(hash_list[[x]]) <- paste0(x, "_", colnames(hash_list[[x]]))
    return(hash_list[[x]])
})
names(hash_list) <- pools

### Remove the unmapped row (creates noise when identifying doublets/singlets)
hash_list <- lapply(hash_list, function(x){
    x[-c(which(rownames(x) == "unmapped")),]
})

### Change the hashtag names for easy interpretation
hash_list6 <- lapply(hash_list[1:6], function(x){
    rownames(x) <- gsub("A0252-TGATGGCCTATTGGG","Brisbane2", rownames(x)) %>%
        gsub("A0253-TTCCGCCTCTCTTTG", "Brisbane3", .) %>%
        gsub("A0254-AGTAAGTTCAGCGTA", "Sydney1", .) %>%
        gsub("A0255-AAGTATCGTTTCGCA", "Sydney2", .) %>%
        gsub("A0256-GGTTGCCAGATGTCA", "Sydney3", .) %>%
        gsub("A0257-TGTCTTTCCTGCCAG", "Melbourne1", .) %>%
        gsub("A0258-CTCCTCTGCAATTAC", "Melbourne2", .)
    return(x)
})

hash_list <- c(hash_list6, hash_list[7:8])

rownames(hash_list[["DRENEA_1"]]) <- gsub("A0251-GTCAACTCTTTAGCG","Brisbane1", rownames(hash_list[["DRENEA_1"]]))
rownames(hash_list[["DRENEA_4"]]) <- gsub("A0251-GTCAACTCTTTAGCG","Brisbane1", rownames(hash_list[["DRENEA_4"]]))
rownames(hash_list[["DRENEA_3"]]) <- gsub("A0251-GTCAACTCTTTAGCG","Melbourne3", rownames(hash_list[["DRENEA_3"]]))
rownames(hash_list[["DRENEA_6"]]) <- gsub("A0251-GTCAACTCTTTAGCG","Melbourne3", rownames(hash_list[["DRENEA_6"]]))

rownames(hash_list[["Village_A_Baseline"]]) <- gsub("A0251-GTCAACTCTTTAGCG","Sydney1", rownames(hash_list[["Village_A_Baseline"]]))
rownames(hash_list[["Village_A_Baseline"]]) <- gsub("A0252-TGATGGCCTATTGGG","Sydney2", rownames(hash_list[["Village_A_Baseline"]]))
rownames(hash_list[["Village_A_Baseline"]]) <- gsub("A0253-TTCCGCCTCTCTTTG","Sydney3", rownames(hash_list[["Village_A_Baseline"]]))
rownames(hash_list[["Village_B_1_week"]]) <- gsub("A0254-AGTAAGTTCAGCGTA","Sydney1", rownames(hash_list[["Village_B_1_week"]]))
rownames(hash_list[["Village_B_1_week"]]) <- gsub("A0255-AAGTATCGTTTCGCA","Sydney2", rownames(hash_list[["Village_B_1_week"]]))
rownames(hash_list[["Village_B_1_week"]]) <- gsub("A0256-GGTTGCCAGATGTCA","Sydney3", rownames(hash_list[["Village_B_1_week"]]))

## Get just the cells that are also in the hash list ##
counts_list <- lapply(names(counts_list), function(x){
    print(x)
    colnames(counts_list[[x]]) <- gsub("-1", "", colnames(counts_list[[x]]))
    counts_list[[x]] <- counts_list[[x]][,which(colnames(counts_list[[x]]) %in% colnames(hash_list[[x]]))]
    return(counts_list[[x]])
})



##### Doublet Detection #####
seurat_list_norm <- lapply(counts_list, function(x){
    CreateSeuratObject(counts = x)
})
names(seurat_list_norm) <- pools

seurat_list_norm <- lapply(seurat_list_norm, function(x){
    tmp <- NormalizeData(x, verbose = TRUE)
    tmp <- FindVariableFeatures(tmp, selection.method = "mean.var.plot")
    tmp <- ScaleData(tmp, features = VariableFeatures(tmp))
    return(tmp)
})
names(seurat_list_norm) <- pools

### Add hashtag data ###
seurat_list_norm <- lapply(names(seurat_list_norm), function(x){
    seurat_list_norm[[x]][["HTO"]] <- CreateAssayObject(counts = hash_list[[x]])
    return(seurat_list_norm[[x]])
})
names(seurat_list_norm) <- pools

seurat_list_norm <- lapply(seurat_list_norm, function(x){
    NormalizeData(x, assay = "HTO", normalization.method = "RC")
})


seurat_list_norm <- lapply(seurat_list_norm, function(x){
    MULTIseqDemux(x, assay = "HTO", autoThresh = TRUE)
})


##### Get pool names from the barcodes #####
seurat_norm@meta.data$Pool <- gsub("_[A-Z]{16}","",rownames(seurat_norm@meta.data)) %>% gsub("-1","",.)

### Fix barcodes so they are all the same (freeze-thaw have -1 at the end)
barcode_key <- data.frame("Original_Barcode" = colnames(seurat_norm), "Updated_Barcode" = gsub("-1", "", colnames(seurat_norm)))

seurat_norm <- RenameCells(seurat_norm, new.names = barcode_key$Updated_Barcode)

### Update some meta.data columns ###
seurat_norm@meta.data$Time <- ifelse((seurat_norm@meta.data$Pool == "DRENEA_1" | 
                                        seurat_norm@meta.data$Pool == "DRENEA_2" |
                                        seurat_norm@meta.data$Pool == "DRENEA_3"), "Uni-Culture", 
                                ifelse((seurat_norm@meta.data$Pool == "DRENEA_4" | 
                                        seurat_norm@meta.data$Pool == "DRENEA_5" |
                                        seurat_norm@meta.data$Pool == "DRENEA_6"),"Village",
                                ifelse((seurat_norm@meta.data$Pool == "Village_A_Baseline"), "Uni-Culture",
                                ifelse((seurat_norm@meta.data$Pool == "Village_B_1_week"), "Village", NA))))

seurat_norm@meta.data$Cryopreservation <- ifelse((seurat_norm@meta.data$Pool == "DRENEA_1" | 
                                        seurat_norm@meta.data$Pool == "DRENEA_2" |
                                        seurat_norm@meta.data$Pool == "DRENEA_3" |
                                        seurat_norm@meta.data$Pool == "DRENEA_4" | 
                                        seurat_norm@meta.data$Pool == "DRENEA_5" |
                                        seurat_norm@meta.data$Pool == "DRENEA_6"),"Fresh", "Cryopreserved")


### Add gene names to feature metadata
genes_list <- lapply(rev(pools), function(x){
    read.delim(paste0("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/",x,"/outs/filtered_feature_bc_matrix/features.tsv.gz"), header = FALSE, col.names = c("ENSG","Gene_ID", "Descriptor"))
})

genes <- do.call(rbind, genes_list)

genes_unique <- rownames(seurat_norm)
print(length(genes_unique))
print(length(unique(genes_unique)))

genes_unique_df <- genes[match(genes_unique, genes$ENSG),]
rownames(genes_unique_df) <- genes_unique_df$ENSG
genes_unique_df$Descriptor <- NULL

seurat_norm[["RNA"]] <- AddMetaData(seurat_norm[["RNA"]], genes_unique_df)



## Add QC metrics
RbGeneList <- read.delim(file = "RibosomalGeneList_GeneID_ENSG.txt") # available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication
MtGeneList <- read.delim(file = "MtGeneList_GeneID_ENSG.txt") # available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication
print("Calculating Mt %")
seurat_norm <- PercentageFeatureSet(seurat_norm, features = MtGeneList$ENSG, col.name = "percent.mt")
print("Calculating Rb %")
RbGeneList <- RbGeneList[which(RbGeneList$ENSG %in% rownames(seurat_norm)),]
seurat_norm <- PercentageFeatureSet(seurat_norm, features = RbGeneList$ENSG, col.name = "percent.rb")


##### Add in Demultiplexing and Doublet Detecting Results #####
demultiplexing <- lapply(pools, function(x){
    read_delim(paste0(dir,"output/CombinedResults/", x, "/CombinedDropletAssignments.tsv"), delim = "\t")
})
names(demultiplexing) <- pools

demultiplexing <- lapply(names(demultiplexing), function(x){
    demultiplexing[[x]] <- as.data.frame(demultiplexing[[x]])
    rownames(demultiplexing[[x]]) <- paste0(x, "_",gsub("-1","",demultiplexing[[x]]$Barcode))
    return(demultiplexing[[x]])
})



### Replace NAs in demuxlet_DropletType and scSplit_Assignment and scSplit_DropletType with unassigned 
demultiplexing <- lapply(demultiplexing, function(x){
    temp <- x
    temp[c("demuxlet_DropletType", "scSplit_Assignment", "scSplit_DropletType")][is.na(temp[c("demuxlet_DropletType", "scSplit_Assignment", "scSplit_DropletType")])] <- "unassigned"
    temp$souporcell_Assignment <- ifelse(temp$souporcell_DropletType == "unassigned", "unassigned", ifelse(temp$souporcell_DropletType == "doublet", "doublet", temp$souporcell_Assignment))
    temp$demuxlet_Assignment <- ifelse(temp$demuxlet_DropletType == "doublet", "doublet", ifelse(temp$demuxlet_DropletType == "unassigned", "unassigned", temp$demuxlet_Assignment))
    return(temp)
})


##### Rename each of the clusters and then combine into one column for identifying common assignment #####
demultiplexing <- lapply(demultiplexing, function(x){
    colnames(x) <- c("Barcode",colnames(x)[2:ncol(x)])
    x$scSplit_Assignment <- paste0("scSplit_", x$scSplit_Assignment)
    x$souporcell_Assignment <- paste0("souporcell_",x$souporcell_Assignment)
    x$freemuxlet_Assignment <- paste0("freemuxlet_", x$freemuxlet_Assignment)
    x$vireo_Assignment <- paste0("vireo_", x$vireo_Assignment)
    x$demuxlet_Assignment <- paste0("demuxlet_",x$demuxlet_Assignment)
    x$combined_assignments <- paste(x$scSplit_Assignment,x$souporcell_Assignment,x$freemuxlet_Assignment,x$vireo_Assignment, x$demuxlet_Assignment, sep = "-")
    return(x)
})
names(demultiplexing) <- pools

### Make a table of the most common combinations ###
joined_assignment_counts <- lapply(demultiplexing, function(x){
    df <- as.data.frame(t(table(x$combined_assignments)))
    df <- df[order(df$Freq, decreasing = TRUE),]
    return(df)
})


### Take just the top three to assign cells to those that were shared as singlets across all to reassign to a common assignment ###
joined_assignment_counts_top <- lapply(joined_assignment_counts, function(x){
    x$Var1 <- NULL
    x <- x[1:3,]
    return(x)
})

joined_assignment_counts_top <- lapply(joined_assignment_counts_top, function(x){
    df <- separate(x,col = Var2, into = c("scSplit","souporcell","freemuxlet","vireo","demuxlet"), sep = "-")
    df$common_assignment <- gsub("demuxlet_","",df$demuxlet)
    return(df)
})

joined_assignment_key <- lapply(joined_assignment_counts_top, function(x){
    df <- pivot_longer(x,cols = c("scSplit","souporcell","freemuxlet","vireo","demuxlet"), names_to = "software")
    df$Freq <- NULL
    return(df)
})
names(joined_assignment_key) <- pools

### Pivot dataframe longer to bind with key ###
demultiplexing_long <- lapply(demultiplexing, function(x){
    temp <- pivot_longer(x, cols = c("freemuxlet_Assignment", "scSplit_Assignment", "souporcell_Assignment", "vireo_Assignment", "demuxlet_Assignment"), names_to = "Software", values_to = "Original_Assignment")
    temp$Software <- gsub("_Assignment","", temp$Software)
    return(temp)
})
names(demultiplexing_long) <- pools

demultiplexing_long <- lapply(names(demultiplexing_long), function(x){
    df <- left_join(demultiplexing_long[[x]], joined_assignment_key[[x]], by = c("Original_Assignment" = "value", "Software" = "software"))
    df$common_assignment <- ifelse(gsub("[a-zA-Z]+_","",df$Original_Assignment) == "unassigned", "unassigned", ifelse(gsub("[a-zA-Z]+_","",df$Original_Assignment) == "doublet", "doublet", df$common_assignment))
    df$Pool <- x
    return(df)
})
names(demultiplexing_long) <- pools




##### demultiplexing_long for droplet type #####
demultiplexing_droptype_long <- lapply(names(demultiplexing), function(x){
    temp <- pivot_longer(demultiplexing[[x]], cols = colnames(demultiplexing[[x]])[grep("DropletType", colnames(demultiplexing[[x]]))], names_to = "Software", values_to = "DropletType")
    temp$Software <- gsub("_DropletType","", temp$Software)
    temp$DropletType <- ifelse(temp$DropletType == "unassigned", "unassigned", ifelse(temp$DropletType == "doublet", "doublet", ifelse(temp$DropletType == "FSA0006", "singlet", ifelse(temp$DropletType == "MBE1006", "singlet", ifelse(temp$DropletType == "TOB0421", "singlet", temp$DropletType)))))
    temp$Pool <- x
    return(temp)
})
names(demultiplexing_droptype_long) <- pools




##### Will use majority to call cell type and individual assignment
### 1. Pivot wider the dataframes
droplettype_wide <- pivot_wider(joined_doublet_df4facet[,c("Barcode", "Software", "DropletType", "Pool")], names_from = "Software", values_from = "DropletType", names_prefix = "DropletType_")
assignment_wide <- pivot_wider(joined_df4facet[,c("Barcode", "Software", "common_assignment", "Pool")], names_from = "Software", values_from = "common_assignment", names_prefix = "Assignment_")

### 2. Add up the number of singlets, doublets, FSA0006, MBE10006, TOB0421 (place sums in their own columns)
for (assignment in c(unique(joined_doublet_df4facet$DropletType))){
    droplettype_wide[,c(assignment)] <- rowSums(droplettype_wide == assignment)
}

for (assignment in c(unique(joined_df4facet$common_assignment))){
    assignment_wide[,c(assignment)] <- rowSums(assignment_wide == assignment)
}

assignment_wide$doublet <- NULL
assignment_wide$unassigned <- NULL

### 3. Combine DropletType and Assignment Results into a single dataframe with barcodes using Barcode, Software and Pool as the join_by argument
droplettype_assignments <- left_join(droplettype_wide, assignment_wide, by = c("Barcode", "Pool"))


### 4. Use if else statements to call singlets and individuals assigned to those cells
droplettype_assignments$Final_Assignment <- ifelse(droplettype_assignments$doublet >= 4, "doublet", 
                                                ifelse(droplettype_assignments$unassigned >= 4, "unassigned", 
                                                ifelse(droplettype_assignments$FSA0006 >= 3, "FSA0006",
                                                ifelse(droplettype_assignments$TOB0421 >= 3, "TOB0421",
                                                ifelse(droplettype_assignments$MBE1006 >= 3, "MBE1006", "unassigned")))))


droplettype_assignments <- as.data.frame(droplettype_assignments)
rownames(droplettype_assignments) <- gsub("-1","",paste0(droplettype_assignments$Pool, "_", droplettype_assignments$Barcode))


### Add final assignments to seurat object ###
seurat_norm <- AddMetaData(seurat_norm, droplettype_assignments)



### Make pre-QC figures ###
seurat_norm <- NormalizeData(seurat_norm, verbose = TRUE)
seurat_norm <- FindVariableFeatures(seurat_norm, selection.method = "mean.var.plot")
seurat_norm <- ScaleData(seurat_norm, features = VariableFeatures(seurat_norm))
seurat_norm <- RunPCA(seurat_norm, features = VariableFeatures(object = seurat_norm))
seurat_norm <- FindNeighbors(seurat_norm, dims = 1:10)
seurat_norm <- FindClusters(seurat_norm, resolution = 0.5)
seurat_norm <- RunUMAP(seurat_norm, dims = 1:10)


### Remove doublets selected ###
Idents(seurat_norm) <- "Final_Assignment"
seurat_norm <- subset(seurat_norm,  idents = c("FSA0006", "TOB0421", "MBE1006"))

seurat_norm@meta.data$Site_rep <- seurat_norm@meta.data$MULTI_ID
Idents(seurat_norm) <- "Site_rep"
seurat_norm <- subset(seurat_norm,  idents = c("Sydney1","Brisbane1", "Melbourne1", "Brisbane2", "Sydney2", "Melbourne2", "Melbourne3", "Sydney3", "Brisbane3"))


##### Remove the outliers #####
### Only filter based on mt% and low number of genes since others are pretty well within expected distributions
seurat_filt <- subset(seurat, subset = percent.mt_mad == "NotOutlier" & nFeature_RNA > 1750) 


seurat_filt <- NormalizeData(seurat_filt, verbose = TRUE)
seurat_filt <- FindVariableFeatures(seurat_filt, selection.method = "mean.var.plot")
seurat_filt <- ScaleData(seurat_filt, features = VariableFeatures(seurat_filt))
seurat_filt <- RunPCA(seurat_filt, features = VariableFeatures(object = seurat_filt))
seurat_filt <- FindNeighbors(seurat_filt, dims = 1:10)
seurat_filt <- FindClusters(seurat_filt, resolution = 0.5)
seurat_filt <- RunUMAP(seurat_filt, dims = 1:10)


##### Compute Cell Cycle Using scran package and compare the results from seurat_filt and scran #####
print("Loading in data for cell cycle determination")
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))

### Run the cell cycle identification ###
assigned <- cyclone(seurat_filt[["RNA"]]@counts, pairs=hs.pairs)  #Note, this takes hours

table(assigned$phases)
write.table(assigned, file = paste0(outdir,"CellCycleProportions.txt"), quote = F, sep = "\t") #Save so that can read in and don't have to wait to recompute again

rownames(CellCycle) <- colnames(seurat_filt)
seurat_filt <- AddMetaData(seurat_filt, CellCycle)


##### Add gene IDs for easy identification downstream #####
GeneConversion <- read_delim("iPSC_Village_Publication/scripts/GeneConversions.tsv", delim = "\t") # available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication

meta <- data.frame("Gene_ID" = GeneConversion$Gene_ID)
rownames(meta) <- GeneConversion$ENSG


seurat_filt[["RNA"]] <- AddMetaData(seurat_filt[["RNA"]], meta)
seurat_filt@meta.data$Location_Time <- gsub(" ", "_", paste0(gsub("\\d","",seurat_filt@meta.data$MULTI_ID), "_", seurat_filt@meta.data$Time))


saveRDS(seurat_filt, paste0(outdir, "village_3_hiPSC_lines_high_quality_cells.rds"))