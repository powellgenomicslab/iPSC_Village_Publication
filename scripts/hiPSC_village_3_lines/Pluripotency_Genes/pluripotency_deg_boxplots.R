##############################################################################
# Author: Drew Neavin
# Date: 2022-10-04
# Description: Pluripotency gene differential expression and violin plots.
##############################################################################

library(data.table)
library(tidyverse)
library(Seurat)
library(ggnewscale)
library(ggpubr)
library(facefuns)
library(facetscales)

##### Set up Directories #####
dir <- '/path/to/output/'
outdir <- paste0(dir,"/Pluripotency_Genes/")
dir.create(outdir, recursive = TRUE)


##### Make funciton to make dataframees for each gene #####
expression_df <- function(seurat_list, gene, meta_data){
	Expression_list <- lapply(seurat_list, function(x){
		data.frame(Counts = x[["SCT"]]@counts[gene,], Normalized = x[["SCT"]]@scale.data[gene,], x@meta.data[,c(meta_data)])
	})
	do.call(rbind, Expression_list)
}


save_figs <- function(plot, basename, width = 17, height = 17, units = "cm"){
    ggsave(plot, filename = paste0(basename,".png"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".pdf"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".eps"), height = height, width = width, units = units)
}



###### Read in data ######
seurat <- readRDS(paste0(dir, "Integrated/village_3_hiPSC_lines_high_quality_cells_integrated.rds"))
pluri_genes <- fread("iPSC_Village_Publication/scripts/hiPSC_village_3_lines/Pluripotency_Genes/pluripotency_genes.tsv") # available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication




##### Create a single column with locations x cryopreservation x line for subsetting ######
seurat@meta.data$Cryopreserved <- ifelse(grepl("Thawed", seurat@meta.data$Location_Time), "Cryopreserved", "Fresh")
seurat@meta.data$Location_Cryopreserved_Line <- paste0(seurat@meta.data$Location, "_", seurat@meta.data$Cryopreservation, "_", seurat@meta.data$Final_Assignment)



##### Subset different locations x cryopreservation x line conditions #####
seurat_list <- lapply(unique(seurat@meta.data$Location_Cryopreserved_Line), function(x){
	subset(seurat, subset = Location_Cryopreserved_Line == x)
})
names(seurat_list) <- unique(seurat@meta.data$Location_Cryopreserved_Line)



seurat_list <- lapply(seurat_list, function(x){
    Idents(x) <- "Village"
    return(x)
})



seurat_list <- lapply(seurat_list, function(x){
    SCTransform(x, vars.to.regress = c("scores.G1", "scores.S", "scores.G2M", "percent.mt", "percent.rb"), return.only.var.genes = FALSE)
})

seurat_list <- lapply(seurat_list, function(x){
    PrepSCTFindMarkers(x)
})




### Try with logistic regression  downsampled to lowest number across all pools###
### Since the smallest condition has 572 cells, downsample to that number for all conditions for equivalent comparison ###
LR_DEGs_sub <- lapply(seurat_list, function(x){
    FindMarkers(x, ident.1 = "Uni-Culture", ident.2 = "Village", latent.vars = "MULTI_classification", test.use = "LR", logfc.threshold = 0, max.cells.per.ident = 572, assay = "RNA")
})



### Update multiple testing for all sites ###
LR_DEGs_sub_nrow_list <- lapply(LR_DEGs_sub, function(x) data.table(nrow(x)))
LR_DEGs_sub_nrow <- sum(do.call(rbind, LR_DEGs_sub_nrow_list)$V1)

LR_DEGs_sub <- lapply(LR_DEGs_sub, function(x){
    x$p_val_adj_updated <- p.adjust(x$p_val, method = "bonferroni", n = LR_DEGs_sub_nrow)
    return(x)
})


LR_DEGs_sub_pluri <- lapply(LR_DEGs_sub, function(x){
    x$ENSG <- rownames(x)
    x <- data.table(x)
    x <- x[pluri_genes, on = c("ENSG")]
    return(x)
})


LR_DEGs_sub_pluri_sig <- lapply(LR_DEGs_sub_pluri, function(x){
    x <- x[!is.na(p_val) & p_val_adj_updated < 0.05]
    return(x)
})



LR_DEGs_sub_pluri_sig_4 <- lapply(LR_DEGs_sub_pluri_sig, function(x){
    x[GeneID %in% c("MYC", "NANOG", "POU5F1", "SOX2")]
})




##### Make volcano of significants #####
LR_DEGs_sub_pluri <- lapply(names(LR_DEGs_sub_pluri), function(x){
    LR_DEGs_sub_pluri[[x]]$Group <- x
    return(LR_DEGs_sub_pluri[[x]])
})


LR_DEGs_sub_pluri_dt <- do.call(rbind, LR_DEGs_sub_pluri)
LR_DEGs_sub_pluri_dt$p_val_adj_updated <- ifelse(is.na(LR_DEGs_sub_pluri_dt$p_val_adj_updated), 1, LR_DEGs_sub_pluri_dt$p_val_adj_updated)
LR_DEGs_sub_pluri_dt$significant <- ifelse(LR_DEGs_sub_pluri_dt$p_val_adj_updated > 0.05, "not significant", ifelse(abs(LR_DEGs_sub_pluri_dt$avg_log2FC) > 1,"significant large effect size", "significant small effect size"))
LR_DEGs_sub_pluri_dt$significant  <- factor(LR_DEGs_sub_pluri_dt$significant, levels = c("significant large effect size", "significant small effect size", "not significant"))



pVolcano <- ggplot(LR_DEGs_sub_pluri_dt[!grepl("Cryopreserved", Group)], aes(avg_log2FC, -log2(p_val_adj_updated), color = significant)) +
    geom_point() +
    theme_classic() +
    geom_vline(xintercept = -1, linetype = "dashed") +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    scale_color_manual(values = c("black", "grey50", "grey80")) +
    ylab("-log2(P Value)") +
    xlab("log2(Fold Change)") +
    ggtitle("Pluripotent Gene Differential Expression\nBetween Uni-Culture and Village") +
    theme(plot.title = element_text(hjust = 0.5))


ggsave(pVolcano, filename = paste0(outdir, "volcano.png"), width = 4.5, height = 3)



pVolcano_cryo <- ggplot(LR_DEGs_sub_pluri_dt[grepl("Sydney", Group)], aes(avg_log2FC, -log2(p_val_adj_updated), color = significant)) +
    geom_point() +
    theme_classic() +
    geom_vline(xintercept = -1, linetype = "dashed") +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    scale_color_manual(values = c("black", "grey50", "grey80")) +
    ylab("-log2(P Value)") +
    xlab("log2(Fold Change)") +
    ggtitle("Pluripotent Gene Differential Expression\nBetween Uni-Culture and Village") +
    theme(plot.title = element_text(hjust = 0.5))


ggsave(pVolcano_cryo, filename = paste0(outdir, "volcano_cryo.png"), width = 4.5, height = 3)




##### Pluri Genes #####
pluri_genes <- data.frame(Gene = c("MYC", "NANOG", "POU5F1", "SOX2"), ENSG = c("ENSG00000136997", "ENSG00000111704", "ENSG00000204531", "ENSG00000181449"))

df_list <- list()
for (gene in pluri_genes$Gene){
	df_list[[gene]] <- expression_df(seurat_list, pluri_genes[which(pluri_genes$Gene == gene),"ENSG"], c("Location", "Time", "Final_Assignment"))
	df_list[[gene]]$Gene <- gene
}
df <- data.table(do.call(rbind, df_list))



scales_y <- list(
  `MYC` = scale_y_continuous(limits = c(-3, 9)),
  `NANOG` = scale_y_continuous(limits = c(-3, 9)),
  `POU5F1` = scale_y_continuous(limits = c(-5, 8)),
  `SOX2` = scale_y_continuous(limits = c(-4, 8))
)

p_scaled_vio <- ggplot() +
					geom_density(data=subset(df[Location != "Cryopreserved"], Time=="Uni-Culture"), aes(y = Normalized, fill=Time, x= -..density.., group = Final_Assignment), alpha = 0.7, lwd = 0) +
					geom_density(data=subset(df[Location != "Cryopreserved"], Time=="Village"), aes(y = Normalized, fill=Time, x= ..density.., group = Final_Assignment), alpha = 0.7, lwd = 0) +
					theme_classic() +
					facet_grid_sc(Gene ~ paste0(Location, Final_Assignment), scales = list(y = scales_y, x = "free")) +
					scale_color_manual(values = village_colors) +
					scale_fill_manual(values = alpha(village_colors, 0.4)) +
					theme(axis.title.x = element_blank()) +
					ylab("Normalized Expression") +
					theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
					geom_text(
						data = pluri_deg_dt[Cryopreservation != "Cryopreserved"],
						aes(x = 0, y = scaled_position,label = Symbol), 
						color = "black",
						size = 3)  +
					theme(axis.title.x=element_blank(),
						axis.text.x=element_blank(),
						axis.ticks.x=element_blank(),
						panel.spacing.x = unit(0, "lines"))

save_figs(p_scaled_vio, paste0(outdir,"pluri_genes/pluri_normalized_vio"), width = 14.25, height = 8.75)




df_cryo <- df[grepl("Site 3", df$Location),]
df_cryo$Final_Assignment <- factor(df_cryo$Final_Assignment, levels = c("FSA0006", "MBE1006", "TOB0421"))

scales_y_cryo <- list(
  `MYC` = scale_y_continuous(limits = c(-3, 10)),
  `NANOG` = scale_y_continuous(limits = c(-3, 10)),
  `POU5F1` = scale_y_continuous(limits = c(-6, 10)),
  `SOX2` = scale_y_continuous(limits = c(-4, 9))
)


p_scaled_cryo_vio <- ggplot() +
					geom_density(data=subset(df_cryo, Time=="Uni-Culture"), aes(y = Normalized, fill=Time, x= -..density.., group = Final_Assignment), alpha = 0.7, lwd = 0) +
					geom_density(data=subset(df_cryo, Time=="Village"), aes(y = Normalized, fill=Time, x= ..density.., group = Final_Assignment), alpha = 0.7, lwd = 0) +
					theme_classic() +
					facet_grid_sc(Gene ~ factor(paste0(Cryopreservation, Final_Assignment), levels = c(paste0("Fresh", levels(df_cryo$Final_Assignment)), paste0("Cryopreserved", levels(df_cryo$Final_Assignment)))), scales = list(y = scales_y_cryo, x = "free")) +
					scale_color_manual(values = village_colors) +
					scale_fill_manual(values = alpha(village_colors, 0.4)) +
					theme(axis.title.x = element_blank()) +
					ylab("Normalized Expression") +
					theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
					geom_text(
						data = pluri_deg_dt_cryo,
						aes(x = 0, y = counts_position,label = Symbol), 
						color = "black",
						size = 3) +
					theme(axis.title.x=element_blank(),
						axis.text.x=element_blank(),
						axis.ticks.x=element_blank(),
						panel.spacing.x = unit(0, "lines"))

save_figs(p_scaled_cryo_vio, paste0(outdir,"pluri_genes/pluri_normalized_cryo_vio"), width = 13.5, height = 8.5)
