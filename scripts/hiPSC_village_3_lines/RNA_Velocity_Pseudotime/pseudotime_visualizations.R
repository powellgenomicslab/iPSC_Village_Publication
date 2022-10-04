library(tidyverse)
library(Seurat)
library(data.table)
library(Nebulosa)
library(viridis)
library(gridExtra)
library(ComplexHeatmap)
library('circlize')


##### Make function for saving figures in different formats #####
save_figs <- function(plot, basename, width = 17, height = 17, units = "cm"){
    ggsave(plot, filename = paste0(basename,".png"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".pdf"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".eps"), height = height, width = width, units = units)
}


##### Setting up Directories
datadir <- "/path/to/output/Integrated/"
outdir <- "/path/to/output/RNA_Velocity_Pseudotime/Visualization/"

dir.create(outdir, recursive = TRUE)

##### Set up Colors #####
cell_line_colors <- c("FSA0006" = "#F79E29", "MBE1006" = "#9B2C99", "TOB0421"= "#35369C")
replicate_colors <- c("1" = "#ACD39E", "2" = "#5AAE61", "3" = "#1B7837")
time_colors <- c("Baseline" = "#b9cee4", "Village Day 4" = "#8a92bb", "Thawed Village Day 0" = "#7d57a0", "Thawed Village Day 7" = "#853786")
village_colors <- c("Uni-Culture" = "#613246", "Village" = "#A286AA")
cycle_colors <- c("G1" = "#4393C3", "S" = "#92C5DE", "G2M" = "#D1E5F0")
site_colors <- c("Brisbane" = "#536CB4", "Sydney" = "#62BD67", "Melbourne" = "#D24F72")
cluster_colors <- c("0" = "#C9D8EA", "1" = "#928CC5", "2" = "#A7C9A9", "3" = "#179085", "4" = "#F79F9F", "5" = "#C35B76", "6" = "#F4C893", "7" = "#F6AA4B")
cryo_colors <- c("Fresh" = "#e7bbc6", "Cryopreserved" = "#d08790")


##### Read in data #####
seurat = readRDS(paste0(datadir,"village_3_hiPSC_lines_high_quality_cells_integrated.rds"))
meta = fread(paste0(outdir,"metadata.csv"), sep = ",")
gene_meta <- fread(paste0(outdir, "var_meta.csv"), sep = ",")


velo_meta <- meta[,c("n_unspliced_counts", "latent_time")]
rownames(velo_meta) <- meta$V1


##### Add metadata time to seurat object
seurat <- AddMetaData(seurat, velo_meta)
DefaultAssay(seurat) <- "SCT"


##### Make plots of velocity genes
#### Conversion dataframe

##### Add gene IDs for easy identification downstream #####
GeneConversion <- fread("iPSC_Village_Publication/scripts/GeneConversions.tsv", sep = "\t")
colnames(GeneConversion) <- c("ENSG_ID", "Gene_ID")

gene_df <- GeneConversion[gene_df, on = c("Gene_ID" = "Gene")]




##### Make latent time and expression plot #####
### Remove cells without latent time
seurat_noNA <- subset(seurat, subset = latent_time >= 0)


##### Make UMAPS #####
## With Village ##
UMAP_village <- DimPlot(seurat, reduction = "umap", group.by = c("Village"), cols = village_colors) + labs(color="Village") + ggtitle(NULL) 
save_figs(UMAP_village, basename = paste0(outdir,"Village_umap"),width = 15, height = 15)

## By Cell Line ##
UMAP_cell_line <- DimPlot(seurat, reduction = "umap", group.by = c("Final_Assignment"), cols = cell_line_colors) + labs(color="Cell Line") + ggtitle(NULL) 
save_figs(UMAP_cell_line, basename = paste0(outdir,"Cell_Line_umap"),width = 15, height = 15)

## By Location ##
UMAP_site <- DimPlot(seurat, reduction = "umap", group.by = c("Location"), cols = alpha(site_colors, 0.7), shuffle = TRUE) + labs(color="Location") + ggtitle(NULL)
UMAP_site[[1]]$layers[[1]]$aes_params$shape = 16
UMAP_site[[1]]$layers[[1]]$aes_params$size = 0.3
save_figs(UMAP_site, basename = paste0(outdir,"Location_umap"),width = 15, height = 15)


## By Cryopreserved ##
UMAP_cryo <- DimPlot(seurat, reduction = "umap", group.by = c("Cryopreserved"), order = c("Cryopreserved", "Fresh"), cols = alpha(cryo_colors, 0.7)) + labs(color="Cryopreserved") + ggtitle(NULL)
UMAP_cryo[[1]]$layers[[1]]$aes_params$shape = 16
UMAP_cryo[[1]]$layers[[1]]$aes_params$size = 0.3
save_figs(UMAP_cryo, basename = paste0(outdir,"Cryopreserved_umap"),width = 15, height = 15)


## By Cell Cycle Phase ##
Idents(seurat) <- "phases"
UMAP_cell_cycle <- DimPlot(seurat, reduction = "umap", group.by = c("phases"), cols = cycle_colors, order = c("G2M", "S", "G1")) + labs(color="Cell Cycle\nPhase") + ggtitle(NULL)
UMAP_cell_cycle[[1]]$layers[[1]]$aes_params$shape = 16
UMAP_cell_cycle[[1]]$layers[[1]]$aes_params$size = 0.3
save_figs(UMAP_cell_cycle, basename = paste0(outdir,"CellCycle_umap"),width = 15, height = 15)


### Plot Pseudotime ###
pLatent <- FeaturePlot(seurat_noNA, feature = "latent_time", pt.size = 0.1, order = TRUE) +
			theme(text = element_text(size=20)) +
			scale_color_viridis(name = "Pseudo-\ntime")

legend_latent <- cowplot::get_legend(pLatent)


pLatent <- FeaturePlot(seurat_noNA, feature = "latent_time", pt.size = 0.1) + 
			scale_color_viridis(alpha = 0.25) + 
			ggtitle("RNA Velocity Pseudotime") + 
			theme(text = element_text(size=20)) +
			xlab("UMAP 1") +
			ylab("UMAP 2") +
			theme(legend.position = "none")

platent_combined  <- grid.arrange(pLatent, legend_latent, ncol = 2, widths = 2:0.5)
ggsave(platent_combined, filename = paste0(outdir, "umap_latent.png"), height = 5, width = 6.5)
ggsave(platent_combined, filename = paste0(outdir, "umap_latent.pdf"), height = 5, width = 6.5)


##### Plot POU5F1 on umap #####
pPOU5F1 <- FeaturePlot(seurat_noNA, feature = GeneConversion$ENSG_ID[which(GeneConversion$Gene_ID == "POU5F1")], pt.size = 0.1) + 
			theme(text = element_text(size=20)) +
			scale_color_viridis(option = "inferno", name = "Expression")

legend_POU5F1 <- cowplot::get_legend(pPOU5F1)

pPOU5F1 <- FeaturePlot(seurat_noNA, feature = GeneConversion$ENSG_ID[which(GeneConversion$Gene_ID == "POU5F1")], pt.size = 0.1) +
			scale_color_viridis(option = "inferno", alpha = 0.25) + 
			ggtitle("POU5F1 Expression") + 
			theme(text = element_text(size=20)) +
			xlab("UMAP 1") +
			ylab("UMAP 2") +
			theme(legend.position = "none")

pPOU5F1_combined  <- grid.arrange(pPOU5F1, legend_POU5F1, ncol = 2, widths = 2:0.5)
ggsave(pPOU5F1_combined, filename = paste0(outdir, "umap_POU5F1.png"), height = 5, width = 6.5, res = 600)
ggsave(pPOU5F1_combined, filename = paste0(outdir, "umap_POU5F1.pdf"), height = 5, width = 6.5)


##### Plot LIX1 on umap #####
pLIX1 <- FeaturePlot(seurat_noNA, feature = GeneConversion$ENSG_ID[which(GeneConversion$Gene_ID == "LIX1")], pt.size = 0.1) + 
			theme(text = element_text(size=20)) +
			scale_color_viridis(option = "inferno", name = "Expression")

legend_LIX1 <- cowplot::get_legend(pLIX1)

pLIX1 <- FeaturePlot(seurat_noNA, feature = GeneConversion$ENSG_ID[which(GeneConversion$Gene_ID == "LIX1")], pt.size = 0.1) +
			scale_color_viridis(option = "inferno", alpha = 0.25) + 
			ggtitle("LIX1 Expression") + 
			theme(text = element_text(size=20)) +
			xlab("UMAP 1") +
			ylab("UMAP 2") +
			theme(legend.position = "none")

pLIX1_combined  <- grid.arrange(pLIX1, legend_LIX1, ncol = 2, widths = 2:0.5)
ggsave(pLIX1_combined, filename = paste0(outdir, "umap_LIX1.png"), height = 5, width = 6.5, res = 600)
ggsave(pLIX1_combined, filename = paste0(outdir, "umap_LIX1.pdf"), height = 5, width = 6.5)


##### Plot PTN on umap #####
pPTN <- FeaturePlot(seurat_noNA, feature = GeneConversion$ENSG_ID[which(GeneConversion$Gene_ID == "PTN")], pt.size = 0.1) + 
			theme(text = element_text(size=20)) +
			scale_color_viridis(option = "inferno", name = "Expression")

legend_PTN <- cowplot::get_legend(pPTN)

pPTN <- FeaturePlot(seurat_noNA, feature = GeneConversion$ENSG_ID[which(GeneConversion$Gene_ID == "PTN")], pt.size = 0.1) +
			scale_color_viridis(option = "inferno", alpha = 0.25) + 
			ggtitle("PTN Expression") + 
			theme(text = element_text(size=20)) +
			xlab("UMAP 1") +
			ylab("UMAP 2") +
			theme(legend.position = "none")

pPTN_combined  <- grid.arrange(pPTN, legend_PTN, ncol = 2, widths = 2:0.5)
ggsave(pPTN_combined, filename = paste0(outdir, "umap_PTN.png"), height = 5, width = 6.5, res = 600)
ggsave(pPTN_combined, filename = paste0(outdir, "umap_PTN.pdf"), height = 5, width = 6.5)



pLocation_Time_latent_bar <- ggplot(seurat_noNA@meta.data, aes(latent_time, fill = Village)) +
	geom_histogram(alpha = 0.75, position="identity") +
	theme_classic() +
	facet_grid(Location ~ Final_Assignment,scales = "free_y") +
	scale_fill_manual(values = village_colors) +
	scale_x_continuous(breaks=c(0,0.5,1))

ggsave(pLocation_Time_latent_bar, filename = paste0(outdir,"faceted_histogram_latent_bar.png"), width = 4.5, height = 5)
ggsave(pLocation_Time_latent_bar, filename = paste0(outdir,"faceted_histogram_latent_bar.pdf"), width = 4.5, height = 5)

