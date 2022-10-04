##### Author: Drew Neavin
##### Date: 4 October, 2022
##### Description: Look at the proportions of each line at each time of Nona's multi-ome experiment



#### Setting up directories #####
print("Set up directories")
dir <- '/path/to/output/hiPSC_village_18_lines/iPSC_maintenance'
outdir <- paste0(outdir, "QC/")



##### Load in libraries #####
library(data.table)
library(Seurat)
library(tidyverse)
library("ggpomological")


##### Set up directories #####
seurat <- readRDS(, paste0(outdir,"seurat_qc_integrated.rds"))

summary <- data.table(prop.table(table(seurat@meta.data[,c("AtLeastHalfSinglet_Individual_Assignment", "Pool")]), margin = 2))


summary$AtLeastHalfSinglet_Individual_Assignment <- factor(summary$AtLeastHalfSinglet_Individual_Assignment, levels = rev(summary[Pool == "Village_P8"]$AtLeastHalfSinglet_Individual_Assignment[order(summary[Pool == "Village_P8"]$N)]))


##### Make proportion plots (area plot) #####
p_stacked_area_filt <- ggplot(village_summary_singlets, aes(x = Day, y = N, fill = factor(AtLeastHalfSinglet_Individual_Assignment), group = AtLeastHalfSinglet_Individual_Assignment)) +
    geom_area(alpha=0.6 , size=0.1, colour="black") +
    theme_classic() +
    xlab("Passage") +
    ylab("Proportion of Cells")

ggsave(p_stacked_area_filt, filename = paste0(outdir,"stacked_area_filtered.png"), width = 4, height = 2)
ggsave(p_stacked_area_filt, filename = paste0(outdir,"stacked_area_filtered.pdf"), width = 4, height = 2)



