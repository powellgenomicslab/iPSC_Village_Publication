##### Author: Drew Neavin
##### Date: 4 October, 2022
##### Reason: Look at the proportions of each line at each time of Nona's multi-ome experiment


##### Load in libraries #####
library(data.table)
library(Seurat)
library(tidyverse)
library("ggpomological")


##### Set up directories #####
metadata <- fread("iPSC_Village_Publication/data/cardiac/cardiac_droplet_annotation.tsv", sep = "\t") # available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication


summary <- data.table(prop.table(table(metadata[,c("Assignment", "Day_updated")]), margin = 2))

summary_singlets <- data.table(prop.table(table(summary[Assignment != "unassigned" & Assignment != "doublet",c("Assignment", "Day_updated")]), margin = 2))

summary_singlets$Assignment <- factor(summary_singlets$Assignment, levels = rev(summary_singlets[Day_updated == 15]$Assignment[order(summary_singlets[Day_updated == 15]$N)]))


##### Make proportion plots (area plot) #####
p_stacked_area <- ggplot(summary_singlets, aes(x = as.numeric(as.character(Pool_ID)), y = N, fill = factor(Assignment), group = Assignment)) +
    geom_area(alpha=0.6 , size=0.1, colour="black") +
    theme_classic() +
    scale_fill_manual(values = c("#f44336", "#e81f63", "#9c27b0", "#673ab7", "#3f51b5", "#2096f3","#2096f3", "#009688", "#4caf50", "#8bc34a", "#cddc39", "#ffeb3b", "#ffc108", "#ff9801", "#ff5723" ,"#795548", "#9e9e9e", "#607d8b")) +
    xlab("Days") +
    ylab("Proportion of Cells")
ggsave(p_stacked_area, filename = paste0(outdir,"stacked_area.png"), width = 4, height = 2)
ggsave(p_stacked_area, filename = paste0(outdir,"stacked_area.pdf"), width = 4, height = 2)


