##### Reason: combine the results for the variance explained by different factors for each gene
##### Author: Drew Neavin
##### Date: 14 March, 2022


##### Load in libraries #####
library(data.table)
library(tidyverse)
library(ggridges)
library(raincloudplots)
library(ggdist)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)
library(RColorBrewer)



##### Set up directories #####
icc_dir <- "/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/Variance/gene_separated/icc/"
icc_interaction_dir <- "/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/Variance/gene_separated/icc_interaction/"
outdir <- "/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/Variance/combined/"

dir.create(outdir, recursive = TRUE)


vars <- c("Line", "Passage", "Line:Passage", "Residual")
selected_vars <- c("Line", "Passage", "Line:Passage", "Residual")
var_colors <- c("#4734a9", "#78c0fe", "#97cf8a", "gray80")

names(var_colors) <- vars

var_colors <- var_colors[selected_vars]



##### Get list of icc files #####
icc_files <- list.files(icc_dir)



##### Read in icc results #####
icc_results_list <- lapply(icc_files, function(x){
    readRDS(paste0(icc_dir,x))
})
names(icc_results_list) <- icc_files


##### Get list of icc interaction files #####
icc_interaction_files <- list.files(icc_interaction_dir, pattern = "_icc.rds")



##### Read in icc results #####
icc_interaction_results_list <- lapply(icc_interaction_files, function(x){
    readRDS(paste0(icc_interaction_dir,x))
})
names(icc_interaction_results_list) <- icc_interaction_files



##### Merge icc results into a single data.table #####
icc_dt <- do.call(rbind, icc_results_list)



##### Update columns for easy plotting and interpretation #####
icc_dt$percent_round <- round(icc_dt$percent)


### Make grp a factor so ordered as desired for plotting ###
icc_dt$grp <- factor(icc_dt$grp, levels= rev(c("Line", "Passage", "Line:Passage", "Residual")))


### Add size of groups ###
group_size  <- data.table(table(icc_dt$grp))
colnames(group_size) <- c("grp", "size")
group_size$grp_size <- paste0(group_size$grp, "\nN = ", group_size$size)

icc_dt <- group_size[icc_dt, on = "grp"]
icc_dt$grp_size <- factor(icc_dt$grp_size, levels = unique(group_size$grp_size))


##### Merge icc_interaction results into a single data.table #####
icc_interaction_dt <- do.call(rbind, icc_interaction_results_list)

icc_interaction_dt$percent_round <- round(icc_interaction_dt$percent)

### Make grp a factor so ordered as desired for plotting ###
icc_interaction_dt$grp <- factor(icc_interaction_dt$grp, levels= rev(c("Line", "Passage", "Line:Passage", "Residual")))


### Add size of groups ###
group_size  <- data.table(table(icc_interaction_dt$grp))
colnames(group_size) <- c("grp", "size")
group_size$grp_size <- paste0(group_size$grp, "\nN = ", group_size$size)

icc_interaction_dt <- group_size[icc_interaction_dt, on = "grp"]
icc_interaction_dt$grp_size <- factor(icc_interaction_dt$grp_size, levels = unique(group_size$grp_size))

### *** Need to add individual effects without interaction in to interaction dt *** ###
icc_interaction_plus_dt <- rbind(icc_interaction_dt, icc_dt[!(gene %in% icc_interaction_dt$gene)])


group_size  <- data.table(table(icc_interaction_plus_dt$grp))
colnames(group_size) <- c("grp", "size")
group_size$grp_size <- paste0(group_size$grp, "\nN = ", group_size$size)

icc_interaction_plus_dt <- group_size[icc_interaction_plus_dt, on = "grp"]
icc_interaction_plus_dt$grp_size <- factor(icc_interaction_plus_dt$grp_size, levels = unique(group_size$grp_size))
icc_interaction_plus_dt$grp_size <- factor(icc_interaction_plus_dt$grp_size, levels = c("Line\nN = 750", "Passage\nN = 750", "Line:Passage\nN = 619", "Residual\nN = 750"))



##### Generate raincloud plot of variance explained #####
pRaincloud_interaction <- ggplot(icc_interaction_plus_dt, aes(x = percent, y = factor(grp_size, levels = levels(grp_size)), fill = factor(grp, levels = rev(vars)))) + 
                geom_density_ridges(size = 0.1,stat = "binline", bins = 100, scale = 0.7, draw_baseline = FALSE, aes(height =..ndensity..)) +
                geom_point(size =1, position = position_nudge(y=-0.09), shape = "|", aes(color = factor(grp, levels = rev(vars)))) +
                coord_cartesian(xlim = c(1.2, NA), clip = "off") +
                theme_classic() +
                theme(axis.title.y=element_blank()) +
                xlab("Percent Variance Explained") +
                scale_y_discrete(expand = c(0.1, 0)) +
                scale_fill_manual(values = var_colors, name = "Variable") +
                scale_color_manual(values = var_colors, name = "Variable") +
                geom_vline(xintercept = 1, lty="11", color = "grey50", size = 0.5)

ggsave(pRaincloud_interaction, filename = paste0(outdir, "variance_explained_raincloud_interaction.png"), height = 2, width = 4)
ggsave(pRaincloud_interaction, filename = paste0(outdir, "variance_explained_raincloud_interaction.pdf"), height = 2, width = 4)





##### Add gene IDs for easy identification downstream #####
GeneConversion <- fread("iPSC_Village_Publication/scripts/GeneConversions.tsv", sep = "\t") # available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication
colnames(GeneConversion) <- c("gene", "Gene_ID")

icc_interaction_plus_dt <- GeneConversion[icc_interaction_plus_dt, on = "gene"]



#### Interrogate variance explained of Y chromosome genes (positive control) and Pluripotency genes
## Read in gtf used as reference and pull just Y chromosome genes from it ##
gtf <- fread("iPSC_Village_Publication/scripts/refdata-cellranger-GRCh38-3.0.0_gene.gtf", sep = "\t", autostart = 6, header = FALSE) # available from github: https://github.com/powellgenomicslab/iPSC_Village_Publication

gtf_genes <- gtf[!(grep("transcript_id", V9))]

gtf_genes$V9 <- gsub("gene_id \"", "",gtf_genes$V9 ) %>%
            gsub("\"; gene_version \"", ";", .) %>%
            gsub("\"; gene_name \"", ";", .) %>%
            gsub("\"; gene_source \"", ";", .) %>%
            gsub("\"; gene_biotype \"", ";", .) %>%
            gsub("\"", "", .)

gtf_genes[, c("gene_id", "gene_version", "gene_name", "gene_source", "gene_biotype") := data.table(str_split_fixed(V9,";", 5))]


Y_chromosome_genes <- gtf_genes[V1 == "Y"]
Y_genelist <- Y_chromosome_genes$gene_id[Y_chromosome_genes$gene_id %in% genes]


### Make stacked bar plots of the variance explained by different factors for these gene groups
## Figure of y chromosome genes ##
icc_y <- icc_interaction_plus_dt[data.table(gene = icc_interaction_plus_dt[grp == "Residual"][gene %in% Y_genelist][order(percent_round)]$gene), on = "gene"]
icc_y$grp <- factor(icc_y$grp, levels = rev(selected_vars))
icc_y$gene <- factor(icc_y$gene, levels = unique(icc_y$gene))
icc_y$Gene_ID <- factor(icc_y$Gene_ID, levels = unique(icc_y$Gene_ID))


bar_proportions_y <- ggplot(icc_y, aes(x = Gene_ID, y = percent, fill = grp)) +
    geom_bar(position="stack", stat="identity", alpha = 0.75) +
    theme_classic() +
    theme(axis.title.x=element_blank(),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(values = var_colors) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle("Variance Explained of\nY Chromosome Genes") +
    ylab("Percent")


ggsave(bar_proportions_y, filename = paste0(outdir, "variance_explained_bar_y_genes.png"), width = 4.5, height = 4.5)
ggsave(bar_proportions_y, filename = paste0(outdir, "variance_explained_bar_y_genes.pdf"), width = 4.5, height = 4.5)



### Plot Pluripotency Genes ###
pluri_genes <- fread("iPSC_Village_Publication/scripts/pluripotency_genes.tsv", sep = "\t") # available from github: https://github.com/powellgenomicslab/iPSC_Village_Publication


pluri_genes <- GeneConversion[pluri_genes, on = "Gene_ID"]


icc_dt_pluri_genes <- icc_interaction_plus_dt[pluri_genes,on = c("gene")]
icc_dt_pluri_genes$grp <- factor(icc_dt_pluri_genes$grp, levels = rev(selected_vars))
icc_dt_pluri_genes <- icc_dt_pluri_genes[data.table(gene = icc_dt_pluri_genes[grp == "Residual"][order(percent)]$gene), on = "gene"]
icc_dt_pluri_genes$Gene_ID <- factor(icc_dt_pluri_genes$Gene_ID, levels = unique(icc_dt_pluri_genes$Gene_ID))
 


pPluri_Genes_Cont <- ggplot() +
						geom_bar(data = icc_dt_pluri_genes, aes(Gene_ID, percent, fill = grp), position = "stack", stat = "identity", alpha = 0.75) +
						theme_classic() +
						scale_fill_manual(values = var_colors) +
						theme(plot.title = element_text(hjust = 0.5),
                            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
						ylab("Percent Gene Expression Variance Explained") +
                        ggtitle("Variance Explained of\nPluripotency Genes") +
						theme(axis.title.x=element_blank())

ggsave(pPluri_Genes_Cont, filename = paste0(outdir, "Pluripotent_Gene_Variable_Contributions.png"), width = 6, height = 4.5)
ggsave(pPluri_Genes_Cont, filename = paste0(outdir, "Pluripotent_Gene_Variable_Contributions.pdf"), width = 6, height = 4.5)


### Four common genes (MYC, NANOG, POU5F1, SOX2) ###
pPluri_Genes_Cont_top <- ggplot() +
						geom_bar(data = icc_dt_pluri_genes[GeneID %in% c("MYC", "NANOG", "POU5F1", "SOX2")], aes(Gene_ID, percent, fill = grp), position = "stack", stat = "identity", alpha = 0.75) +
						theme_classic() +
						scale_fill_manual(values = var_colors) +
						theme(plot.title = element_text(hjust = 0.5),
                            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
						ylab("Percent Gene Expression Variance Explained") +
                        ggtitle("Variance of\nPluripotency Genes Explained") +
						theme(axis.title.x=element_blank())

ggsave(pPluri_Genes_Cont_top, filename = paste0(outdir, "Pluripotent_Gene_Variable_Contributions_four.png"), width = 3, height = 4.5)
ggsave(pPluri_Genes_Cont_top, filename = paste0(outdir, "Pluripotent_Gene_Variable_Contributions_four.pdf"), width = 3, height = 4.5)
