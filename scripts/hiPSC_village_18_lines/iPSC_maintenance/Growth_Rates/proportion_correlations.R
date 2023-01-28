library(data.table)
library(tidyverse)
library(ggplot2)
library(Seurat)


##### Set up directories #####
indir <- "/path/to/iPSC_Village_Publication/data/village_18_hiPSC_lines/"
outdir <- "/path/to/iPSC_Village_Publication/output/hiPSC_village_18_lines/iPSC_maintenance/Growth_Rates/"
rate_dir <- "/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/Growth_Rates/"

dir.create(outdir, recursive = TRUE)



##### Read in proportions #####
diff_props <- fread(paste0(indir, "cardiac_diff_prop_lines.tsv"), sep = "\t") ### Table provided on Github/Zenodo
multi_passage_props <- fread(paste0(indir, "multi_passage_prop_lines.tsv"), sep = "\t") ### Table provided on Github/Zenodo



##### Calculate the ratio of proportions to baseline "proportion #####
diff_ratios_dt <- unique(diff_props[,c("Assignment", "Day")])
diff_ratios_dt$ratio <- as.numeric(NA)

multi_passage_ratios_dt <-  unique(multi_passage_props[,c("Assignment", "Day")])
multi_passage_ratios_dt$ratio <- as.numeric(NA)

for (sample in unique(diff_props$Assignment)){
    for (day in unique(diff_props$Day)){
        diff_ratios_dt[Assignment == sample & Day == day]$ratio <- diff_props[Assignment == sample & Day == day]$N/diff_props[Assignment == sample & Day == 0]$N
    }
    for (day in unique(multi_passage_ratios_dt$Day)){
        multi_passage_ratios_dt[Assignment == sample & Day == day]$ratio <- multi_passage_props[Assignment == sample & Day == day]$N/diff_props[Assignment == sample & Day == 0]$N
    }
}

multi_passage_ratios_dt <- rbind(diff_ratios_dt[Day ==0], multi_passage_ratios_dt)




##### Read in results #####
rate_results.2 <- lapply(unique(diff_props$Assignment), function(sample){
    fread(paste0(rate_dir, sample, "_ratrack.fit.csv"))
})
names(rate_results.2) <- unique(diff_props$Assignment)


rate_results_2.2 <- lapply(rate_results.2, function(x){
    x[model_index == 1]
})



##### Calculate the rate in the middle of the growth curve #####
rate_results_2.2 <- lapply(rate_results_2.2, function(x){
    x$name <- gsub(" ", "", x$name)
    return(x)
})

rate_results_2.2_middle <- lapply(rate_results_2.2, function(x){
    tmp <- data.table(name = unique(x$name), Assignment = gsub("_rep[1-3]","",unique(x$name)), rate = as.numeric(NA))
    for (line in tmp$name){
        tmp[name == line, c("rate")] <- (x[rate_position == 1 & name == line]$rate_mean + x[rate_position == 0 & name == line]$rate_mean)/2
    }
    return(tmp)
})


multi_passage_ratios_dt_2 <- do.call(rbind, rate_results_2.2_middle)




##### Just proportion of the village at that time #####
multi_passage_proportion_combined_dt <- multi_passage_props[multi_passage_ratios_dt_2, on = "Assignment", allow.cartesian=TRUE]
colnames(multi_passage_proportion_combined_dt) <- gsub("N", "Proportion", colnames(multi_passage_proportion_combined_dt))
multi_passage_proportion_combined_dt <- starting_prop[multi_passage_proportion_combined_dt, on = "Assignment"]



multi_passage_proportion_cor <- data.table(Day = unique(multi_passage_proportion_combined_dt$Day), rho = as.numeric(NA), P = as.numeric(NA))

for (day in unique(multi_passage_proportion_combined_dt$Day)){
    print(day)
    test <- cor.test(multi_passage_proportion_combined_dt[Day == day]$Proportion, multi_passage_proportion_combined_dt[Day == day]$rate, method = "spearman")
    multi_passage_proportion_cor[Day == day, c("rho")] <- test$estimate
    multi_passage_proportion_cor[Day == day, c("P")] <- test$p.value
}

multi_passage_proportion_cor$FDR <- p.adjust(multi_passage_proportion_cor$P, method="BH")



p_multi_passage_prop_mean2 <- ggplot(multi_passage_proportion_combined_dt[Day != 0], aes(Proportion, rate, color = Assignment, size = N)) +
        geom_point(alpha = 0.8) +
        facet_wrap(vars(paste0("Passage ", Day)), nrow = 1, scales = "free_x") +
        scale_color_manual(values = colors) +
        theme_classic() +
        geom_smooth() +
        geom_smooth(aes(Proportion, rate), color = "black", span = 5) +
        xlab("Proportion of Village") +
        ylab("Growth Rate") +
        theme(axis.text.x = element_text(size = 12, angle = 90, hjust = .5, vjust = .5),
                axis.text.y = element_text(size = 12),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                strip.text.x = element_text(size = 20))

ggsave(p_multi_passage_prop_mean2, filename = paste0(outdir, "multi_passage_prop_rate_comparison2.png"), width =8, height = 3)
ggsave(p_multi_passage_prop_mean2, filename = paste0(outdir, "multi_passage_prop_rate_comparison2.pdf"), width =8, height = 3)



diff_proportion_combined_dt <- diff_props[multi_passage_ratios_dt_2, on = "Assignment", allow.cartesian=TRUE]
colnames(diff_proportion_combined_dt) <- gsub("N", "Proportion", colnames(diff_proportion_combined_dt))
diff_proportion_combined_dt <- starting_prop[diff_proportion_combined_dt, on = "Assignment"]


diff_proportion_cor_list <- list()
diff_proportion_cor <- data.table(Day = unique(diff_proportion_combined_dt$Day), rho = as.numeric(NA), P = as.numeric(NA))

for (day in unique(diff_proportion_combined_dt$Day)){
    print(day)
    test <- cor.test(diff_proportion_combined_dt[Day == day]$Proportion, diff_proportion_combined_dt[Day == day]$rate, method = "spearman")
    diff_proportion_cor[Day == day, c("rho")] <- test$estimate
    diff_proportion_cor[Day == day, c("P")] <- test$p.value
}

diff_proportion_cor$FDR <- p.adjust(diff_proportion_cor$P, method="BH")



p_diff_prop_mean2_d0 <- ggplot(diff_proportion_combined_dt, aes(Proportion, rate, color = Assignment, size = N)) +
        geom_point(alpha = 0.8) +
        facet_wrap(vars(factor(paste0("Day ",Day), levels = paste0("Day ",c(0:5,7,15)))), nrow = 1, scales = "free_x") +
        scale_color_manual(values = colors) +
        theme_classic() +
        geom_smooth() +
        geom_smooth(aes(Proportion, rate), color = "black", span = 5) +
        xlab("Proportion of Village") +
        ylab("Growth Rate") +
        theme(axis.text.x = element_text(size = 12, angle = 90, hjust = .5, vjust = .5),
                axis.text.y = element_text(size = 12),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                strip.text.x = element_text(size = 20))

ggsave(p_diff_prop_mean2_d0, filename = paste0(outdir, "dif_prop_rate_comparison2_w_d0.png"), width = 15, height = 3)
ggsave(p_diff_prop_mean2_d0, filename = paste0(outdir, "dif_prop_rate_comparison2_w_d0.pdf"), width = 15, height = 3)


