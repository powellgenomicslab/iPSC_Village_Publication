import pandas as pd


genes_file = "/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/QC/time-integrated_filtered_seurat_1pct_expressing_genes.tsv" ## Produced with iPSC_Village_Publication/scripts/hiPSC_village_18_lines/iPSC_maintenance/Quality_Control_and_Processing/hiPSC_village_18_lines_QC.R
genes = pd.read_csv(genes_file, sep = "\t")



rule all:
    input:
        expand("/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/Variance/gene_separated/fit_models/{gene}_fitted_models.rds", gene = genes.Gene),


rule partition_variance_integratedSCT:
    input:
        seurat =  ancient("/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/QC/time-integrated_filtered_seurat_1pct_expressing.rds") ## Produced with iPSC_Village_Publication/scripts/hiPSC_village_18_lines/iPSC_maintenance/Quality_Control_and_Processing/hiPSC_village_18_lines_QC.R
    output:
        "/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/Variance/gene_separated/icc/{gene}_icc.rds",
        "/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/Variance/gene_separated/fit_models/{gene}_fitted_models.rds"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 8,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 8
    threads: 4
    params:
        script = "iPSC_Village_Publication/scripts/hiPSC_village_18_lines/iPSC_maintenance/Variance_Explained/variance_partition_multipassage.R",
        out_icc="/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/Variance/gene_separated/icc/",
        out_model="/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/Variance/gene_separated/fit_models/",
        out_resids="/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/Variance/gene_separated/residuals4qtl/",
        out_icc_interaction = "/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/Variance/gene_separated/icc_interaction/",
        out_model_interaction = "/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/Variance/gene_separated/icc_interaction/"
    log:
    shell:
        """
        mkdir -p {params.out_icc_interaction}
        mkdir -p {params.out_model_interaction}
        mkdir -p {params.out_resids}
        
        /directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/baseR402/bin/Rscript {params.script} {params.out_icc_interaction} {params.out_icc} {params.out_model_interaction} {params.out_model} {params.out_resids} {wildcards.gene}
        """