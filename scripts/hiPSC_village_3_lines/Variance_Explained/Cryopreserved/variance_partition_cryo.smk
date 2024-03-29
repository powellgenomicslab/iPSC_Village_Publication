import pandas as pd


genes_file = "/path/to/output/variance_partition_cryo/seurat_integrated_cryopreserved_1pct_expressing_genes.tsv"
genes = pd.read_csv(genes_file, sep = "\t")


rule all:
    input:
        expand("/path/to/output/variance_partition_cryo/gene_separated/fit_models/{gene}_fitted_models.rds", gene = genes.Gene),
        expand("/path/to/output/variance_partition_cryo/uniculture/gene_separated/icc/{gene}_icc.rds", gene = genes.Gene),
        expand("/path/to/output/variance_partition_cryo/village/gene_separated/icc/{gene}_icc.rds", gene = genes.Gene)


rule partition_variance:
    input:
        seurat = "/path/to/output/variance_partition_cryo/seurat_integrated_cryopreserved_1pct_expressing_genes.rds"
    output:
        "/path/to/output/variance_partition_cryo/gene_separated/icc/{gene}_icc.rds",
        "/path/to/output/variance_partition_cryo/gene_separated/fit_models/{gene}_fitted_models.rds"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 16,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 16
    threads: 4
    params:
        script = "iPSC_Village_Publication/scripts/hiPSC_village_3_lines/Variance_Explained/Cryopreserved/variance_partition_post_review_cryo.R", # This script is available from the github: https://github.com/powellgenomicslab/iPSC_Village_Publication
        out_icc="/path/to/output/variance_partition_cryo/gene_separated/icc/",
        out_model="/path/to/output/variance_partition_cryo/gene_separated/fit_models/",
        out_resids="/path/to/output/variance_partition_cryo/gene_separated/residuals4qtl/",
        out_icc_interaction = "/path/to/output/variance_partition_cryo/gene_separated/icc_interaction/",
        out_model_interaction = "/path/to/output/variance_partition_cryo/gene_separated/icc_interaction/"
    log:
    shell:
        """
        mkdir -p {params.out_icc_interaction}
        mkdir -p {params.out_model_interaction}
        mkdir -p {params.out_resids}
        
        Rscript {params.script} {params.out_icc_interaction} {params.out_icc} {params.out_model_interaction} {params.out_model} {params.out_resids} {wildcards.gene}
        """
        
