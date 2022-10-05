import pandas as pd


genes_file = "/path/to/output/eQTL_Comparison/variance_partition_cryo/seurat_integrated_cryopreserved_1pct_expressing_genes.tsv"
genes = pd.read_csv(genes_file, sep = "\t")



rule all:
    input:
        expand("/path/to/output/eQTL_Comparison/variance_partition_cryo/uniculture/gene_separated/icc/{gene}_icc.rds", gene = genes.Gene),
        expand("/path/to/output/eQTL_Comparison/variance_partition_cryo/village/gene_separated/icc/{gene}_icc.rds", gene = genes.Gene)



rule partition_variance_uniculture:
    input:
        seurat = "/path/to/output/eQTL_Comparison/variance_partition_cryo/seurat_integrated_cryopreserved_1pct_expressing_genes.rds"
    output:
        "/path/to/output/eQTL_Comparison/variance_partition_cryo/uniculture/gene_separated/icc/{gene}_icc.rds",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 16,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 16
    threads: 4
    params:
        script = "iPSC_Village_Publication/scripts/hiPSC_village_3_lines/eQTL_Comparison/Cryopreserved/variance_partition_cryo_uni_culture.R", # This script is available from the github: https://github.com/powellgenomicslab/iPSC_Village_Publication
        out_icc="/path/to/output/eQTL_Comparison/variance_partition_cryo/uniculture/gene_separated/icc/",
        out_model="/path/to/output/eQTL_Comparison/variance_partition_cryo/uniculture/gene_separated/fit_models/",
        out_resids="/path/to/output/eQTL_Comparison/variance_partition_cryo/uniculture/gene_separated/residuals4qtl/",
        out_icc_interaction = "/path/to/output/eQTL_Comparison/variance_partition_cryo/uniculture/gene_separated/icc_interaction/",
        out_model_interaction = "/path/to/output/eQTL_Comparison/variance_partition_cryo/uniculture/gene_separated/icc_interaction/"
    log:
    shell:
        """
        mkdir -p {params.out_icc_interaction}
        mkdir -p {params.out_model_interaction}
        mkdir -p {params.out_resids}
        
        Rscript {params.script} {params.out_icc_interaction} {params.out_icc} {params.out_model_interaction} {params.out_model} {params.out_resids} {wildcards.gene}
        """


rule partition_variance_village:
    input:
        seurat = "/path/to/output/variance_partition_fresh/seurat_integrated_all_times_clustered_1pct_expressing.rds"
    output:
        "/path/to/output/eQTL_Comparison/variance_partition_cryo/village/gene_separated/icc/{gene}_icc.rds",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 16,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 16
    threads: 4
    params:
        script = "iPSC_Village_Publication/scripts/hiPSC_village_3_lines/eQTL_Comparison/Cryopreserved/variance_partition_cryo_village.R", # This script is available from the github: https://github.com/powellgenomicslab/iPSC_Village_Publication
        out_icc = "/path/to/output/eQTL_Comparison/variance_partition_cryo/village/gene_separated/icc/",
        out_model = "/path/to/output/eQTL_Comparison/variance_partition_cryo/village/gene_separated/fit_models/",
        out_resids = "/path/to/output/eQTL_Comparison/variance_partition_cryo/village/gene_separated/residuals4qtl/",
        out_icc_interaction = "/path/to/output/eQTL_Comparison/variance_partition_cryo/village/gene_separated/icc_interaction/",
        out_model_interaction = "/path/to/output/eQTL_Comparison/variance_partition_cryo/village/gene_separated/icc_interaction/"
    log:
    shell:
        """
        mkdir -p {params.out_icc_interaction}
        mkdir -p {params.out_model_interaction}
        mkdir -p {params.out_resids}
        
        Rscript {params.script} {params.out_icc_interaction} {params.out_icc} {params.out_model_interaction} {params.out_model} {params.out_resids} {wildcards.gene}
        """
