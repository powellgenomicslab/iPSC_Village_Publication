import pandas as pd


genes_file = "/path/to/output/RNA_Velocity_Pseudotime/Variance_Explained/village_3_hiPSC_lines_high_quality_cells_integrated_1pct_expressing.tsv" ## Produced with "prepare_pseudotime.R" script from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication
genes = pd.read_csv(genes_file, sep = "\t")



rule all:
    input:
        expand("/path/to/output/RNA_Velocity_Pseudotime/Variance_Explained/gene_separated/icc/{gene}_icc.rds", gene = genes.Gene)


rule partition_variance:
    input:
        seurat = "/path/to/output/RNA_Velocity_Pseudotime/Variance_Explained/village_3_hiPSC_lines_high_quality_cells_integrated_1pct_expressing.rds" ## Produced with "prepare_pseudotime.R" script from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication
    output:
        "/path/to/output/RNA_Velocity_Pseudotime/Variance_Explained/gene_separated/icc/{gene}_icc.rds"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 64,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 64
    threads: 1
    params:
        script = "iPSC_Village_Publication/scripts/hiPSC_village_3_lines/Variance_Explained/RNA_Velocity_Pseudotime/pseudotime_effect.R", ## Available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication
        out_icc="/path/to/output/RNA_Velocity_Pseudotime/Variance_Explained/gene_separated/icc/",
        out_icc_interaction = "/path/to/output/RNA_Velocity_Pseudotime/Variance_Explained/gene_separated/icc_interaction/",
        out_model_interaction = "/path/to/output/RNA_Velocity_Pseudotime/Variance_Explained/gene_separated/icc_interaction/",
        out_plot = "/path/to/output/RNA_Velocity_Pseudotime/Variance_Explained/gene_separated/plots/",
        out_effects = "/path/to/output/RNA_Velocity_Pseudotime/Variance_Explained/gene_separated/effect_betas/"
    log:
    shell:
        """
        mkdir -p {params.out_icc_interaction}
        mkdir -p {params.out_model_interaction}
        mkdir -p {params.out_plot}
        mkdir -p {params.out_effects}
        
        /directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/baseR402/bin/Rscript {params.script} {params.out_icc_interaction} {params.out_icc} {params.out_plot} {wildcards.gene} {params.out_effects}
        """


