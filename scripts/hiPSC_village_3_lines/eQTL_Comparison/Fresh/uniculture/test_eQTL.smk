import pandas as pd


genes_file = "/path/to/output/eQTL_Comparison/Fresh/uni_village/uniculture/deboever_imputed_overlapping_filtered_header_pruned_snp_gene.tsv"
genes = pd.read_csv(genes_file, sep = "\t")


rule all:
    input:
        expand("/path/to/output/eQTL_Comparison/Fresh/uni_village/uniculture/gene_separate/beds/{gene}_deboever_eQTL_results.bed", gene = genes.gene_id),


rule eqtl:
    input:
        bed = "/path/to/output/eQTL_Comparison/deboever_imputed_overlapping_filtered_header_pruned.bed"
    output:
        "/path/to/output/eQTL_Comparison/Fresh/uni_village/uniculture/gene_separate/beds/{gene}_deboever_eQTL_results.bed"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 4,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 4
    threads: 4
    params:
        script = "iPSC_Village_Publication/scripts/hiPSC_village_3_lines/eQTL_Comparison/Fresh/uniculture/test_eQTL.R",
        outdir="/path/to/output/eQTL_Comparison/Fresh/uni_village/uniculture/gene_separate/",
    log:
    shell:
        """        
        Rscript {params.script} {wildcards.gene} {params.outdir} {input.bed}
        """

