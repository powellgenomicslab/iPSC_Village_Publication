import pandas as pd


genes_file = "/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/eQTL_Comparison/deboever_imputed_overlapping_filtered_header_pruned_snp_gene.tsv"

genes = pd.read_csv(genes_file, sep = "\t")


rule all:
    input:
        expand("/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/eQTL_Comparison/deboever/gene_separate/beds/{gene}_deboever_eQTL_results.bed", gene = genes.gene_id),


rule eqtl:
    input:
        bed = "/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/eQTL_Comparison/deboever_imputed_overlapping_filtered_header_pruned.bed"
    output:
        "/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/eQTL_Comparison/deboever/gene_separate/beds/{gene}_deboever_eQTL_results.bed"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 4,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 4
    threads: 4
    params:
        script = "/path/to/scripts/hiPSC_village_18_lines/iPSC_maintenance/eQTL_Comparison/test_eQTL.R",
        outdir="/path/to/output/hiPSC_village_18_lines/iPSC_maintenance/eQTL_Comparison/deboever/gene_separate/",
        snp=lambda wildcards: genes.ID[genes.gene_id == wildcards.gene],
    log:
    shell:
        """        
        Rscript {params.script} {wildcards.gene} {params.outdir} {input.bed}
        """
