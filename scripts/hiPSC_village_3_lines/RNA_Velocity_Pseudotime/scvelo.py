##### run with conda environment "scvelo"

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import os
import anndata as ad
import matplotlib as plt
import numba
import umap
import itertools
from combat.pycombat import pycombat
import math
from scipy.sparse import csr_matrix
import skmisc
from plotnine import *

datadir = "/path/to/output/RNA_Velocity_Pseudotime/preprocess" ## Produced with "Pull_Seurat_Data.R" script from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication
outdir="/path/to/output/RNA_Velocity_Pseudotime/scvelo/"
os.makedirs(outdir)
os.chdir(outdir)

# scvelo settings
scv.logging.print_version()
# Running scvelo 0.2.3 (python 3.8.5) on 2021-05-31 10:45.
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.set_figure_params('scvelo', figsize=(40,40))  # for beautified visualization

# scanpy settings
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
# scanpy==1.7.2 anndata==0.7.6 numpy==1.20.3 scipy==1.6.3 pandas==1.2.4 scikit-learn==0.24.2
sc.settings.set_figure_params(dpi=200)


### Load in sample data info
samples = pd.read_table('iPSC_Village_Publication/scripts/hiPSC_village_3_lines/Quality_Control_and_Processing/Demultiplexing_Doublet_Detecting/Sample_meta.tsv', delim_whitespace=True, header=0) # available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication


### Read in data
adata = list()

for i in range(0,8):
	adata.append(scv.read('/path/to/cellranger/gene_expression/' + samples['Pool'][i] + '/velocyto/' + samples['Pool'][i] + '.loom'], cache=True)) ## Produced by iPSC_Village_Publication/scripts/hiPSC_village_3_lines/RNA_Velocity_Pseudotime/velocyto.sh available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication


### Update droplet names to match seurat names
for i in range(0,8):
	adata[i].obs_names = [row.replace("x", "") for row in [row.replace(":", "_") for row in adata[i].obs_names]]
	adata[i].var_names_make_unique()


#### Combine the data together
adata_concat = ad.concat(adata)


### Get metadata
meta = pd.read_csv(datadir + "metadata.csv")
# ## PHATE
## Seurat
umap_cord = pd.read_csv(datadir + "cell_embeddings.csv")


### Filter anndata object for droplet ids
sample_obs = pd.read_csv(datadir + "cellID_obs.csv")
adata_concat = adata_concat[np.isin(adata_concat.obs.index,sample_obs["x"])]


### Order the metadata droplets in same orderas adata
adata_barcodes = pd.DataFrame(adata_concat.obs.index)
adata_barcodes = adata_barcodes.rename(columns = {0:'Cell ID'})

meta = meta.rename(columns = {'Unnamed: 0':'Cell ID'})
meta_ordered = adata_barcodes.merge(meta, on = "Cell ID")


### Order the umap droplets in same orderas adata
umap = umap_cord.rename(columns = {'Unnamed: 0':'Cell ID'})
umap_ordered = adata_barcodes.merge(umap, on = "Cell ID")


### Update metadata columns
meta_ordered['Site'] = meta_ordered['MULTI_ID'].str.replace('\d','')

## Add cryopreserved to site for necessary syndey samples
meta_ordered['Site'][["Thawed" in i for i in meta_ordered['Time']]] = meta_ordered['Site']+ "cryopreserved"
meta_ordered['Time'][["Village Day 4" in i for i in meta_ordered['Time']]] = "Village"
meta_ordered['Time'][["Thawed Village Day 7" in i for i in meta_ordered['Time']]] = "Village"
meta_ordered['Time'][["Thawed Village Day 0" in i for i in meta_ordered['Time']]] = "Baseline"

meta_ordered['Site_Individual_Time'] = meta_ordered['Site'] + "-" + meta_ordered['Final_Assignment'] + "-" + meta_ordered['Time']



### Add metadata to anndata object
adata_concat.obs['clusters'] = meta_ordered['integrated_snn_res.0.28'].values
adata_concat.obs['individual'] = meta_ordered['Final_Assignment'].values
adata_concat.obs['cell_cycle'] = meta_ordered['phases'].values
adata_concat.obs['site_rep'] = meta_ordered['Site_rep'].values
adata_concat.obs['Location'] = meta_ordered['Site'].values
adata_concat.obs['Time'] = meta_ordered['Time'].values
adata_concat.obs['Site_Individual_Time'] = meta_ordered['Site_Individual_Time'].values
adata_concat.obs['Seurat_UMAP_1'] = umap_ordered['UMAP_1'].values
adata_concat.obs['Seurat_UMAP_2'] = umap_ordered['UMAP_2'].values



### Filter
# Also filter cells with too few unspliced counts
adata_concat.obs['n_unspliced_counts'] = adata_concat.layers['unspliced'].sum(axis=1).A1
adata_concat = adata_concat[adata_concat.obs['n_unspliced_counts'] >= 1000]

# Filter genes
sc.pp.filter_genes(adata_concat, min_cells=20)

# Filter genes that don't have unspliced counts detected in at least 10 cells
adata_concat.var['n_unspliced_counts'] = np.count_nonzero(adata_concat.layers['unspliced'].toarray(), axis=0)
adata_concat = adata_concat[:,adata_concat.var['n_unspliced_counts']>10]
scv.pp.filter_and_normalize(adata_concat)
sc.pp.highly_variable_genes(adata_concat, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = "Site_Individual_Time")


numerics = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']


Mtmp = pd.DataFrame.sparse.from_spmatrix(adata_concat.layers['spliced'] + adata_concat.layers['unspliced'])
M = Mtmp.T
Rtmp = pd.DataFrame(adata_concat.layers['spliced']/(adata_concat.layers['spliced'] + adata_concat.layers['unspliced']))
R = Rtmp.T



##### Batch Normalize the data #####
Mb = pycombat(M, adata_concat.obs['Site_Individual_Time'])



##### Convert back to S and U matrices
Sb = Mb * R
Sb = Sb.replace(np.nan, 0).T
Ub = Mb * (1 - R)
Ub = Ub.replace(np.nan, 0).T



##### Add the updated dataframes to the adata
adata_concat.layers['spliced'] = csr_matrix(Sb)
adata_concat.layers['unspliced'] = csr_matrix(Ub)
adata_concat.layers['matrix'] = csr_matrix(Mb.T)


### Get abundances of spliced to unspliced
scv.utils.show_proportions(adata_concat)


##### RNA Velocity #####
scv.pp.moments(adata_concat)
scv.tl.velocity(adata_concat, mode = "stochastic")
scv.tl.velocity_graph(adata_concat)
scv.tl.velocity_pseudotime(adata_concat)



### scvelo umap
## UMAP + clustering
# Identify highly-variable genes.
sc.pl.highly_variable_genes(adata_concat, save= "variable_genes.png")
# Scale each gene to unit variance. Clip values exceeding standard deviation 10.
sc.pp.scale(adata_concat, max_value=1)
# Reduce the dimensionality of the data by running principal component analysis (PCA), which reveals the main axes of variation and denoises the data.
sc.tl.pca(adata_concat, svd_solver='arpack', n_comps = 50)
# compute the neighborhood graph of cells using the PCA representation of the data matrix. 
# You might simply use default values here. For the sake of reproducing Seurat’s results, let’s take the following values.
sc.pp.neighbors(adata_concat, n_neighbors=10)
# Embedding the neighborhood graph - reccomend embedding the graph in 2 dimensions using UMAP (McInnes et al., 2018).
#  It is potentially more faithful to the global connectivity of the manifold than tSNE, i.e., it better preservers trajectories.
#  In some ocassions, you might still observe disconnected clusters and similar connectivity violations. They can usually be remedied by running:
sc.tl.umap(adata_concat)
# Clustering the neighborhood graph - recommend the Leiden graph-clustering method (community detection based on optimizing modularity) by Traag *et al.* (2018).
# Note that Leiden clustering directly clusters the neighborhood graph of cells, which we already computed in the previous section.
sc.tl.leiden(adata_concat)
sc.pl.umap(adata_concat, color=['clusters'], save= "umap_seurat_clusters.pdf")
sc.pl.umap(adata_concat, color=['clusters'], save= "umap_seurat_clusters.png")
sc.pl.umap(adata_concat, color=['Time'], save= "umap_time.png")
scv.pl.velocity_embedding_stream(adata_concat, basis='umap', color=['leiden'], save= "stream.png", size=10)
# Creates a plot AND adds embedded velocity vectors 'velocity_umap' to the obsm slot
scv.pl.velocity_embedding(adata_concat, basis='umap', arrow_length=10, arrow_size=50, dpi=600, color=['clusters'], save= "arrow.png")
scv.pl.velocity_embedding(adata_concat, basis='umap', arrow_length=5, arrow_size=5, dpi=600, color=['Site_Individual_Time'], save= "arrow_batch.png")
# Just a different representation, nothing new is added to adata




scv.tl.recover_dynamics(adata_concat, fit_basal_transcription=True, max_iter=20)
scv.tl.velocity(adata_concat, mode='dynamical')
scv.tl.velocity_graph(adata_concat)
scv.tl.recover_latent_time(adata_concat)
adata_concat.write(filename = outdir +  "adata_concat.h5ad")


scv.pl.scatter(adata_concat, color='latent_time', fontsize=24, size=100, color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0,1], save= "latent_time.png")
top_gene = adata_concat.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.heatmap(adata_concat, var_names = top_gene[0:100], sortby='latent_time', col_color='clusters', yticklabels=True, n_convolve=100, figsize=(40, 20), save= "heatmap.png")
scv.pl.heatmap(adata_concat, var_names = top_gene[0:100], sortby='velocity_pseudotime', col_color='clusters', yticklabels=True, n_convolve=100, figsize=(40, 20), save= "pseudotime_heatmap.png")
scv.pl.heatmap(adata_concat, var_names = top_gene[0:100], sortby='latent_time', col_color='Site_Individual_Time', yticklabels=True, n_convolve=100, figsize=(40, 20), save= "line_location_time_heatmap.png")
scv.pl.scatter(adata_concat, color='clusters', fontsize=24, size=10, color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0,1], save= "seurat_clusters.png")
scv.pl.scatter(adata_concat, color='velocity_pseudotime', fontsize=24, size=10, color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0,1], save= "velocity_pseudotime_umap.png")
scv.pl.velocity_embedding_stream(adata_concat, basis='umap', color=['latent_time'], save= "stream_latent_colors.png", size=10)
scv.pl.velocity_embedding_stream(adata_concat, basis='umap', color=['clusters'], save= "stream_seurat_clusters.png", size=10)

plot = (ggplot(adata_concat.obs, aes(x='latent_time', color='Time', fill='Time')) + geom_density(alpha=0.1))
plot.save(filename = "density_plot.png")

plot2 = (ggplot(adata_concat.obs, aes(x='latent_time', color='Time', fill='Time')) + geom_density(alpha=0.1) + facet_grid('Location ~ individual'))
plot2.save(filename = "density_plot_facet.png")

plot3 = (ggplot(adata_concat.obs, aes(x='Seurat_UMAP_1', y = 'Seurat_UMAP_2', color='latent_time')) + geom_point(size = 0.5, alpha=0.5) + theme_classic())
plot3.save(filename = "umap_seurat_lateny.png", dpi = 600, width = 10, height = 10)


scv.pl.scatter(adata_concat, basis=top_gene[:15], ncols=5, frameon=False, color=['clusters'], save= "top_genes.png", size=20, alpha = 0.01)


adata_concat.obs['UMAP_1'] = pd.DataFrame(adata_concat.obsm['X_umap'])[0].values
adata_concat.obs['UMAP_2'] = pd.DataFrame(adata_concat.obsm['X_umap'])[1].values
adata_concat.obs.to_csv(outdir + "metadata.csv")
adata_concat.var.to_csv(outdir + "var_meta.csv")


##### Experimenting with figures #####
adata_concat = scv.read(filename = outdir +  "adata_concat.h5ad")


### Read in latent genes for heatmap figure ###
genes4latent = pd.read_table('iPSC_Village_Publication/scripts/hiPSC_village_3_lines/RNA_Velocity_Pseudotime/latent_heatmap_markers.tsv', delim_whitespace=True, header=0, sep = '\t') # available from iPSC_Village_Publication github: https://github.com/powellgenomicslab/iPSC_Village_Publication


### Make heatmap of desired genes ###
scv.pl.heatmap(adata_concat, var_names = genes4latent.Gene.values, sortby='latent_time', col_color='latent_time', color_map = 'inferno', yticklabels=True, n_convolve=100, figsize=(10, 5), row_cluster=False, save= "latent_genes_heatmap.pdf")


