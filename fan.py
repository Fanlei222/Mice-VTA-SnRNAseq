import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import re
import matplotlib.gridspec as gridspec
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import os
from scipy import stats
import scvelo as scv 
import cellrank as cr  
import anndata as ad 
scv.settings.verbosity = 3 
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False) 
cr.settings.verbosity = 2

##输入文件
adata = sc.read_h5ad("cluster4_splice.h5ad")
pca =pd.read_csv("cluster4_subcluster_pca.csv")
adata.obsm["X_pca"] = pca.loc[:,pca.columns != 'Unnamed: 0'].to_numpy()
umap = pd.read_csv("cluster4_subcluster_umap.csv")
adata.obsm["X_umap"] = umap.loc[:,umap.columns != 'Unnamed: 0'].to_numpy()
cluster = pd.read_csv("cluster4_subcluster_metadata.csv")
adata.obs['cluster'] = cluster.loc[:,'cell_005'].to_numpy()

##绘图
sc.pl.umap(adata,color=['cluster','group'],save='_umap_group.pdf')

sc.tl.embedding_density(adata,basis='umap', groupby='group')
sc.pl.embedding_density(adata,basis='umap', key='umap_density_group',save='_umap_density.pdf')

###Add-vs-Pre
addvspre = adata[adata.obs.group != "Post-addiction"]
addvspre.obs['group'].unique()

sc.pl.umap(addvspre,color=['cluster','group'])
sc.pl.umap(addvspre,color=['Addicting/pre_density_diff'])

###RNA速率分析
scv.pl.proportions(addvspre,groupby="cluster",save=False)

scv.pp.filter_and_normalize(addvspre, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(addvspre, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(addvspre)
scv.tl.velocity(addvspre, mode='dynamical')
scv.tl.velocity_graph(addvspre)

scv.pl.velocity_embedding_stream(addvspre, basis='umap',color='cluster',save='_add_vs_pre_rna_scv.pdf')
scv.pl.velocity_embedding_stream(addvspre, basis='umap',color='Addicting/pre_density_diff',save='_add_vs_pre_rna_scv_density.pdf')

##Post-vs-Pre
postvspre = adata[adata.obs.group != "Addicting"]
postvspre.obs['group'].unique()

sc.pl.umap(postvspre,color=['cluster','group'])
sc.pl.umap(postvspre,color=['post/pre_density_diff'])

##RNA速率分析
scv.pl.proportions(postvspre,groupby="cluster",save=False)

scv.pp.filter_and_normalize(postvspre, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(postvspre, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(postvspre)
scv.tl.velocity(postvspre, mode='dynamical')
scv.tl.velocity_graph(postvspre)

scv.pl.velocity_embedding_stream(postvspre, basis='umap',color='cluster',save='_post_vs_pre_rna_scv.pdf')
scv.pl.velocity_embedding_stream(postvspre, basis='umap',color='post/pre_density_diff',save='_post_vs_pre_rna_scv_density.pdf')