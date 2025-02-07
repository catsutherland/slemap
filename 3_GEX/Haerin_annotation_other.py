# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sci
import scvi
import copy
import scanpy as sc
import re

# %config InlineBackend.print_figure_kwargs={'facecolor' : "w"}
# %config InlineBackend.figure_format='retina'
# %matplotlib inline  

pd.set_option('display.max_columns', None)
sc.set_figure_params(figsize=(5, 5),dpi=200)

# %% [markdown]
# ## Read in data with B/T/Other annotation and check annotations

# %%
SLEmap = sc.read("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/hj10/Results/2_YASCP_final/6_annotation/For_annotation.h5ad")
protein_SLEmap = sc.read("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/hj10/Results/2_YASCP_final/6_annotation/For_annotation_protein.h5ad")

# %%
SLEmap

# %%
sc.pl.umap(
    SLEmap,
    color=["Azimuth:predicted.celltype.l1","Azimuth:predicted.celltype.l2","VDJ","annotation_BTother"],
    legend_loc='on data',
    legend_fontsize='xx-small',
    ncols=2,
)

# %% [markdown]
# ## Subset to only the 'Other' cells and recluster

# %%
SLEmap_other = SLEmap[SLEmap.obs["annotation_BTother"] == "Other"]
protein_SLEmap_other = protein_SLEmap[protein_SLEmap.obs["annotation_BTother"] == "Other"]

# %%
sc.pp.neighbors(SLEmap_other, use_rep="X_totalVI",n_neighbors=15)
sc.tl.umap(SLEmap_other)
sc.tl.leiden(SLEmap_other, key_added="leiden_totalVI",resolution=1)

# %%
sc.pl.umap(
    SLEmap_other,
    color=["leiden_totalVI"],
    legend_loc='on data',
    legend_fontsize='xx-small',
    ncols=2,
)

# %%
sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(
    SLEmap_other,
    color=["Azimuth:predicted.celltype.l1","Azimuth:predicted.celltype.l2","VDJ","leiden_totalVI"],
    legend_fontsize='xx-small',
    ncols=2,
)

# %%
sc.pl.umap(
    SLEmap_other,
    color=['condition', 'recruitment_centre','poolID',"leiden_totalVI"],
    legend_fontsize='xx-small',
    ncols=2,
)

# %%
sc.pl.umap(
    SLEmap_other,
    color=['log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_gene_group__mito_transcript', 'log1p_total_counts_gene_group__ribo_protein'],
    legend_fontsize='xx-small',
    ncols=2,
)

# %% [markdown]
# ## Annotation using RNA-seq gene markers
# three ways of plotting: dot plots, violin plots, umaps

# %%
# to use gene symbolsl instead of ensembl ids 
SLEmap_other.var_names = SLEmap_other.var["gene_symbols"]

# %%
# level1 cell type markers
sc.settings.set_figure_params(figsize=(4,3))
marker_genes_dict = {
    'T-cell': ['CD3D'],
    'B-cell': ['CD79A', 'MS4A1','CD19'],
    'NK': ['GNLY', 'NKG7'],
    'Myeloid':['HLA-DRA','CST3'],
    'Monocytes': ['ITGAM'],
    'Dendritic': ['CD74','HLA-DPA1','FLT3'],
    'ILC':['IL7R','KIT'],
    'HSPC':['CD34']
}
sc.pl.dotplot(SLEmap_other, marker_genes_dict, 'leiden_totalVI', dendrogram=False,use_raw=False)

# %%
# level2 cell type markers
sc.settings.set_figure_params(figsize=(4,3))
marker_genes_dict = {
    'NK_CD56Bright': ['NCAM1','SELL','GZMK'],
    'NK_proliferating': ['MKI67'],
    'cDC': ['BST2'],
    'pDC':['IL3RA'],
    'Classical_Monocytes': ['CD14'],
    'NonClassical_Monocytes': ['FCGR3A'],
}
sc.pl.dotplot(SLEmap_other, marker_genes_dict, 'leiden_totalVI', dendrogram=False,use_raw=False)

# %%
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(SLEmap_other, ['MKI67'], groupby='leiden_totalVI', use_raw=False)

# %%
sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(SLEmap_other, use_raw=False,color=['MKI67'])

# %% [markdown]
# ## Annotation using CITE-seq surface protein markers¶

# %%
protein_SLEmap_other.obs["leiden_totalVI"] = SLEmap_other.obs["leiden_totalVI"] 
protein_SLEmap_other.obsm["X_umap"] = SLEmap_other.obsm["X_umap"] 
protein_SLEmap_other.uns["leiden_totalVI_colors"] = SLEmap_other.uns["leiden_totalVI_colors"] 

# %%
sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(
    protein_SLEmap_other,
    color=['total_ADTlevels', 'log1p_total_ADTlevels','leiden_totalVI'],
    legend_loc='on data',
    legend_fontsize='xx-small',
    ncols=3,
)

# %%
#list of all CITE-seq markers
protein_SLEmap_other.var.index.tolist()

# %%
sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(protein_SLEmap_other, use_raw=False,color=['anti-human_CD3','anti-human_CD4','anti-human_CD8','anti-human_CD19','anti-human_CD123','anti-human_CD45RA','anti-human_CD11c','anti-human_CD16','anti-human_CD14','anti-human_CD56','anti-human_CD161','anti-human_CD141_(Thrombomodulin)'])

# %%
##Surface protein levels
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(protein_SLEmap_other, ['anti-human_CD11c'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_other, ['anti-human_CD16'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_other, ['anti-human_CD56'], groupby='leiden_totalVI', use_raw=False)#CD56Bright NK 
sc.pl.violin(protein_SLEmap_other, ['anti-human_CD14'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_other, ['anti-human_HLA-DR'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_other, ['anti-human_CD123'], groupby='leiden_totalVI', use_raw=False)

# %%
sc.settings.set_figure_params(figsize=(4,3))
marker_genes_dict = {
    'NK_CD56Bright': ['anti-human_CD56','anti-human_CD62L'],
    'Myeloid':['anti-human_HLA-DR'],
    'Dendritic':['anti-human_CD127_(IL-7R_)', 'anti-human_CD11c', 'anti-human_CD141_(Thrombomodulin)'],
    'cDC':['anti-human_CD1c'],
    'pDC':['anti-human_CD123','anti-human_CD45RA'],
    'Monocytes': ['anti-human_CD11b'],
    'Classical_Monocytes':['anti-human_CD14'],
    'Non_Classical_Monocytes':['anti-human_CD16'],
}
sc.pl.dotplot(protein_SLEmap_other, marker_genes_dict, 'leiden_totalVI', dendrogram=False,use_raw=False)

# %% [markdown]
# ## Add new annotation to adata

# %%
old_to_new = {
'0':'NK',
'1':'NK',
'2':'NK',
'3':'Mono',
'4':'NK',
'5':'Mono',
'6':'Mono',
'7':'NK',
'8':'Mono',
'9':'DC',
'10':'HSPC',
'11':'Mono',
'12':'DC',
'13':'ILC',
'14':'NK',
'15':'Mono',
'16':'Mono',
'17':'NK',
}
SLEmap_other.obs['manual_annotation_l1'] = (
SLEmap_other.obs['leiden_totalVI']
.map(old_to_new)
.astype('category')
)

# %%
old_to_new = {
    '0':'CD56Dim_NK',
    '1':'CD56Dim_NK',
    '2':'CD56Dim_NK',
    '3':'Classical_Mono',
    '4':'CD56Bright_NK',
    '5':'Classical_Mono',
    '6':'Classical_Mono',
    '7':'CD56Dim_NK',
    '8':'Non_classical_Mono',
    '9':'cDC',
    '10':'HSPC',
    '11':'Classical_Mono',
    '12':'pDC',
    '13':'ILC',
    '14':'CD56Dim_NK',
    '15':'Classical_Mono',
    '16':'Classical_Mono',
    '17':'Proliferating_NK',
}
SLEmap_other.obs['manual_annotation_l2'] = (
SLEmap_other.obs['leiden_totalVI']
.map(old_to_new)
.astype('category')
)

# %%
# check obs for annotation
SLEmap_other.obs

# %%
sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(SLEmap_other, use_raw=False,color=['manual_annotation_l1','manual_annotation_l2'])

# %%
SLEmap_other.obs['manual_annotation_l1'].value_counts()

# %%
SLEmap_other.obs['manual_annotation_l2'].value_counts()

# %%
SLEmap_other.write("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/hj10/Results/2_YASCP_final/6_annotation/manualannotation_Other.h5ad")
