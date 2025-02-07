# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: groups/team282/SLEmap_py/5
#     language: python
#     name: chl0ag9uaehhss9zb2z0cgfjay9ncm91chmvdgvhbti4mi9ttevtyxbfchkvnqo.
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
import seaborn as sns

# %config InlineBackend.print_figure_kwargs={'facecolor' : "w"}
# %config InlineBackend.figure_format='retina'
# %matplotlib inline  

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 100)

sc.set_figure_params(figsize=(5, 5),dpi=200)


# %% [markdown]
# # Read in data with B/T/Other annotation and check annotations

# %%
SLEmap = sc.read("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/hj10/Results/2_YASCP_final/6_annotation/For_annotation.h5ad")
protein_SLEmap = sc.read("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/hj10/Results/2_YASCP_final/6_annotation/For_annotation_protein.h5ad")

# %%
SLEmap.obs.externalID.nunique()

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
# # Subset to only the 'B' cells and recluster

# %%
SLEmap_B = SLEmap[SLEmap.obs["annotation_BTother"] == "B"]
protein_SLEmap_B = protein_SLEmap[protein_SLEmap.obs["annotation_BTother"] == "B"]

# %%
annot_info = pd.read_csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/cs54/final_data/2_BCR_TCR/1_IgBLAST_assignment/bcr_df_heavy_single_bcr.csv")

annot_info

# %%
obs_tmp = SLEmap_B.obs.copy()

obs_tmp["cell_barcode"] = obs_tmp.index.str.split("-SLE").str[0]

obs_tmp = obs_tmp.reset_index().merge(how = 'left', right = annot_info, left_on = ["cell_barcode", "poolID"], right_on = ["barcode", "pool"]).set_index('index')

obs_tmp = obs_tmp.drop(["cell_barcode", "barcode", "pool"], axis = 1)

# %%
SLEmap_B.obs = obs_tmp

# %%
sc.pp.neighbors(SLEmap_B, use_rep="X_totalVI",n_neighbors=15)
sc.tl.umap(SLEmap_B)

# %%
sc.tl.leiden(SLEmap_B, key_added="leiden_totalVI_0.3",resolution=0.3)
sc.tl.leiden(SLEmap_B, key_added="leiden_totalVI_0.4",resolution=0.4)
sc.tl.leiden(SLEmap_B, key_added="leiden_totalVI_0.5",resolution=0.5)
sc.tl.leiden(SLEmap_B, key_added="leiden_totalVI_0.6",resolution=0.6)

# %%
sc.tl.leiden(SLEmap_B, key_added="leiden_totalVI_0.7",resolution=0.7)
sc.tl.leiden(SLEmap_B, key_added="leiden_totalVI_0.8",resolution=0.8)
sc.tl.leiden(SLEmap_B, key_added="leiden_totalVI_0.9",resolution=0.9)
sc.tl.leiden(SLEmap_B, key_added="leiden_totalVI_1",resolution=1)


# %%
sc.tl.leiden(SLEmap_B, key_added="leiden_totalVI",resolution=1)


# %%
sc.set_figure_params(figsize=(5, 5),dpi=150)

sc.pl.umap(
    SLEmap_B,
    color=["leiden_totalVI_0.3", "leiden_totalVI_0.4", "leiden_totalVI_0.5", "leiden_totalVI_0.6",
          "leiden_totalVI_0.7", "leiden_totalVI_0.8", "leiden_totalVI_0.9", "leiden_totalVI_1"],
    legend_loc='on data',
    legend_fontsize='xx-small',
    ncols=2,
)

# %%
c_call_props = SLEmap_B.obs.groupby("leiden_totalVI").c_call.value_counts(normalize = True).reset_index()
c_call_props.query("leiden_totalVI == '15'")

# %%
sc.set_figure_params(figsize=(5, 5),dpi=150)

sc.pl.umap(
    SLEmap_B,
    color=["c_call", "leiden_totalVI"],
    legend_fontsize='xx-small',
    ncols=2,
)

# %%
sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(SLEmap_B, use_raw=False,color=['CD27', "CXCR5", "FOXP3", "ITGAX", "TBX21"])

# %%

# %%
sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(
    SLEmap_B,
    color=["Azimuth:predicted.celltype.l1","Azimuth:predicted.celltype.l2","VDJ","leiden_totalVI", 
           "Celltypist:Immune_All_Low:majority_voting"],
    legend_fontsize='xx-small',
    ncols=2,
)


# %%
def cluster_small_multiples(adata, clust_key, size=1, frameon=False, legend_loc=None, **kwargs):
    tmp = adata.copy()

    for i,clust in enumerate(adata.obs[clust_key].cat.categories):
        tmp.obs[clust] = adata.obs[clust_key].isin([clust]).astype('category')
        tmp.uns[clust+'_colors'] = ['#d3d3d3', adata.uns[clust_key+'_colors'][i]]

    sc.pl.umap(tmp, groups=tmp.obs[clust].cat.categories[1:].values, 
               color=adata.obs[clust_key].cat.categories.tolist(), size=size, 
               frameon=frameon, legend_loc=legend_loc, **kwargs)

cluster_small_multiples(SLEmap_B, "Celltypist:Immune_All_Low:majority_voting")

# %%
cluster_small_multiples(SLEmap_B, "Azimuth:predicted.celltype.l2")

# %%
cluster_6_cells = SLEmap_B.obs.query("(recruitment_centre == 'Kings') & (leiden_totalVI == '6')").externalID.value_counts().reset_index()

# %%
SLEmap_B.obs.leiden_totalVI.value_counts()

# %%
SLEmap_B.obs.groupby("leiden_totalVI").recruitment_centre.value_counts(normalize = True)

# %%
total_B_cells = SLEmap_B.obs.query("(recruitment_centre == 'Kings')").externalID.value_counts().reset_index()

# %%
total_B_cells = total_B_cells.rename({"count":"total_count"}, axis =1)

# %%
cluster_6_kings_total = pd.merge(total_B_cells,cluster_6_cells )

# %%
sc.set_figure_params(figsize=(3,3),dpi=150)

sns.scatterplot(data = cluster_6_kings_total, x = "total_count", y = "count")

# %%
sc.pl.umap(
    SLEmap_B,
    color=['condition', 'recruitment_centre','poolID',"leiden_totalVI"],
    legend_fontsize='xx-small',
    ncols=2,
)

# %%
sc.pl.umap(
    SLEmap_B,
    color=['log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_gene_group__mito_transcript', 'log1p_total_counts_gene_group__ribo_protein'],
    legend_fontsize='xx-small',
    ncols=2,
)

# %% [markdown]
# # Annotation using RNA-seq gene markers
# three ways of plotting: dot plots, violin plots, umaps

# %%
# to use gene symbolsl instead of ensembl ids 
SLEmap_B.var_names = SLEmap_B.var["gene_symbols"]

# %%
SLEmap_B_raw = SLEmap_B.raw.to_adata().copy()
SLEmap_B_raw.var_names = SLEmap_B_raw.var["gene_symbols"]

# %%
sc.pp.normalize_total(SLEmap_B_raw)
sc.pp.log1p(SLEmap_B_raw)

# %%
SLEmap_B_raw

# %% [markdown]
# # Marker gene identification

# %%
sc.pp.highly_variable_genes(SLEmap_B, n_top_genes=2000, batch_key="poolID")

# %%
# sc.tl.rank_genes_groups(SLEmap_B, groupby="leiden_totalVI", method="wilcoxon", use_raw = False)

# %% [markdown]
# # Key B cell markers

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
sc.pl.dotplot(SLEmap_B, marker_genes_dict, 'leiden_totalVI', dendrogram=False,use_raw=False)

# %%
B_markers = ["CD27", "CD38", "CD24", "CR2", "FAS", "CD86", "ITGAX", "TBX21", 
               "SLAMF7", "CXCR5",  "SDC1", "IL10", "IL12A", "EBI3", "FCRL4", "FCRL5",
               "IL4R", "MKI67", "JCHAIN", "TCL1A", "CLEC2B"]

# %%
sc.pl.dotplot(SLEmap_B, B_markers, 'leiden_totalVI', dendrogram=False,use_raw=False)

# %%
sc.pl.dotplot(SLEmap_B_raw, B_markers, 'leiden_totalVI', dendrogram=False,use_raw=False)

# %%
B_markers = ["CD27", "CD38", "CD24", "CR2", "FAS", "CD86", "ITGAX", "TBX21", 
               "SLAMF7", "CXCR5",  "SDC1", "IL10", "IL12A", "EBI3", "FCRL4", "FCRL5",
               "IL4R", "MKI67", "TCL1A", "CLEC2B", "CD19", "MS4A1"]

sc.pl.dotplot(SLEmap_B, B_markers, 'leiden_totalVI', dendrogram=False,use_raw=False)

# %%
sc.settings.set_figure_params(figsize=(5,5))
isotype_markers = ["IGHA1", "IGHA2","IGHD", "IGHE", "IGHG1", "IGHG2", "IGHG3", "IGHM"]

sc.pl.dotplot(SLEmap_B_raw, isotype_markers, 'leiden_totalVI', dendrogram=False)


# %%
sc.set_figure_params(figsize=(3,3),dpi=150)

sc.pl.umap(
    SLEmap_B,
    color=["leiden_totalVI"],
    legend_loc='on data',
    legend_fontsize='xx-small',
    ncols=2,
)

# %%
sc.set_figure_params(figsize=(3,3),dpi=150)

sc.pl.umap(
    SLEmap_B,
    color=["Celltypist:Immune_All_Low:majority_voting"],
    legend_fontsize='xx-small',
    ncols=2,
)

# %%
sc.pl.stacked_violin(SLEmap_B, B_markers, 'leiden_totalVI', dendrogram=False,use_raw=False)

# %%
isotype_markers = ["IGHA1", "IGHA2","IGHD", "IGHE", "IGHG1", "IGHG2", "IGHG3", "IGHM"]

sc.pl.dotplot(SLEmap_B_raw, isotype_markers, 'leiden_totalVI', dendrogram=False)

# %%
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(SLEmap_B, ['CD27'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(SLEmap_B, ['CR2'], groupby='leiden_totalVI', use_raw=False)

# %%
sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(SLEmap_B, use_raw=False,color=['CD27', "CXCR5", "FOXP3", "ITGAX", "TBX21"])

# %% [markdown]
# # Annotation using CITE-seq surface protein markers

# %%
protein_SLEmap_B.obs["leiden_totalVI"] = SLEmap_B.obs["leiden_totalVI"] 
protein_SLEmap_B.obsm["X_umap"] = SLEmap_B.obsm["X_umap"] 
protein_SLEmap_B.uns["leiden_totalVI_colors"] = SLEmap_B.uns["leiden_totalVI_colors"] 

# %%
sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(
    protein_SLEmap_B,
    color=['total_ADTlevels', 'log1p_total_ADTlevels','leiden_totalVI'],
    legend_loc='on data',
    legend_fontsize='xx-small',
    ncols=3,
)

# %%
#list of all CITE-seq markers
protein_SLEmap_B.var.index.tolist()

# %%
sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(protein_SLEmap_B, use_raw=False,color=['anti-human_IgM','anti-human_IgD', 'leiden_totalVI'])

# %%
sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(protein_SLEmap_B, use_raw=False,color=["anti-human_CD185_(CXCR5)", "anti-human_CD124_(IL-4R_)", 
                                                  "anti-human_Ig_light_chain_K", "anti-human_Ig_light_chain_G", 
                                                  "anti-human_CD38", "anti-human_CD11c", "anti-human_CD27"])


# %%
sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(protein_SLEmap_B, use_raw=False,color=['anti-human_CD3','anti-human_CD4','anti-human_CD8','anti-human_CD19','anti-human_CD123','anti-human_CD45RA','anti-human_CD11c','anti-human_CD16','anti-human_CD14','anti-human_CD56','anti-human_CD161','anti-human_CD141_(Thrombomodulin)'])

# %%
##Surface protein levels
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(protein_SLEmap_B, ['anti-human_CD11c'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_B, ['anti-human_CD124_(IL-4R_)'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_B, ['anti-human_Ig_light_chain_K'], groupby='leiden_totalVI', use_raw=False)#CD56Bright NK 
sc.pl.violin(protein_SLEmap_B, ['anti-human_Ig_light_chain_G'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_B, ['anti-human_CD38'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_B, ['anti-human_CD27'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_B, ['anti-human_IgM'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_B, ['anti-human_IgD'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_B, ['anti-human_CD19'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_B, ['anti-human_CD20'], groupby='leiden_totalVI', use_raw=False)

# %%
# sc.settings.set_figure_params(figsize=(4,3))
# marker_genes_dict = {
#     'NK_CD56Bright': ['anti-human_CD56','anti-human_CD62L'],
#     'Myeloid':['anti-human_HLA-DR'],
#     'Dendritic':['anti-human_CD127_(IL-7R_)', 'anti-human_CD11c', 'anti-human_CD141_(Thrombomodulin)'],
#     'cDC':['anti-human_CD1c'],
#     'pDC':['anti-human_CD123','anti-human_CD45RA'],
#     'Monocytes': ['anti-human_CD11b'],
#     'Classical_Monocytes':['anti-human_CD14'],
#     'Non_Classical_Monocytes':['anti-human_CD16'],
# }
# sc.pl.dotplot(protein_SLEmap_B, marker_genes_dict, 'leiden_totalVI', dendrogram=False,use_raw=False)

# %% [markdown]
# # Add new annotation to adata

# %%
old_to_new = {
'0':'Naive',
'1':'Naive',
'2':'Naive',
'3':'Naive',
'4':'Memory',
'5':'Memory',
'6':'Naive',
'7':'Memory',
'8':'Memory',
'9':'Naive',
'10':'Naive',
'11':'Naive',
'12':'Naive',
'13':'Naive',
'14':'ASC',
'15':'Naive'
}
SLEmap_B.obs['manual_annotation_l1'] = (
SLEmap_B.obs['leiden_totalVI']
.map(old_to_new)
.astype('category')
)

# %%
old_to_new = {
'0':'Naive',
'1':'Naive',
'2':'Naive',
'3':'Naive',
'4':'Unswitched memory',
'5':'Switched memory',
'6':'Naive',
'7':'Switched memory',
'8':'ABC',
'9':'Naive',
'10':'Naive',
'11':'Naive',
'12':'Naive',
'13':'Naive',
'14':'ASC',
'15':'Naive'
}
SLEmap_B.obs['manual_annotation_l2'] = (
SLEmap_B.obs['leiden_totalVI']
.map(old_to_new)
.astype('category')
)

# %%
SLEmap_B.obs.groupby("manual_annotation_l1").externalID.value_counts().reset_index().tail(n = 30)

# %%
SLEmap_B.obs.manual_annotation_l1.value_counts()

# %%
SLEmap_B.obs.externalID.nunique()

# %%
SLEmap_B.obs.manual_annotation_l2.value_counts()

# %%
sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(SLEmap_B, use_raw=False,color=['manual_annotation_l1','manual_annotation_l2'])

# %%
SLEmap_B.obs['manual_annotation_l1'].value_counts()

# %%
SLEmap_B.obs['manual_annotation_l2'].value_counts()

# %%
SLEmap_B.write("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/hj10/Results/2_YASCP_final/6_annotation/manualannotation_Other.h5ad")
