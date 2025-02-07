# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: scanpy
#     language: python
#     name: scanpy
# ---

# %%
import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42



# %%
freeze = sc.read("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/hj10/Results/2_YASCP_final/5_QC_afterHLA_afterClusterremoval/DataFreeze_tmp.h5ad")

# %%
freeze

# %%
annot_info = pd.read_csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/cs54/final_data/2_BCR_TCR/1_IgBLAST_assignment/has_tcr_bcr_barcodes.csv")

# %%
annot_info

# %%
obs_tmp = freeze.obs.copy()

# %%
obs_tmp["cell_barcode"] = obs_tmp.index.str.split("-SLE").str[0]

# %%
obs_tmp = obs_tmp.reset_index().merge(how = 'left', right = annot_info[["barcode", "pool", "VDJ"]], left_on = ["cell_barcode", "poolID"], right_on = ["barcode", "pool"]).set_index('index')

# %%
obs_tmp = obs_tmp.drop(["cell_barcode", "barcode", "pool"], axis = 1)

# %%
freeze.obs = obs_tmp

# %%
freeze.write("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/hj10/Results/2_YASCP_final/5_QC_afterHLA_afterClusterremoval/DataFreeze_tmp_vdj.h5ad")

# %%
freeze

# %% [markdown]
# # Plot UMAP

# %%
sc.set_figure_params(figsize=(3,3), dpi=500, dpi_save=500)
sc.settings.vector_friendly = True
sc.settings.figdir = "/nfs/users/nfs_c/cs54/OTAR2064_cs54/final_data/Figures/umaps"

ax = sc.pl.umap(freeze, 
           color=["Azimuth:predicted.celltype.l1"], 
           legend_fontsize='xx-small',
           frameon=False,title=[''],
#            palette = ["red","#709AE1FF"],
           size=0.1,
           save = "_azimuth_l1")

# %%
sc.set_figure_params(figsize=(3,3), dpi=500, dpi_save=500)
sc.settings.vector_friendly = True
sc.settings.figdir = "/nfs/users/nfs_c/cs54/OTAR2064_cs54/final_data/Figures/umaps"

ax = sc.pl.umap(freeze, 
           color=["VDJ"], 
           legend_fontsize='xx-small',
           frameon=False,title=[''],
           palette = ["red","#709AE1FF"],
           size=0.1,
           save = "_vdj_info_freeze")

# %%
