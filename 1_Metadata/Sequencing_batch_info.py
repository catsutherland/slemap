# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: general-venv
#     language: python
#     name: general-venv
# ---

# %%
import pandas as pd
import numpy as np

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 200)

# %%
metadata

# %%
sequencing_info = pd.read_csv("/nfs/users/nfs_c/cs54/OTAR2064_cs54/final_data/1_Metadata/slemap_wh_22_01_2025.txt", sep = "\t")

# %%
sequencing_info["PoolID"] = sequencing_info["supplier_name"]

# %%
sequencing_info["Flow_cell"] = sequencing_info["id_run"]

# %%
sequencing_info["Lane"] = sequencing_info["id_run"].astype(str) + "_" + sequencing_info["position"].astype(str)

# %%
len(sequencing_info)

# %%
sequencing_info.PoolID.value_counts()

# %%
sequencing_info.groupby("id_pool_lims").id_run.value_counts()

# %%
sequencing_info.groupby("id_pool_lims").id_run.value_counts()

# %%
sequencing_info.query("id_pool_lims == 'NT1815119G'")

# %%
sequencing_info.query("id_pool_lims == 'NT1815119G'").supplier_name.unique()

# %%
sequencing_info.query("id_pool_lims == 'NT1815119G'").supplier_name.unique()

# %%
sequencing_info.groupby(["PoolID", "Lane"]).tag_index.nunique().reset_index()

# %%
sequencing_info.query("PoolID == 'LDP69'")


# %%
sequencing_info[["PoolID", "Flow_cell", "id_pool_lims"]].to_csv("/nfs/users/nfs_c/cs54/OTAR2064_cs54/final_data/1_Metadata/sequencing_info.txt", index = False, sep = "\t")

# %%
metadata = pd.read_csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/metadata/SLEmap_metadata_191224.csv")

# %%
slemap_pools = metadata.poolID.unique()

# %%
sequencing_info_per_10x_pool = sequencing_info.query("PoolID in @slemap_pools")


# %%
sequencing_info_per_10x_pool.groupby("id_pool_lims").id_run.value_counts()

# %%
sequencing_info_per_10x_pool.groupby("id_pool_lims").Lane.value_counts()

# %%
sequencing_info_per_10x_pool.groupby("id_pool_lims").supplier_name.nunique()

# %%
sequencing_info_per_10x_pool.query("id_pool_lims == 'NT1815119G'")

# %%
sequencing_info_per_10x_pool.groupby(["id_pool_lims", "Lane"]).tag_index.nunique().reset_index()

# %%
sequencing_info_per_10x_pool.groupby(["PoolID", "Lane"]).tag_index.nunique().reset_index()

# %%
sequencing_info_per_10x_pool = sequencing_info_per_10x_pool.drop_duplicates(subset = "PoolID")


# %%
sequencing_info_per_10x_pool

# %%
sequencing_info_per_10x_pool.groupby(["PoolID", "Lane"]).tag_index.nunique().reset_index()

# %%
sequencing_info_per_10x_pool[["PoolID", "Flow_cell", "id_pool_lims"]].to_csv("/nfs/users/nfs_c/cs54/OTAR2064_cs54/final_data/1_Metadata/sequencing_info_per_slemap_10x_pool.txt", index = False, sep = "\t")

# %%
metadata = pd.merge(metadata, sequencing_info[["PoolID", "Flow_cell", "id_pool_lims"]], left_on = "poolID", right_on = "PoolID").drop_duplicates()

# %%
metadata.groupby("scSeq_batch").id_pool_lims.unique()

# %%

# %%
used_pools = metadata.id_pool_lims.unique()

# %%
sequencing_info.query("id_pool_lims in @used_pools").groupby("id_pool_lims").Lane.unique()


# %%
sequencing_info.query("id_pool_lims in @used_pools").groupby("PoolID").id_pool_lims.unique()

# %%
sequencing_info.groupby("PoolID").id_pool_lims.unique()

# %%
metadata.groupby("scSeq_batch").Flow_cell.unique()

# %% [markdown]
#
