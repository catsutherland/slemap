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

# %%
metadata = pd.read_csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/metadata/SLEmap_metadata_191224.csv")

# %% [markdown]
# # Correct date processed for EQTL118P

# %%
metadata["date_processed"] = metadata["date_processed"].replace("07/12/2022", "07/11/2022")

# %%
len(metadata)

# %% [markdown]
# # Add sequencing info
#
#

# %%
sequencing_info = sequencing_info = pd.read_csv("/nfs/users/nfs_c/cs54/OTAR2064_cs54/final_data/1_Metadata/sequencing_info_per_slemap_10x_pool.txt", sep = "\t")
sequencing_info["poolID"] = sequencing_info["PoolID"]


# %%
sequencing_info.poolID

# %%
metadata = pd.merge(metadata, sequencing_info[["poolID", "Flow_cell", "id_pool_lims"]])

# %%
len(metadata)

# %%
metadata.to_csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/metadata/SLEmap_metadata_040225.csv", index = False)

# %%
metadata
