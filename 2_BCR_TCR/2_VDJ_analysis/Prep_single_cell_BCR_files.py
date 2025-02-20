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

# %% [markdown]
# # Import packages and setup

# %%
# %matplotlib inline
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 400)


# %% [markdown]
# # Define custom functions

# %%
def mask_d_germline(germline,v_germline_end, d_sequence_end, d_sequence_start, np1_length, chain):
    if chain == 'Light': #no D gene to mask in the light chains
        return germline
    if d_sequence_start != d_sequence_start:
        return np.nan
    elif d_sequence_end != d_sequence_end:
        return np.nan
    else:
        d_length = int(d_sequence_end - d_sequence_start + 1)
        d_mask = "N" * d_length
        mask_germline = germline[:v_germline_end + np1_length] + d_mask + germline[v_germline_end + np1_length + d_length:]
        return mask_germline

def get_mut_count_v(sequence, germline, v_germline_end):
    mutations = 0
    for i, base in enumerate(sequence[:v_germline_end]):
        if base != "." and base != "-" and base != "N" and germline[i] != "N" and germline[i] != "-" and germline[i] !=  "." :
            if base != germline[i]:
                mutations += 1
    return mutations


def get_mut_freq(sequence, germline_d_mask):
    if germline_d_mask != germline_d_mask:
        return np.nan
    if len(sequence) != len(germline_d_mask):
        return(np.nan)
    mut_count = 0
    length = 0
    for i, base in enumerate(sequence):
        if base != "." and base != "-" and base != "N" and germline_d_mask[i] != "N" and germline_d_mask[i] != "-" and germline_d_mask[i] !=  "." :
            length += 1
            if base != germline_d_mask[i]:
                mut_count += 1
    return(mut_count / length)
  
def get_mut_freq_v(sequence, germline, v_germline_end):
    mutations = 0
    length = 0
    for i, base in enumerate(sequence[:v_germline_end]):
        if base != "." and base != "-" and base != "N" and germline[i] != "N" and germline[i] != "-" and germline[i] !=  "." :
            length += 1
            if base != germline[i]:
                mutations += 1
    if length == 0:
        return np.nan
    else:
        return mutations/length
    

