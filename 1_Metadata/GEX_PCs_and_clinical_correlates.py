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
from scipy.stats import false_discovery_control, pearsonr, spearmanr


# %% [markdown]
# # Clinical Data

# %%
clinical_all = pd.read_csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/phenotype/combined_data_tidied_sc_used_25_11_24.csv")

# %%
clinical_all.head()

# %%
continuous_list = ["Age_at_sampling", "Years_from_diagnosis", "Albumin", "Creatinine", "eGFR", "uPCR", "SLEDAI_2K",
                'SLEDAI_2KG','cSLEDAI_2K','Haemoglobin']

binary_list = ['Recruit_centre','Lupus_status','IFN_status','LN',
                       'Thyroid_Disease_Hashimotos','Graves_Disease','Type_I_Diabetes',
                       'Rheumatoid_Arthritis','MMF','Pred','HCQ','AZA','MTX','TAC','ASP','CYC_ever',
                       'BEL_ever','RTX_ever','OFA_ever','dsDNA','dsDNA_Ab_ever','ANA_ever','C3_status',
                       'C4_status', "Lymphopenia_ever", "Lymphopenia_now"]

# %%
replace_dict = {"Yes":1, "No":0, "ACTIVE":1, "INACTIVE":0, "IFN-high":1, "IFN-low":0, "IFN-neg":0, 
                "Awaited":np.nan, "awaited":np.nan, "Positive":1, "Negative":0, "NORMAL":0, "LOW":1,
                "Normal":0, "Low":1,"Kings":0, "Imperial":1, "YEs":1}

# %% [markdown]
# # Load metadata to get sample names

# %%
metadata = pd.read_csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/metadata/SLEmap_metadata_191224.csv")

# %%
metadata["experiment_id_yascp"] = metadata.experiment_id.str.replace("_L", "-L")
metadata["experiment_id_tensor_qtl"] = metadata.experiment_id.str.rsplit("_", n = 1).str[0]

# %% [markdown]
# # PCA for all cells

# %%
PCA_all = pd.read_csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/hj10/Results/2_YASCP_final/5_QC_afterHLA_afterClusterremoval/PCA_level0_n281_withvdj.csv")

# %%
len(PCA_all)

# %%
PCA_all = PCA_all.rename({"Unnamed: 0":"Sample"}, axis = 1)

# %% [markdown]
# # Add info

# %%
PCA_all = pd.merge(metadata[["experiment_id_yascp", "externalID"]], PCA_all, right_on = "Sample", left_on = "experiment_id_yascp")

# %%
PCA_all = pd.merge(PCA_all, clinical_all)

# %%
PCA_all.head()

# %% [markdown]
# # Plot

# %%
plt.figure(figsize=(5,5))
sns.scatterplot(data=PCA_all, x="Dim.1", y="Dim.2")
plt.xlabel("PC1")
plt.ylabel("PC2");

# %%
plt.figure(figsize=(5,5))
sns.scatterplot(data=PCA_all, x="Dim.1", y="Dim.2", hue = "Recruit_centre")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5));

# %%
plot_x = "Dim.1"
plot_y = "Dim.2"
label_x = f"PC{plot_x.split('.')[1]}"
label_y = f"PC{plot_y.split('.')[1]}"
hue_value = "LN"

plt.figure(figsize=(5,5))
sns.scatterplot(data=PCA_all, x=plot_x, y=plot_y, hue = hue_value)
plt.xlabel(label_x)
plt.ylabel(label_y)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5));

# %%
plot_x = "Dim.1"
plot_y = "Dim.2"
label_x = f"PC{plot_x.split('.')[1]}"
label_y = f"PC{plot_y.split('.')[1]}"
hue_value = "Site"

plt.figure(figsize=(5,5))
sns.scatterplot(data=PCA_all, x=plot_x, y=plot_y, hue = hue_value)
plt.xlabel(label_x)
plt.ylabel(label_y)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5));

# %%
plot_x = "Dim.4"
plot_y = "Dim.5"
label_x = f"PC{plot_x.split('.')[1]}"
label_y = f"PC{plot_y.split('.')[1]}"
hue_value = "LN"

plt.figure(figsize=(5,5))
sns.scatterplot(data=PCA_all, x=plot_x, y=plot_y, hue = hue_value)
plt.xlabel(label_x)
plt.ylabel(label_y)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5));

# %%
plot_x = "Dim.4"
plot_y = "Dim.5"
label_x = f"PC{plot_x.split('.')[1]}"
label_y = f"PC{plot_y.split('.')[1]}"
hue_value = "Age_at_sampling"

plt.figure(figsize=(5,5))
sns.scatterplot(data=PCA_all, x=plot_x, y=plot_y, hue = hue_value)
plt.xlabel(label_x)
plt.ylabel(label_y)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5));

# %%
binary_variables = PCA_all.copy()
binary_variables = binary_variables[binary_list]

binary_variables = binary_variables.map(lambda x: replace_dict.get(x, x))

clinvar_data = pd.concat([PCA_all[continuous_list], binary_variables], axis = 1)

clinvars = {"continuous":continuous_list, "categorical":binary_list }

pcs = PCA_all.iloc[:, 3:13].to_numpy()

signif_threshold = 0.05

# Loop through each clinical variable type (categorical or continuous)
# Construct correlation matrix between PCs and clinical variables.
for clinvar_type, clinvar_list in clinvars.items():
    corr_x = pcs.shape[1]
    corr_y = len(clinvar_list)
    corr_matrix = np.zeros((corr_x, corr_y))
    pvalue_matrix = np.zeros((corr_x, corr_y))
    for pc_idx in range(corr_x):
        for clinvar_idx, clinvar in enumerate(clinvar_list):
            x = clinvar_data[clinvar].to_numpy()
            y = pcs[:, pc_idx]
            # Use Spearman's R if categorical variable
            if clinvar_type == "categorical":
                # Compute correlation
                s = spearmanr(x, y, nan_policy="omit")
            # Use Pearson's R if continuous variable
            if clinvar_type == "continuous":
                # Ignore NaNs
                nans = np.logical_or(np.isnan(x), np.isnan(y))
                s = pearsonr(x[~nans], y[~nans])
            r = s.statistic
            p = s.pvalue
            # Store correlation and p-value in matrices
            corr_matrix[pc_idx, clinvar_idx] = r
            pvalue_matrix[pc_idx, clinvar_idx] = p
    # Adjust p-values using Benjamini-Yekutieli because some variables
    # may be correlated to each other.
    adjpvalue_matrix = false_discovery_control(pvalue_matrix, method="by")
    # Construct array to use as annotation to indicate which correlations
    # are significant
    signif_matrix = adjpvalue_matrix < signif_threshold
    # U+2022 is a bullet.  It is centred, so indicates cells in a clearer way
    # than asterisk, especially if the font is Helvetica.
    signif_matrix = np.where(signif_matrix, "\u2022", " ")
    # Draw heatmap to represent correlation matrix
    fig, ax = plt.subplots(figsize=(7, 7))
    sns.heatmap(
        corr_matrix,
        annot=signif_matrix,
        fmt="",
        cmap="vlag",
        center=0,
        xticklabels=clinvar_list,
        # Off-by-one: first PC is PC1, not PC0
        yticklabels=list(range(1, corr_x + 1)),
        ax=ax,
    )
    ax.set_xlabel("Clinical variable")
    ax.set_ylabel("Principal component");
    # ax.set_title(f"{txtype}, melioid")
    # Write heatmap to file
    # filename = "pca_correl_" + txtype + "_" + clinvar_type + ".pdf"
    # filepath = "../reports/" + filename
    # fig.savefig(filepath)


# %% [markdown]
# # Different normalisation methods

# %%
source_folder = "/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/hj10/Results/3_eQTL_prep/PCA"

file_list = [x for x in os.listdir(f'{source_folder}')]
file_list

# %%
all_pcs_all

# %%
all_pcs = pd.read_csv(f"{source_folder}/{file}")
    
all_pcs = all_pcs.rename({"Unnamed: 0":"Sample"}, axis = 1)

all_pcs["Sample"] = all_pcs["Sample"].str.replace("_clusterAllcelltypes", "")

all_pcs["Sample"] = all_pcs["Sample"].str.replace("_L", "-L")

all_pcs


# %%
all_pcs

# %%
metadata[["experiment_id_yascp", "externalID"]]

# %%
source_folder = "/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/hj10/Results/3_eQTL_prep/PCA"

file_list = [x for x in os.listdir(f'{source_folder}')]
file_list

norm_list = [x.replace("PCA_level0_n281_", "") for x in file_list]
norm_list = [x.replace("_withvdj", "") for x in norm_list]

norm_pcs = {}

for norm,file in zip(norm_list, file_list):
    all_pcs = pd.read_csv(f"{source_folder}/{file}")
    
    all_pcs = all_pcs.rename({"Unnamed: 0":"Sample"}, axis = 1)
    
    all_pcs["Sample"] = all_pcs["Sample"].str.replace("_clusterAllcelltypes", "")

    all_pcs["Sample"] = all_pcs["Sample"].str.replace("_L", "-L")

    # Add info

    all_pcs = pd.merge(metadata[["experiment_id_yascp", "externalID"]], all_pcs, right_on = "Sample", left_on = "experiment_id_yascp")

    all_pcs_all = pd.merge(all_pcs, clinical_all)

    binary_variables = all_pcs_all.copy()
    binary_variables = binary_variables[binary_list]

    binary_variables = binary_variables.map(lambda x: replace_dict.get(x, x))

    # Remove categorical variables where all non-NaN values are the same

    binary_variables = binary_variables.loc[:, binary_variables.nunique(dropna=True) > 1]
    
    temp_binary_list = binary_variables.columns
    
    # Concatenate continuous variables with the cleaned binary variables

    clinvar_data = pd.concat([all_pcs_all[continuous_list], binary_variables], axis = 1)

    # Dictionary of continuous and categorical variables

    clinvars = {"continuous":continuous_list, "categorical":temp_binary_list }

    # Extract principal components
    
    pcs_df =  all_pcs_all.iloc[:, 7:17]
    pcs_df = pcs_df.astype(float)
    pcs = pcs_df.to_numpy()

    signif_threshold = 0.05

    # Loop through each clinical variable type (categorical or continuous)
    # Construct correlation matrix between PCs and clinical variables.
    for clinvar_type, clinvar_list in clinvars.items():
        corr_x = pcs.shape[1]
        corr_y = len(clinvar_list)
        corr_matrix = np.zeros((corr_x, corr_y))
        pvalue_matrix = np.zeros((corr_x, corr_y))
        for pc_idx in range(corr_x):
            for clinvar_idx, clinvar in enumerate(clinvar_list):
                x = clinvar_data[clinvar].to_numpy()
                y = pcs[:, pc_idx]
                # Use Spearman's R if categorical variable
                if clinvar_type == "categorical":
                    # Compute correlation
                    s = spearmanr(x, y, nan_policy="omit")
                # Use Pearson's R if continuous variable
                if clinvar_type == "continuous":
                    # Ignore NaNs
                    nans = np.logical_or(np.isnan(x), np.isnan(y))
                    s = pearsonr(x[~nans], y[~nans])
                r = s.statistic
                p = s.pvalue
                # Store correlation and p-value in matrices
                corr_matrix[pc_idx, clinvar_idx] = r
                pvalue_matrix[pc_idx, clinvar_idx] = p
        # Adjust p-values using Benjamini-Yekutieli because some variables
        # may be correlated to each other.
        adjpvalue_matrix = false_discovery_control(pvalue_matrix, method="by")
        # Construct array to use as annotation to indicate which correlations
        # are significant
        signif_matrix = adjpvalue_matrix < signif_threshold
        # U+2022 is a bullet.  It is centred, so indicates cells in a clearer way
        # than asterisk, especially if the font is Helvetica.
        signif_matrix = np.where(signif_matrix, "\u2022", " ")
        # Draw heatmap to represent correlation matrix
        fig, ax = plt.subplots(figsize=(7, 7))
        sns.heatmap(
            corr_matrix,
            annot=signif_matrix,
            fmt="",
            cmap="vlag",
            center=0,
            xticklabels=clinvar_list,
            # Off-by-one: first PC is PC1, not PC0
            yticklabels=list(range(1, corr_x + 1)),
            ax=ax,
        )
        ax.set_xlabel("Clinical variable")
        ax.set_ylabel("Principal component");
        ax.set_title(f"{norm}")
        plt.show();
        # Write heatmap to file
        # filename = "pca_correl_" + txtype + "_" + clinvar_type + ".pdf"
        # filepath = "../reports/" + filename
        # fig.savefig(filepath)

    
    
    
    
    

# %% [markdown]
# # Plots per cell type (tensorQTL PCs)

# %%
source_folder = "/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/mo11/8.eqtl_run/TensorQTL_eQTLS"

folder_list = [x for x in os.listdir(f'{source_folder}')]
folder_list

celltype_list = [x.replace("dMean__", "") for x in folder_list]

celltype_pcs = {}

for celltype,folder in zip(celltype_list, folder_list):
    file_path = f"{source_folder}/{folder}/20pcs/base_output/base/Covariates.tsv"
    # Some celltypes don't have 20 PCs (need to check why)
    if os.path.exists(file_path):
        tensor_pcs = pd.read_csv(file_path, sep = "\t")
    else:
        continue
    tensor_pcs = tensor_pcs.T

    tensor_pcs.columns = tensor_pcs.iloc[0]

    # Drop the first row since it's now been set as the column names
    tensor_pcs = tensor_pcs.drop("Unnamed: 0")

    tensor_pcs = tensor_pcs.reset_index().rename({"index":"Sample"}, axis = 1)

    tensor_pcs = pd.merge(metadata[["experiment_id_tensor_qtl", "externalID"]], tensor_pcs, right_on = "Sample", left_on = "experiment_id_tensor_qtl")

    tensor_pcs_all = pd.merge(tensor_pcs, clinical_all)

    binary_variables = tensor_pcs_all.copy()
    binary_variables = binary_variables[binary_list]

    binary_variables = binary_variables.map(lambda x: replace_dict.get(x, x))

    # Remove categorical variables where all non-NaN values are the same

    binary_variables = binary_variables.loc[:, binary_variables.nunique(dropna=True) > 1]
    
    temp_binary_list = binary_variables.columns
    
    # Concatenate continuous variables with the cleaned binary variables

    clinvar_data = pd.concat([tensor_pcs_all[continuous_list], binary_variables], axis = 1)

    # Dictionary of continuous and categorical variables

    clinvars = {"continuous":continuous_list, "categorical":temp_binary_list }

    # Extract principal components
    
    pcs_df =  tensor_pcs_all.iloc[:, 7:17]
    pcs_df = pcs_df.astype(float)
    pcs = pcs_df.to_numpy()

    signif_threshold = 0.05

    # Loop through each clinical variable type (categorical or continuous)
    # Construct correlation matrix between PCs and clinical variables.
    for clinvar_type, clinvar_list in clinvars.items():
        corr_x = pcs.shape[1]
        corr_y = len(clinvar_list)
        corr_matrix = np.zeros((corr_x, corr_y))
        pvalue_matrix = np.zeros((corr_x, corr_y))
        for pc_idx in range(corr_x):
            for clinvar_idx, clinvar in enumerate(clinvar_list):
                x = clinvar_data[clinvar].to_numpy()
                y = pcs[:, pc_idx]
                # Use Spearman's R if categorical variable
                if clinvar_type == "categorical":
                    # Compute correlation
                    s = spearmanr(x, y, nan_policy="omit")
                # Use Pearson's R if continuous variable
                if clinvar_type == "continuous":
                    # Ignore NaNs
                    nans = np.logical_or(np.isnan(x), np.isnan(y))
                    s = pearsonr(x[~nans], y[~nans])
                r = s.statistic
                p = s.pvalue
                # Store correlation and p-value in matrices
                corr_matrix[pc_idx, clinvar_idx] = r
                pvalue_matrix[pc_idx, clinvar_idx] = p
        # Adjust p-values using Benjamini-Yekutieli because some variables
        # may be correlated to each other.
        adjpvalue_matrix = false_discovery_control(pvalue_matrix, method="by")
        # Construct array to use as annotation to indicate which correlations
        # are significant
        signif_matrix = adjpvalue_matrix < signif_threshold
        # U+2022 is a bullet.  It is centred, so indicates cells in a clearer way
        # than asterisk, especially if the font is Helvetica.
        signif_matrix = np.where(signif_matrix, "\u2022", " ")
        # Draw heatmap to represent correlation matrix
        fig, ax = plt.subplots(figsize=(7, 7))
        sns.heatmap(
            corr_matrix,
            annot=signif_matrix,
            fmt="",
            cmap="vlag",
            center=0,
            xticklabels=clinvar_list,
            # Off-by-one: first PC is PC1, not PC0
            yticklabels=list(range(1, corr_x + 1)),
            ax=ax,
        )
        ax.set_xlabel("Clinical variable")
        ax.set_ylabel("Principal component");
        ax.set_title(f"{celltype}")
        plt.show();
        # Write heatmap to file
        # filename = "pca_correl_" + txtype + "_" + clinvar_type + ".pdf"
        # filepath = "../reports/" + filename
        # fig.savefig(filepath)

    
    
    
    
    

# %%
