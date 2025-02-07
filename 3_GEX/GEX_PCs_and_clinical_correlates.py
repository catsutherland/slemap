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


# %%
figures_folder = "/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/cs54/final_data/Figures/PCA"

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
metadata = pd.read_csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/metadata/SLEmap_metadata_040225.csv")

# %%

# %%
metadata["experiment_id_yascp"] = metadata.experiment_id.str.replace("_L", "-L")
metadata["experiment_id_tensor_qtl"] = metadata.experiment_id.str.rsplit("_", n = 1).str[0]

# %%
sequencing_info = pd.read_csv("/nfs/users/nfs_c/cs54/OTAR2064_cs54/final_data/1_Metadata/sequencing_info_per_slemap_10x_pool.txt", sep = "\t")

# %%
metadata = pd.merge(metadata, sequencing_info, left_on = "poolID", right_on = "PoolID")

# %%
metadata["date_processed"].nunique()

# %%
metadata["Flow_cell"] = metadata["Flow_cell"].astype(str)
metadata["date_processed"] = metadata["date_processed"].replace("07/12/2022", "07/11/2022")

metadata['date_processed_dt'] = pd.to_datetime(metadata.date_processed, format = '%d/%m/%Y')
metadata['date_num'] = (metadata['date_processed_dt'] - metadata['date_processed_dt'].min()).dt.days


# %%
len(metadata['date_num'].unique())

# %% [markdown]
# # PCA for all cells

# %% [markdown]
# # PCA all cells - different normalisation methods

# %%
source_folder = "/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/hj10/Results/3_eQTL_prep/PCA"

file_list = [x for x in os.listdir(f'{source_folder}')]
file_list

# %%
all_pcs_all

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

    all_pcs = pd.merge(metadata[["experiment_id_yascp", "externalID", "date_num",
                                 "date_processed", "poolID", "Flow_cell", "id_pool_lims"]], 
                       all_pcs, right_on = "Sample", left_on = "experiment_id_yascp")

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
    
    for plot_x, plot_y in zip(["Dim.1", "Dim.3"], ["Dim.2", "Dim.4"]):
        label_x = f"PC{plot_x.split('.')[1]}"
        label_y = f"PC{plot_y.split('.')[1]}"

        plt.figure(figsize=(5,5))
        sns.scatterplot(data=all_pcs_all, x=plot_x, y=plot_y)
        plt.xlabel(label_x)
        plt.ylabel(label_y)
        plt.gca().set_aspect('equal')

        filename = f"pca_{norm}_PCS_{plot_x}_{plot_y}.pdf"
        plt.savefig(f"{figures_folder}/{filename}", bbox_inches = "tight");
        
    for hue_variable in ["Recruit_centre", "Site", "Ethnicity_cerner", "SLEDAI_2K", 
                         "Age_at_sampling", "date_processed","date_num",
                         "poolID", "Flow_cell", "id_pool_lims"]:
        plot_x = "Dim.1"
        plot_y = "Dim.2"
        label_x = f"PC{plot_x.split('.')[1]}"
        label_y = f"PC{plot_x.split('.')[1]}"
        plt.figure(figsize=(5,5))
        sns.scatterplot(data=all_pcs_all, x=plot_x, y=plot_y, hue = hue_variable)
        plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title = hue_variable)
        plt.xlabel(label_x)
        plt.ylabel(label_y)
        plt.title(f"{norm}")
        plt.gca().set_aspect('equal')

        filename = f"pca_{norm}_PCS_{plot_x}_{plot_y}_{hue_variable}.pdf"
        plt.savefig(f"{figures_folder}/{filename}", bbox_inches = "tight");
    

    # Extract principal components
    
    pcs_df =  all_pcs_all.iloc[:, 8:18]
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
        plt.xticks(rotation=50, ha='right')
        plt.show();
        # Write heatmap to file
        filename = f"pca_correl_{norm}_{clinvar_type}.pdf"
        fig.savefig(f"{figures_folder}/{filename}", bbox_inches = "tight")
    


# %%
norm_list

# %%
source_folder = "/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/hj10/Results/3_eQTL_prep/PCA"

file_list = [x for x in os.listdir(f'{source_folder}')]
file_list

norm_list = [x.replace("PCA_level0_n281_", "") for x in file_list]
norm_list = [x.replace("_withvdj", "") for x in norm_list]

norm_list = ["meanPS_simplenorm.csv"]

norm_pcs = {}

for norm,file in zip(norm_list, file_list):
    all_pcs = pd.read_csv(f"{source_folder}/{file}")
    
    all_pcs = all_pcs.rename({"Unnamed: 0":"Sample"}, axis = 1)
    
    all_pcs["Sample"] = all_pcs["Sample"].str.replace("_clusterAllcelltypes", "")

    all_pcs["Sample"] = all_pcs["Sample"].str.replace("_L", "-L")

    # Add info

    all_pcs = pd.merge(metadata[["experiment_id_yascp", "externalID", "date_num",
                                 "date_processed", "poolID", "Flow_cell", "id_pool_lims"]], 
                       all_pcs, right_on = "Sample", left_on = "experiment_id_yascp")

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
    
    for plot_x, plot_y in zip(["Dim.1", "Dim.3"], ["Dim.2", "Dim.4"]):
        label_x = f"PC{plot_x.split('.')[1]}"
        label_y = f"PC{plot_y.split('.')[1]}"

        plt.figure(figsize=(5,5))
        sns.scatterplot(data=all_pcs_all, x=plot_x, y=plot_y)
        x_min, x_max = all_pcs_all[plot_x].min(), all_pcs_all[plot_x].max()
        y_min, y_max = all_pcs_all[plot_y].min(), all_pcs_all[plot_y].max()

        # Set the axis limits for all subplots
        plt.xlim(1.2 * min(x_min, y_min), 1.2 * max(x_max, y_max)) 
        plt.ylim(1.2 * min(x_min, y_min), 1.2 * max(x_max, y_max))
        plt.gca().set_aspect('equal')
      
        
        plt.xlabel(label_x)
        plt.ylabel(label_y)
        filename = f"pca_{norm}_PCS_{plot_x}_{plot_y}.pdf"
        plt.savefig(f"{figures_folder}/{filename}", bbox_inches = "tight");
        
    for hue_variable in ["date_processed", "Flow_cell", "id_pool_lims"]:
        plot_x = "Dim.1"
        plot_y = "Dim.2"
        label_x = f"PC{plot_x.split('.')[1]}"
        label_y = f"PC{plot_x.split('.')[1]}"
        plt.figure(figsize=(14,14))
        g = sns.relplot(data=all_pcs_all, x=plot_x, y=plot_y, hue=hue_variable, 
                        col=hue_variable, col_wrap = 6, s=100)
        x_min, x_max = all_pcs_all[plot_x].min(), all_pcs_all[plot_x].max()
        y_min, y_max = all_pcs_all[plot_y].min(), all_pcs_all[plot_y].max()

        for ax in g.axes.flat:
            ax.set_aspect('equal')
            ax.set_xlim(1.2 * min(x_min, y_min), 1.2 * max(x_max, y_max))
            ax.set_ylim(1.2 * min(x_min, y_min), 1.2 * max(x_max, y_max))

#         plt.xlabel(label_x)
#         plt.ylabel(label_y)
#         plt.title(f"{norm}")
        filename_1 = f"pca_{norm}_PCS_{plot_x}_{plot_y}_{hue_variable}_split.pdf"
        filename_2 = f"pca_{norm}_PCS_{plot_x}_{plot_y}_{hue_variable}_split.png"
        plt.savefig(f"{figures_folder}/{filename_1}", bbox_inches = "tight")
        plt.savefig(f"{figures_folder}/{filename_2}", bbox_inches = "tight");
        
                
    for hue_variable in ["poolID"]:
        plot_x = "Dim.1"
        plot_y = "Dim.2"
        label_x = f"PC{plot_x.split('.')[1]}"
        label_y = f"PC{plot_x.split('.')[1]}"
        plt.figure(figsize=(7,7))
        g = sns.relplot(data=all_pcs_all, x=plot_x, y=plot_y, hue=hue_variable, 
                        col=hue_variable, col_wrap = 6, s= 200)
        x_min, x_max = all_pcs_all[plot_x].min(), all_pcs_all[plot_x].max()
        y_min, y_max = all_pcs_all[plot_y].min(), all_pcs_all[plot_y].max()

        for ax in g.axes.flat:
            ax.set_aspect('equal')
            ax.set_xlim(1.2 * min(x_min, y_min), 1.2 * max(x_max, y_max))
            ax.set_ylim(1.2 * min(x_min, y_min), 1.2 * max(x_max, y_max))

#         plt.xlabel(label_x)
#         plt.ylabel(label_y)
#         plt.title(f"{norm}")
        filename_1 = f"pca_{norm}_PCS_{plot_x}_{plot_y}_{hue_variable}_split.pdf"
        filename_2 = f"pca_{norm}_PCS_{plot_x}_{plot_y}_{hue_variable}_split.png"
        plt.savefig(f"{figures_folder}/{filename_1}", bbox_inches = "tight")
        plt.savefig(f"{figures_folder}/{filename_2}", bbox_inches = "tight");



# %%
g = sns.relplot(data=all_pcs_all, x=plot_x, y=plot_y, hue=hue_variable, col=hue_variable, col_wrap = 4)
for ax in g.axes.flat:
    ax.set_aspect('equal')

# Display the plot
plt.show()


# %% [markdown]
# # Plots per cell type (tensorQTL PCs)

# %%
storage_biologics = pd.read_csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/metadata/old/storage_biologics.csv")

# %%
storage_biologics

# %%
source_folder = "/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/mo11/8.eqtl_run/TensorQTL_eQTLS__1mbTSS"

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
    
    tensor_pcs_all = pd.merge(tensor_pcs_all, storage_biologics)


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
    
    for plot_x, plot_y in zip(["Phenotype PC1", "Phenotype PC3"], 
                              ["Phenotype PC2", "Phenotype PC4"]):
        label_x = f"{plot_x.split(' ')[1]}"
        label_y = f"{plot_y.split(' ')[1]}"

        plt.figure(figsize=(5,5))
        sns.scatterplot(data=tensor_pcs_all, x=plot_x, y=plot_y)
        plt.xlabel(label_x)
        plt.ylabel(label_y)
        filename = f"pca_tensorqtl_{celltype}_PCS_{label_x}_{label_y}.pdf"
        plt.savefig(f"{figures_folder}/{filename}", bbox_inches = "tight");
        
        
    for hue_variable in ["Recruit_centre", "Site", "Ethnicity_cerner", "SLEDAI_2K", "storage_biologics"]:
        plot_x = "Phenotype PC1"
        plot_y = "Phenotype PC2"
        label_x = f"{plot_x.split(' ')[1]}"
        label_y = f"{plot_y.split(' ')[1]}"
        plt.figure(figsize=(5,5))
        sns.scatterplot(data=tensor_pcs_all, x=plot_x, y=plot_y, hue = hue_variable)
        plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title = hue_variable)
        plt.xlabel(label_x)
        plt.ylabel(label_y)
        filename = f"pca_tensorqtl_{celltype}_PCS_1_2_{hue_variable}.pdf"
        plt.savefig(f"{figures_folder}/{filename}", bbox_inches = "tight");

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
        plt.xticks(rotation=50, ha='right')
        plt.show();
        # Write heatmap to file
        filename = f"pca_correl_tensorqtl_{celltype}_{clinvar_type}.pdf"
        fig.savefig(f"{figures_folder}/{filename}", bbox_inches = "tight")
    
    
    
    
    

# %%
