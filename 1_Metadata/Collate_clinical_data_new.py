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
# # Set-up

# %%
# %matplotlib inline
import pandas as pd
import numpy as np

pd.set_option('display.max_columns', None)


# %% [markdown]
# # Guy's data - first batch

# %%
guys = pd.read_csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/phenotype/raw/08_03_24 GUYS HOPITAL-SANGER PHENOTYPE.csv", header = 1)

guys = guys.drop("Patient Information", axis =1)

guys = guys.T

guys.columns = guys.iloc[0]

guys = guys.dropna(axis = 0, how = 'all')

guys = guys.drop(index = ["Patient_ID", "Unnamed: 2", "L1"])

guys = guys.reset_index(names = "Patient_ID")

guys["Site"] = "Guys"

guys["Recruit_centre"] = "Kings"

drop_patients = ["EQTL0166p", "EQTL0180p", "EQTL0123p", "EQTL0110p", "EQTL0105p"] 

# 0166p and 0180p duplicated in other spreadsheet, 
# 123p, 110p and 034p are multiple samples from same individ, keep 034p (because 
# sample collection date in metadata matches clinical info for this sample)
# 105p and 134p are dups, keep 134p

guys = guys[~guys["Patient_ID"].isin(drop_patients)]

# EQTL0156p also in guys_new, need to check because some values differ 
# (but guys_new appears to have more appropriate test dates)
guys = guys.query("Patient_ID != 'EQTL0156p'") 

# %% [markdown]
# # Guy's data - second batch

# %%
guys_new = pd.read_csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/phenotype/raw/14_03_24 GUYS HOPITAL-SANGER PHENOTYPE TEMPLATE UPDATED.csv", header = 1)

guys_new = guys_new.drop("Patient Information", axis =1)

guys_new = guys_new.T

guys_new.columns = guys_new.iloc[0]

guys_new = guys_new.dropna(axis = 0, how = 'all')

guys_new = guys_new.drop(index = ["Patient_ID", "Unnamed: 2", "L1"])

guys_new = guys_new.reset_index(names = "Patient_ID")

guys_new["Site"] = "Guys"

guys_new["Recruit_centre"] = "Kings"


# %% [markdown]
# # King's data

# %%
kings = pd.read_csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/phenotype/raw/21_02_2024 KINGS HOPITAL-SANGER PHENOTYPE.csv", header = 1)

kings = kings.drop("Patient Information", axis =1)

kings = kings.T

kings.columns = kings.iloc[0]

kings = kings.dropna(axis = 0, how = 'all')

kings = kings.drop(index = ["Patient_ID", "Unnamed: 2", "L1"])

kings = kings.reset_index(names = "Patient_ID")

kings["Site"] = "Kings"

kings["Recruit_centre"] = "Kings"

# %% [markdown]
# # Combine Guy's and King's

# %%
kings_both = pd.concat([guys, guys_new, kings])
kings_both.SLEDAI_2K = kings_both.SLEDAI_2K.astype("Int64")
kings_both = kings_both.rename({"C4 _level":"C4_level"}, axis = 1)
kings_both = kings_both.rename({"Year of Birth":"Year_of_Birth"}, axis = 1)


# %%
kings_both = kings_both.replace("Unspecified", np.nan).infer_objects(copy=False)


# %%
kings_both = kings_both.rename({"Date of sampling":"Sample_date"}, axis = 1)
kings_both["Year_of_Birth"] = kings_both["Year_of_Birth"].astype(int)
kings_both["Sample_year"] = kings_both["Sample_date"].str.split("/").str[2]
kings_both["Sample_year"] = kings_both["Sample_year"].astype(int)
kings_both["Age_at_sampling"] = kings_both["Sample_year"] - kings_both["Year_of_Birth"]

kings_both["Diagnosis_year"] = kings_both["Diagnosis_date"]
kings_both["Diagnosis_year"] = kings_both["Diagnosis_year"].replace({"Sep-19":"2019", "Aug-21":"2021",
                                                                       "Jan-23":"2023"})
kings_both["Diagnosis_year"] = kings_both["Diagnosis_year"].astype(float)


# %% [markdown]
# # Output transposed King's

# %%
kings_both_t = kings_both.T 

# %%
kings_both_t.columns = kings_both_t.iloc[0]
kings_both_t = kings_both_t[1:]

# %%
kings_both_t.to_csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/phenotype/raw/Kings_Guys_combined.csv", na_rep='NA')

# %% [markdown]
# # Imperial data

# %%
# imperial = pd.read_excel("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/phenotype/raw/OTAR_Phenotyping_ICL_300424.xlsx")

# imperial.to_csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/phenotype/raw/OTAR_Phenotyping_ICL_300424.csv", index = False)
#L116 and S79 are dups, keep S79

# %%
imperial = pd.read_csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/phenotype/raw/OTAR_Phenotyping_ICL_300424.csv")

# %%
imperial["Site"] = "Imperial"

imperial["Recruit_centre"] = "Imperial"

imperial = imperial.rename({"Patient_N":"Patient_ID"}, axis = 1)

imperial["Patient_ID"] = imperial["Patient_ID"].str.replace(" ", "")

# %%
imperial = imperial.replace({"Not Done":np.nan, "ND":np.nan, "NaT":np.nan, "NK":np.nan}).infer_objects(copy=False)


# %%
imperial["Sample_year"] = imperial["Sample_date"].str.split("-").str[0]
imperial["Sample_year"] = imperial["Sample_year"].astype(int)
imperial["Diagnosis_year"] = imperial["Diagnosis_date"].str.split("-").str[0]
imperial["Diagnosis_year"] = imperial["Diagnosis_year"].astype(float)
imperial["Year_of_Birth"] = imperial["Sample_year"] - imperial["Age_at_sampling"]

# %%
imperial

# %% [markdown]
# # Combine Imperial and Guy's/King's

# %%
both_sites = pd.concat([kings_both, imperial], join = "outer").copy()

# %%
both_sites = both_sites.rename({" B_cell_% ":"B_cell_%"}, axis = 1)

both_sites.columns = both_sites.columns.str.replace(" ", "_")


# %%
both_sites.Recruit_centre.unique()

# %%
both_sites["Years_from_diagnosis"] = both_sites["Sample_year"] - both_sites["Diagnosis_year"]

# %%
both_sites["Age_at_diagnosis"] = both_sites["Diagnosis_year"] - both_sites["Year_of_Birth"] 

# %% [markdown]
# # Tidy up IDs to match metadata 
# Some metadata sample IDs have a trailing "r" which the clinical data doesn't have

# %%
both_sites.Patient_ID = both_sites.Patient_ID.str.replace("p", "P")

patient_ids = both_sites["Patient_ID"]

tidy_ids = list()

for x in patient_ids:
    if len(x) == 9 and x[4] == '0':
        new_x = x[:4] + x[5:]
        tidy_ids.append(new_x)
    else:
        tidy_ids.append(x)
        
both_sites.insert(loc=1, column='externalID', value=tidy_ids)

# %% [markdown]
# # Make case and NA usage consistent 

# %%
both_sites = both_sites.reset_index(drop = True)
both_sites = both_sites.rename({"C4 _level":"C4_level"}, axis = 1)

with pd.option_context('future.no_silent_downcasting', True):
    both_sites = both_sites.replace({"Not Done":np.nan, "ND":np.nan, "NaT":np.nan, "NK":np.nan}).infer_objects(copy=False)
    both_sites = both_sites.replace({"POSITIVE":"Positive", "NEGATIVE":"Negative"}).infer_objects(copy=False)
    both_sites = both_sites.replace({"YES":"Yes", "NO":"No"}).infer_objects(copy=False)

# %% [markdown]
# # Sort dtypes

# %%
both_sites.C3_level = both_sites.C3_level.astype(float)
both_sites.C4_level = both_sites.C4_level.astype(float)
both_sites.uPCR = both_sites.uPCR.str.replace("<", "") # need to decide what to do with < values
both_sites.uPCR = both_sites.uPCR.astype(float)
both_sites.Albumin = both_sites.Albumin.astype(float)
both_sites.Creatinine = both_sites.Creatinine.astype(float)
both_sites.eGFR = both_sites.eGFR.replace(">90", 90)
both_sites.eGFR = both_sites.eGFR.astype(float)
both_sites.Lymphocyte_count = both_sites.Lymphocyte_count.astype(float)
both_sites.SLEDAI_2K = both_sites.SLEDAI_2K.astype("float")
both_sites.SLEDAI_2K_Clinical = both_sites.SLEDAI_2K_Clinical.astype("float")


# %% [markdown]
# # Condense lupus nephritis classes

# %%
both_sites["LN_Class"].unique()

# %%
both_sites["LN_Class_tidy"] = both_sites["LN_Class"].str.replace("Class","")
both_sites["LN_Class_tidy"] = both_sites["LN_Class_tidy"].str.replace(" ","")
both_sites["LN_Class_tidy"] = both_sites["LN_Class_tidy"].replace({"III/V":"III+V", "IV/V":"IV+V", 
                                                         "III/IV":"III+IV", "IV+III":"III+IV"})

# %%
both_sites["LN_Class_tidy"].unique()

# %% [markdown]
# # Adjust eGFR values

# %%
#account for the fact that Imperial does not calculate eGFR > 90
both_sites["eGFR_adjust"] = np.where((both_sites['eGFR'] > 90), 90, both_sites['eGFR'] )



# %% [markdown]
# # Calculate clinical SLEDAI

# %%
def calculate_clinical_sledai(df):
    new_sledai = df["SLEDAI_2K"]
    if df["C3_status"] == "LOW" or df["C4_status"] == "LOW":
        new_sledai = new_sledai - 2
    if df["dsDNA_SLEDAI_status"] == "HIGH":
        new_sledai = new_sledai - 2
    if new_sledai < 0 :
        return 0
    else:
        return new_sledai

    

# %%
both_sites["cSLEDAI_2K"] = both_sites.apply(calculate_clinical_sledai, axis = 1)


# %%
list(both_sites.columns)

# %% [markdown]
# # Re-order columns

# %%
new_order = ['Patient_ID',
 'externalID',
 'Sample_ID',
 'Site',
 'Recruit_centre',
 'Date_of_consent',
 'Sample_date',
 'Sample_year',
 'Year_of_Birth',
 'Age_at_sampling',
 'Ethnicity_cerner',
 'Diagnosis_date',
 'Diagnosis_year',
 'Years_from_diagnosis',
 'Age_at_diagnosis',
 'Sex',
 'Lupus_status',
 'IFN_score',
 'IFN_status',
 'SLEDAI_2K',
 'SLEDAI_2K_Comments',
 'SLEDAI_2KG',
 'SLEDAI_2K_Clinical',
 'cSLEDAI_2K',
 'LN',
 'LN_date',
 'LN_Class',
 'LN_Class_tidy',
 'Kidney_biopsies',
 'Organ_CVS',
 'Organ_Arterial',
 'Organ_Venous',
 'Organ_Lung_Fibrosis',
 'Multiple_Sclerosis',
 'Thyroid_Disease_Hashimotos',
 'Graves_Disease',
 'Type_I_Diabetes',
 'Rheumatoid_Arthritis',
 'Any_other_comorbidities',
 'MMF',
 'MMF_dose',
 'Pred',
 'Pred_dose',
 'HCQ',
 'HCQ_dose',
 'AZA',
 'AZA_dose',
 'MTX',
 'MTX_dose',
 'TAC',
 'TAC_dose',
 'ASP',
 'ASP_dose',
 'CYC_ever',
 'CYC_dose',
 'CYC_date',
 'CYC_cycles',
 'BEL_ever',
 'BEL_dose',
 'BEL_date',
 'RTX_ever',
 'RTX_dose',
 'RTX_date',
 'RTX_cycles',
 'OFA_ever',
 'OFA_dose',
 'OFA_date',
 'OFA_cycles',
 'dsDNA',
 'dsDNA_titre',
 'dsDNA_units',
 'dsDNA_Ab_ever',
 'dsDNA_date',
 'ANA',
 'ANA_ever',
 'Ribosomal_P',
 'Ro_52',
 'Ro_60',
 'La',
 'Centromere',
 'Sm',
 'Sm_RNP',
 'RNP_68',
 'RNP_A',
 'Scl70',
 'Jo_1',
 'MPO',
 'PR3',
 'LAC_status',
 'ACA_status',
 'ACA_IgG',
 'ACA_IgM',
 'B2G_status',
 'B2G_IgG',
 'B2G_IgM',
 'Vitamin_D',
 'Vitamin_D_date',
 'TSH',
 'TSH_date',
 'C3_level',
 'C4_level',
 'Comp_date',
 'uPCR',
 'uPCR_units',
 'uPCR_date',
 'Albumin',
 'Creatinine',
 'eGFR',
 'eGFR_adjust',
 'Haemoglobin',
 'Lymphopenia_now',
 'Lymphopenia_ever',
 'Lymphocyte_count',
 'Neutrophil_count',
 'Monocyte_count',
 'FBC_Biochem_date',
 'B_cell_n',
 'B_cell_%',
 'CD4_n',
 'CD4_%',
 'CD8_n',
 'CD8_%',
 'Lymph_subsets_date',
 'IgA',
 'IgG',
 'IgM',
 'Ig_date',
 'C3_status',
 'C4_status',
 'dsDNA_SLEDAI_status',
 'ACR_Malar',
 'ACR_Discoid',
 'ACR_Photosensitivity',
 'ACR_oral_ulcers',
 'ACR_Arthritis',
 'ACR_serositis',
 'ACR_seizures/psychosis',
 'ACR_Haem',
 'ACR_Immunologic',
 'KBX_2',
 'KBX_2_LN_Class',
 'KBX_3',
 'KBX_3_LN_Class',
 'KBX_4',
 'KBX_4_LN_Class',
 'KBX_5',
 'KBX_5_LN_Class',
 'KBX_6',
 'KBX_6_LN_Class',
 'KBX_7',
 'KBX_7_LN_Class']

# %%
both_sites = both_sites[new_order]

# %%
both_sites.to_csv(
    "/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/phenotype/combined_data_tidied_30_10_24.csv", index = False)

# %% [markdown]
# # Check that calculated clinical SLEDAI matches that from Imperial

# %%
both_sites.query(
    "Recruit_centre == 'Imperial' and SLEDAI_2K_Clinical != cSLEDAI_2K")[["SLEDAI_2K","SLEDAI_2K_Clinical", "cSLEDAI_2K", "Patient_ID", "SLEDAI_2K_Comments"]]

# %%
both_sites.Ethnicity_cerner.value_counts()

# %% [markdown]
# # Drop individuals not used for single cell sequencing

# %%
metadata = pd.read_csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/metadata/SLEmap_metadata_191224.csv")

# %%
metadata.SC_included.value_counts()

# %%
metadata_sc_use = metadata.query("SC_included == 'use'").copy()

metadata_sc_use["externalID_stripr"] = metadata_sc_use["externalID"].str.replace("r","")

metadata_sc_use_samples = list(metadata_sc_use.externalID_stripr)

# %%
both_sites_sc_use = both_sites.query("externalID in @metadata_sc_use_samples")

# %%
both_sites_sc_use.to_csv(
    "/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/phenotype/combined_data_tidied_sc_used_25_11_24.csv", index = False)

# %% [markdown]
# # Healthy info

# %%
metadata_sc_use_healthy = metadata_sc_use.query("condition == 'HC'").copy()

# %%
metadata_sc_use_healthy

# %% [markdown]
# # King's samples with no C3/C4

# %%
c3_c4_dsdna = both_sites.copy()

# %%
c3_c4_dsdna_kings = c3_c4_dsdna.query("Recruit_centre == 'Kings'")

# %%
c3_c4_dsdna_kings = c3_c4_dsdna_kings[["Patient_ID", "Site", "Sample_date", "C3_level", "C4_level", "dsDNA"]]

# %%
c3_c4_dsdna_kings_missing = c3_c4_dsdna_kings[c3_c4_dsdna_kings.isnull().any(axis=1)]

# %%
c3_c4_dsdna_kings_missing = c3_c4_dsdna_kings_missing.fillna("No")

# %%
c3_c4_dsdna_kings_missing.to_excel("/nfs/users/nfs_c/cs54/OTAR2064_cs54/final_data/1_Metadata/Kings_samples_missing_C3_C4_dsDNA.xlsx")  

# %%
both_sites

# %%
both_sites.Lupus_status.value_counts(dropna=False)

# %%
