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
# # Set up

# %%
# %matplotlib inline
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats
import pingouin as pg
pd.set_option('display.max_columns', None)


# %%
site_palette = {"Kings":"#FC766AFF", "Imperial":"#5B84B1FF"}
active_palette = {"ACTIVE":"#FED439","INACTIVE":"#8A9197"}
ln_palette = {"No":"#D2AF81","Yes":"#FD7446"}

# %%
big_palette = ["#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2","#197EC0","#F05C3B","#46732E","#71D0F5",
"#370335","#075149","#C80813","#91331F", "#1A9993","#FD8CC1"]

# %% [markdown]
# # Read in tidied clinical data 
# Only samples that are used in single cell data (excludes duplicates and unidentifiable samples but not QC outliers)

# %%
both_sites = pd.read_csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/phenotype/combined_data_tidied_sc_used_25_11_24.csv")

# %%
len(both_sites)

# %%
both_sites.Recruit_centre.value_counts()

# %%
both_sites.Site.value_counts()

# %%
figures_folder = "/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/cs54/pilot_run/1_Metadata/clinical_plots"

# %% [markdown]
# # Age

# %%
plt.figure(figsize=(4, 4))

g = sns.stripplot(data = both_sites, x = "Recruit_centre", y = "Age_at_sampling", 
              color = "black", jitter = True, dodge = True, alpha = 0.5)
sns.boxplot(data=both_sites,
        x="Recruit_centre",
        y="Age_at_sampling",
        hue = "Recruit_centre",
        palette = site_palette,
        showfliers = False)
g.yaxis.get_major_locator().set_params(integer=True)

plt.ylim(0,85)
plt.xlabel("Site")
plt.ylabel("Age");

plt.savefig(f"{figures_folder}/age_site_boxplot.png", bbox_inches = "tight")

# %%
both_sites.Age_at_sampling.max()

# %%
both_sites.Age_at_sampling.min()

# %%
both_sites.Age_at_sampling.median()

# %%
both_sites.query("Recruit_centre == 'Imperial'").Age_at_sampling.describe()

# %%
both_sites.query("Recruit_centre == 'Kings'").Age_at_sampling.describe()

# %%
for test_result in ["Age_at_sampling"]:
    x = list(both_sites.dropna(subset = test_result).query("Recruit_centre == 'Imperial'")[test_result])
    y = list(both_sites.dropna(subset = test_result).query("Recruit_centre == 'Kings'")[test_result])
    print(test_result)
    print(pg.normality(dv = test_result,
             group = "Recruit_centre",
             data = both_sites))
    print(pg.homoscedasticity(dv = test_result,
                    group = "Recruit_centre",
                    method = "levene",
                    data = both_sites.dropna(subset = test_result)))
#     print(pg.ttest(x, y))

# %%
p_values = []
for test_result in ["Age_at_sampling"]:
    x = list(both_sites.dropna(subset = test_result).query("Recruit_centre == 'Imperial'")[test_result])
    y = list(both_sites.dropna(subset = test_result).query("Recruit_centre == 'Kings'")[test_result])
    print(test_result)
    pg_out = pg.mwu(x, y)
    print(pg_out)
    print(scipy.stats.mannwhitneyu(x,y))
    p_values.append(pg_out["p-val"].iloc[0])

# %%
sns.scatterplot(data = both_sites, x = "Age_at_sampling", y = "SLEDAI_2K", alpha = 0.5)

# %%
sns.scatterplot(data = both_sites, x = "Age_at_sampling", y = "Years_from_diagnosis", alpha = 0.5)

# %%
sns.scatterplot(data = both_sites, x = "Age_at_sampling", y = "Age_at_diagnosis", alpha = 0.5)

# %%
sns.histplot(both_sites["Age_at_diagnosis"], binwidth = 5)

# %%
both_sites.query("Years_from_diagnosis < 1")

# %%
g = sns.histplot(data = both_sites, x = "Age_at_diagnosis", hue= "Recruit_centre", 
               stat = "percent", common_norm = False, palette = site_palette, legend = True, binwidth = 5)

g.get_legend().set_title("Site")

plt.ylabel("Percentage of patients (per site)")
plt.xlabel("Age at diagnosis");

# plt.savefig(f"{figures_folder}/sledai_2k_site_barplot.png", bbox_inches = "tight")

# %%
g = sns.histplot(data = both_sites, x = "Years_from_diagnosis", hue= "Recruit_centre", 
               stat = "percent", common_norm = False, palette = site_palette, legend = True, binwidth = 5)

g.get_legend().set_title("Site")

plt.ylabel("Percentage of patients (per site)")
plt.xlabel("Years from diagnosis");

# plt.savefig(f"{figures_folder}/sledai_2k_site_barplot.png", bbox_inches = "tight")

# %%
both_sites.query("Age_at_diagnosis < 10")

# %% [markdown]
# # Plot general test results

# %%
plt.figure(figsize=(9, 9))
# plt.subplots_adjust(hspace=0.5)

for n, test_result in enumerate(['C3_level', 'C4_level', 'uPCR', 'Albumin', 'Lymphocyte_count',
       'eGFR', 'eGFR_adjust', 'Creatinine']):
    if n > 1:
        ax = plt.subplot(4, 3, n + 1, sharex =  plt.subplot(4, 3, n))
    else:
        ax = plt.subplot(4, 3, n + 1)
    
    sns.boxplot(data = both_sites, x = "Recruit_centre", y = test_result, 
            color = "white", showfliers = False, ax = ax)
    sns.stripplot(data = both_sites, x = "Recruit_centre", y = test_result, 
            color = "black", ax = ax, alpha = 0.3)
    plt.ylim(None, both_sites[test_result].max() + (0.2*both_sites[test_result].max()))
    plt.xlabel("Recruitment Centre")
    
    if test_result == "Creatinine":
        ax = plt.subplot(4, 3, n + 2, sharex =  plt.subplot(4, 3, n + 1))
        
        sns.boxplot(data = both_sites.query("Patient_ID != 'EQTL0125p'"), x = "Recruit_centre", y = test_result, 
            color = "white", showfliers = False, ax = ax)
        sns.stripplot(data = both_sites.query("Patient_ID != 'EQTL0125p'"), x = "Recruit_centre", y = test_result, 
                    color = "black", ax = ax, alpha = 0.3)
    plt.xlabel("Recruitment Centre")

plt.tight_layout()

# plt.savefig("Figures/fbc_and_biochem_boxes.png")

# %%
plt.figure(figsize=(12, 4))
plt.subplots_adjust(hspace=0.1)

for n, test_result in enumerate(['Haemoglobin', 'Neutrophil_count', 'Monocyte_count', 'Vitamin_D']):
    ax = plt.subplot(1, 4, n + 1)
    
    sns.boxplot(data = both_sites, x = "Recruit_centre", y = test_result, 
            color = "white", showfliers = False, ax = ax)
    sns.stripplot(data = both_sites, x = "Recruit_centre", y = test_result, 
            color = "black", ax = ax,alpha = 0.3)
plt.tight_layout()


#         plt.show()

# %%
for test_result in ['C3_level', 'C4_level', 'uPCR', 'Albumin', 'Lymphocyte_count',
       'eGFR', 'eGFR_adjust', 'Creatinine']:
    x = list(both_sites.dropna(subset = test_result).query("Recruit_centre == 'Imperial'")[test_result])
    y = list(both_sites.dropna(subset = test_result).query("Recruit_centre == 'Kings'")[test_result])
    print(test_result)
    print(pg.normality(dv = test_result,
             group = "Recruit_centre",
             data = both_sites))
    print(pg.homoscedasticity(dv = test_result,
                    group = "Recruit_centre",
                    method = "levene",
                    data = both_sites.dropna(subset = test_result)))
#     print(pg.ttest(x, y))
    

# %%
p_values = []
for test_result in ['C3_level', 'C4_level', 'Albumin', 'Lymphocyte_count',
       'eGFR_adjust', 'Creatinine']:
    x = list(both_sites.dropna(subset = test_result).query("Recruit_centre == 'Imperial'")[test_result])
    y = list(both_sites.dropna(subset = test_result).query("Recruit_centre == 'Kings'")[test_result])
    print(test_result)
    pg_out = pg.mwu(x, y)
    print(pg_out)
    print(scipy.stats.mannwhitneyu(x,y))
    p_values.append(pg_out["p-val"].iloc[0])

# %%
scipy.stats.false_discovery_control(p_values)

# %%
reject, pvals_corr = pg.multicomp(p_values, method='bonf')
print(reject, pvals_corr)

# %% [markdown]
# # Look at differences in Albumin

# %%
test_result = 'Albumin'

g = sns.boxplot(data = both_sites, x = "Site", y = test_result, 
        showfliers = False, color = "white")
sns.stripplot(data = both_sites, x = "Site", y = test_result, 
        color = "black", dodge = True)
handles, labels = g.get_legend_handles_labels()
# plt.legend(handles = handles[0:2], title = "Lupus nephritis")


# %%
print(pg.normality(dv = "Albumin",
         group = "Site",
         data = both_sites))
print(pg.homoscedasticity(dv = "Albumin",
                group = "Site",
                method = "levene",
                data = both_sites.dropna(subset = "Albumin")))

# %%
test_result = 'Albumin'

g = sns.boxplot(data = both_sites, x = "Site", y = test_result, 
        palette = "Set1", showfliers = False, hue = "LN")
sns.stripplot(data = both_sites, x = "Site", y = test_result, 
        color = "black", hue = "LN", dodge = True)
handles, labels = g.get_legend_handles_labels()
plt.legend(handles = handles[0:2], title = "Lupus nephritis")


# %%
test_result = 'Albumin'

g = sns.boxplot(data = both_sites, x = "Recruit_centre", y = test_result, 
        palette = "Set1", showfliers = False, hue = "LN")
sns.stripplot(data = both_sites, x = "Recruit_centre", y = test_result, 
        color = "black", hue = "LN", dodge = True)
handles, labels = g.get_legend_handles_labels()
plt.legend(handles = handles[0:2], title = "Lupus nephritis")


# %% [markdown]
# # Autoantibodies

# %%
plt.figure(figsize=(21, 12))
plt.subplots_adjust(hspace=0.5)

for n, feature in enumerate(['dsDNA_Ab_ever', 'ANA', 'ANA_ever', 'Ribosomal_P', 'Ro_52', 'Ro_60',
       'La', 'Centromere', 'Sm', 'Sm_RNP', 'RNP_68', 'RNP_A', 'Scl70', 'Jo_1',
       'MPO', 'PR3', 'LAC_status', 'ACA_status', 'B2G_status']):

    feature_df = both_sites.groupby("Recruit_centre")[feature].value_counts(normalize = True, dropna = False).reset_index()
    feature_df = feature_df.replace(np.nan, "NA")
    
    ax = plt.subplot(4, 5, n + 1)

    sns.barplot(data = feature_df, x = feature, y = "proportion", legend = False,
            hue = "Recruit_centre", order = ["Positive", "Negative", "NA"], palette = site_palette, ax = ax)
    plt.ylabel("Proportion of patients")
    plt.xlabel("Test result")
    plt.title(f"{feature}")


# %% [markdown]
# # Lupus nephritis

# %%
plt.figure(figsize=(6, 6))

feature_df = both_sites.groupby("Recruit_centre")["LN"].value_counts(normalize = True, dropna = False).reset_index()
feature_df = feature_df.replace(np.nan, "NA")

sns.barplot(data = feature_df, x = "LN", y = "proportion", 
        hue = "Recruit_centre", palette = site_palette)
plt.ylabel("Proportion of patients")
plt.xlabel("Lupus Nephritis Present")

plt.legend(title = "Site")


plt.savefig(f"{figures_folder}/nephritis_site_barplot.png", bbox_inches = "tight")

# %%
both_sites.LN.value_counts()

# %%
both_sites.groupby("Recruit_centre").LN.value_counts()

# %%
plt.figure(figsize=(12, 6))
feature_df = both_sites.query("LN == 'Yes'").copy()
feature_df = feature_df.groupby("Recruit_centre")["LN_Class"].value_counts(normalize = True, dropna = False).reset_index()
feature_df = feature_df.replace(np.nan, "NA")


sns.barplot(data = feature_df, x = "LN_Class", y = "proportion", 
        hue = "Recruit_centre", palette = site_palette)
plt.ylabel("Proportion of patients with lupus nephritis")
plt.xlabel("Class")
plt.xticks(rotation = 30, ha = "right")
plt.title("Lupus Nephritis Class");

# %% [markdown]
# # Medications

# %%
plt.figure(figsize=(14, 8))
plt.subplots_adjust(hspace=0.5)

for n, feature in enumerate(['MMF', 'Pred', 'HCQ', 'AZA', 'ASP']):

    feature_df = both_sites.groupby("Recruit_centre")[feature].value_counts(normalize = True, dropna = False).reset_index()
    feature_df = feature_df.replace(np.nan, "NA")

    ax = plt.subplot(2, 3, n + 1)
    sns.barplot(data = feature_df, x = feature, y = "proportion", 
            hue = "Recruit_centre", palette = site_palette, order = ["Yes", "No"], ax = ax, legend = True)
    plt.ylabel("Proportion of patients")
    plt.xlabel("Drug prescribed")
    plt.ylim(None, 1)
    plt.title(f"{feature}")
    plt.legend(title = "Site")


plt.savefig(f"{figures_folder}/drugs_site_barplot.png", bbox_inches = "tight")

# %% [markdown]
# # SLEDAI Scores

# %%
g = sns.histplot(data = both_sites, x = "SLEDAI_2K", hue= "Recruit_centre", 
               stat = "percent", common_norm = False, palette = site_palette, legend = True, binwidth = 1)

g.get_legend().set_title("Site")

plt.ylabel("Percentage of patients with score (per site)")
plt.xlabel("SLEDAI Score");

plt.savefig(f"{figures_folder}/sledai_2k_site_barplot.png", bbox_inches = "tight")

# %%
g = sns.histplot(data = both_sites, x = "SLEDAI_2K", hue= "Recruit_centre", 
               stat = "percent", common_norm = False, palette = site_palette, legend = True, binwidth = 2,
                hue_order = ["Imperial", "Kings"] )

g.get_legend().set_title("Recruitment Centre")

plt.ylabel("Percentage of patients with score (per site)")
plt.xlabel("SLEDAI Score");

# plt.savefig("SLEDAI_score_by_centre_histogram_with_new_data_Mar24.png", bbox_inches = "tight")

# %%
plt.figure(figsize=(4, 4))


g = sns.stripplot(data = both_sites, y = "SLEDAI_2K", x = "Recruit_centre", 
              color = "black", jitter = True, order = ["Imperial", "Kings"])
sns.boxplot(data=both_sites,
        x="Recruit_centre",
        y="SLEDAI_2K",
        palette = site_palette,
        showfliers = False,
        order = ["Imperial", "Kings"])
g.yaxis.get_major_locator().set_params(integer=True)


plt.xlabel("Site")
plt.ylabel("SLEDAI_2K Score");

plt.savefig(f"{figures_folder}/sledai_site_boxplot.png", bbox_inches = "tight")

# plt.savefig("Figures/SLEDAI_score_by_centre_boxplot.png", bbox_inches = "tight")

# %%
both_sites

# %%
plt.figure(figsize=(4, 4))


g = sns.stripplot(data = both_sites, y = "SLEDAI_2K", x = "Site", 
              color = "black", jitter = True, order = ["Imperial", "Kings", "Guys"])
sns.boxplot(data=both_sites,
        x="Site",
        y="SLEDAI_2K",
#         palette = site_palette,
        showfliers = False,
        order = ["Imperial", "Kings", "Guys"])
g.yaxis.get_major_locator().set_params(integer=True)


plt.xlabel("Site")
plt.ylabel("SLEDAI_2K Score");

plt.savefig(f"{figures_folder}/sledai_discrete_site_boxplot.png", bbox_inches = "tight")

# plt.savefig("Figures/SLEDAI_score_by_centre_boxplot.png", bbox_inches = "tight")

# %%
both_sites

# %%
plt.figure(figsize=(4, 4))


g = sns.stripplot(data = both_sites, y = "cSLEDAI_2K", x = "Site", 
              color = "black", jitter = True, order = ["Imperial", "Kings", "Guys"])
sns.boxplot(data=both_sites,
        x="Site",
        y="cSLEDAI_2K",
#         palette = site_palette,
        showfliers = False,
        order = ["Imperial", "Kings", "Guys"])
g.yaxis.get_major_locator().set_params(integer=True)


plt.xlabel("Site")
plt.ylabel("cSLEDAI_2K Score");

plt.savefig(f"{figures_folder}/clinical_sledai_discrete_site_boxplot.png", bbox_inches = "tight")

# plt.savefig("Figures/SLEDAI_score_by_centre_boxplot.png", bbox_inches = "tight")

# %%
both_sites.query("Recruit_centre == 'Imperial'").SLEDAI_2K.describe()

# %%
both_sites.query("Recruit_centre == 'Kings'").SLEDAI_2K.describe()

# %%
both_sites.SLEDAI_2K.describe()

# %%
plt.figure(figsize=(4, 4))


g = sns.stripplot(data = both_sites, x = "LN", y = "SLEDAI_2K", hue = "Recruit_centre", 
              color = "black", jitter = True, dodge = True)
sns.boxplot(data=both_sites,
        x="LN",
        y="SLEDAI_2K",
        palette = site_palette,
        showfliers = False,
        hue =  "Recruit_centre")
g.yaxis.get_major_locator().set_params(integer=True)
handles, labels = g.get_legend_handles_labels()
plt.legend(handles = handles[2:4], title = "Recruitment centre", loc='center left', bbox_to_anchor=(1.0, 0.5));

plt.xlabel("Lupus Nephritis")
plt.ylabel("SLEDAI_2K Score");

# plt.savefig("Figures/SLEDAI_score_by_centre_boxplot.png", bbox_inches = "tight")

# %%
plt.figure(figsize=(4, 4))


g = sns.stripplot(data = both_sites, x = "Lupus_status", y = "SLEDAI_2K",
              color = "black", jitter = True, dodge = True)
sns.boxplot(data=both_sites,
        x="Lupus_status",
        y="SLEDAI_2K",
        palette = active_palette,
        showfliers = False)
g.yaxis.get_major_locator().set_params(integer=True)
handles, labels = g.get_legend_handles_labels()
# plt.legend(handles = handles[2:4], title = "Recruitment centre", loc='center left', bbox_to_anchor=(1.0, 0.5));

plt.xlabel("Disease_activity")
plt.ylabel("SLEDAI_2K Score");

# plt.savefig("Figures/SLEDAI_score_by_centre_boxplot.png", bbox_inches = "tight")

# %%
plt.figure(figsize=(4, 4))


g = sns.stripplot(data = both_sites, x = "LN", y = "SLEDAI_2K", 
              color = "black", jitter = True, dodge = True)
sns.boxplot(data=both_sites,
        x="LN",
        y="SLEDAI_2K",
        palette = ln_palette,
        showfliers = False)
g.yaxis.get_major_locator().set_params(integer=True)


plt.xlabel("Lupus Nephritis")
plt.ylabel("SLEDAI_2K Score");

# plt.savefig("Figures/SLEDAI_score_by_centre_boxplot.png", bbox_inches = "tight")

# %% [markdown]
# # Clinical SLEDAI scores

# %%
both_sites.C3_status.value_counts(dropna = False)

both_sites.query("Recruit_centre == 'Kings'").C3_status.value_counts(dropna = False)

both_sites.query("Recruit_centre == 'Kings'").C4_status.value_counts(dropna = False)

both_sites.query("Recruit_centre == 'Kings'").dsDNA.value_counts(dropna = False)

# %%
plt.figure(figsize=(4, 4))


g = sns.stripplot(data = both_sites, y = "cSLEDAI_2K", x = "Recruit_centre", 
              color = "black", jitter = True, order = ["Imperial", "Kings"], alpha = 0.2)
sns.boxplot(data=both_sites,
        x="Recruit_centre",
        y="cSLEDAI_2K",
        palette = site_palette,
        showfliers = False,
        hue = "Recruit_centre",
        order = ["Imperial", "Kings"])
g.yaxis.get_major_locator().set_params(integer=True)


plt.xlabel("Recruitment Centre")
plt.ylabel("Clinical SLEDAI_2K Score");


plt.savefig(f"{figures_folder}/csledai_site_boxplot.png", bbox_inches = "tight")


# %%
plt.figure(figsize=(4, 4))


g = sns.stripplot(data = both_sites, y = "cSLEDAI_2K", x = "Recruit_centre", 
              color = "black", jitter = True, order = ["Imperial", "Kings"], alpha = 0.2)
sns.violinplot(data=both_sites,
        x="Recruit_centre",
        y="cSLEDAI_2K",
        palette = site_palette,
        hue = "Recruit_centre",
        order = ["Imperial", "Kings"],
              inner= "quart")
g.yaxis.get_major_locator().set_params(integer=True)


plt.xlabel("Recruitment Centre")
plt.ylabel("Clinical SLEDAI_2K Score");


plt.savefig(f"{figures_folder}/csledai_site_violinplot.png", bbox_inches = "tight")


# %%
g = sns.histplot(data = both_sites, x = "cSLEDAI_2K", hue= "Recruit_centre", 
               stat = "percent", common_norm = False, palette = site_palette, legend = True, binwidth = 2,
                hue_order = ["Imperial", "Kings"] )

g.get_legend().set_title("Recruitment Centre")

plt.ylabel("Percentage of patients with score (per site)")
plt.xlabel("Clinical SLEDAI Score");

# plt.savefig("SLEDAI_score_by_centre_histogram_with_new_data_Mar24.png", bbox_inches = "tight")

# %%
g = sns.histplot(data = both_sites, x = "cSLEDAI_2K", hue= "Recruit_centre", 
               stat = "percent", common_norm = False, palette = site_palette, legend = True, binwidth = 1,
                hue_order = ["Imperial", "Kings"] )

g.get_legend().set_title("Site")

plt.ylabel("Percentage of patients with score (per site)")
plt.xlabel("Clinical SLEDAI Score");

plt.savefig(f"{figures_folder}/clsedai_site_barplot.png", bbox_inches = "tight")



# %%
both_sites.query("Recruit_centre == 'Imperial'").cSLEDAI_2K.describe()

# %%
both_sites.query("Recruit_centre == 'Kings'").cSLEDAI_2K.describe()

# %%
both_sites.cSLEDAI_2K.describe()

# %%

# %%
print(pg.normality(dv = "SLEDAI_2K",
         group = "Recruit_centre",
         data = both_sites))
print(pg.homoscedasticity(dv = "SLEDAI_2K",
                group = "Recruit_centre",
                method = "levene",
                data = both_sites.dropna(subset = "SLEDAI_2K")))

# %%
x = list(both_sites.dropna(subset = "SLEDAI_2K").query("Recruit_centre == 'Imperial'")["SLEDAI_2K"])
y = list(both_sites.dropna(subset = "SLEDAI_2K").query("Recruit_centre == 'Kings'")["SLEDAI_2K"])

print(scipy.stats.mannwhitneyu(x,y))


# %%
print(pg.normality(dv = "cSLEDAI_2K",
         group = "Recruit_centre",
         data = both_sites))
print(pg.homoscedasticity(dv = "cSLEDAI_2K",
                group = "Recruit_centre",
                method = "levene",
                data = both_sites.dropna(subset = "cSLEDAI_2K")))

# %%
x = list(both_sites.dropna(subset = "cSLEDAI_2K").query("Recruit_centre == 'Imperial'")["cSLEDAI_2K"])
y = list(both_sites.dropna(subset = "cSLEDAI_2K").query("Recruit_centre == 'Kings'")["cSLEDAI_2K"])

print(scipy.stats.mannwhitneyu(x,y))


# %%
both_sites

# %% [markdown]
# # Ancestry distributions

# %%
metadata = pd.read_csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/metadata/SLEmap_metadata_291124.csv")

# %%
metadata["externalID_stripr"] = metadata["externalID"].str.replace("r","")
metadata = metadata.drop("externalID", axis =1)

# %%
both_ancestry = pd.merge(both_sites, metadata, left_on = "externalID", right_on = "externalID_stripr")
both_ancestry = both_ancestry.sort_values(by = "ancestry.l1")
both_ancestry["ancestry_plot"] = both_ancestry["ancestry.l1"].str.replace(" ancestry","")

ancestry_order = list(both_ancestry["ancestry.l1"].unique())

# %%
set(list(both_sites["externalID"])) - set(list(metadata["externalID_stripr"]))

# %%
len(metadata)

# %%
len(both_ancestry)

# %%
ancestry_props = both_ancestry.groupby("Recruit_centre")["ancestry_plot"].value_counts(normalize = True).reset_index()

# %%
ancestry_props

# %%
ancestry_props_clin = both_ancestry.groupby("Recruit_centre")["Ethnicity_cerner"].value_counts(normalize = True).reset_index()

# %%
ancestry_props_clin

# %%
both_ancestry.query("ancestry_plot== 'East Asian'")[["Ethnicity_cerner","ancestry.l1"]]

# %%
both_ancestry.query("ancestry_plot== 'Mixed'")[["Ethnicity_cerner","ancestry.l1"]]

# %%
plt.figure(figsize=(6, 6))

sns.barplot(data = ancestry_props, x = "ancestry_plot", y = "proportion", 
        hue = "Recruit_centre", palette = site_palette)
plt.ylabel("Proportion of patients")
plt.xlabel("Self-reported ancestry")
plt.xticks(rotation = 45, ha = "right");

# %%
g = sns.stripplot(data = both_ancestry, y = "SLEDAI_2K", x = "ancestry.l1", jitter = True, hue = "Recruit_centre",
             dodge = True, palette = ["black","black"], legend = False, order = ancestry_order, 
                  hue_order = ["Imperial", "Kings"])
sns.boxplot(data=both_ancestry,
        x="ancestry.l1",
        y="SLEDAI_2K",
        showfliers = False,
        hue = "Recruit_centre",
        palette = site_palette,
        order = ancestry_order,
           hue_order = ["Imperial", "Kings"])

plt.xlabel("Self-reported ancestry")
plt.ylabel("SLEDAI Score")
plt.xticks(rotation = 45, ha = "right")
plt.ylim(None, 19) 
g.yaxis.get_major_locator().set_params(integer=True)

plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title = "Recruitment Centre");

# plt.savefig("SLEDAI_score_by_ancestry.png", bbox_inches = "tight")

# %%
both_ancestry_props = both_ancestry.groupby("Recruit_centre")["ancestry.l1"].value_counts(normalize = True,dropna = False).reset_index()

both_ancestry_props["percentage"] = both_ancestry_props["proportion"] * 100

# %%
g = sns.barplot(data = both_ancestry_props, x = "ancestry.l1", y = "percentage", hue = "Recruit_centre",
             dodge = True, legend = True, palette = site_palette, order = ancestry_order)


plt.xlabel("Self-reported ancestry")
plt.ylabel("Percentage of samples (per site)")
plt.xticks(rotation = 45, ha = "right")
g.yaxis.get_major_locator().set_params(integer=True)

plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title = "Recruitment Centre");

# %%
g = sns.stripplot(data = both_ancestry, y = "Age_at_sampling", x = "ancestry.l1", jitter = True, hue = "Recruit_centre",
             dodge = True, palette = ["black","black"], legend = False, order = ancestry_order, 
                  hue_order = ["Imperial", "Kings"])
sns.boxplot(data=both_ancestry,
        x="ancestry.l1",
        y="Age_at_sampling",
        showfliers = False,
        hue = "Recruit_centre",
        palette = site_palette,
        order = ancestry_order,
        hue_order = ["Imperial", "Kings"])

plt.xlabel("Self-reported ancestry")
plt.ylabel("Age_at_sampling")
plt.xticks(rotation = 45, ha = "right")
g.yaxis.get_major_locator().set_params(integer=True)

plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title = "Recruitment Centre");

# plt.savefig("SLEDAI_score_by_ancestry.png", bbox_inches = "tight")

# %% [markdown]
# # Genetic ancestry (Wanseon assignment)

# %%
wan_ancestry = pd.read_csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/wl2/wgs/01_3.variant_GT_QCed/01.populations/1.3.SLEmap_estimated_ancestry.txt", sep = "\t")

# %%
wan_ancestry.predicted_group.unique()

# %%
len(both_ancestry)

# %%
both_ancestry = pd.merge(both_ancestry, wan_ancestry[["id", "predicted_group" ]], left_on = "WGS_ID", right_on = "id")

# %%
wan_ancestry_props = both_ancestry.groupby("Recruit_centre").predicted_group.value_counts(normalize = True).reset_index()

# %%
wan_ancestry_props

# %%
plt.figure(figsize=(6, 6))

sns.barplot(data = wan_ancestry_props, x = "predicted_group", y = "proportion", 
        hue = "Recruit_centre", palette = site_palette)
plt.ylabel("Proportion of patients")
plt.xlabel("Predicted genetic ancestry")
plt.legend(title = "Site")

plt.savefig(f"{figures_folder}/wgs_site_barplot.png", bbox_inches = "tight")

# %%

# %%
g = sns.stripplot(data = both_ancestry, y = "Age_at_diagnosis", x = "predicted_group", 
                  jitter = True, hue = "Recruit_centre",
                  dodge = True, palette = ["black","black"], 
                  legend = False,
                  hue_order = ["Imperial", "Kings"])
sns.boxplot(data=both_ancestry,
        x="predicted_group",
        y="Age_at_diagnosis",
        showfliers = False,
        hue = "Recruit_centre",
        palette = site_palette,
        hue_order = ["Imperial", "Kings"])

plt.xlabel("Predicted ancestry")
plt.ylabel("Age_at_diagnosis")
plt.xticks(rotation = 45, ha = "right")
g.yaxis.get_major_locator().set_params(integer=True)

plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title = "Recruitment Centre");

# plt.savefig("SLEDAI_score_by_ancestry.png", bbox_inches = "tight")

# %%
g = sns.stripplot(data = both_ancestry, y = "Age_at_diagnosis", x = "Recruit_centre", 
                  jitter = True, hue = "predicted_group",
                  dodge = True, palette = ["black","black"], 
                  legend = False)
sns.boxplot(data=both_ancestry,
        x="Recruit_centre",
        y="Age_at_diagnosis",
        showfliers = False,
        hue = "predicted_group")

plt.xlabel("Predicted ancestry")
plt.ylabel("Age_at_diagnosis")
plt.xticks(rotation = 45, ha = "right")
g.yaxis.get_major_locator().set_params(integer=True)

plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title = "Recruitment Centre");

# plt.savefig("SLEDAI_score_by_ancestry.png", bbox_inches = "tight")

# %%
g = sns.stripplot(data = both_ancestry, y = "Age_at_sampling", x = "Recruit_centre", 
                  jitter = True, hue = "predicted_group",
                  dodge = True, palette = ["black","black"], 
                  legend = False)
sns.boxplot(data=both_ancestry,
        x="Recruit_centre",
        y="Age_at_sampling",
        showfliers = False,
        hue = "predicted_group")

plt.xlabel("Predicted ancestry")
plt.ylabel("Age_at_sampling")
g.yaxis.get_major_locator().set_params(integer=True)

plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title = "Recruitment Centre");

# plt.savefig("SLEDAI_score_by_ancestry.png", bbox_inches = "tight")

# %%
g = sns.stripplot(data = both_ancestry, y = "Age_at_sampling", x = "predicted_group", 
                  jitter = True, hue = "Recruit_centre",
                  dodge = True, palette = ["black","black"], 
                  legend = False,
                  hue_order = ["Imperial", "Kings"])
sns.boxplot(data=both_ancestry,
        x="predicted_group",
        y="Age_at_sampling",
        showfliers = False,
        hue = "Recruit_centre",
        palette = site_palette,
        hue_order = ["Imperial", "Kings"])

plt.xlabel("Predicted ancestry")
plt.ylabel("Age_at_sampling")
plt.xticks(rotation = 45, ha = "right")
g.yaxis.get_major_locator().set_params(integer=True)

plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title = "Recruitment Centre");

# plt.savefig("SLEDAI_score_by_ancestry.png", bbox_inches = "tight")

# %% [markdown]
# # Check that we have all the clinical data for sequenced samples

# %%
sequenced_kings = list(metadata.query("condition == 'SLE' and recruitment_centre != 'Imperial'")["externalID"].str.rstrip("r"))

data_kings = list(both_sites.query("Recruit_centre == 'Kings'")["externalID"])

# %%
set(sequenced_kings) - set(new_data_kings)

# %%
sequenced_imperial = list(metadata.query("condition == 'SLE' and recruitment_centre == 'Imperial'")["externalID"].str.rstrip("r"))

data_imperial = list(both_sites.query("Recruit_centre == 'Imperial'")["externalID"])

# %%
set(sequenced_imperial) - set(data_imperial)

# %%
