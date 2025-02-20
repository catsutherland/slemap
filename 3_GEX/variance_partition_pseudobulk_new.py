# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: groups/hgi/ubeR/2024-04-21
#     language: R
#     name: ugbir0kvc29mdhbhy2svz3jvdxbzl2hnas91ymvslziwmjqtmdqtmjek
# ---

# %%
library(variancePartition)
library(tidyr)
library(dplyr)
library("biomaRt")
library(edgeR)
library(DGEobj.utils)
library(gridGraphics)
library(glue)

# %%
sessionInfo()

# %%
figures_folder = "/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/cs54/final_data/Figures/variancePartition"

# %%
all_pseudo <- readRDS("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/hj10/Results/3_eQTL_prep/pseudobulk/meanPB_level0.RDS")

# %%
clinical_all <- read.csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/phenotype/combined_data_tidied_sc_used_25_11_24.csv")

# %%
nrow(clinical_all)

# %%
metadata <- read.csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/metadata/SLEmap_metadata_040225.csv")

# %%
metadata$externalID <- sub("r$", "", metadata$externalID)

# %%
clinical_all <- merge(metadata[, c("externalID","experiment_id", "poolID", "Flow_cell",
                                   "id_pool_lims", "date_processed")],
                      clinical_all, by = c("externalID"))

# %%
nrow(clinical_all)

# %%
clinical_all$experiment_id <- gsub("_", "-", clinical_all$experiment_id)

# %%
rownames(clinical_all) <- clinical_all$experiment_id

# %%
clinical_all <- clinical_all[colnames(all_pseudo), ]

# %%
all_pseudo <- as.data.frame(as.matrix(all_pseudo))

# %%
clinical_all$Flow_cell<-as.character(clinical_all$Flow_cell)
clinical_all$date_processed<-as.character(clinical_all$date_processed)


# %%
form <- ~ (1|externalID) + (1|poolID)+ (1|id_pool_lims) + 
            (1|Flow_cell) + (1|Recruit_centre) + Age_at_sampling + (1|date_processed)

# Compute Canonical Correlation Analysis (CCA)
# between all pairs of variables
# returns absolute correlation value
C <- canCorPairs(form, clinical_all)

# Plot correlation matrix
# between all pairs of variables
plotCorrMatrix(C)

# %%
form1 <- ~ (1|poolID) + (1|date_processed) + (1|id_pool_lims) + (1|Flow_cell) + (1|Recruit_centre) + Age_at_sampling + SLEDAI_2K
form2 <- ~ (1|poolID) + (1|Recruit_centre) + Age_at_sampling + SLEDAI_2K
form3 <- ~ (1|date_processed) + (1|Recruit_centre) + Age_at_sampling + SLEDAI_2K
form4 <- ~ (1|id_pool_lims) + (1|Recruit_centre) + Age_at_sampling + SLEDAI_2K
form5 <- ~ (1|Flow_cell) + (1|Recruit_centre) + Age_at_sampling + SLEDAI_2K


# %%
forms <- c(form1, form2, form3, form4, form5)

# %%
for (i in seq_along(forms)) {
   varPart <- fitExtractVarPartModel(all_pseudo, forms[[i]], clinical_all)
    vp <- sortCols( varPart )
    p1 <- plotPercentBars( vp[1:50,] )
    p2 <- plotVarPart( vp) 
    print(p1)
    print(p2)
    ggsave(glue("{figures_folder}/variance_bar_pseudol0_form{i}.pdf"), p1)
    ggsave(glue("{figures_folder}/variance_violin_pseudol0_form{i}.pdf"), p2)
}

# %%

# %%
pseudo_l1 <- readRDS("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/hj10/Results/3_eQTL_prep/pseudobulk/meanPB_level1.RDS")

# %%
# Define the suffixes you want to add
suffixes <- c("_B", "_CD4 T", "_CD8 T", "_Mono", "_NK", "_other T", "_other", "_DC")

# Create an empty list to store the expanded rows
expanded_rows <- list()

# Loop through each row and duplicate it with new row names
for (row in rownames(clinical_all)) {
  # Get the base row name (e.g., "SLE-map13436155-LDP41")
  base_name <- row
  
  # Create new rows with suffixes
  for (suffix in suffixes) {
    new_rowname <- paste0(base_name, suffix)
    
    # Add the row to the expanded_rows list
    expanded_rows[[new_rowname]] <- clinical_all[row, ]
  }
}

# Convert the list to a data frame
expanded_clinical_all <- do.call(rbind, expanded_rows)


# %%
expanded_clinical_all$celltype <- sub(".*_(.*)", "\\1", rownames(expanded_clinical_all))

# %%
expanded_clinical_all <- expanded_clinical_all[colnames(pseudo_l1), ]

# %%
pseudo_l1 <- as.data.frame(as.matrix(pseudo_l1))

# %%
form1 <- ~ (1|celltype) + (1|poolID) + (1|date_processed) + (1|id_pool_lims) + (1|Flow_cell) + (1|Recruit_centre) + Age_at_sampling + SLEDAI_2K
form2 <- ~ (1|celltype) + (1|poolID) + (1|Recruit_centre) + Age_at_sampling + SLEDAI_2K
form3 <- ~ (1|celltype) + (1|date_processed) + (1|Recruit_centre) + Age_at_sampling + SLEDAI_2K
form4 <- ~ (1|celltype) + (1|id_pool_lims) + (1|Recruit_centre) + Age_at_sampling + SLEDAI_2K
form5 <- ~ (1|celltype) + (1|Flow_cell) + (1|Recruit_centre) + Age_at_sampling + SLEDAI_2K

forms <- c(form1, form2, form3, form4, form5)


for (i in seq_along(forms)) {
    varPart <- fitExtractVarPartModel(pseudo_l1, forms[[i]], expanded_clinical_all)
    vp <- sortCols( varPart )
    p1 <- plotPercentBars( vp[1:30,] )
    p2 <- plotVarPart( vp)
    print(p1)
    print(p2)
    ggsave(glue("{figures_folder}/variance_bar_pseudol1_form{i}.pdf"), p1)
    ggsave(glue("{figures_folder}/variance_violin_pseudol1_form{i}.pdf"), p2)
}

# %%
form1 <- ~ (1|poolID) + (1|date_processed) + (1|id_pool_lims) + (1|Flow_cell) + (1|Recruit_centre) + Age_at_sampling + SLEDAI_2K
form2 <- ~ (1|poolID) + (1|Recruit_centre) + Age_at_sampling + SLEDAI_2K
form3 <- ~ (1|date_processed) + (1|Recruit_centre) + Age_at_sampling + SLEDAI_2K
form4 <- ~ (1|id_pool_lims) + (1|Recruit_centre) + Age_at_sampling + SLEDAI_2K
form5 <- ~ (1|Flow_cell) + (1|Recruit_centre) + Age_at_sampling + SLEDAI_2K

forms <- c(form1, form2, form3, form4, form5)


for (i in seq_along(forms)) {
    for (celltype in c("_B", "_CD4 T", "_CD8 T", "_other T", "_NK", "_DC", "_Mono", "_other")) {
        print(gsub("_", "", celltype))
        pseudo_celltype <- pseudo_l1[, grep(celltype, names(pseudo_l1))]
        clinical_celltype <- expanded_clinical_all[colnames(pseudo_celltype), ]
        varPart <- fitExtractVarPartModel(pseudo_celltype, forms[[i]], clinical_celltype)
        vp <- sortCols( varPart )
        p1 <- plotPercentBars( vp[1:30,] )
        p2 <- plotVarPart( vp)
        print(p1)
        print(p2)
        ggsave(glue("{figures_folder}/variance_violin{celltype}_form{i}.pdf"), p2)
        ggsave(glue("{figures_folder}/variance_bar{celltype}_form{i}.pdf"), p1)
    }
   } 

# %%
for (celltype in c("_B", "_CD4 T", "_CD8 T", "_other T", "_NK", "_DC", "_Mono", "_other")) {
    print(gsub("_", "", celltype))
    pseudo_celltype <- pseudo_l1[, grep(celltype, names(pseudo_l1))]
    clinical_celltype <- expanded_clinical_all[colnames(pseudo_celltype), ]
    form <- ~ (1|poolID) + (1|Recruit_centre) + Age_at_sampling + SLEDAI_2K
    varPart <- fitExtractVarPartModel(pseudo_celltype, form, clinical_celltype)
    vp <- sortCols( varPart )
    p1 <- plotPercentBars( vp[1:10,] )
    p2 <- plotVarPart( vp)
    print(p1)
    print(p2)
    ggsave(glue("{figures_folder}/variance_violin{celltype}.pdf"), p2)
    ggsave(glue("{figures_folder}/variance_bar{celltype}.pdf"), p1)

}

# %%
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <-  rownames(all_pseudo)
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
              values = genes, mart= mart)
