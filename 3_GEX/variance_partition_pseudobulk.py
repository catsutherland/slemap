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
pseudo_haerin<- readRDS("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/hj10/Results/3_eQTL_prep/pseudobulk/meanPB_level0.RDS")

# %%
pseudo_haerin

# %%
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <-  row.names(all_pseudo)
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
              values = genes, mart= mart)

# %%
gene_IDs_matched <- gene_IDs[!gene_IDs$hgnc_symbol=="", ]

# %%
nrow(gene_IDs_matched)

# %%
length(unique(gene_IDs$hgnc_symbol))

# %%
length(unique(gene_IDs$ensembl_gene_id))

# %%
clinical_all <- read.csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/phenotype/combined_data_tidied_sc_used_25_11_24.csv")

# %%
metadata <- read.csv("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/metadata/SLEmap_metadata_191224.csv")

# %%
sequencing_info <- read.csv("/nfs/users/nfs_c/cs54/OTAR2064_cs54/final_data/1_Metadata/sequencing_info_per_slemap_10x_pool.txt", sep = "\t")

# %%
metadata <- merge(metadata, sequencing_info, by.x = "poolID", by.y = "PoolID")

# %%
metadata$date_processed

# %%
clinical_all <- merge(metadata[, c("externalID","experiment_id", "poolID", 
                                   "Flow_cell", "id_pool_lims", "date_processed")],
                      clinical_all, by = c("externalID"))

# %%
clinical_all$experiment_id <- gsub("_", "-", clinical_all$experiment_id)

# %%
rownames(clinical_all) <- clinical_all$experiment_id

# %%
clinical_all <- clinical_all[colnames(all_pseudo), ]

# %%
all_pseudo <- as.data.frame(as.matrix(all_pseudo))

# %%
pseudo_symbols <- merge(all_pseudo, gene_IDs_matched, by.x = 0, by.y = "ensembl_gene_id")

# %%
length(unique(pseudo_symbols$hgnc_symbol))

# %%
pseudo_symbols <- pseudo_symbols[!pseudo_symbols$hgnc_symbol=='GOLGA8M', ]

# %%
# Replace the row names of df1 with the corresponding gene names
rownames(pseudo_symbols) <- pseudo_symbols$hgnc_symbol

# %%
pseudo_symbols <- subset(pseudo_symbols, select = -c(Row.names, hgnc_symbol))


# %%
clinical_all$Flow_cell<-as.character(clinical_all$Flow_cell)


# %%
form <- ~ (1|externalID) + (1|poolID)+ (1|id_pool_lims) + (1|Flow_cell) + (1|Recruit_centre) + Age_at_sampling

# Compute Canonical Correlation Analysis (CCA)
# between all pairs of variables
# returns absolute correlation value
C <- canCorPairs(form, clinical_all)

# Plot correlation matrix
# between all pairs of variables
plotCorrMatrix(C)

# %%
length(unique(clinical_all$Flow_cell))

# %%
#### 1. Run with just first timepoint and Healthy - for age, ancestry, sex etc. ####

form <- ~ (1|poolID) + (1|id_pool_lims) + (1|Flow_cell) + (1|Recruit_centre) + Age_at_sampling
varPart <- fitExtractVarPartModel(all_pseudo, form, clinical_all)
vp <- sortCols( varPart )
plotPercentBars( vp[1:10,] )
plotVarPart( vp)

# %%
#### 1. Run with just first timepoint and Healthy - for age, ancestry, sex etc. ####

form <- ~ (1|id_pool_lims) + (1|Flow_cell) + (1|Recruit_centre) + Age_at_sampling
varPart <- fitExtractVarPartModel(all_pseudo, form, clinical_all)
vp <- sortCols( varPart )
plotPercentBars( vp[1:10,] )
plotVarPart( vp)

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
length(colnames(pseudo_l1))

# %%
expanded_clinical_all <- expanded_clinical_all[colnames(pseudo_l1), ]

# %%
nrow(expanded_clinical_all)

# %%
pseudo_l1 <- as.data.frame(as.matrix(pseudo_l1))

# %%
expanded_clinical_all

# %%
expanded_clinical_all

# %%
metadata

# %%
form <- ~ (1|externalID) + (1|celltype) + (1|poolID) + (1|Flow_cell) + (1|id_pool_lims) + (1|Recruit_centre) + Age_at_sampling
varPart <- fitExtractVarPartModel(pseudo_l1, form, expanded_clinical_all)
vp <- sortCols( varPart )
plotPercentBars( vp[1:10,] )
plotVarPart( vp)

# %%
for (celltype in c("_B", "_CD4 T", "_CD8 T", "_other T", "_NK", "_DC", "_Mono", "_other")) {
    print(gsub("_", "", celltype))
    pseudo_celltype <- pseudo_l1[, grep(celltype, names(pseudo_l1))]
    clinical_celltype <- expanded_clinical_all[colnames(pseudo_celltype), ]
    form <- ~ (1|poolID) + (1|Recruit_centre) + Age_at_sampling
    varPart <- fitExtractVarPartModel(pseudo_celltype, form, clinical_celltype)
    vp <- sortCols( varPart )
    p1 <- plotPercentBars( vp[1:10,] )
    p2 <- plotVarPart( vp)
    print(p1)
    print(p2)
    ggsave(glue("{figures_folder}/variance_violin{celltype}.pdf", plot = p2))
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
    ggsave(glue("{figures_folder}/variance_violin{celltype}.pdf", plot = p2))
}

# %%
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <-  rownames(all_pseudo)
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
              values = genes, mart= mart)
