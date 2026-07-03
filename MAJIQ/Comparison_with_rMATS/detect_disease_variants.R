#!/bin/R

####---Author: Maria Kontili
####----Date : 22/06/2026
####Subject: Find if common-rmats/majiq detected genes overlap with Muscle-Disease-Gene-Table (Bonne,Rivier)



library(tidyverse)
library(dplyr)
library(purrr)


epissav_dir <- "~/MYO_WORKSPACE/PROJECTS/EPIS_SAV/"
majiq_wsp <- paste0(epissav_dir,"MAJIQ/")
work_dir <- paste0(majiq_wsp,"Scripts/R_dowstream/")

common_genes <- readRDS("Comparison_with_rMATS/Common_genes_found_rmats_majiq_filtered.rds")


###
### Filter by Muscle-disease-Genes
###

gene_tables_dir <- paste0(epissav_dir,"GeneTables_Bonne-Rivier/")

##> Updated Merged_GeneTables because diseases were omitted first time!
mrg_genetables_collapsed <- readr::read_delim(file =
                                                paste0(gene_tables_dir,
                                                       "Merged_GeneTables_collapsed_diseases_03-06-26.txt"),
                                              delim = "\t",
                                              col_names = TRUE,
                                              na = "NA")

genes_found_muscledisease_idx <- map(common_genes, ~which((.x %in% mrg_genetables_collapsed$GeneSymbol)))


genes_found_muscledisease <- map2(.x = genes_found_muscledisease_idx,
                                    .y = common_genes,
                                    ~(.y[.x]) )


print(genes_found_muscledisease)

bind_rows(genes_found_muscledisease)
# > genes_found_muscledisease
# $C002I3K
# [1] "ALDH3A2"
# 
# $C002I3L
# [1] "NEB"
# 
# $C002I3M
# [1] "TOR1AIP1"  "TPM1"      "SVIL"      "TPM2"      "DTNA"      "ITGA7"     "TPM3"      "NEXN"      "MACF1"     "ANK2"      "FLNC"      "MICU1"    
# [13] "BIN1"      "ABCD3"     "DDHD1"     "SYNE1"     "COL6A3"    "DNM2"      "CAPN3"     "FHOD3"     "PMP22"     "KIDINS220" "DMD"       "PFKM"     
# 
# $C002I3N
# character(0)
# 
# $C002I3P
# [1] "TTN"
# 
# $C002I3Q
# character(0)
# 
# $C002I3R
# [1] "DMD"     "ALDH3A2"
# 
# $C002I3S
# [1] "DPAGT1"  "ALDH3A2"
# 
# $C002I3U
# character(0)
# 
# $C002I3W
# [1] "TOR1AIP1" "DMD"      "ALDH3A2" 
# 
# $C002I3Y
# character(0)
# 
# $C002I3Z
# [1] "TNNT3"   "ALDH3A2"
# 
# $C002I44
# [1] "NEB"     "ALDH3A2"
# 
# $C002I45
# [1] "DMD"   "POMT1"
# 
# $C002I47
# character(0)
# 
# $C002I4F
# [1] "DMD"
# 
# $C002I4G
# [1] "TTN"     "ALDH3A2"

saveRDS(genes_found_muscledisease,"Comparison_with_rMATS/Genes_common_btw_rmats-majiq_and_MuscleGeneTable.rds")
