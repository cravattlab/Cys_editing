

rm(list=ls())

set.seed(123)

library("data.table")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("gtools")
library("beeswarm")
library("RColorBrewer")
library("scales")

load("Rdat/ABE_annot.Rdat")
load("Rdat/evoCDA_annot.Rdat")
load("Rdat/combined_exp_data.Rdat")
load("Rdat/combined_con_data.Rdat")
load("Rdat/df_scout_20230411.Rdat")
load("Rdat/df_denature_20220813.Rdat")
load("Rdat/unique_ID_Ev87_all.Rdat")
load("Rdat/AF_RSA.Rdat")
load("Rdat/ceres_used.Rdat")
load("Rdat/HUGO_dep_list.Rdat")
load("Rdat/dep_summary_list.Rdat")


# source major functions
source("code/Cys_functions.R")

#Remove sgRNAs that are close to slicing site or causes stop codon.
data_list_non_splice_dist3 <- extract_non_splice_data(con_data, exp_data,
                                                      cell_line_sel=c('PC14', 'KMS26'),
                                                      ABE_annot, evoCDA_annot,
                                                      splice_dist='splice_dist3')

LFC_input_pan <- get_LFC_input(data_list_non_splice_dist3, cell_names='PC|KMS')
LFC_input_PC <- get_LFC_input(data_list_non_splice_dist3, cell_names='PC')
LFC_input_KMS <- get_LFC_input(data_list_non_splice_dist3, cell_names='KMS')

# Get the null distribution based on non-targeting sgRNAs.
null_list_pan <- get_null_list(LFC_input_pan, n_sim=1e5)
null_list_PC <- get_null_list(LFC_input_PC, n_sim=1e5)
null_list_KMS <- get_null_list(LFC_input_KMS, n_sim=1e5)

save(null_list_pan, file='null_list_pan.Rdat')
save(null_list_PC, file='null_list_PC.Rdat')
save(null_list_KMS, file='null_list_KMS.Rdat')


Cys_dropout <- run_Cys_pipeline(LFC_input_pan, null_list_pan,
                                  LFC_input_PC, null_list_PC,
                                  LFC_input_KMS, null_list_KMS,
                                  HUGO_dep_list, AF_RSA,
                                  df_scout, df_denature)

save(Cys_dropout, file='Cys_dropout.Rdat')


null_res_list <- simulate_NT_pos_rate(LFC_input_pan, null_list_pan,
                                 LFC_input_PC, null_list_PC,
                                 LFC_input_KMS, null_list_KMS,
                                 HUGO_dep_list, AF_RSA, df_scout,
                                 df_denature, nsim=10)

save(null_res_list, file='null_res_list.Rdat')













