

rm(list=ls())

library('mgcv')
library("scales")
library("gplots")
library('gtools')
library("RColorBrewer")
library("pheatmap")
library("data.table")
library("bio3d")
library("Biostrings")

load("Rdat/Fig1_ABE_annot.Rdat")
load("Rdat/Fig1_evoCDA_annot.Rdat")
load("Rdat/Fig1_exp_data.Rdat")
load("Rdat/Fig1_con_data.Rdat")
load("Rdat/Fig1_combined_seq_df.Rdat")

source("code/Fig1_functions.R")

## Remove sgRNAs that are close to splicing site or causes stop codon.
data_list_non_splice_dist3 <- extract_non_splice_data(con_data, exp_data,
                                                      cell_line_sel=c('PC14', 'KMS26'),
                                                      ABE_annot, evoCDA_annot,
                                                      splice_dist='splice_dist3')

LFC_input_both <- get_LFC_input(data_list_non_splice_dist3, cell_names='PC|KMS')
LFC_input_PC <- get_LFC_input(data_list_non_splice_dist3, cell_names='PC')
LFC_input_KMS <- get_LFC_input(data_list_non_splice_dist3, cell_names='KMS')

res_PC_EGFR <- get_saturated_result(HUGO_sel='EGFR', LFC_input=LFC_input_PC, Fig1_combined_seq_df)
res_KMS_EGFR <- get_saturated_result(HUGO_sel='EGFR', LFC_input=LFC_input_KMS, Fig1_combined_seq_df)
res_KMS_XPO1 <- get_saturated_result(HUGO_sel='XPO1', LFC_input=LFC_input_KMS, Fig1_combined_seq_df)

save(res_PC_EGFR, file='res_PC_EGFR.Rdat')
save(res_KMS_EGFR, file='res_KMS_EGFR.Rdat')
save(res_KMS_XPO1, file='res_KMS_XPO1.Rdat')


###### EGFR in PC14  PDF 7*6
plot_residue_waterfall(res_PC_EGFR$unique_ID, LFC=res_PC_EGFR$LFC_AC, p=res_PC_EGFR$p_AC, 
                       LFC_cut=-0.5, pcut=-1.3)
abline(v=797, col='red')

plot_residue_waterfall(res_PC_EGFR$unique_ID, LFC=res_PC_EGFR$LFC_A, p=res_PC_EGFR$p_A, 
                       LFC_cut=-0.5, pcut=-1.3)
abline(v=797, col='red')

plot_residue_waterfall(res_PC_EGFR$unique_ID, LFC=res_PC_EGFR$LFC_C, p=res_PC_EGFR$p_C, 
                       LFC_cut=-0.5, pcut=-1.3)
abline(v=797, col='red')

#fwrite(res_PC_EGFR, file='res_PC_EGFR.csv')


###### EGFR in KMS26  PDF 7*6
plot_residue_waterfall(res_KMS_EGFR$unique_ID, LFC=res_KMS_EGFR$LFC_AC, p=res_KMS_EGFR$p_AC,
                       LFC_cut=-0.5, pcut=-1.3)
abline(v=797, col='red')

plot_residue_waterfall(res_KMS_EGFR$unique_ID, LFC=res_KMS_EGFR$LFC_A, p=res_KMS_EGFR$p_A,
                       LFC_cut=-0.5, pcut=-1.3)
abline(v=797, col='red')

plot_residue_waterfall(res_KMS_EGFR$unique_ID, LFC=res_KMS_EGFR$LFC_C, p=res_KMS_EGFR$p_C,
                       LFC_cut=-0.5, pcut=-1.3)
abline(v=797, col='red')


###### XPO1  PDF 6.5*6
plot_residue_waterfall(res_KMS_XPO1$unique_ID, LFC=res_KMS_XPO1$LFC_AC, p=res_KMS_XPO1$p_AC,
                       LFC_cut=-0.5, pcut=-1.3)
abline(v=528, col='red')

#fwrite(res_KMS_XPO1, file='res_KMS_XPO1.csv')


####################################################
# Show that only PC14 has selective Cys dropouts but not KMS26.
# PDF 4*6
PC14_vs_KMS26 <- compare_2_Cys(LFC_input1=res_PC_EGFR$LFC_AC, LFC_input2=res_KMS_EGFR$LFC_AC, 
                               p_input1=res_PC_EGFR$p_AC, p_input2=res_KMS_EGFR$p_AC, 
                               unique_ID=res_PC_EGFR$unique_ID, LFC_cut=-0.5)

####################################################
# coef_vec is a 1 column matrix, each row corresponds to an ensembl_aa.
coef_vec_EGFR <- get_coef_vec_from_p(res_PC_EGFR$p_AC, cutoff=-1.3)
coef_vec_XPO1 <- get_coef_vec_from_p(res_KMS_XPO1$p_AC, cutoff=-1.3)

HUGO <- 'EGFR'
PDB_ID <- "6JXT" #chain A EGFR 696-1022 WT in complex with AZD9291
get_clean_PDB(PDB_ID=PDB_ID, HUGO=HUGO, combined_seq_df=Fig1_combined_seq_df, 
              chain='A', coef_vec=coef_vec_EGFR, output_folder='data')

HUGO <- 'XPO1'
PDB_ID <- "6TVO"  #chain A, Human CRM1 (XPO1) -RanGTP in complex with Leptomycin B
get_clean_PDB(PDB_ID=PDB_ID, HUGO=HUGO, combined_seq_df=Fig1_combined_seq_df, 
              chain='A', coef_vec=coef_vec_XPO1, output_folder='data')

 

