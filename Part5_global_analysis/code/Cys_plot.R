
rm(list=ls())
set.seed(1234)

library("scales")
library("gplots")
library("RColorBrewer")
library("pheatmap")
library("data.table")
library("beeswarm")
library("gtools")


load("Rdat/unique_ID_Ev87_all.Rdat")
load("Rdat/ABE_annot.Rdat")
load("Rdat/evoCDA_annot.Rdat")
load("Rdat/combined_exp_data.Rdat")
load("Rdat/combined_con_data.Rdat")
load("Rdat/df_scout_20230411.Rdat")
load("Rdat/df_denature_20220813.Rdat")
load("Rdat/ceres_used.Rdat")
load("Rdat/AF_RSA.Rdat")
load("Rdat/Cys_conserve_score.Rdat")
load("Rdat/Cys_dropout.Rdat")
load("Rdat/null_res_list.Rdat")
load("Rdat/dep_summary_list.Rdat")


source("code/Cys_functions.R")

###############################
# Define those with significant dropouts.
Cys_res <- filter_dropout_data(Cys_dropout, dep_score=ceres_used, sel_filter=T)
Cys_res_NF <- filter_dropout_data(Cys_dropout, dep_score=ceres_used, sel_filter=F)
# Add EGFR_C797 ligandability manually.
Cys_res[Cys_res$unique_ID=='EGFR_797', 'KB_engage'] <- 100


## PDF 3.5*6
QC_sel_hits(Cys_res, Cys_res_NF, cell_sel='PC14', dep_score=ceres_used)
QC_sel_hits(Cys_res, Cys_res_NF, cell_sel='KMS26', dep_score=ceres_used)


# Compare the hit rate in the actual data vs the simulated non-targeting control data.
# PDF 3.5*6
compare_hit_rate(Cys_res, null_res_list, dep_score=ceres_used)

######
# Show common essentials are essential in both cell lines. 
# save PDF 5.5*6
plot_common_sig(Cys_res, Cys_sel=Cys_res$A_sig_sel, x='PC_mean_LFC_A', y='KMS_mean_LFC_A', sig_col="#66A61E")
plot_common_sig(Cys_res, Cys_sel=Cys_res$C_sig_sel, x='PC_mean_LFC_C', y='KMS_mean_LFC_C', sig_col="#E6AB02")

#########
# PDF 5.5*6
# Compare ABE vs CBE in common essential genes where we have sgRNA.
CE_sel <- Cys_res$cell=='both' & (Cys_res$A_guide_num>0) & (Cys_res$C_guide_num>0)
plot(Cys_res$mean_LFC_A[CE_sel], Cys_res$mean_LFC_C[CE_sel], pch=21, 
     bg=alpha('black', 0.3), col=NULL)
cor.test(Cys_res$mean_LFC_A[CE_sel], Cys_res$mean_LFC_C[CE_sel])



### Plot dep in PC14 and dep in KMS26 and common essentials. PDF 5.5 * 6
plot_dep_class(dep_score=ceres_used, all_covered_dep=unique(Cys_dropout$HUGO), dep_cutoff=-0.4)


### Pie chart for the dependency type.
S_dep <- unique(Cys_res$HUGO[Cys_res$cell!='both'])
CE_dep <- unique(Cys_res$HUGO[Cys_res$cell=='both'])
dep_class_vec <- c(length(S_dep), length(CE_dep))
names(dep_class_vec) <- c('Selective', 'Common Essential')
make_part_pie(dep_class_vec, brewer_col_name="Greys")


## Show sgRNA per Cys for ABE/evoCDA library 
## One for all Cys in dependencies.
## The other for all lig Cys in dependencies.
dep_sel_1 <- unique_ID_Ev87_all$HUGO%in%unique(Cys_res$HUGO)
unique_ID_1 <- unique_ID_Ev87_all[dep_sel_1, c('HUGO', 'unique_ID_Ev87')]
dep_sel_2 <- (df_scout$HUGO%in%unique(Cys_res$HUGO))&(df_scout$KB_engage>50)
unique_ID_2 <- df_scout[dep_sel_2, c('HUGO', 'unique_ID_Ev87')]
## PDF 3*6
plot_stack_bar_sgRNA(unique_ID_input=unique_ID_1, Cys_res)
plot_stack_bar_sgRNA(unique_ID_input=unique_ID_2, Cys_res)


## Scatter plot to show ligandable Cys with dropouts.
## save as PDF 6*6.5
plot_lig_essential_Cys(Cys_res)


## Use heatmap to show KB02/05 liganding.
## PDF 3*5.8
HUGO_sel <- 'CHUK'
data_1 <- extract_data_for_plotting(Cys_res, HUGO_sel, non_NA_name='KB_engage')
plot_lig_heatmap(data_1)
plot_dropout_lollipop(data_1)
data_1

## PDF 3*4.2
HUGO_sel <- 'ADSL'
data_1 <- extract_data_for_plotting(Cys_res, HUGO_sel, non_NA_name='KB_engage')
plot_lig_heatmap(data_1)
plot_dropout_lollipop(data_1)
data_1

## PDF 3*5.2
HUGO_sel <- 'METAP1'
data_1 <- extract_data_for_plotting(Cys_res, HUGO_sel, non_NA_name='KB_engage')
plot_lig_heatmap(data_1)
plot_dropout_lollipop(data_1)
data_1

## PDF 3*5.2
HUGO_sel <- 'DNMT1'
data_1 <- extract_data_for_plotting(Cys_res, HUGO_sel, non_NA_name='KB_engage')
data_2 <- data_1[data_1$target_Cys>1200, ]
plot_lig_heatmap(data_2)
plot_dropout_lollipop(data_2)
data_2


# Use stacked barplots side by side to show essential Cys tend to be structural.
# PDF 4*6
plot_stack_bar_IADTB(Cys_res, df_denature, df_scout)


## Use a scatter plot to show top hits that are reactive / unreactive
## save as PDF 6*6.5
plot_reactive_essential_Cys(Cys_res)


# Examples of essential Cys with different reactivities (solvent accessibility)
# PDF 3*4
HUGO_sel <- 'RBM22'
data_1 <- extract_data_for_plotting(Cys_res, HUGO_sel, non_NA_name='denature_LFC')
plot_denature_heatmap(data_1)
plot_dropout_lollipop(data_1)

# PDF 3*3
HUGO_sel <- 'UTP15'
data_1 <- extract_data_for_plotting(Cys_res, HUGO_sel, non_NA_name='denature_LFC')
plot_dropout_lollipop(data_1)

# PDF 3*3
HUGO_sel <- 'RAN'
data_1 <- extract_data_for_plotting(Cys_res, HUGO_sel, non_NA_name='denature_LFC')
plot_dropout_lollipop(data_1)

# PDF 3*4
HUGO_sel <- 'EEF1A1'
data_1 <- extract_data_for_plotting(Cys_res, HUGO_sel, non_NA_name='denature_LFC')
plot_denature_heatmap(data_1)
plot_dropout_lollipop(data_1)

# PDF 3*4.7
HUGO_sel <- 'UTP18'
data_1 <- extract_data_for_plotting(Cys_res, HUGO_sel, non_NA_name='denature_LFC')
plot_denature_heatmap(data_1)
plot_dropout_lollipop(data_1)

# PDF 3*4
HUGO_sel <- 'RRM1'
data_1 <- extract_data_for_plotting(Cys_res, HUGO_sel, non_NA_name='denature_LFC')
plot_denature_heatmap(data_1)
plot_dropout_lollipop(data_1)

# PDF 3*4.7
HUGO_sel <- 'BRIP1'
data_1 <- extract_data_for_plotting(Cys_res, HUGO_sel, non_NA_name='denature_LFC')
plot_denature_heatmap(data_1)
plot_dropout_lollipop(data_1)

# Show AF2 RSA vs LFC changes upon denaturation.
# Save PDF 5.5*6
df_denature <- add_AF_RSA(df_denature, AF_RSA)
plot_RSA_denature_scatter(df_denature)


# Cys conservation score analysis.
# Save PDF 5.5*6
Cys_res$conserve_1 <- Cys_conserve_score[Cys_res$unique_ID, 'score1']
plot_conserve_scatter(Cys_res, col_name='p_A')
plot_conserve_scatter(Cys_res, col_name='p_C')




