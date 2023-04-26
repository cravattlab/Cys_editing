
rm(list=ls())
set.seed(1234)
library("xlsx")
library("scales")
library("gplots")
library("RColorBrewer")
library("pheatmap")
library("data.table")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("gtools")


load("Rdat/panel_exp_list.Rdat")
load("Rdat/panel_ABE_annot.Rdat")
load("Rdat/panel_evoCDA_annot.Rdat")
load("Rdat/df_scout_20230411.Rdat")
load("Rdat/ceres_used.Rdat")
load("Rdat/expanded_target_Cys.Rdat")
source("code/panel_functions.R")


HUGO_panel <- unique(expanded_target_Cys$HUGO)

cell_panel <- c("22RV1", "DLD1", "GSS", "KMS26", "KMS34", "MCC142", 
                "MM1S", "PANC1005", "PC14", "SNU216", "SUDHL5", "UACC257")

ceres_panel <- clean_dep_score(dep_score=ceres_used, cell_panel, HUGO_panel)


cor_res_1 <- get_combined_cor(data_list=exp_list, 
                              dep_score_panel=ceres_panel, 
                              HUGO_panel=HUGO_panel, 
                              ABE_annot=ABE_annot, 
                              evoCDA_annot=evoCDA_annot, 
                              expanded_target_Cys=expanded_target_Cys)

filtering_metrics <- get_filtering_metrics(cor_res_1, 
                                           data_list=exp_list, 
                                           dep_score_panel=ceres_panel, 
                                           ABE_annot=ABE_annot, 
                                           evoCDA_annot=evoCDA_annot,
                                           splice_type='splice_dist1')

cor_res_2 <- filter_cor_res(cor_res_1,
                            filtering_metrics=filtering_metrics,
                            BE_cutoff=-0.5, 
                            KO_cutoff=-0.5)

cor_res_3 <- remove_failed(cor_res_2)

cor_res_4 <- get_FDR(cor_res_3)

cor_res_5 <- filter_lig(cor_res_4, df_scout) 

# Make scatter plots for BE and KO
cor_res_A_sig <- cor_res_5[cor_res_5$FDR_A>1 & !is.na(cor_res_5$FDR_A), ]
cor_res_C_sig <- cor_res_5[cor_res_5$FDR_C>1 & !is.na(cor_res_5$FDR_C), ]

plot_BE_KO_bulk(cor_res=cor_res_A_sig, exp_list=exp_list, dep_score_panel=ceres_panel,
                lib='ABE', output_folder='plot')
plot_BE_KO_bulk(cor_res=cor_res_C_sig, exp_list=exp_list, dep_score_panel=ceres_panel,
                lib='evoCDA', output_folder='plot')


