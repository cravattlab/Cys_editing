
clean_dep_score <- function(dep_score, cell_panel, HUGO_panel){
  # Some genes have newer names in dep_score.
  # HUGO_panel[!HUGO_panel%in%colnames(dep_score)]
  # If needed, update these gene symbols.
  colnames(dep_score) <- gsub(c("CIP2A"), c("KIAA1524"), colnames(dep_score))
  colnames(dep_score) <- gsub(c("WASHC5"), c("KIAA0196"), colnames(dep_score))
  dep_score_sel <- dep_score[cell_panel, HUGO_panel]
  return(dep_score_sel)
}


# Calculate the correlations between KO data and gene-associated sgRNA base editing data.
get_BE_KO_cor <- function(sgRNA_data, sgRNA_sel, ref_vec){
  
  col_sel <- colnames(sgRNA_data) %in% names(ref_vec)
  data_sel <- t(sgRNA_data[sgRNA_sel, col_sel, drop=F])
  n <- dim(data_sel)[2]
  p_vec <- rep(0, n)
  cor_vec <- rep(0, n)

  for (i in 1:n){
    p_temp <- cor.test(data_sel[, i], ref_vec[rownames(data_sel)], method='pearson', alternative='greater')$p.value
    p_vec[i] <- round(-log10(p_temp), 1)
    cor_temp <- cor.test(data_sel[, i], ref_vec[rownames(data_sel)], method='pearson', alternative='greater')$estimate
    cor_vec[i] <- round(cor_temp, 2)
  }

  cor_res <- list()
  cor_res$cor <- cor_vec
  cor_res$p <- p_vec
  return(cor_res)
}

get_res_per_HUGO <- function(data_list, dep_score_panel, HUGO_current, ABE_annot, evoCDA_annot, 
                             expanded_target_Cys){
  
  all_sgRNA <- rownames(ABE_annot)
  annot_HUGO <- ABE_annot$HUGO
  sgRNA_sel <- all_sgRNA[annot_HUGO == HUGO_current]
  
  # Get the KO dependency values in order. repeat for 2 batches.
  ref_vec <- dep_score_panel[, HUGO_current]
  
  cor_res_A <- get_BE_KO_cor(data_list$data_A, sgRNA_sel, ref_vec)
  cor_res_C <- get_BE_KO_cor(data_list$data_C, sgRNA_sel, ref_vec)

  aa_channge_A <- (ABE_annot$aa_change)[annot_HUGO == HUGO_current]
  aa_channge_C <- (evoCDA_annot$aa_change)[annot_HUGO == HUGO_current]
  
  group_A <- assign_Cys_region(aa_channge_A, expanded_target_Cys, HUGO_current)
  group_C <- assign_Cys_region(aa_channge_C, expanded_target_Cys, HUGO_current)
  
  group <- group_A
  group[is.na(group)] <- group_C[is.na(group)]
  
  result <- data.frame(HUGO=HUGO_current, 
                       aa_channge_A=aa_channge_A, aa_channge_C=aa_channge_C, 
                       group_A=group_A, group_C=group_C, group=group,
                       cor_A=cor_res_A$cor, cor_C=cor_res_C$cor, 
                       p_A=cor_res_A$p, p_C=cor_res_C$p)

  # Add sgRNA without PAM as rownames
  rownames(result) <- sgRNA_sel
  return(result)
}

assign_Cys_region <- function(raw_annot, expanded_target_Cys, HUGO_current){
  all_letters <- toupper(letters)
  n <- length(raw_annot)
  group_vec <- rep(NA, n)
  
  for (i in 1:n){
    mutation_split <- sapply(raw_annot[i], strsplit, split='_')[[1]]
    mutation_split <- mutation_split[!(mutation_split%in%all_letters)]
    mutation_split <- mutation_split[!is.na(mutation_split)]
    if (length(mutation_split)>0){
      ID_temp <- paste0(HUGO_current, '_', mutation_split)
      match_res <- match(ID_temp, expanded_target_Cys$ID)
      match_res <- match_res[!is.na(match_res)][1]
      group_vec[i] <- expanded_target_Cys$center[match_res]
    }
  }
  
  return(group_vec)
} 

# Get the combined df for different BE-KO correlations.
get_combined_cor <- function(data_list, dep_score_panel, HUGO_panel, ABE_annot, evoCDA_annot, 
                             expanded_target_Cys){

  cor_combined <- data.frame()
  for(HUGO_current in HUGO_panel){
    print(HUGO_current)
    res_temp <- get_res_per_HUGO(data_list, dep_score_panel, HUGO_current, ABE_annot, evoCDA_annot, 
                                 expanded_target_Cys)
    cor_combined <- rbind(cor_combined, res_temp)
  }
  return(cor_combined)
}


get_filtering_metrics <- function(cor_res, data_list, dep_score_panel,
                                  ABE_annot, evoCDA_annot, splice_type='splice_dist1'){
  
  sgRNA_vec <- rownames(cor_res)
  max_A <- apply(data_list$data_A_ave[sgRNA_vec, ], 1, max)
  min_A <- apply(data_list$data_A_ave[sgRNA_vec, ], 1, min)
  
  max_C <- apply(data_list$data_C_ave[sgRNA_vec, ], 1, max)
  min_C <- apply(data_list$data_C_ave[sgRNA_vec, ], 1, min)
  
  max_KO <- apply(dep_score_panel[, cor_res$HUGO], 2, max)
  min_KO <- apply(dep_score_panel[, cor_res$HUGO], 2, min)
  
  # Either to stop codon or from stop codon.
  A_stop <- ABE_annot[sgRNA_vec, 'stop']
  C_stop <- evoCDA_annot[sgRNA_vec, 'stop']
  
  splice_A <- ABE_annot[sgRNA_vec, splice_type]
  splice_C <- evoCDA_annot[sgRNA_vec, splice_type]
  
  filtering_df <- data.frame(max_A, min_A, max_C, min_C, max_KO, min_KO, A_stop, C_stop, splice_A, splice_C)
  
  return(filtering_df)
}

# We only care about BE-KO correlations for the sgRNAs that led to dropouts.
# Remove sgRNAs that don't cause dropouts.
# Remove sgRNAs that may lead to stop codons or target splicing sites.
# Remove sgRNAs whose associated genes don't have a dependent cell line in the panel.
# Remove common essential genes from correlation analysis.
filter_cor_res <- function(cor_res, filtering_metrics, cor_cutoff=0.5, 
                           BE_cutoff=-0.5, KO_cutoff=-0.5){
  filtering_metrics <- filtering_metrics[rownames(cor_res), ]
  range_sel <- (filtering_metrics$max_KO - filtering_metrics$min_KO) > 0.4
  KO_sel <- range_sel & (filtering_metrics$min_KO <= KO_cutoff) & (filtering_metrics$max_KO>=-0.2)

  A_sel <- KO_sel & (filtering_metrics$min_A <= BE_cutoff) & 
    !(filtering_metrics$A_stop) & 
    !(filtering_metrics$splice_A) & 
    cor_res$cor_A >= cor_cutoff
  C_sel <- KO_sel & (filtering_metrics$min_C <= BE_cutoff) &
    !(filtering_metrics$C_stop) &
    !(filtering_metrics$splice_C) &
    cor_res$cor_C >= cor_cutoff

  cor_res$A_sel <- A_sel
  cor_res$C_sel <- C_sel
  
  cor_res <- cbind(cor_res, filtering_metrics)
  return(cor_res)
}

plot_BE_KO_bulk <- function(cor_res, exp_list, dep_score_panel, lib, output_folder){
  
  folder_path <- file.path(output_folder, lib)
  dir.create(folder_path)
  n <- dim(cor_res)[1]
  if (lib=='ABE'){
    BE_data <- exp_list$data_A
  } else {
    BE_data <- exp_list$data_C
  }
  
  for (i in 1:n){
    sg <- rownames(cor_res)[i]
    group <- cor_res$group[i]
    HUGO <- cor_res$HUGO[i]
    x <- BE_data[sg, ]
    y <- dep_score_panel[names(x), HUGO]
    d_x <- (max(x)-min(x))/20
    d_y <- (max(y)-min(y))/20
    rho <- round(cor.test(x, y, alternative=c("greater"))$estimate, 2)
    p <- round(-log10(cor.test(x, y, alternative=c("greater"))$p.value), 3)
    title <- sprintf("%s, r: %s, p: %s, lib: %s, \n sgRNA: %s", group, rho, p, lib, sg)
    x1 <- matrix(x, ncol=2)[, 1]
    x2 <- matrix(x, ncol=2)[, 2]
    y1 <- matrix(y, ncol=2)[, 1]
    x_ave <- (x1+x2)/2
    
    col_set1 <- alpha(brewer.pal(3, name = 'Set1'), 0.7)
    
    # Save PDF files in the same folder.
    pdf(file=file.path(folder_path, paste0(title, '.pdf')), width=4.5, height=5)
    plot(x1, y1, pch=21, bg=col_set1[1], cex=1.5, 
         xlab='Log2 fold change (base editing)',
         ylab='Gene-level dependency',
         xlim=c(min(x)-d_x, max(x)+d_x), 
         ylim=c(min(y)-d_y, max(y)+d_y),
         main=title, font.main=1)
    points(x2, y1, pch=21, bg=col_set1[2], cex=1.5)
    segments(x1, y1, x2, y1, col='black')
    text(x_ave, y1+d_y, names(y)[1:(length(y)/2)])
    dev.off()
  }
}

remove_failed <- function(cor_res){
  cor_res$cor_A[!cor_res$A_sel] <- NA
  cor_res$cor_C[!cor_res$C_sel] <- NA
  cor_res$p_A[!cor_res$A_sel] <- NA
  cor_res$p_C[!cor_res$C_sel] <- NA
  return(cor_res)
}

get_FDR <- function(cor_res){

  unique_HUGO <- unique(cor_res$HUGO)
  cor_res$FDR_A <- NA
  cor_res$FDR_C <- NA
  for (i in 1:length(unique_HUGO)){
    HUGO_sel <- cor_res$HUGO==unique_HUGO[i]
    cor_res_temp <- cor_res[HUGO_sel, ]
    cor_res$FDR_A[HUGO_sel] <- -log10(p.adjust(10^-cor_res_temp$p_A, method='fdr'))
    cor_res$FDR_C[HUGO_sel] <- -log10(p.adjust(10^-cor_res_temp$p_C, method='fdr'))
  }
  
  return(cor_res)
}


get_hit_metrics <- function(cor_res, FDR_cutoff=1){
  A_hit <- tapply(cor_res$FDR_A>=FDR_cutoff, cor_res$group, sum, na.rm=T)
  C_hit <- tapply(cor_res$FDR_C>=FDR_cutoff, cor_res$group, sum, na.rm=T)
  total_hit <- A_hit + C_hit
  HUGO_vec <- cor_res$HUGO[match(names(total_hit), cor_res$group)]
  hit_df <- data.frame(HUGO_vec, group=names(total_hit), A_hit, C_hit, total_hit)
  hit_df <- hit_df[order(-hit_df$total_hit), ]
  return(hit_df)
}


filter_lig <- function(cor_res, df_scout){
  lig_sel <- df_scout$KB_engage>50 & !is.na(df_scout$KB_engage)
  lig_HUGO <- unique(c(df_scout$HUGO[lig_sel], 'EGFR'))
  lig_sites <- df_scout$unique_ID_Ev87[lig_sel]
  lig_sites <- unique(c(gsub('_', '_C', lig_sites), "EGFR_C797"))
  clean_res <- cor_res[cor_res$group%in%lig_sites, ]
  return(clean_res)
}


barplot_top_hits <- function(hit_df, top_N=10){
  col_vec <- brewer.pal(n=6, name='Set2')
  
  hit_df_1 <- hit_df[order(-hit_df$total_hit), ][1:top_N, ]
  hit_df_2 <- t(hit_df_1[, c('A_hit','C_hit','total_hit')])
  hit_df_3 <- hit_df_2[1:2, order(hit_df_2['total_hit', ])]
  
  barplot(hit_df_3, horiz=T, las=1,
          ylab=NULL, cex.names=1.2,
          xlab='Number of base editing events that correlate with KO',
          col=alpha(col_vec[5:6], 0.7))
  abline(v=0)
}


plot_dep_overview <- function(dep_score_panel){
  breaksList <- seq(-0.8, 0, by = 0.01)
  col_palette <- c(rev(brewer.pal(n=7, name="Blues")), 'white')
  color <- colorRampPalette(col_palette)(length(breaksList))
  pheatmap(dep_score_panel, color=color, border='white', breaks=breaksList, fontsize=12)
}



