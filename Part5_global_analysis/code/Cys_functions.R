
convert_HUGO_to_name <- function(HUGO_vec) {
  gene_names <- mapIds(org.Hs.eg.db,
                       keys=HUGO_vec,
                       column="GENENAME",
                       keytype="SYMBOL",
                       multiVals="first")
  gene_names <- gsub(',',' ', gene_names)
  
  return(gene_names)
}



extract_non_splice_data <- function(con_data, exp_data, cell_line_sel, ABE_annot, evoCDA_annot, splice_dist){
  con_data <- con_data[, grep(paste(cell_line_sel, collapse='|'), colnames(con_data))]
  exp_data <- exp_data[, grep(paste(cell_line_sel, collapse='|'), colnames(exp_data))]
  splicing_stop_sel_A <- ABE_annot[, splice_dist]|ABE_annot$stop
  splicing_stop_sel_C <- evoCDA_annot[, splice_dist]|evoCDA_annot$stop
  exp_data[splicing_stop_sel_A, grepl('libA', colnames(exp_data))] <- NA
  exp_data[splicing_stop_sel_C, grepl('libC', colnames(exp_data))] <- NA
  
  keep <- (!splicing_stop_sel_A) | (!splicing_stop_sel_C)
  data_list <- list()
  data_list$con_data <- con_data
  data_list$exp_data <- exp_data[keep, ]
  data_list$ABE_annot <- ABE_annot[keep, ]
  data_list$evoCDA_annot <- evoCDA_annot[keep, ]
  
  return(data_list)
}


get_LFC_input <- function(data_list, cell_names='PC|KMS'){
  con_data <- data_list$con_data
  exp_data <- data_list$exp_data
  ABE_annot <- data_list$ABE_annot
  evoCDA_annot <- data_list$evoCDA_annot
  
  # get the distribution from reference guides
  con_A_sel <- grepl('libA', colnames(con_data)) & grepl(cell_names, colnames(con_data))
  con_C_sel <- grepl('libC', colnames(con_data)) & grepl(cell_names, colnames(con_data))
  exp_A_sel <- grepl('libA', colnames(exp_data)) & grepl(cell_names, colnames(exp_data))
  exp_C_sel <- grepl('libC', colnames(exp_data)) & grepl(cell_names, colnames(exp_data))
  
  NT_A <- apply(con_data[, con_A_sel, drop=F], 1, mean, na.rm=T)
  NT_C <- apply(con_data[, con_C_sel, drop=F], 1, mean, na.rm=T)
  NT_AC <- c(NT_A, NT_C)
  
  contrast_A <- apply(exp_data[, exp_A_sel, drop=F], 1, mean, na.rm=T)
  contrast_C <- apply(exp_data[, exp_C_sel, drop=F], 1, mean, na.rm=T)
  
  ABE_ID <- paste(ABE_annot$HUGO, ABE_annot$target_aa, sep='_')
  evoCDA_ID <- paste(evoCDA_annot$HUGO, evoCDA_annot$target_aa, sep='_')
  unique_ID <- unique(c(ABE_ID, evoCDA_ID))
  # Remove unwanted Cys (included those that target two genes)
  unique_ID <- mixedsort(unique_ID[!grepl('_NA', unique_ID)])
  
  LFC_list <- list()
  LFC_list$ABE_ID <- ABE_ID
  LFC_list$evoCDA_ID <- evoCDA_ID
  LFC_list$unique_ID <- unique_ID
  
  LFC_list$contrast_A <- contrast_A
  LFC_list$contrast_C <- contrast_C
  LFC_list$NT_A <- NT_A
  LFC_list$NT_C <- NT_C
  LFC_list$NT_AC <- NT_AC
  
  LFC_list$ABE_annot <- ABE_annot
  LFC_list$evoCDA_annot <- evoCDA_annot
  
  return(LFC_list)
}


get_null_list <- function(LFC_list, n_sim=1e5, sg_n_max=6){
  
  NT_A <- LFC_list$NT_A
  NT_C <- LFC_list$NT_C
  
  n <- length(NT_A)
  null_list <- list()

  for(sg_n in 1:sg_n_max){
    print(sg_n)
    df <- data.frame(mean_LFC_A=NA, mean_LFC_C=NA, A_guide_num=rep(sg_n, n_sim), 
                     C_guide_num=rep(sg_n, n_sim))
    for(i in 1:n_sim){
      A_lfc_sel <- NT_A[sample(1:n, sg_n, replace=T)]
      C_lfc_sel <- NT_C[sample(1:n, sg_n, replace=T)]
      df$mean_LFC_A[i] <- round(mean(A_lfc_sel), 2)
      df$mean_LFC_C[i] <- round(mean(C_lfc_sel), 2)
    }
    null_list[[sg_n]] <- df
  }
  return(null_list)
}


sim_LFC_input <- function(LFC_list, sim_ID){

  LFC_list$contrast_A <- LFC_list$NT_A[sim_ID]
  LFC_list$contrast_C <- LFC_list$NT_C[sim_ID]
  
  return(LFC_list)
}


get_essential_Cys <- function(LFC_list, null_list){
  
  ABE_ID <- LFC_list$ABE_ID
  evoCDA_ID <- LFC_list$evoCDA_ID
  unique_ID <- LFC_list$unique_ID
  
  contrast_A <- LFC_list$contrast_A
  contrast_C <- LFC_list$contrast_C
  
  A_sgRNA_LFC <- list()
  C_sgRNA_LFC <- list()
  df <- data.frame(unique_ID=unique_ID, mean_LFC_A=NA, mean_LFC_C=NA, mean_LFC_AC=NA, 
                   A_guide_num=NA, C_guide_num=NA, p_A=NA, p_C=NA, p_AC=NA)
  
  for(i in 1:length(unique_ID)){
    current_ID <- unique_ID[i]
    sel <- (ABE_ID%in%current_ID) | (evoCDA_ID%in%current_ID)
    A_lfc_sel <- contrast_A[sel]
    C_lfc_sel <- contrast_C[sel]
    A_lfc_sel <- A_lfc_sel[!is.na(A_lfc_sel)]
    C_lfc_sel <- C_lfc_sel[!is.na(C_lfc_sel)]
    target_lfc <- c(A_lfc_sel, C_lfc_sel)

    A_sgRNA_LFC[[i]] <- round(A_lfc_sel[order(A_lfc_sel)], 2)
    C_sgRNA_LFC[[i]] <- round(C_lfc_sel[order(C_lfc_sel)], 2)

    df$A_guide_num[i] <- length(A_lfc_sel)
    df$C_guide_num[i] <- length(C_lfc_sel)
    
    df$mean_LFC_A[i] <- round(mean(A_lfc_sel), 2)
    df$mean_LFC_C[i] <- round(mean(C_lfc_sel), 2)
    df$mean_LFC_AC[i] <- round((mean(target_lfc)), 2)

    if (length(A_lfc_sel)>0){
      p_A <- mean(null_list[[length(A_lfc_sel)]]$mean_LFC_A <= mean(A_lfc_sel))+1e-6
      df$p_A[i] <- round(-log10(p_A), 2)
    }
    if (length(C_lfc_sel)>0){
      p_C <- mean(null_list[[length(C_lfc_sel)]]$mean_LFC_C <= mean(C_lfc_sel))+1e-6
      df$p_C[i] <- round(-log10(p_C), 2)
    }
  }

  df$A_sgRNA_LFC <- A_sgRNA_LFC
  df$C_sgRNA_LFC <- C_sgRNA_LFC

  return(df)
}


clean_Cys_df <- function(Cys_dropout){
  Cys_dropout$HUGO <- sapply(strsplit(Cys_dropout$unique_ID, '_'), "[[", 1)
  Cys_dropout$target_Cys <- sapply(strsplit(Cys_dropout$unique_ID, '_'), "[[", 2)
  Cys_dropout$full_name <- convert_HUGO_to_name(Cys_dropout$HUGO)
  Cys_dropout$mean_LFC_A[is.na(Cys_dropout$mean_LFC_A)] <- 0
  Cys_dropout$mean_LFC_C[is.na(Cys_dropout$mean_LFC_C)] <- 0
  rownames(Cys_dropout) <- Cys_dropout$unique_ID
  return(Cys_dropout)
}



add_proteomic_data <- function(Cys_dropout_df, df_scout, df_denature){
  ID_match_1 <- match(Cys_dropout_df$unique_ID, df_scout$unique_ID_Ev87)
  Cys_dropout_df$KB_engage <- df_scout$KB_engage[ID_match_1]

  ID_match_2 <- match(Cys_dropout_df$unique_ID, df_denature$unique_ID)
  col_to_keep <- !(colnames(df_denature)%in%c('HUGO', 'unique_ID'))
  Cys_dropout_df <- data.frame(Cys_dropout_df, df_denature[ID_match_2, col_to_keep])
  
  return(Cys_dropout_df)
}



plot_dep_class <- function(dep_score, all_covered_dep, dep_cutoff=-0.4){
  
  HUGO_sel <- colnames(dep_score)%in%all_covered_dep
  x <- dep_score['KMS26', HUGO_sel]
  y <- dep_score['PC14', HUGO_sel]
  
  col_vec <- brewer.pal(n=4, name='Set1')
  pt_col <- rep('grey', sum(HUGO_sel))
  pt_col[x<=dep_cutoff & y<=dep_cutoff] <- col_vec[4] # purple
  pt_col[x<=dep_cutoff & y>dep_cutoff] <- col_vec[1] # red
  pt_col[x>dep_cutoff & y<=dep_cutoff] <- col_vec[2] # blue
  
  plot(x, y, cex=1.5, pch=21, bg=alpha(pt_col, 0.5), col=alpha('black', 0.5),
       xlab='KO score in KMS26', ylab='KO score in PC14')
  abline(v=-0.4, col='black', lty=2)
  abline(h=-0.4, col='black', lty=2)
  
  print(table(pt_col))
}


plot_stack_bar_sgRNA <- function(unique_ID_input, input_data){
  ID_match <- match(unique_ID_input$unique_ID_Ev87, input_data$unique_ID)
  A_guide_num <- input_data$A_guide_num[ID_match]
  C_guide_num <- input_data$C_guide_num[ID_match]
  A_guide_num[is.na(A_guide_num)] <- 0
  C_guide_num[is.na(C_guide_num)] <- 0

  libA_table <- table(A_guide_num)
  libC_table <- table(C_guide_num)
  libAC_table <- data.frame(libA_sg=as.numeric(libA_table/sum(libA_table)), 
                            libC_sg=as.numeric(libC_table/sum(libC_table)))
  
  data <- libAC_table[1:4, ]
  summed_part <- apply(libAC_table[5:dim(libAC_table)[1], ], 2, sum)
  data <- rbind(data, summed_part)
  rownames(data) <- c('sg=0', 'sg=1', 'sg=2', 'sg=3', 'sg>=4')
  
  grey_vec <- brewer.pal(n=dim(data)[1], name='Greys')
  barplot(as.matrix(data), col=grey_vec, border="black")
  
  print(round(data*100, 1))
}



plot_lig_essential_Cys <- function(Cys_res){
  
  lig_sel <- (Cys_res$KB_engage > 50) & !is.na(Cys_res$KB_engage)
  data_1 <- Cys_res[lig_sel, ]
  
  x <- data_1$mean_LFC_A
  y <- data_1$mean_LFC_C
  x[x>0] <- 0 # Only focus on the dropouts.
  y[y>0] <- 0 
  sig_sel <- data_1$sig_sel
  
  ID <- gsub("_", " C", data_1$unique_ID)
  
  col_palette <- alpha(brewer.pal(n=4, name='Set1'), 0.7)
  # Use different colors to indicate whether it is PC/KMS dependency.
  col_vec <- rep(col_palette[4], length(ID))
  col_vec[data_1$cell=='KMS'] <- col_palette[1]
  col_vec[data_1$cell=='PC'] <- col_palette[2]
  
  xmin=min(min(x)*1.1, -1)
  ymin=min(min(y)*1.1, -1)
  plot(x[!sig_sel], y[!sig_sel], cex=1.5, pch=21, bg=alpha('grey', 0.3), xlim=c(xmin, 0.05), ylim=c(ymin, 0.05), 
       col=NA, xlab='Library ABE8e L2FC', ylab='Library evoCDA L2FC')
  
  points(x[sig_sel], y[sig_sel], cex=1.5, pch=21, bg=alpha(col_vec[sig_sel], 0.7), xlim=c(xmin, 0.05), ylim=c(ymin, 0.05), 
         col='black')
  
  text(x[sig_sel], y[sig_sel], ID[sig_sel])
  abline(v=0, col='black', lty=2)
  abline(h=0, col='black', lty=2)
  
}


extract_data_for_plotting <- function(data_1, HUGO_sel, non_NA_name='KB_engage'){
  non_NA_sel <- !is.na(data_1[, non_NA_name])
  HUGO_sel <- data_1$HUGO==HUGO_sel
  data_2 <- data_1[non_NA_sel & HUGO_sel, ]
  data_2$target_Cys <- as.numeric(data_2$target_Cys)
  rownames(data_2) <- paste0('C', data_2$target_Cys)
  return(data_2)
}


plot_lig_heatmap <- function(data_1){
  lig_breaksList <- c(0, 50, 80, 100)
  lig_color <- c(alpha('grey', 0.4),  "#6BAED6", "#084594")
  pheatmap(data_1[, 'KB_engage', drop=F], color=lig_color, breaks=lig_breaksList, 
           cluster_rows=F, cluster_col=F, border_color="white", legend=F)
}


# Scale denature color to the range between -1.6 to 1.6 (log2 fold change=3)
assign_denature_col <- function(data_1){
  num_vec <- data_1$denature_LFC
  col_palette <- c("#B21F2C", "#F58F83", 'white' ,"#74BBE7", "#2066AC")
  RdBu <- colorRampPalette(col_palette)(321)
  
  names(RdBu) <- as.character(seq(-160, 160, by=1))
  integer_ID <- round(num_vec, 2)*100
  integer_ID[integer_ID>=160] <- 160
  integer_ID[integer_ID<=-160] <- -160
  
  col_vec <- RdBu[as.character(integer_ID)]
  col_vec[is.na(col_vec)] <- 'grey'
  col_vec[data_1$IADTB_access=='SS'] <- "#E7CF43"
  data_1$denature_col <- col_vec
  return(data_1)
}


plot_denature_heatmap <- function(data_1){
  temp <- matrix(1:dim(data_1)[1], ncol=1)
  rownames(temp) <- rownames(data_1)
  pheatmap(temp, color=data_1$denature_col, cluster_rows=F, cluster_col=F, border_color="white", legend=F)
}


plot_dropout_lollipop <- function(data_1){
  
  n_Cys_sel <- dim(data_1)[1]
  col_palette <- alpha(brewer.pal(n=6, name='Set2'), 0.8)
  y_cord <- -(1:n_Cys_sel)
  
  AC_max <- apply(data_1[, c('FDR_A', 'FDR_C')], 1, max, na.rm=T)
  ## Use heatmap to show dropouts among different residues.
  plot(data_1[, 'FDR_A'], y_cord, cex=0, pch=21, bg=NULL, yaxt='n', bty ="n", 
       xlim=c(0, max(AC_max)*1.1), xlab=NULL, ylab=NULL)
  segments(0, y_cord, AC_max, y_cord, col='dark grey', lty=1, lw=5)
  
  points(data_1[, 'FDR_A'], y_cord, cex=2, pch=21, bg=col_palette[5], yaxt='n', bty ="n")
  points(data_1[, 'FDR_C'], y_cord, cex=2, pch=21, bg=col_palette[6], yaxt='n', bty ="n")
  
  abline(v=0)
  
}


plot_stack_bar_IADTB <- function(Cys_res, df_denature, df_scout){
  #1) All quantified cysteines
  #2) All edited quantified cysteines
  #3) All essential quantified cysteine
  #4) All ligandable cysteines.
  scout_ID_match <- match(df_scout$unique_ID_Ev87, df_denature$unique_ID)
  df_scout$IADTB_access <- df_denature$IADTB_access[scout_ID_match]
  
  sig_sel <- Cys_res$sig_sel
  
  table_1 <- table(df_denature$IADTB_access)
  table_2 <- table(Cys_res$IADTB_access)
  table_3 <- table(Cys_res$IADTB_access[sig_sel])
  table_4 <- table(df_scout$IADTB_access[df_scout$KB_engage>50])
  
  data_1 <- cbind(table_1, table_2, table_3, table_4)
  NS <- data_1['part_buried', ] + data_1['part_free', ]
  data_2 <- rbind(data_1, NS)
  data_3 <- data_2[c('buried','NS','free'), ]
  data_4 <- t(t(data_3)/colSums(data_3))
  # dark blue, grey, dark red
  col_vec <- alpha(c("#2066AC", "grey", "#B21F2C"), 0.8)
  barplot(as.matrix(data_4), col=col_vec, border="white", ylab='Proportions')
  print(data_3)
  print(round(data_4, 3)*100)
}


plot_reactive_essential_Cys <- function(Cys_res){
  
  denature_sel <- Cys_res$IADTB_access%in%c('free', 'buried')
  data_1 <- Cys_res[denature_sel, ]
  
  x <- data_1$FDR_A
  y <- data_1$FDR_C
  sig_sel <- data_1$sig_sel
  
  ID <- gsub("_", " C", data_1$unique_ID)
  col_palette <- alpha(c("#74BBE7", "#F58F83"), 0.7)
  
  col_vec <- rep('grey', length(data_1$IADTB_access))
  col_vec[data_1$IADTB_access== 'buried'] <- col_palette[1] #Light Blue
  col_vec[data_1$IADTB_access== 'free'] <- col_palette[2] #Light Red
  
  # Use different shapes to indicate whether it is PC/KMS dependency.
  # pch=21, circle; pch=22, rectangle; pch=24, triangle.
  pch_vec <- rep(21, length(ID))
  pch_vec[data_1$cell=='PC'] <- 22
  pch_vec[data_1$cell=='KMS'] <- 24
  
  xmax=max(max(x)*1.1, 6)
  ymax=max(max(y)*1.1, 6)
  
  plot(x[!sig_sel], y[!sig_sel], cex=1, pch=pch_vec[!sig_sel], bg=alpha('grey', 0.3), xlim=c(0, xmax), ylim=c(0, ymax), 
       col=NA, xlab='Library ABE8e (-log10 FDR)', ylab='Library evoCDA (-log10 FDR)')
  points(x[sig_sel], y[sig_sel], cex=1.5, pch=pch_vec[sig_sel], bg=col_vec[sig_sel], xlim=c(0, xmax), ylim=c(0, ymax), 
         col='black')
  label_sel <- (data_1$p_A>1.3 & data_1$FDR_A>1 & data_1$mean_LFC_A<=-0.6) &
               (data_1$p_C>1.3 & data_1$FDR_C>1 & data_1$mean_LFC_C<=-0.6)

  text(x[label_sel], y[label_sel], ID[label_sel], cex=0.5)
  
}


make_part_pie <- function(input_vec, brewer_col_name="Greys"){
  table_label <- paste0(names(input_vec), ' N=', input_vec)
  if (length(input_vec)==2){
    col_vec <- brewer.pal(n=3, name=brewer_col_name)[c(1,3)]
  } else {
    col_vec <- brewer.pal(n=length(input_vec), name=brewer_col_name)
  }
  pie(input_vec, labels=table_label, border="white", col=alpha(col_vec, 0.9))
  print(round(input_vec/sum(input_vec)*100, 1))
}

# Add AlphaFold2 relative solvent area information.
add_AF_RSA <- function(data_input, AF_RSA){
  AF_RSA$unique_ID <- paste0(AF_RSA$symbol, '_', AF_RSA$residue_number)
  AF_RSA$RSA[AF_RSA$pLDDT<70] <- NA
  ID_match <- match(data_input$unique_ID, AF_RSA$unique_ID)
  data_input$RSA <- AF_RSA$RSA[ID_match]
  return(data_input)
}


plot_RSA_denature_scatter <- function(df_denature){
  x <- df_denature$RSA
  y <- df_denature$denature_LFC
  non_NA <- !is.na(x) & !is.na(y)
  col_vec <- c("white", brewer.pal(n=9, 'Greys'))
  
  smoothScatter(x[non_NA], y[non_NA], main=NULL, colramp=colorRampPalette(col_vec),
                xlab='Relative % solvent accessbility (AlphaFold2)', 
                ylab='IADTB labeling upon denaturing (L2FC)')
  ss_fit <- smooth.spline(x[non_NA], y[non_NA], cv=T, lambda=1e-2)
  lines(ss_fit, lty=1, col=2, lwd=3)
  print(cor.test(x, y))
}


plot_conserve_scatter <- function(Cys_res, col_name='p_A'){
  x <- Cys_res$conserve_1
  y <- Cys_res[, col_name]
  non_NA <- !is.na(x) & !is.na(y)
  col_vec <- c("white", brewer.pal(n=9, 'Greens'))
  
  smoothScatter(x[non_NA], y[non_NA], main=NULL, colramp=colorRampPalette(col_vec),
                xlab='Cys conservation score', 
                ylab='Cys dropout significance (−log10 p)')
  ss_fit <- smooth.spline(x[non_NA], y[non_NA], cv=T, lambda=1e-2)
  lines(ss_fit, lty=1, col=2, lwd=3)
  print(cor.test(x, y))
}


filter_dropout_data <- function(Cys_dropout, dep_score, sel_filter=T){
  
  # assign FDR for Cys within the same protein.
  Cys_dropout <- assign_FDR(Cys_dropout)
  
  # Use cell line level data to further filter the data.
  filter_both_A <- (Cys_dropout$PC_mean_LFC_A >-0.3) | (Cys_dropout$KMS_mean_LFC_A >-0.3)
  filter_both_C <- (Cys_dropout$PC_mean_LFC_C >-0.3) | (Cys_dropout$KMS_mean_LFC_C >-0.3)
  
  Cys_dropout$FDR_A[filter_both_A & Cys_dropout$cell=='both'] <- 0
  Cys_dropout$FDR_C[filter_both_C & Cys_dropout$cell=='both'] <- 0
  
  if (sel_filter){
    
    filter_PC_A <- (Cys_dropout$PC_mean_LFC_A - Cys_dropout$KMS_mean_LFC_A) >-0.3
    filter_PC_C <- (Cys_dropout$PC_mean_LFC_C - Cys_dropout$KMS_mean_LFC_C) >-0.3
    
    filter_KMS_A <- (Cys_dropout$KMS_mean_LFC_A - Cys_dropout$PC_mean_LFC_A) >-0.3
    filter_KMS_C <- (Cys_dropout$KMS_mean_LFC_C - Cys_dropout$PC_mean_LFC_C) >-0.3
    
    PC_sel_HUGO <- names(which((dep_score['PC14', ]-dep_score['KMS26', ])<=-0.8))
    KMS_sel_HUGO <- names(which((dep_score['KMS26', ]-dep_score['PC14', ])<=-0.8))
    
    Cys_dropout$FDR_A[filter_PC_A & (Cys_dropout$HUGO%in%PC_sel_HUGO) & Cys_dropout$cell=='PC'] <- 0
    Cys_dropout$FDR_C[filter_PC_C & (Cys_dropout$HUGO%in%PC_sel_HUGO) & Cys_dropout$cell=='PC'] <- 0
    Cys_dropout$FDR_A[filter_KMS_A & (Cys_dropout$HUGO%in%KMS_sel_HUGO) & Cys_dropout$cell=='KMS'] <- 0
    Cys_dropout$FDR_C[filter_KMS_C & (Cys_dropout$HUGO%in%KMS_sel_HUGO) & Cys_dropout$cell=='KMS'] <- 0
  }

  # Define Cys with significant dropouts.
  Cys_dropout$A_sig_sel <- (Cys_dropout$p_A>1.3 & Cys_dropout$FDR_A>1 & Cys_dropout$mean_LFC_A<=-0.6)
  Cys_dropout$C_sig_sel <- (Cys_dropout$p_C>1.3 & Cys_dropout$FDR_C>1 & Cys_dropout$mean_LFC_C<=-0.6)
  Cys_dropout$sig_sel <- (Cys_dropout$A_sig_sel) | (Cys_dropout$C_sig_sel)
  
  return(Cys_dropout)
}


assign_FDR <- function(Cys_dropout){
  
  Cys_dropout$FDR_A <- NA
  Cys_dropout$FDR_C <- NA
  unique_HUGO <- unique(Cys_dropout$HUGO)
  n <- length(unique_HUGO)
  
  for (i in 1:n){
    HUGO_temp <- unique_HUGO[i]
    row_sel <- which(Cys_dropout$HUGO == HUGO_temp)
    df_temp <- Cys_dropout[row_sel, ]
    Cys_dropout$FDR_A[row_sel] <- round(-log10(p.adjust(10^-df_temp$p_A, method='fdr')), 2)
    Cys_dropout$FDR_C[row_sel] <- round(-log10(p.adjust(10^-df_temp$p_C, method='fdr')), 2)
  }
  
  Cys_dropout$FDR_A[is.na(Cys_dropout$FDR_A)] <- 0
  Cys_dropout$FDR_C[is.na(Cys_dropout$FDR_C)] <- 0
  
  return(Cys_dropout)
}


compare_hit_rate <- function(Cys_dropout, null_res_list, dep_score){
  
  null_rate_vec <- rep(NA, length(null_res_list))
  
  for (i in 1:length(null_res_list)){
    null_res <- filter_dropout_data(null_res_list[[i]], dep_score)
    null_rate_vec[i] <- mean(null_res$sig_sel)
  }
  
  x <- c(mean(Cys_dropout$sig_sel), null_rate_vec)
  g <- c('1', rep('2', length(null_rate_vec)))
  x_mean <- c(x[1], mean(null_rate_vec))
  print(x_mean)
  
  beeswarm(x ~ g, cex=1.5, pch=21, bg=c("red", 'grey'), ylim=c(0, 0.13))
  abline(h=0)
  points(c(1, 2), x_mean, pch='_', cex=3)
  
}

plot_common_sig <- function(Cys_res, Cys_sel, x, y, sig_col){
  non_sig <- Cys_res$unique_ID[Cys_res$cell=='both' & !Cys_sel]
  sig <- Cys_res$unique_ID[Cys_res$cell=='both' & Cys_sel]
  
  plot(Cys_res[non_sig, x], Cys_res[non_sig, y], xlab=x, ylab=y,
       xlim=c(-3, 1.5), ylim=c(-3, 1.5), pch=21, col=NA, bg=alpha("grey", 0.5))
  points(Cys_res[sig, x], Cys_res[sig, y], 
         pch=21, col=alpha('black', 0.5), bg=alpha(sig_col, 0.5))
  abline(v=0, lty=2)
  abline(h=0, lty=2)
}

# Make beeswarm plots with connected lines.
plot_matched_diff <- function(v1, v2, v1_v2_col){
  v1_v2 <- c(v1, v2)
  v_group <- c(rep(1, length(v1)), rep(2, length(v2)))
  beeswarm(v1_v2 ~ v_group, pch=19, col=alpha("dark grey", 0.8), cex=1.5,
           pwcol=v1_v2_col, spacing=0.2)
  segments(rep(1, length(v1)), v1, rep(2, length(v2)), v2)
}

# For selective dependencies that are separated by >0.8 CERES scores,
# we can check whether the Cys dropout also shows strong selectivity.
# Each line connects the same Cys. The grey ones have Cys dropout selectivity < 0.3
# The grey ones in the old data were removed. 
QC_sel_hits <- function(Cys_res, Cys_res_NF, cell_sel, dep_score){
  
  if (cell_sel == 'PC14'){
    HUGO_sel <- names(which((dep_score['PC14', ]-dep_score['KMS26', ])<=-0.8))
    HUGO_cell_sel <- Cys_res$HUGO%in%HUGO_sel & Cys_res$cell=='PC'
  }
  
  if (cell_sel == 'KMS26'){
    HUGO_sel <- names(which((dep_score['KMS26', ]-dep_score['PC14', ])<=-0.8))
    HUGO_cell_sel <- Cys_res$HUGO%in%HUGO_sel & Cys_res$cell=='KMS'
  }
  
  A_ID_sel_NF <- Cys_res$unique_ID[Cys_res_NF$A_sig_sel & HUGO_cell_sel]
  C_ID_sel_NF <- Cys_res$unique_ID[Cys_res_NF$C_sig_sel & HUGO_cell_sel]
  
  A_ID_sel <- Cys_res$unique_ID[Cys_res$A_sig_sel & HUGO_cell_sel]
  C_ID_sel <- Cys_res$unique_ID[Cys_res$C_sig_sel & HUGO_cell_sel]
  
  v1 <- c(Cys_res[A_ID_sel_NF, 'PC_mean_LFC_A'], Cys_res[C_ID_sel_NF, 'PC_mean_LFC_C'])
  v2 <- c(Cys_res[A_ID_sel_NF, 'KMS_mean_LFC_A'], Cys_res[C_ID_sel_NF, 'KMS_mean_LFC_C'])
  names(v1) <- c(A_ID_sel_NF, C_ID_sel_NF)
  names(v2) <- c(A_ID_sel_NF, C_ID_sel_NF)
  
  print(t.test(v1-v2))
  v1_col <- rep('black', length(v1))
  v1_col[c(A_ID_sel_NF%in%A_ID_sel, C_ID_sel_NF%in%C_ID_sel)] <- 'red'
  v1_v2_col <- alpha(rep(v1_col, 2), 0.5)
  
  plot_matched_diff(v1, v2, v1_v2_col)
}


run_Cys_pipeline <- function(LFC_input_pan, null_list_pan,
                             LFC_input_PC, null_list_PC,
                             LFC_input_KMS, null_list_KMS,
                             HUGO_dep_list, AF_RSA,
                             df_scout, df_denature){
  
  HUGO_PC_dep <- HUGO_dep_list$HUGO_PC_dep
  HUGO_KMS_dep <- HUGO_dep_list$HUGO_KMS_dep
  HUGO_pan_dep <- HUGO_dep_list$HUGO_pan_dep
  
  Cys_dropout_pan_1 <- get_essential_Cys(LFC_list=LFC_input_pan, null_list=null_list_pan)
  Cys_dropout_PC_1 <- get_essential_Cys(LFC_list=LFC_input_PC, null_list=null_list_PC)
  Cys_dropout_KMS_1 <- get_essential_Cys(LFC_list=LFC_input_KMS, null_list=null_list_KMS)
  
  Cys_dropout_pan_2 <- clean_Cys_df(Cys_dropout_pan_1)
  Cys_dropout_PC_2 <- clean_Cys_df(Cys_dropout_PC_1)
  Cys_dropout_KMS_2 <- clean_Cys_df(Cys_dropout_KMS_1)
  
  ## Focus on the dependencies only in PC9 (=PC14) or KMS26 or both.
  Cys_dropout_pan_3 <- Cys_dropout_pan_2[Cys_dropout_pan_2$HUGO%in%HUGO_pan_dep, ]
  Cys_dropout_PC_3 <- Cys_dropout_PC_2[Cys_dropout_PC_2$HUGO%in%HUGO_PC_dep, ]
  Cys_dropout_KMS_3 <- Cys_dropout_KMS_2[Cys_dropout_KMS_2$HUGO%in%HUGO_KMS_dep, ]
  
  ## Add proteomic data.
  Cys_dropout_pan_4 <- add_proteomic_data(Cys_dropout_pan_3, df_scout, df_denature)
  Cys_dropout_PC_4 <- add_proteomic_data(Cys_dropout_PC_3, df_scout, df_denature)
  Cys_dropout_KMS_4 <- add_proteomic_data(Cys_dropout_KMS_3, df_scout, df_denature)
  
  # Combine data from different cell lines.
  Cys_dropout_4 <- rbind(Cys_dropout_pan_4, Cys_dropout_PC_4, Cys_dropout_KMS_4)
  Cys_dropout_4$cell <- c(rep('both', dim(Cys_dropout_pan_4)[1]), 
                          rep('PC', dim(Cys_dropout_PC_4)[1]), 
                          rep('KMS', dim(Cys_dropout_KMS_4)[1]))
  
  # Add color based on the denaturing proteomics data and add AlphaFold 2 RSA.
  Cys_dropout_5 <- assign_denature_col(Cys_dropout_4)
  Cys_dropout_5 <- add_AF_RSA(Cys_dropout_5, AF_RSA)
  
  # Add cell line dropout data.
  Cys_dropout_5$PC_mean_LFC_A <- Cys_dropout_PC_2[Cys_dropout_5$unique_ID, 'mean_LFC_A']
  Cys_dropout_5$PC_mean_LFC_C <- Cys_dropout_PC_2[Cys_dropout_5$unique_ID, 'mean_LFC_C']
  Cys_dropout_5$KMS_mean_LFC_A <- Cys_dropout_KMS_2[Cys_dropout_5$unique_ID, 'mean_LFC_A']
  Cys_dropout_5$KMS_mean_LFC_C <- Cys_dropout_KMS_2[Cys_dropout_5$unique_ID, 'mean_LFC_C']
  
  return(Cys_dropout_5)
}

# Note, these sampled null data come from non-targeting sgRNAs.
# They map to the same Cys profiled in this study except the LFC is from sampling.
simulate_NT_pos_rate <- function(LFC_input_pan, null_list_pan,
                                 LFC_input_PC, null_list_PC,
                                 LFC_input_KMS, null_list_KMS, 
                                 HUGO_dep_list, AF_RSA, df_scout, 
                                 df_denature, nsim=10){
  null_res <- list()
  
  for (i in 1:nsim){
    print(i)
    # Use the same sim_ID for all to keep the data structure.
    sim_ID <- sample(1:length(LFC_input_pan$NT_A), length(LFC_input_pan$contrast_A), replace=T)
    
    LFC_input_pan_sim <- sim_LFC_input(LFC_input_pan, sim_ID)
    LFC_input_PC_sim <- sim_LFC_input(LFC_input_PC, sim_ID)
    LFC_input_KMS_sim <- sim_LFC_input(LFC_input_KMS, sim_ID)
    
    null_res[[i]] <- run_Cys_pipeline(LFC_input_pan_sim, null_list_pan,
                                      LFC_input_PC_sim, null_list_PC,
                                      LFC_input_KMS_sim, null_list_KMS,
                                      HUGO_dep_list, AF_RSA, df_scout, 
                                      df_denature)
    
  }
  return(null_res)
}


