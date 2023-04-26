

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
  unique_ID <- mixedsort(unique_ID[unique_ID!='NA_NA'])
  
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


get_score_matrix <- function(HUGO_sel, annot_input, contrast_input, unique_ID){
  
  annot <- annot_input[annot_input$HUGO==HUGO_sel, ]
  contrast <- contrast_input[annot_input$HUGO==HUGO_sel]
  sgRNA_n <- dim(annot)[1]
  m <- matrix(NA, nrow=sgRNA_n, ncol=length(unique_ID))
  rownames(m) <- rownames(annot)
  colnames(m) <- unique_ID
  
  for (i in 1:sgRNA_n){
    aa_change_temp <- annot$aa_change[i]
    if (!is.na(aa_change_temp)){
      aa_m_temp <- matrix(unlist(strsplit(aa_change_temp, split='_')), ncol=3, byrow=T)
      m[i, as.integer(aa_m_temp[, 2])] <- contrast[i]
    }
  }
  return(m)
}


get_flat_score <- function(m){
  non_NA_count <- apply(!is.na(m), 2, sum)
  ID <- rep(colnames(m), non_NA_count)
  score <- as.vector(m[!is.na(m)])
  flat_score_df <- data.frame(ID=ID, score=score)
  return(flat_score_df)
}


compute_stat_flat <- function(flat_A, flat_C, NT_A, NT_C, n_sim=1e5){

  LFC_A <- tapply(flat_A$score, flat_A$ID, mean)
  LFC_C <- tapply(flat_C$score, flat_C$ID, mean)
  LFC_AC <- tapply(c(flat_A$score, flat_C$score), c(flat_A$ID, flat_C$ID), mean)
  
  A_outlier <- 0
  C_outlier <- 0
  AC_outlier <- 0
  
  for (i in 1:n_sim){
    sim_A <- as.numeric(NT_A[sample(1:length(NT_A), dim(flat_A)[1], replace=T)])
    sim_C <- as.numeric(NT_C[sample(1:length(NT_C), dim(flat_C)[1], replace=T)])
    sim_AC <- c(sim_A, sim_C)
    A_outlier <- A_outlier + (tapply(sim_A, flat_A$ID, mean) < LFC_A)
    C_outlier <- C_outlier + (tapply(sim_C, flat_C$ID, mean) < LFC_C)
    AC_outlier <- AC_outlier + (tapply(sim_AC, c(flat_A$ID, flat_C$ID), mean) < LFC_AC)
  }

  p_A <- A_outlier/n_sim + 1e-6
  p_C <- C_outlier/n_sim + 1e-6
  p_AC <- AC_outlier/n_sim + 1e-6
  
  res_list <- list()
  res_list$LFC_A <- LFC_A
  res_list$LFC_C <- LFC_C
  res_list$LFC_AC <- LFC_AC
 
  res_list$p_A <- log10(p_A)
  res_list$p_C <- log10(p_C)
  res_list$p_AC <- log10(p_AC)
  
  return(res_list)
}

get_LFC_list <- function(flat_data, unique_ID){
  LFC_list <- list()
  n <- length(unique_ID)
  for (i in 1:n){
    row_sel <- flat_data$ID==unique_ID[i]
    LFC_list[[i]] <- sort(round(as.vector(flat_data$score[row_sel]), 2))
  }
  return(LFC_list)
}


get_saturated_result <- function(HUGO_sel, LFC_input, combined_seq_df) {
  
  protein_seq <- combined_seq_df[which(combined_seq_df$HUGO==HUGO_sel), "Ev87_seq"]
  aa_seq <- strsplit(protein_seq, split='')[[1]]
  unique_ID <- paste0(HUGO_sel, '_', aa_seq, seq(1:length(aa_seq)))
  
  m_A <- get_score_matrix(HUGO_sel, annot_input=LFC_input$ABE_annot, 
                                 contrast_input=LFC_input$contrast_A, unique_ID)
  m_C <- get_score_matrix(HUGO_sel, annot_input=LFC_input$evoCDA_annot, 
                                 contrast_input=LFC_input$contrast_C, unique_ID)
  
  flat_A <- get_flat_score(m_A)
  flat_C <- get_flat_score(m_C)
  
  res_list <- compute_stat_flat(flat_A, flat_C, LFC_input$NT_A, LFC_input$NT_C, n_sim=1e5)
  
  result <- data.frame(HUGO=HUGO_sel, unique_ID=unique_ID,
                       p_A=round(res_list$p_A[unique_ID], 2), 
                       p_C=round(res_list$p_C[unique_ID], 2),
                       p_AC=round(res_list$p_AC[unique_ID], 2),
                       LFC_A=round(res_list$LFC_A[unique_ID], 2), 
                       LFC_C=round(res_list$LFC_C[unique_ID], 2),
                       LFC_AC=round(res_list$LFC_AC[unique_ID], 2))
  
  result$A_sgRNA_LFC <- get_LFC_list(flat_A, unique_ID)
  result$C_sgRNA_LFC <- get_LFC_list(flat_C, unique_ID)
  
  result$A_guide_num <- as.integer(lapply(result$A_sgRNA_LFC, length))
  result$C_guide_num <- as.integer(lapply(result$C_sgRNA_LFC, length))
  
  # Only focus on residues with at least 2 sgRNA data points.
  result$p_A[result$A_guide_num <= 1] <- 0
  result$p_C[result$C_guide_num <= 1] <- 0
  result$p_AC[(result$A_guide_num + result$C_guide_num) <= 1] <- 0

  return(result)
}


plot_residue_waterfall <- function(unique_ID, LFC, p, LFC_cut=-0.4, pcut=-1.3, title=NA){
  
  x <- (1:length(unique_ID))
  y <- p
  non_sig <- (LFC>LFC_cut) | is.na(LFC) | p>pcut | is.na(p)
  y[non_sig] <- 0
  
  ymin <- min(y, na.rm=T)
  sel_Cys_pos <- grep('_C', unique_ID)
  
  plot(x, y, cex=0, pch=21, ylim=c(ymin*1.1, 0), xlab='Residue', ylab='Dropout log10 p value')
  segments(x, 0, x, y, col=alpha('#555555', 0.5), lty=1, lw=5)
  points(x, y, cex=1, pch=21, col=NA, bg='dark grey', xlab=NULL, ylab=NULL)
  points(sel_Cys_pos, y[sel_Cys_pos], cex=1.5, pch=21, bg="#FFD92FB2", yaxt='n', bty ="n")
  
}



compare_2_Cys <- function(LFC_input1, LFC_input2, p_input1, p_input2, unique_ID, 
                          LFC_cut=-0.4, pcut=-1.3, title=NA){

  sel_Cys_pos <- grep('_C', unique_ID)
  LFC1 <- LFC_input1[sel_Cys_pos]
  LFC2 <- LFC_input2[sel_Cys_pos]
  p1 <- p_input1[sel_Cys_pos]
  p2 <- p_input2[sel_Cys_pos]
  
  non_sig_1 <- (LFC1>LFC_cut) | is.na(LFC1) | p1>pcut | is.na(p1)
  non_sig_2 <- (LFC2>LFC_cut) | is.na(LFC2) | p2>pcut | is.na(p2)
  
  p1[non_sig_1] <- 0
  p2[non_sig_2] <- 0
  
  col_palette <- alpha(brewer.pal(n=4, name='Set1'), 0.7)
  max_dropout <- min(c(p1, p2), na.rm=T)
  
  plot(p1-0.02, 1:length(p1), cex=1.5, xlim=c(max_dropout*1.1, 0.5), col=NA,
       bg=col_palette[1], pch=21, xlab='Dropout p value', ylab='Residue')
  points(p2+0.02, 1:length(p2), cex=1.5, col=NA,
         bg=col_palette[2], pch=22)
  abline(v=0, lw=2)
  
  df <- data.frame(sel_Cys_pos, LFC1, p1, LFC2, p2)
  return(df)
}


# This is to assign ID one-by-one and skip the '-' gaps (ensembl: continuous)
assign_ensembl_ID <- function(matched_seq, st) {
  matched_ID <- rep(NA, length(matched_seq))
  for(i in 1:length(matched_seq)){
    s <- matched_seq[i]
    if (s != '-'){
      matched_ID[i] <- st
      st <- st + 1
    }
  }
  return(matched_ID)
}

# PDB ID might be discontinuous, e.g., aa 988 then aa 1004
assign_pdb_ID <- function(matched_seq, rel_st, pdb_sel) {
  all_IDs <- names(pdbseq(pdb_sel))
  all_IDs <- all_IDs[rel_st:length(all_IDs)]
  matched_ID <- rep(NA, length(matched_seq))
  for(i in 1:length(matched_seq)){
    s <- matched_seq[i]
    if (s != '-'){
      matched_ID[i] <- all_IDs[1]
      all_IDs <- all_IDs[-1]
    }
  }
  return(as.integer(matched_ID))
}



get_clean_PDB <- function(PDB_ID, HUGO, combined_seq_df, chain='A', coef_vec=NA, output_folder='~/Desktop'){
  pdb_data_w_H2O <- read.pdb(PDB_ID)
  # Remove water molecules
  pdb_data <- trim.pdb(pdb_data_w_H2O, atom.select(pdb_data_w_H2O, "water", inverse=T))
  
  chain_index <- atom.select(pdb_data, type="ATOM", chain=chain) 
  pdb_data$atom$b <- 0
  pdb_sel <- trim.pdb(pdb_data, chain_index)
  pdb_seq_raw <- as.vector(pdbseq(pdb_sel))
  pdb_seq <- paste(pdb_seq_raw, collapse='')
  ensembl_seq <- combined_seq_df[combined_seq_df$HUGO==HUGO, "Ev87_seq"]
  
  # Use local alignment to find overlapping regions.
  localAlign <- pairwiseAlignment(pattern=pdb_seq, subject=ensembl_seq, type = "local")
  
  # Note the relative aa position in ensembl is the same as the absolute position (start from 1)
  # But for PDB, the relative position 1 might differ from the absolute position (protein fragment).
  ensembl_st <- start(subject(localAlign))
  ensembl_ed <- end(subject(localAlign))
  pdb_rel_st <- start(pattern(localAlign))
  pdb_rel_ed <- end(pattern(localAlign))
  
  pdb_matched_seq <- unlist(strsplit(as.character(pattern(localAlign)), ''))
  ensembl_matched_seq <- unlist(strsplit(as.character(subject(localAlign)), ''))
  
  pdb_matched_ID <- assign_pdb_ID(matched_seq=pdb_matched_seq, rel_st=pdb_rel_st, pdb_sel=pdb_sel)
  ensembl_matched_ID <- assign_ensembl_ID(matched_seq=ensembl_matched_seq, st=ensembl_st)
  
  # Label the residues that I want to discard as -1
  tmp_resno <- pdb_data$atom$resno
  tmp_resno[chain_index$atom] <- -1
  
  for(n in 1:length(pdb_matched_ID)){
    pdb_aa <- pdb_matched_ID[n]
    ensembl_aa <- ensembl_matched_ID[n]
    if (!is.na(pdb_aa) & !is.na(ensembl_aa)){
      tmp_aa_index <- atom.select(pdb_data, type="ATOM", chain=chain, resno=pdb_aa)
      tmp_resno[tmp_aa_index$atom] <- ensembl_aa
      # If there is no coef_vec, we just updated the amino acid positions.
      if (!is.na(sum(coef_vec))){
        pdb_data$atom$b[tmp_aa_index$atom] <- coef_vec[ensembl_aa, 1] 
      }
    }
  }
  pdb_data$atom$resno <- tmp_resno
  trim_pdb_data <- trim.pdb(pdb_data, atom.select(pdb_data, resno=-1, inverse=T))
  
  output_name <- paste0(HUGO, "_", PDB_ID, "_cleaned.pdb")
  write.pdb(trim_pdb_data, file=file.path(output_folder, output_name))
  coverage_metric <- sum(pdb_matched_ID&ensembl_matched_ID, na.rm=T)/length(pdb_seq_raw)
  print(paste0("Percentage PDB sequence covered: ", 100*round(coverage_metric, 2)))
}
  

get_coef_vec_from_p <- function(p_vec, cutoff=-1.3){
  p_vec[is.na(p_vec) | p_vec>cutoff] <- 0
  coef_vec <- matrix(p_vec, ncol=1)
  return(coef_vec)
}


