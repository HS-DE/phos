# ============================================================
# Clean PhosR kinase activity workflow for mouse phosphoproteomics
#
# Purpose:
#   This script replaces the messy exploratory KSEA流程.R logic with a
#   function-based PhosR workflow. It does NOT depend on temporary variables
#   typed in the R console.
#
# Required upstream objects / files:
#   1) phos_df: site-level phosphosite table, recommended from site_mat2
#      Required columns: Protein.Id, Genes, abs_pos, abs_residue, sample columns
#   2) fasta_map: named protein sequence vector/list, generated from mouse FASTA
#   3) sample_info: data.frame with at least two columns:
#      - sample: sample column names matching phos_df
#      - group : biological group / condition for each sample
#   4) sample_cols: character vector of quantitative sample columns in phos_df
#
# Notes:
#   - KSEAapp was originally designed around human PSP/NetworKIN usage.
#     For mouse data, this script uses PhosR::PhosphoSite.mouse and
#     PhosR::kinaseSubstrateScore().
#   - If there are no biological replicates, the script only reports logFC-like
#     changes and does not pretend to calculate valid P values / FDR.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(PhosR)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(limma)
  library(pheatmap)
})

# ------------------------------------------------------------
# Helper 1. Safely fetch protein sequence from fasta_map
# ------------------------------------------------------------
get_protein_sequence <- function(protein_id, fasta_map) {
  protein_id <- as.character(protein_id)

  fetch_one <- function(x) {
    if (is.null(x)) return(NA_character_)
    if (is.list(x)) return(as.character(x[[1]]))
    as.character(x)
  }

  seq <- fetch_one(fasta_map[[protein_id]])
  if (!is.na(seq)) return(seq)

  # If isoform is missing in FASTA, try the canonical accession.
  # Example: P12345-2 -> P12345
  if (str_detect(protein_id, "-")) {
    protein_id_base <- sub("-.*$", "", protein_id)
    seq2 <- fetch_one(fasta_map[[protein_id_base]])
    if (!is.na(seq2)) return(seq2)
  }

  NA_character_
}

# ------------------------------------------------------------
# Helper 2. Generate 15-aa flanking sequence around phosphosite
# ------------------------------------------------------------
get_flank15 <- function(prot_seq, pos, win = 7, pad = "X") {
  # prot_seq: full protein sequence
  # pos     : 1-based absolute phosphosite position
  # output  : length-15 sequence: 7 aa upstream + center site + 7 aa downstream
  if (is.na(prot_seq) || is.na(pos)) return(NA_character_)

  pos <- as.integer(pos)
  L <- nchar(prot_seq)
  if (pos < 1 || pos > L) return(NA_character_)

  st <- pos - win
  ed <- pos + win

  left_pad  <- if (st < 1) paste(rep(pad, 1 - st), collapse = "") else ""
  right_pad <- if (ed > L) paste(rep(pad, ed - L), collapse = "") else ""

  core <- substr(prot_seq, max(1, st), min(L, ed))
  seq15 <- paste0(left_pad, core, right_pad)

  if (nchar(seq15) != 2 * win + 1) return(NA_character_)
  seq15
}

# ------------------------------------------------------------
# Helper 3. Compute pairwise logFC / limma results
# ------------------------------------------------------------
run_pairwise_site_change <- function(mat_norm, grp_vec, out_dir) {
  # mat_norm: rows = phosphosites, columns = samples
  # grp_vec : group vector named by sample columns

  stopifnot(length(grp_vec) == ncol(mat_norm))
  names(grp_vec) <- colnames(mat_norm)

  grp_fac <- factor(grp_vec)
  old_lvls <- levels(grp_fac)
  new_lvls <- make.names(old_lvls)
  levels(grp_fac) <- new_lvls

  group_name_map <- data.frame(group_raw = old_lvls, group_safe = new_lvls)
  write.csv(group_name_map, file.path(out_dir, "05_group_name_map.csv"), row.names = FALSE)

  tab_n <- table(grp_fac)
  has_residual_df <- ncol(mat_norm) > length(levels(grp_fac))

  grp_levels <- levels(grp_fac)
  pair_mat <- combn(grp_levels, 2)
  contrast_formula <- apply(pair_mat, 2, function(x) paste0(x[2], "-", x[1]))
  contrast_names <- apply(pair_mat, 2, function(x) paste0(x[2], "_vs_", x[1]))

  de_list <- setNames(vector("list", length(contrast_names)), contrast_names)

  if (has_residual_df) {
    # With residual degrees of freedom, limma can estimate variance.
    design <- model.matrix(~ 0 + grp_fac)
    colnames(design) <- levels(grp_fac)

    fit <- limma::lmFit(mat_norm, design)
    cont_mat <- limma::makeContrasts(contrasts = contrast_formula, levels = design)
    colnames(cont_mat) <- contrast_names

    fit2 <- limma::contrasts.fit(fit, cont_mat)
    fit2 <- limma::eBayes(fit2)

    for (cn in contrast_names) {
      tt <- limma::topTable(fit2, coef = cn, number = Inf, sort.by = "P")
      tt$site_label <- rownames(tt)
      tt$contrast <- cn
      de_list[[cn]] <- tt
      write.csv(tt, file.path(out_dir, paste0("DE_", cn, ".csv")), row.names = FALSE)
    }
  } else {
    # No replicate / no residual df: report logFC only.
    # P.Value and adj.P.Val are intentionally set to NA.
    message("No residual degrees of freedom detected. Only logFC-like changes will be exported.")

    for (k in seq_along(contrast_names)) {
      cn <- contrast_names[k]
      g1 <- pair_mat[1, k]
      g2 <- pair_mat[2, k]

      col1 <- which(grp_fac == g1)[1]
      col2 <- which(grp_fac == g2)[1]
      logFC <- mat_norm[, col2] - mat_norm[, col1]

      tt <- data.frame(
        logFC = as.numeric(logFC),
        AveExpr = rowMeans(mat_norm, na.rm = TRUE),
        t = NA_real_,
        P.Value = NA_real_,
        adj.P.Val = NA_real_,
        B = NA_real_,
        site_label = rownames(mat_norm),
        contrast = cn,
        stringsAsFactors = FALSE
      )

      de_list[[cn]] <- tt
      write.csv(tt, file.path(out_dir, paste0("DE_", cn, ".csv")), row.names = FALSE)
    }
  }

  saveRDS(de_list, file.path(out_dir, "05_DE_list.rds"))
  de_list
}

# ------------------------------------------------------------
# Main function
# ------------------------------------------------------------
run_phosr_kinase_activity <- function(phos_df,
                                      fasta_map,
                                      sample_info,
                                      sample_cols,
                                      out_dir = "Results_2",
                                      control_group = NULL,
                                      keep_only_single_site = TRUE,
                                      ratio_method = c("row_median", "control_mean"),
                                      do_phosr_impute = TRUE,
                                      topN = 5000) {
  ratio_method <- match.arg(ratio_method)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ----------------------------
  # Step 1. Input checks
  # ----------------------------
  required_cols <- c("Protein.Id", "Genes", "abs_pos", "abs_residue", sample_cols)
  missing_cols <- setdiff(required_cols, names(phos_df))
  if (length(missing_cols) > 0) {
    stop("phos_df is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!all(c("sample", "group") %in% names(sample_info))) {
    stop("sample_info must contain columns: sample and group")
  }

  # Explicitly generate grp_vec here. This is the key fix compared with the
  # original KSEA流程.R, which depended on a temporary console variable.
  grp_vec <- sample_info$group[match(sample_cols, sample_info$sample)]
  if (anyNA(grp_vec)) {
    missing_samples <- sample_cols[is.na(grp_vec)]
    stop("These sample columns are not found in sample_info$sample: ",
         paste(missing_samples, collapse = ", "))
  }
  names(grp_vec) <- sample_cols

  if (ratio_method == "control_mean") {
    if (is.null(control_group)) {
      stop("control_group must be provided when ratio_method = 'control_mean'.")
    }
    if (!control_group %in% grp_vec) {
      stop("control_group is not found in sample_info$group: ", control_group)
    }
  }

  # ----------------------------
  # Step 2. Clean site-level table
  # ----------------------------
  phos_df2 <- phos_df %>%
    mutate(
      GeneSymbol = toupper(Genes),
      Residue = as.character(abs_residue),
      Site = as.integer(abs_pos)
    ) %>%
    filter(!is.na(Site), !is.na(Residue), Residue %in% c("S", "T", "Y"))

  if (keep_only_single_site && "Residue.Both" %in% names(phos_df2)) {
    phos_df2 <- phos_df2 %>%
      mutate(n_site_in_row = str_count(Residue.Both, ";") + 1) %>%
      filter(n_site_in_row == 1) %>%
      select(-n_site_in_row)
  }

  phos_df2 <- phos_df2 %>%
    mutate(
      site_label = paste0(GeneSymbol, ";", Residue, Site, ";"),
      prot_seq = map_chr(Protein.Id, ~ get_protein_sequence(.x, fasta_map)),
      Sequence15 = map2_chr(prot_seq, Site, ~ get_flank15(.x, .y, win = 7))
    ) %>%
    select(-prot_seq) %>%
    filter(!is.na(Sequence15))

  write.csv(phos_df2, file.path(out_dir, "00_phos_df_for_phosr.csv"), row.names = FALSE)

  # ----------------------------
  # Step 3. Build quantification matrix
  # ----------------------------
  quant_mat <- phos_df2 %>%
    select(all_of(sample_cols)) %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()

  rownames(quant_mat) <- phos_df2$site_label
  colnames(quant_mat) <- sample_cols

  quant_log2 <- log2(quant_mat + 1)

  if (ratio_method == "row_median") {
    # Each phosphosite is centered by its median across all samples.
    row_ref <- apply(quant_log2, 1, median, na.rm = TRUE)
    quant_ratio <- sweep(quant_log2, 1, row_ref, "-")
  } else {
    # Each phosphosite is centered by the mean of the control group.
    ctl_cols <- grp_vec == control_group
    ctl_mean <- rowMeans(quant_log2[, ctl_cols, drop = FALSE], na.rm = TRUE)
    quant_ratio <- sweep(quant_log2, 1, ctl_mean, "-")
  }

  write.csv(
    data.frame(site_label = rownames(quant_ratio), quant_ratio, check.names = FALSE),
    file = file.path(out_dir, "01_quant_log2ratio_raw.csv"),
    row.names = FALSE
  )

  # ----------------------------
  # Step 4. Build PhosphoExperiment object
  # ----------------------------
  coldata <- S4Vectors::DataFrame(
    group = factor(grp_vec[colnames(quant_ratio)]),
    row.names = colnames(quant_ratio)
  )

  ppe <- PhosR::PhosphoExperiment(
    assays = list(Quantification = quant_ratio),
    colData = coldata,
    UniprotID = phos_df2$Protein.Id,
    GeneSymbol = phos_df2$GeneSymbol,
    Site = phos_df2$Site,
    Residue = phos_df2$Residue,
    Sequence = phos_df2$Sequence15
  )

  saveRDS(ppe, file.path(out_dir, "02_ppe_raw.rds"))

  # ----------------------------
  # Step 5. Filtering / imputation / scaling
  # ----------------------------
  grps <- as.character(SummarizedExperiment::colData(ppe)$group)
  names(grps) <- colnames(ppe)

  ppe_filt <- PhosR::selectGrps(ppe, grps, percent = 0.7, n = 1, assay = "Quantification")
  mat_filt <- SummarizedExperiment::assay(ppe_filt, "Quantification")

  # scImpute is most meaningful when at least one group has replicates.
  # If every group has only one sample, skip scImpute instead of forcing it.
  if (do_phosr_impute && any(table(grps) >= 2)) {
    set.seed(123)
    mat_sc <- PhosR::scImpute(mat_filt, percent = 0.7, grps)
  } else {
    message("Skipping scImpute because no group has replicates or do_phosr_impute = FALSE.")
    mat_sc <- mat_filt
  }

  set.seed(123)
  mat_t <- PhosR::tImpute(mat_sc)
  mat_scaled <- PhosR::medianScaling(mat_t, scale = FALSE)

  SummarizedExperiment::assay(ppe_filt, "imputed") <- mat_t
  SummarizedExperiment::assay(ppe_filt, "scaled") <- mat_scaled

  saveRDS(ppe_filt, file.path(out_dir, "03_ppe_scaled.rds"))
  write.csv(
    data.frame(site_label = rownames(mat_scaled), mat_scaled, check.names = FALSE),
    file = file.path(out_dir, "03_quant_scaled.csv"),
    row.names = FALSE
  )

  mat_norm <- mat_scaled

  # ----------------------------
  # Step 6. Pairwise site-level changes
  # ----------------------------
  de_list <- run_pairwise_site_change(mat_norm, grp_vec[colnames(mat_norm)], out_dir)

  # ----------------------------
  # Step 7. PhosR kinase activity analysis
  # ----------------------------
  data("PhosphoSitePlus", package = "PhosR")

  row_sd <- apply(mat_norm, 1, sd, na.rm = TRUE)
  min_non_na <- max(2, floor(0.5 * ncol(mat_norm)))
  ok_keep <- rowSums(!is.na(mat_norm)) >= min_non_na

  row_sd2 <- row_sd
  row_sd2[!ok_keep] <- NA_real_

  idx <- order(row_sd2, decreasing = TRUE, na.last = NA)
  idx <- idx[seq_len(min(topN, length(idx)))]
  mat_reg <- mat_norm[idx, , drop = FALSE]

  mat_std <- PhosR::standardise(mat_reg)

  seq_all <- PhosR::Sequence(ppe_filt)
  names(seq_all) <- rownames(SummarizedExperiment::assay(ppe_filt, "Quantification"))
  seq_reg <- seq_all[rownames(mat_reg)]

  keep_seq <- !is.na(seq_reg) & nzchar(seq_reg)
  mat_std2 <- mat_std[keep_seq, , drop = FALSE]
  seq_reg2 <- seq_reg[keep_seq]

  KSR <- PhosR::kinaseSubstrateScore(
    substrate.list = PhosphoSite.mouse,
    mat = mat_std2,
    seqs = seq_reg2,
    numMotif = 5,
    numSub = 1,
    species = "mouse",
    verbose = TRUE
  )

  saveRDS(KSR, file.path(out_dir, "06_KSR_kinaseSubstrateScore.rds"))

  ks_act <- KSR[["ksActivityMatrix"]]
  write.csv(
    data.frame(Kinase = rownames(ks_act), ks_act, check.names = FALSE),
    file = file.path(out_dir, "06_kinase_activity_matrix.csv"),
    row.names = FALSE
  )

  pdf(file.path(out_dir, "06_kinaseSubstrateHeatmap.pdf"), width = 10, height = 8)
  PhosR::kinaseSubstrateHeatmap(KSR)
  dev.off()

  png(file.path(out_dir, "06_kinaseSubstrateHeatmap.png"), width = 10, height = 8, units = "in", res = 600)
  PhosR::kinaseSubstrateHeatmap(KSR)
  dev.off()

  pdf(file.path(out_dir, "06_kinase_activity_heatmap_pheatmap.pdf"), width = 10, height = 10)
  pheatmap::pheatmap(ks_act, fontsize_row = 6, fontsize_col = 8)
  dev.off()

  png(file.path(out_dir, "06_kinase_activity_heatmap_pheatmap.png"), width = 10, height = 10, units = "in", res = 600)
  pheatmap::pheatmap(ks_act, fontsize_row = 6, fontsize_col = 8)
  dev.off()

  # ----------------------------
  # Step 8. Kinase activity delta between samples/groups
  # ----------------------------
  ks_act_mat <- as.matrix(ks_act)
  pairs <- combn(colnames(ks_act_mat), 2, simplify = FALSE)

  delta_mat <- do.call(
    cbind,
    lapply(pairs, function(p) {
      ks_act_mat[, p[2]] - ks_act_mat[, p[1]]
    })
  )
  colnames(delta_mat) <- vapply(pairs, function(p) paste0(p[2], "_vs_", p[1]), character(1))
  rownames(delta_mat) <- rownames(ks_act_mat)

  write.csv(
    data.frame(Kinase = rownames(delta_mat), delta_mat, check.names = FALSE),
    file = file.path(out_dir, "06_kinase_activity_delta_all_contrasts.csv"),
    row.names = FALSE
  )

  pdf(file.path(out_dir, "06_kinase_activity_delta_heatmap.pdf"), width = 10, height = 10)
  pheatmap::pheatmap(delta_mat, fontsize_row = 6, fontsize_col = 8)
  dev.off()

  png(file.path(out_dir, "06_kinase_activity_delta_heatmap.png"), width = 10, height = 10, units = "in", res = 300)
  pheatmap::pheatmap(delta_mat, fontsize_row = 6, fontsize_col = 8)
  dev.off()

  # ----------------------------
  # Return key objects for interactive inspection
  # ----------------------------
  invisible(list(
    phos_df_for_phosr = phos_df2,
    quant_ratio = quant_ratio,
    ppe = ppe,
    ppe_filt = ppe_filt,
    mat_norm = mat_norm,
    de_list = de_list,
    KSR = KSR,
    kinase_activity = ks_act,
    kinase_activity_delta = delta_mat
  ))
}

# ============================================================
# Example usage
# ============================================================
# result <- run_phosr_kinase_activity(
#   phos_df = site_mat2,
#   fasta_map = fasta_map,
#   sample_info = sample_info,
#   sample_cols = sample_cols,
#   out_dir = "../Results/PhosR_kinase_activity",
#   ratio_method = "row_median",
#   keep_only_single_site = TRUE,
#   topN = 5000
# )
