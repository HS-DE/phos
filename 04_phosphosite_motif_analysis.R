# ============================================================
# 04. Phosphosite-centered motif analysis
#
# Purpose:
#   Extract flanking amino acid sequences around phosphorylation sites
#   and draw sequence logos for pS / pT / pY sites.
#
# Main input:
#   site_mat2   : output from 01_build_phosphosite_matrix_from_diann.R
#   fasta_file  : UniProt protein FASTA used in step 01
#
# Main output:
#   01_all_phosphosite_flanking_sequences.csv
#   02_motif_summary_by_residue.csv
#   frequency_matrix/all_pS_frequency_matrix.csv
#   frequency_matrix/all_pT_frequency_matrix.csv
#   frequency_matrix/all_pY_frequency_matrix.csv
#   all_sites/all_pS_logo.pdf/png
#   all_sites/all_pT_logo.pdf/png
#   all_sites/all_pY_logo.pdf/png
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggseqlogo)
})

# ------------------------------------------------------------
# 0. Small utilities
# ------------------------------------------------------------

.safe_dir_create <- function(path) {
  if (!is.null(path) && nzchar(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

.read_table_if_path <- function(x) {
  if (is.character(x) && length(x) == 1 && file.exists(x)) {
    utils::read.csv(x, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    x
  }
}

.safe_slug <- function(x) {
  x <- as.character(x)
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[/\\\\:;*?\"<>|]", "_", x)
  x <- gsub("\\+", "plus", x)
  x <- gsub("[^[:alnum:]_.-]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# ------------------------------------------------------------
# 1. Read UniProt FASTA as named sequence vector
# ------------------------------------------------------------
# 这个函数和 01 脚本里的 read_fasta_as_map() 逻辑一致：
# 把 FASTA 读成：
#   names(fasta_map) = UniProt accession
#   fasta_map[[accession]] = protein sequence

read_fasta_as_map <- function(fasta_file) {
  x <- readLines(fasta_file, warn = FALSE)
  
  header_idx <- which(startsWith(x, ">"))
  
  if (length(header_idx) == 0) {
    stop("No FASTA headers found in: ", fasta_file)
  }
  
  end_idx <- c(header_idx[-1] - 1, length(x))
  headers <- x[header_idx]
  
  seqs <- vapply(seq_along(header_idx), function(i) {
    st <- header_idx[i]
    ed <- end_idx[i]
    paste0(x[(st + 1):ed], collapse = "")
  }, character(1))
  
  # UniProt header 通常长这样：
  # >sp|Q3UH06|RREB1_MOUSE Ras-responsive element-binding protein 1 OS=Mus musculus ...
  # 所以这里优先取第二段 Q3UH06 作为 Protein.Id
  acc <- vapply(headers, function(h) {
    h <- sub("^>", "", h)
    parts <- strsplit(h, "\\|")[[1]]
    
    if (length(parts) >= 2) {
      parts[2]
    } else {
      strsplit(h, " ")[[1]][1]
    }
  }, character(1))
  
  stats::setNames(seqs, acc)
}

# ------------------------------------------------------------
# 2. Find protein sequence by Protein.Id
# ------------------------------------------------------------
# 这里做一个兼容：
# 如果 Protein.Id 是 isoform，例如 P12345-2，
# 而 FASTA 里没有 P12345-2，就尝试找 P12345。

.get_protein_sequence <- function(protein_id, fasta_map) {
  protein_id <- as.character(protein_id)
  
  if (protein_id %in% names(fasta_map)) {
    return(as.character(fasta_map[[protein_id]]))
  }
  
  if (grepl("-", protein_id)) {
    protein_id_base <- sub("-.*$", "", protein_id)
    
    if (protein_id_base %in% names(fasta_map)) {
      return(as.character(fasta_map[[protein_id_base]]))
    }
  }
  
  NA_character_
}

# ------------------------------------------------------------
# 3. Extract one flanking motif sequence
# ------------------------------------------------------------
# 输入：
#   protein_id   : 蛋白 accession
#   abs_pos      : 磷酸化位点在蛋白上的绝对位置，例如 673
#   abs_residue  : 中心残基 S/T/Y
#   fasta_map    : FASTA 序列表
#   flank        : 左右各取几个 aa，默认 6
#
# 输出：
#   motif_seq    : 长度 2*flank+1 的序列，例如 13 aa
#   valid_motif  : 是否成功提取
#   fail_reason  : 如果失败，原因是什么

.extract_one_motif <- function(protein_id, abs_pos, abs_residue, fasta_map, flank = 6) {
  prot_seq <- .get_protein_sequence(protein_id, fasta_map)
  
  if (is.na(prot_seq) || !nzchar(prot_seq)) {
    return(tibble(
      protein_length = NA_integer_,
      motif_start = NA_integer_,
      motif_end = NA_integer_,
      center_from_fasta = NA_character_,
      motif_seq = NA_character_,
      valid_motif = FALSE,
      fail_reason = "protein_not_found_in_fasta"
    ))
  }
  
  protein_length <- nchar(prot_seq)
  abs_pos <- suppressWarnings(as.integer(abs_pos))
  
  if (is.na(abs_pos) || abs_pos < 1 || abs_pos > protein_length) {
    return(tibble(
      protein_length = protein_length,
      motif_start = NA_integer_,
      motif_end = NA_integer_,
      center_from_fasta = NA_character_,
      motif_seq = NA_character_,
      valid_motif = FALSE,
      fail_reason = "abs_pos_out_of_range"
    ))
  }
  
  motif_start <- abs_pos - flank
  motif_end <- abs_pos + flank
  
  # 第一版先要求必须能取满完整窗口。
  # 例如 flank=6 时，必须完整取到 13 aa。
  # 靠近蛋白 N/C 端的位点先过滤掉。
  if (motif_start < 1 || motif_end > protein_length) {
    center_from_fasta <- substr(prot_seq, abs_pos, abs_pos)
    
    return(tibble(
      protein_length = protein_length,
      motif_start = motif_start,
      motif_end = motif_end,
      center_from_fasta = center_from_fasta,
      motif_seq = NA_character_,
      valid_motif = FALSE,
      fail_reason = "incomplete_flanking_window"
    ))
  }
  
  center_from_fasta <- substr(prot_seq, abs_pos, abs_pos)
  
  if (!identical(center_from_fasta, as.character(abs_residue))) {
    return(tibble(
      protein_length = protein_length,
      motif_start = motif_start,
      motif_end = motif_end,
      center_from_fasta = center_from_fasta,
      motif_seq = NA_character_,
      valid_motif = FALSE,
      fail_reason = "center_residue_mismatch"
    ))
  }
  
  motif_seq <- substr(prot_seq, motif_start, motif_end)
  
  tibble(
    protein_length = protein_length,
    motif_start = motif_start,
    motif_end = motif_end,
    center_from_fasta = center_from_fasta,
    motif_seq = motif_seq,
    valid_motif = TRUE,
    fail_reason = NA_character_
  )
}

# ------------------------------------------------------------
# 4. Build motif table from site_mat2
# ------------------------------------------------------------

build_phosphosite_motif_table <- function(site_mat2,
                                          fasta_file,
                                          flank = 6,
                                          residues = c("S", "T", "Y"),
                                          keep_unique_site = TRUE) {
  site_mat2 <- .read_table_if_path(site_mat2)
  site_mat2 <- as.data.frame(site_mat2, check.names = FALSE, stringsAsFactors = FALSE)
  
  required_cols <- c("Protein.Id", "Genes", "site_id", "abs_pos", "abs_residue")
  
  missing_cols <- setdiff(required_cols, colnames(site_mat2))
  if (length(missing_cols) > 0) {
    stop("site_mat2 缺少这些列: ", paste(missing_cols, collapse = ", "))
  }
  
  fasta_map <- read_fasta_as_map(fasta_file)
  
  site_tbl <- site_mat2 %>%
    mutate(
      Protein.Id = as.character(Protein.Id),
      Genes = as.character(Genes),
      site_id = as.character(site_id),
      abs_pos = suppressWarnings(as.integer(abs_pos)),
      abs_residue = toupper(as.character(abs_residue)),
      Residue.Both = if ("Residue.Both" %in% colnames(.)) as.character(Residue.Both) else NA_character_,
      n_site_in_peptide = ifelse(
        !is.na(Residue.Both) & Residue.Both != "",
        stringr::str_count(Residue.Both, fixed(";")) + 1L,
        NA_integer_
      )
    ) %>%
    filter(abs_residue %in% residues)
  
  # 这里是一个重要选择：
  # 对全局 motif 来说，同一个 phosphosite 不应该因为被多个肽段检测到而重复计数。
  # 所以默认按 Protein.Id + site_id + abs_pos + abs_residue 去重。
  if (isTRUE(keep_unique_site)) {
    site_tbl <- site_tbl %>%
      arrange(Protein.Id, abs_pos, site_id) %>%
      distinct(Protein.Id, site_id, abs_pos, abs_residue, .keep_all = TRUE)
  }
  
  motif_info <- lapply(seq_len(nrow(site_tbl)), function(i) {
    .extract_one_motif(
      protein_id = site_tbl$Protein.Id[i],
      abs_pos = site_tbl$abs_pos[i],
      abs_residue = site_tbl$abs_residue[i],
      fasta_map = fasta_map,
      flank = flank
    )
  }) %>%
    bind_rows()
  
  motif_tbl <- bind_cols(site_tbl, motif_info) %>%
    mutate(
      flank = flank,
      motif_width = 2 * flank + 1,
      center_position_in_motif = flank + 1,
      center_residue = abs_residue
    )
  
  motif_tbl
}

# ------------------------------------------------------------
# 5. Make amino acid position frequency matrix
# ------------------------------------------------------------
# 输入一组等长 motif 序列，输出每个位置每种氨基酸的频率矩阵。

make_aa_frequency_matrix <- function(seqs, flank = 6) {
  seqs <- as.character(seqs)
  seqs <- seqs[!is.na(seqs) & nzchar(seqs)]
  
  expected_width <- 2 * flank + 1
  seqs <- seqs[nchar(seqs) == expected_width]
  
  if (length(seqs) == 0) {
    return(data.frame())
  }
  
  aa_levels <- c(
    "A", "R", "N", "D", "C",
    "Q", "E", "G", "H", "I",
    "L", "K", "M", "F", "P",
    "S", "T", "W", "Y", "V"
  )
  
  pos_labels <- as.character(-flank:flank)
  pos_labels[pos_labels == "0"] <- "p0"
  
  char_mat <- do.call(
    rbind,
    strsplit(seqs, split = "")
  )
  
  freq_list <- lapply(seq_len(ncol(char_mat)), function(j) {
    tab <- table(factor(char_mat[, j], levels = aa_levels))
    as.numeric(tab) / sum(tab)
  })
  
  freq_mat <- do.call(cbind, freq_list)
  rownames(freq_mat) <- aa_levels
  colnames(freq_mat) <- pos_labels
  
  as.data.frame(freq_mat, check.names = FALSE) %>%
    rownames_to_column("AA")
}

# ------------------------------------------------------------
# 6. Plot sequence logo for one residue type
# ------------------------------------------------------------

plot_one_logo <- function(seqs,
                          title,
                          out_prefix,
                          flank = 6,
                          min_sequences = 5,
                          width = 8,
                          height = 4) {
  seqs <- as.character(seqs)
  seqs <- seqs[!is.na(seqs) & nzchar(seqs)]
  
  expected_width <- 2 * flank + 1
  seqs <- seqs[nchar(seqs) == expected_width]
  
  if (length(seqs) < min_sequences) {
    warning("序列数量少于 min_sequences=", min_sequences, "，跳过 logo: ", title)
    return(character(0))
  }
  
  p <- ggseqlogo::ggseqlogo(
    seqs,
    method = "prob",
    seq_type = "aa"
  ) +
    ggplot2::labs(
      title = title,
      x = paste0("Position relative to phosphosite, flank = ±", flank),
      y = "Probability"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  
  f_pdf <- paste0(out_prefix, ".pdf")
  f_png <- paste0(out_prefix, ".png")
  
  ggplot2::ggsave(f_pdf, p, width = width, height = height)
  ggplot2::ggsave(f_png, p, width = width, height = height, dpi = 300)
  
  c(f_pdf, f_png)
}

# ------------------------------------------------------------
# 7. Main wrapper
# ------------------------------------------------------------

run_phosphosite_motif_analysis <- function(site_mat2,
                                           fasta_file,
                                           out_dir = "Results/phosphosite_motif",
                                           flank = 6,
                                           residues = c("S", "T", "Y"),
                                           keep_unique_site = TRUE,
                                           make_logo = TRUE,
                                           min_sequences_for_logo = 5) {
  .safe_dir_create(out_dir)
  
  all_dir <- file.path(out_dir, "all_sites")
  freq_dir <- file.path(out_dir, "frequency_matrix")
  
  .safe_dir_create(all_dir)
  .safe_dir_create(freq_dir)
  
  # ----------------------------------------------------------
  # Step 1: 从 site_mat2 + FASTA 提取 ±flank aa motif 序列
  # ----------------------------------------------------------
  motif_tbl <- build_phosphosite_motif_table(
    site_mat2 = site_mat2,
    fasta_file = fasta_file,
    flank = flank,
    residues = residues,
    keep_unique_site = keep_unique_site
  )
  
  utils::write.csv(
    motif_tbl,
    file.path(out_dir, "01_all_phosphosite_flanking_sequences.csv"),
    row.names = FALSE,
    fileEncoding = "GBK"
  )
  
  # ----------------------------------------------------------
  # Step 2: 统计每种中心残基的位点数量、成功提取数量、失败原因
  # ----------------------------------------------------------
  summary_by_residue <- motif_tbl %>%
    group_by(abs_residue) %>%
    summarise(
      n_sites = n(),
      n_valid_motif = sum(valid_motif, na.rm = TRUE),
      n_invalid_motif = sum(!valid_motif, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(abs_residue)
  
  utils::write.csv(
    summary_by_residue,
    file.path(out_dir, "02_motif_summary_by_residue.csv"),
    row.names = FALSE,
    fileEncoding = "GBK"
  )
  
  fail_summary <- motif_tbl %>%
    filter(!valid_motif) %>%
    count(abs_residue, fail_reason, name = "n") %>%
    arrange(abs_residue, desc(n))
  
  utils::write.csv(
    fail_summary,
    file.path(out_dir, "03_motif_failed_reason_summary.csv"),
    row.names = FALSE,
    fileEncoding = "GBK"
  )
  
  # ----------------------------------------------------------
  # Step 3: 对 pS / pT / pY 分别输出序列、频率矩阵、logo
  # ----------------------------------------------------------
  plot_files <- character(0)
  freq_list <- list()
  
  for (res in residues) {
    res_name <- paste0("p", res)
    
    seqs_res <- motif_tbl %>%
      filter(valid_motif, abs_residue == res) %>%
      pull(motif_seq)
    
    # 保存序列
    seq_file <- file.path(all_dir, paste0("all_", res_name, "_motif_sequences.txt"))
    writeLines(seqs_res, con = seq_file)
    
    # 保存频率矩阵
    freq_res <- make_aa_frequency_matrix(seqs_res, flank = flank)
    freq_list[[res_name]] <- freq_res
    
    utils::write.csv(
      freq_res,
      file.path(freq_dir, paste0("all_", res_name, "_frequency_matrix.csv")),
      row.names = FALSE,
      fileEncoding = "GBK"
    )
    
    # 画 logo
    if (isTRUE(make_logo)) {
      out_prefix <- file.path(all_dir, paste0("all_", res_name, "_logo"))
      
      files_i <- plot_one_logo(
        seqs = seqs_res,
        title = paste0("All ", res_name, " phosphosite motifs"),
        out_prefix = out_prefix,
        flank = flank,
        min_sequences = min_sequences_for_logo
      )
      
      plot_files <- c(plot_files, files_i)
    }
  }
  
  utils::write.csv(
    data.frame(plot_file = unique(plot_files)),
    file.path(out_dir, "04_motif_plot_files.csv"),
    row.names = FALSE,
    fileEncoding = "GBK"
  )
  
  invisible(list(
    motif_table = motif_tbl,
    summary_by_residue = summary_by_residue,
    fail_summary = fail_summary,
    frequency_matrix = freq_list,
    plot_files = unique(plot_files)
  ))
}

# ============================================================
# Example usage
# ============================================================
# source("04_phosphosite_motif_analysis.R")
#
# fasta_file <- "C:/Work/SH/Pub_database/添加GO注释/data/UP000000589_17155_Reviewed_20230517_uniprot-download_true_format_fasta_includeIsoform_true_query__28Mus_-2023.05.17-06.21.58.02.fasta"
#
# motif <- run_phosphosite_motif_analysis(
#   site_mat2 = prep$site_mat2,
#   fasta_file = fasta_file,
#   out_dir = "./demo/Results/phosphosite_motif",
#   flank = 6,
#   keep_unique_site = TRUE,
#   make_logo = TRUE
# )



# ------------------------------------------------------------
# 8. Comparison-specific up/down trend motif analysis
# ------------------------------------------------------------
# Purpose:
#   Use de$all_results logFC to split phosphosites into
#   up-trend and down-trend groups for each comparison,
#   then draw pS / pT / pY sequence logos separately.
#
# Important:
#   For FC-only comparisons, these should be interpreted as
#   up-trend / down-trend phosphosites, not statistically significant sites.
# ------------------------------------------------------------
run_comparison_motif_analysis <- function(motif_table,
                                          de_results,
                                          out_dir = "Results/phosphosite_motif/by_comparison",
                                          comparison_col = "comparison",
                                          logfc_col = "logFC",
                                          logfc_thresh = 1,
                                          residues = c("S", "T", "Y"),
                                          min_sequences_for_logo = 5,
                                          make_logo = TRUE,
                                          feature_id_cols = c("Protein.Id", "Genes", "Residue.Both", "Modified.Sequence")) {
  .safe_dir_create(out_dir)
  
  motif_table <- as.data.frame(motif_table, check.names = FALSE, stringsAsFactors = FALSE)
  de_results <- as.data.frame(de_results, check.names = FALSE, stringsAsFactors = FALSE)
  
  # ----------------------------------------------------------
  # 1. 检查必要列
  # ----------------------------------------------------------
  required_motif_cols <- c(
    feature_id_cols,
    "site_id", "abs_pos", "abs_residue",
    "motif_seq", "valid_motif"
  )
  
  missing_motif_cols <- setdiff(required_motif_cols, colnames(motif_table))
  
  if (length(missing_motif_cols) > 0) {
    stop("motif_table 缺少这些列: ", paste(missing_motif_cols, collapse = ", "))
  }
  
  required_de_cols <- c("feature_id", comparison_col, logfc_col)
  missing_de_cols <- setdiff(required_de_cols, colnames(de_results))
  
  if (length(missing_de_cols) > 0) {
    stop("de_results 缺少这些列: ", paste(missing_de_cols, collapse = ", "))
  }
  
  # ----------------------------------------------------------
  # 2. 给 motif_table 构造和 de$all_results 一样的 feature_id
  # ----------------------------------------------------------
  motif_for_join <- motif_table %>%
    mutate(
      feature_id_raw = do.call(paste, c(across(all_of(feature_id_cols)), sep = "|")),
      feature_id = make.unique(feature_id_raw),
      site_id = as.character(site_id),
      abs_residue = toupper(as.character(abs_residue)),
      motif_seq = as.character(motif_seq)
    ) %>%
    filter(valid_motif, abs_residue %in% residues) %>%
    select(
      feature_id,
      feature_id_raw,
      site_id,
      Protein.Id,
      Genes,
      Residue.Both,
      Modified.Sequence,
      abs_pos,
      abs_residue,
      motif_seq,
      center_from_fasta,
      motif_start,
      motif_end
    )
  
  # ----------------------------------------------------------
  # 3. 用 feature_id 合并 DE/logFC 和 motif_seq
  # ----------------------------------------------------------
  de2 <- de_results %>%
    mutate(
      feature_id = as.character(feature_id),
      comparison = as.character(.data[[comparison_col]]),
      logFC = suppressWarnings(as.numeric(.data[[logfc_col]]))
    ) %>%
    filter(!is.na(comparison), comparison != "", !is.na(logFC)) %>%
    left_join(motif_for_join, by = "feature_id", suffix = c("", ".motif"))
  
  utils::write.csv(
    de2,
    file.path(out_dir, "01_de_results_with_motif_sequence.csv"),
    row.names = FALSE,
    fileEncoding = "GBK"
  )
  
  # ----------------------------------------------------------
  # 4. 匹配情况汇总
  # ----------------------------------------------------------
  match_summary <- de2 %>%
    group_by(comparison) %>%
    summarise(
      n_de_rows = n(),
      n_rows_with_valid_motif = sum(!is.na(motif_seq) & motif_seq != ""),
      n_unique_sites_with_valid_motif = n_distinct(site_id[!is.na(motif_seq) & motif_seq != ""]),
      .groups = "drop"
    ) %>%
    arrange(comparison)
  
  utils::write.csv(
    match_summary,
    file.path(out_dir, "02_comparison_motif_match_summary.csv"),
    row.names = FALSE,
    fileEncoding = "GBK"
  )
  
  # ----------------------------------------------------------
  # 5. 按 comparison + up/down 趋势生成序列、频率矩阵和 logo
  # ----------------------------------------------------------
  comparisons <- unique(de2$comparison)
  comparisons <- comparisons[!is.na(comparisons) & comparisons != ""]
  
  all_trend_tables <- list()
  plot_files <- character(0)
  
  for (cn in comparisons) {
    cn_slug <- .safe_slug(cn)
    comp_dir <- file.path(out_dir, cn_slug)
    .safe_dir_create(comp_dir)
    
    comp_df <- de2 %>%
      filter(comparison == cn, !is.na(motif_seq), motif_seq != "") %>%
      mutate(
        trend_direction = case_when(
          logFC >= logfc_thresh ~ "up_trend",
          logFC <= -logfc_thresh ~ "down_trend",
          TRUE ~ "not_selected"
        )
      )
    
    # 一个 site 可能来自多个 feature，这里按 site_id 合并，避免重复计数。
    trend_df <- comp_df %>%
      filter(trend_direction %in% c("up_trend", "down_trend")) %>%
      group_by(comparison, trend_direction, site_id) %>%
      summarise(
        Protein.Id = dplyr::first(Protein.Id),
        Genes = dplyr::first(Genes),
        abs_pos = dplyr::first(abs_pos),
        abs_residue = dplyr::first(abs_residue),
        motif_seq = dplyr::first(motif_seq),
        logFC = median(logFC, na.rm = TRUE),
        n_features_collapsed = n(),
        .groups = "drop"
      ) %>%
      arrange(trend_direction, abs_residue, desc(abs(logFC)))
    
    all_trend_tables[[cn]] <- trend_df
    
    utils::write.csv(
      trend_df,
      file.path(comp_dir, paste0(cn_slug, "_up_down_trend_motif_sequences.csv")),
      row.names = FALSE,
      fileEncoding = "GBK"
    )
    
    trend_summary <- trend_df %>%
      count(trend_direction, abs_residue, name = "n_sites") %>%
      arrange(trend_direction, abs_residue)
    
    utils::write.csv(
      trend_summary,
      file.path(comp_dir, paste0(cn_slug, "_up_down_trend_motif_summary.csv")),
      row.names = FALSE,
      fileEncoding = "GBK"
    )
    
    for (direction_i in c("up_trend", "down_trend")) {
      for (res in residues) {
        res_name <- paste0("p", res)
        
        seqs_i <- trend_df %>%
          filter(trend_direction == direction_i, abs_residue == res) %>%
          pull(motif_seq)
        
        seq_file <- file.path(
          comp_dir,
          paste0(cn_slug, "_", direction_i, "_", res_name, "_motif_sequences.txt")
        )
        
        writeLines(seqs_i, con = seq_file)
        
        if (length(seqs_i) > 0) {
          flank_i <- (nchar(seqs_i[1]) - 1) / 2
        } else {
          flank_i <- 6
        }
        
        freq_i <- make_aa_frequency_matrix(
          seqs = seqs_i,
          flank = flank_i
        )
        
        freq_file <- file.path(
          comp_dir,
          paste0(cn_slug, "_", direction_i, "_", res_name, "_frequency_matrix.csv")
        )
        
        utils::write.csv(
          freq_i,
          freq_file,
          row.names = FALSE,
          fileEncoding = "GBK"
        )
        
        if (isTRUE(make_logo)) {
          out_prefix <- file.path(
            comp_dir,
            paste0(cn_slug, "_", direction_i, "_", res_name, "_logo")
          )
          
          files_i <- plot_one_logo(
            seqs = seqs_i,
            title = paste0(cn, " | ", direction_i, " | ", res_name),
            out_prefix = out_prefix,
            flank = flank_i,
            min_sequences = min_sequences_for_logo
          )
          
          plot_files <- c(plot_files, files_i)
        }
      }
    }
  }
  
  all_trend_long <- bind_rows(all_trend_tables)
  
  utils::write.csv(
    all_trend_long,
    file.path(out_dir, "03_all_comparison_up_down_trend_motif_sequences.csv"),
    row.names = FALSE,
    fileEncoding = "GBK"
  )
  
  utils::write.csv(
    data.frame(plot_file = unique(plot_files)),
    file.path(out_dir, "04_comparison_motif_plot_files.csv"),
    row.names = FALSE,
    fileEncoding = "GBK"
  )
  
  invisible(list(
    de_results_with_motif = de2,
    match_summary = match_summary,
    trend_motif_table = all_trend_long,
    plot_files = unique(plot_files)
  ))
}