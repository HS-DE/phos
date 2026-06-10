# ============================================================
# 04b. R-only phosphosite motif enrichment analysis
#
# Purpose:
#   不依赖 MoMo / MEME Suite。
#   直接在 R 里基于 phosphosite-centered 13-mer 序列做 motif enrichment。
#
# Main outputs:
#   1. 预测 Motif 对应磷酸化修饰位点数量图 Top20
#   2. 预测保守 Motif 富集统计图 Top20
#
# Core idea:
#   foreground = 目标 phosphosite motif sequences
#   background = 背景 phosphosite motif sequences 或 shuffled sequences
#
#   对每个 position-AA 或 position-AA 组合计算：
#     - foreground matched count
#     - background matched count
#     - fold increase
#     - Fisher exact test p value
#     - BH adjusted q value
#
# Inputs:
#   motif$motif_table
#   motif_comp$trend_motif_table
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(readr)
  library(ggplot2)
})

# ------------------------------------------------------------
# 0. Small utilities
# ------------------------------------------------------------

.safe_dir_create <- function(path) {
  if (!is.null(path) && nzchar(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
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

.get_valid_motif_table <- function(x,
                                   flank = 6,
                                   require_valid_motif = TRUE) {
  x <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
  
  required_cols <- c("site_id", "motif_seq", "abs_residue")
  missing_cols <- setdiff(required_cols, colnames(x))
  
  if (length(missing_cols) > 0) {
    stop("输入表缺少这些列: ", paste(missing_cols, collapse = ", "))
  }
  
  expected_width <- 2 * flank + 1
  
  out <- x %>%
    mutate(
      motif_seq = toupper(as.character(motif_seq)),
      abs_residue = toupper(as.character(abs_residue))
    ) %>%
    filter(
      !is.na(site_id),
      !is.na(motif_seq),
      motif_seq != "",
      nchar(motif_seq) == expected_width,
      abs_residue %in% c("S", "T", "Y"),
      grepl("^[ARNDCQEGHILKMFPSTWYV]+$", motif_seq)
    )
  
  if (require_valid_motif && "valid_motif" %in% colnames(out)) {
    out <- out %>% filter(valid_motif)
  }
  
  # 一个 phosphosite 只保留一次，避免同一 site 被多个肽段重复计数。
  out %>%
    distinct(site_id, .keep_all = TRUE)
}

.seqs_to_matrix <- function(seqs) {
  seqs <- as.character(seqs)
  mat <- do.call(rbind, strsplit(seqs, split = "", fixed = TRUE))
  colnames(mat) <- seq_len(ncol(mat))
  mat
}

.make_relative_positions <- function(flank = 6) {
  rel_pos <- -flank:flank
  tibble(
    index = seq_along(rel_pos),
    rel_pos = rel_pos
  )
}

.make_motif_label <- function(rel_pos_vec,
                              aa_vec,
                              center_residue,
                              flank = 6) {
  pos_map <- .make_relative_positions(flank)
  tokens <- rep(".", nrow(pos_map))
  
  center_i <- which(pos_map$rel_pos == 0)
  tokens[center_i] <- paste0(center_residue, "_P")
  
  for (i in seq_along(rel_pos_vec)) {
    idx <- pos_map$index[pos_map$rel_pos == rel_pos_vec[i]]
    tokens[idx] <- aa_vec[i]
  }
  
  paste0(tokens, collapse = "")
}

.match_conditions <- function(mat, idx_vec, aa_vec) {
  hit <- rep(TRUE, nrow(mat))
  
  for (i in seq_along(idx_vec)) {
    hit <- hit & mat[, idx_vec[i]] == aa_vec[i]
  }
  
  hit
}

.calc_one_motif <- function(fg_mat,
                            bg_mat,
                            rel_pos_vec,
                            aa_vec,
                            center_residue,
                            set_name,
                            flank = 6,
                            pseudocount = 0.5) {
  pos_map <- .make_relative_positions(flank)
  
  idx_vec <- pos_map$index[match(rel_pos_vec, pos_map$rel_pos)]
  
  fg_hit <- .match_conditions(fg_mat, idx_vec, aa_vec)
  bg_hit <- .match_conditions(bg_mat, idx_vec, aa_vec)
  
  n_fg <- sum(fg_hit)
  n_bg <- sum(bg_hit)
  
  fg_total <- nrow(fg_mat)
  bg_total <- nrow(bg_mat)
  
  fg_freq <- n_fg / fg_total
  bg_freq <- n_bg / bg_total
  
  # 用 pseudocount 避免 background 为 0 时 fold increase = Inf
  fold_increase <- ((n_fg + pseudocount) / (fg_total + 2 * pseudocount)) /
    ((n_bg + pseudocount) / (bg_total + 2 * pseudocount))
  
  fisher_mat <- matrix(
    c(
      n_fg, fg_total - n_fg,
      n_bg, bg_total - n_bg
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  p_value <- tryCatch(
    fisher.test(fisher_mat, alternative = "greater")$p.value,
    error = function(e) NA_real_
  )
  
  tibble(
    set_name = set_name,
    center_residue = center_residue,
    motif_order = length(rel_pos_vec),
    motif_label = .make_motif_label(
      rel_pos_vec = rel_pos_vec,
      aa_vec = aa_vec,
      center_residue = center_residue,
      flank = flank
    ),
    motif_positions = paste(rel_pos_vec, collapse = ";"),
    motif_aas = paste(aa_vec, collapse = ";"),
    n_foreground_matches = n_fg,
    n_background_matches = n_bg,
    foreground_size = fg_total,
    background_size = bg_total,
    foreground_frequency = fg_freq,
    background_frequency = bg_freq,
    fold_increase = fold_increase,
    p_value = p_value
  )
}

.make_shuffle_background <- function(fg_tbl,
                                     flank = 6,
                                     n_shuffle = 5) {
  center_idx <- flank + 1
  seqs <- fg_tbl$motif_seq
  
  out <- vector("character", length(seqs) * n_shuffle)
  out_site_id <- vector("character", length(seqs) * n_shuffle)
  out_residue <- vector("character", length(seqs) * n_shuffle)
  
  k <- 1
  
  for (i in seq_along(seqs)) {
    chars <- strsplit(seqs[i], split = "", fixed = TRUE)[[1]]
    flank_idx <- setdiff(seq_along(chars), center_idx)
    
    for (j in seq_len(n_shuffle)) {
      chars2 <- chars
      chars2[flank_idx] <- sample(chars2[flank_idx], length(flank_idx), replace = FALSE)
      
      out[k] <- paste0(chars2, collapse = "")
      out_site_id[k] <- paste0(fg_tbl$site_id[i], "_shuffle", j)
      out_residue[k] <- fg_tbl$abs_residue[i]
      
      k <- k + 1
    }
  }
  
  tibble(
    site_id = out_site_id,
    motif_seq = out,
    abs_residue = out_residue
  )
}

# ------------------------------------------------------------
# 1. Build motif enrichment table for one foreground set
# ------------------------------------------------------------

run_r_motif_enrichment_one_set <- function(foreground_tbl,
                                           background_tbl = NULL,
                                           set_name,
                                           out_dir,
                                           flank = 6,
                                           background_mode = c("shuffle", "provided"),
                                           n_shuffle = 5,
                                           max_order = 2,
                                           top_seed = 30,
                                           min_count = 5,
                                           min_fold = 1.2,
                                           seed_p_cutoff = 0.05,
                                           top_n = 20,
                                           residues = c("S", "T", "Y")) {
  background_mode <- match.arg(background_mode)
  
  .safe_dir_create(out_dir)
  
  table_dir <- file.path(out_dir, "tables", .safe_slug(set_name))
  plot_dir <- file.path(out_dir, "plots", .safe_slug(set_name))
  
  .safe_dir_create(table_dir)
  .safe_dir_create(plot_dir)
  
  fg_tbl <- .get_valid_motif_table(foreground_tbl, flank = flank)
  
  if (nrow(fg_tbl) == 0) {
    warning("foreground 为空: ", set_name)
    return(list(
      set_name = set_name,
      status = "empty_foreground",
      enrichment = data.frame(),
      plot_files = character(0)
    ))
  }
  
  if (background_mode == "shuffle") {
    bg_tbl <- .make_shuffle_background(
      fg_tbl = fg_tbl,
      flank = flank,
      n_shuffle = n_shuffle
    )
  } else {
    if (is.null(background_tbl)) {
      stop("background_mode = 'provided' 时必须提供 background_tbl")
    }
    
    bg_tbl <- .get_valid_motif_table(background_tbl, flank = flank)
  }
  
  pos_map <- .make_relative_positions(flank)
  pos_map <- pos_map %>% filter(rel_pos != 0)
  
  aa_vec_all <- c(
    "A", "R", "N", "D", "C",
    "Q", "E", "G", "H", "I",
    "L", "K", "M", "F", "P",
    "S", "T", "W", "Y", "V"
  )
  
  all_res <- list()
  
  for (res in residues) {
    fg_res <- fg_tbl %>% filter(abs_residue == res)
    bg_res <- bg_tbl %>% filter(abs_residue == res)
    
    if (nrow(fg_res) < min_count || nrow(bg_res) < min_count) {
      next
    }
    
    fg_mat <- .seqs_to_matrix(fg_res$motif_seq)
    bg_mat <- .seqs_to_matrix(bg_res$motif_seq)
    
    # --------------------------------------------------------
    # A. 单位置 motif: position-AA
    # --------------------------------------------------------
    single_stats <- list()
    m <- 1
    
    for (rp in pos_map$rel_pos) {
      for (aa in aa_vec_all) {
        single_stats[[m]] <- .calc_one_motif(
          fg_mat = fg_mat,
          bg_mat = bg_mat,
          rel_pos_vec = rp,
          aa_vec = aa,
          center_residue = res,
          set_name = set_name,
          flank = flank
        )
        m <- m + 1
      }
    }
    
    single_stats <- bind_rows(single_stats)
    
    # --------------------------------------------------------
    # B. 选择 seed features，用于构造二阶 / 三阶组合 motif
    # --------------------------------------------------------
    seed_df <- single_stats %>%
      filter(
        n_foreground_matches >= min_count,
        fold_increase >= min_fold,
        !is.na(p_value)
      ) %>%
      mutate(seed_rank_score = -log10(p_value + 1e-300) * log2(fold_increase + 1e-6)) %>%
      arrange(desc(seed_rank_score), p_value, desc(fold_increase)) %>%
      slice_head(n = top_seed)
    
    # 如果阈值过严导致没有 seed，就退而求其次取最靠前的单点特征。
    if (nrow(seed_df) == 0) {
      seed_df <- single_stats %>%
        filter(n_foreground_matches >= min_count, fold_increase > 1) %>%
        arrange(p_value, desc(fold_increase)) %>%
        slice_head(n = top_seed)
    }
    
    combo_stats <- list()
    
    # --------------------------------------------------------
    # C. 二阶 / 三阶 motif
    # --------------------------------------------------------
    if (max_order >= 2 && nrow(seed_df) >= 2) {
      for (ord in 2:max_order) {
        if (nrow(seed_df) < ord) {
          next
        }
        
        cb <- combn(seq_len(nrow(seed_df)), ord, simplify = FALSE)
        
        for (idx in cb) {
          rp_vec <- as.integer(strsplit(seed_df$motif_positions[idx], ";") %>% unlist())
          aa_vec <- strsplit(seed_df$motif_aas[idx], ";") %>% unlist()
          
          # 不允许同一个相对位置出现多个限制条件。
          if (length(unique(rp_vec)) < length(rp_vec)) {
            next
          }
          
          # 按相对位置排序，让 motif label 稳定。
          ord_idx <- order(rp_vec)
          rp_vec <- rp_vec[ord_idx]
          aa_vec <- aa_vec[ord_idx]
          
          one <- .calc_one_motif(
            fg_mat = fg_mat,
            bg_mat = bg_mat,
            rel_pos_vec = rp_vec,
            aa_vec = aa_vec,
            center_residue = res,
            set_name = set_name,
            flank = flank
          )
          
          if (one$n_foreground_matches >= min_count && one$fold_increase > 1) {
            combo_stats[[length(combo_stats) + 1]] <- one
          }
        }
      }
    }
    
    combo_stats <- bind_rows(combo_stats)
    
    res_tbl <- bind_rows(single_stats, combo_stats) %>%
      distinct(motif_label, .keep_all = TRUE) %>%
      mutate(
        q_value = p.adjust(p_value, method = "BH")
      ) %>%
      arrange(q_value, p_value, desc(fold_increase), desc(n_foreground_matches))
    
    all_res[[res]] <- res_tbl
  }
  
  enrichment <- bind_rows(all_res)
  
  if (nrow(enrichment) == 0) {
    warning("没有得到 motif enrichment 结果: ", set_name)
    
    return(list(
      set_name = set_name,
      status = "empty_enrichment",
      enrichment = data.frame(),
      plot_files = character(0)
    ))
  }
  
  enrich_file <- file.path(table_dir, paste0(.safe_slug(set_name), "_motif_enrichment_all.csv"))
  
  utils::write.csv(
    enrichment,
    enrich_file,
    row.names = FALSE,
    fileEncoding = "GBK"
  )
  
  plot_files <- plot_r_motif_enrichment_top20(
    enrichment = enrichment,
    set_name = set_name,
    plot_dir = plot_dir,
    table_dir = table_dir,
    top_n = top_n,
    min_count = min_count
  )
  
  list(
    set_name = set_name,
    status = "ok",
    enrichment_file = enrich_file,
    enrichment = enrichment,
    plot_files = plot_files
  )
}

# ------------------------------------------------------------
# 2. Plot Top20 horizontal barplots
# ------------------------------------------------------------

.plot_one_r_motif_bar <- function(plot_df,
                                  value_col,
                                  xlab,
                                  title,
                                  out_prefix,
                                  top_n = 20) {
  if (nrow(plot_df) == 0) {
    return(character(0))
  }
  
  plot_df <- plot_df %>%
    filter(!is.na(.data[[value_col]]), is.finite(.data[[value_col]])) %>%
    arrange(desc(.data[[value_col]]), q_value, p_value) %>%
    slice_head(n = top_n) %>%
    mutate(
      motif_label_plot = factor(motif_label, levels = rev(motif_label))
    )
  
  if (nrow(plot_df) == 0) {
    return(character(0))
  }
  
  p <- ggplot(plot_df, aes(x = motif_label_plot, y = .data[[value_col]])) +
    geom_col(width = 0.75) +
    coord_flip() +
    labs(
      title = title,
      x = NULL,
      y = xlab
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 8)
    )
  
  f_pdf <- paste0(out_prefix, ".pdf")
  f_png <- paste0(out_prefix, ".png")
  
  ggsave(f_pdf, p, width = 7.8, height = max(5, 0.28 * nrow(plot_df) + 1.5))
  ggsave(f_png, p, width = 7.8, height = max(5, 0.28 * nrow(plot_df) + 1.5), dpi = 300)
  
  c(f_pdf, f_png)
}

plot_r_motif_enrichment_top20 <- function(enrichment,
                                          set_name,
                                          plot_dir,
                                          table_dir,
                                          top_n = 20,
                                          min_count = 5) {
  .safe_dir_create(plot_dir)
  .safe_dir_create(table_dir)
  
  set_slug <- .safe_slug(set_name)
  
  # 对应：预测 Motif 对应磷酸化修饰位点数量图
  count_df <- enrichment %>%
    filter(
      n_foreground_matches >= min_count,
      fold_increase > 1
    ) %>%
    arrange(desc(n_foreground_matches), q_value, p_value)
  
  count_top <- count_df %>% slice_head(n = top_n)
  
  utils::write.csv(
    count_top,
    file.path(table_dir, paste0(set_slug, "_top", top_n, "_motif_site_count.csv")),
    row.names = FALSE,
    fileEncoding = "GBK"
  )
  
  # 对应：预测保守 Motif 富集统计图
  # 为避免极少数低 count motif 因 background 很低导致 fold 虚高，
  # 这里仍然要求 n_foreground_matches >= min_count。
  fold_df <- enrichment %>%
    filter(
      n_foreground_matches >= min_count,
      fold_increase > 1
    ) %>%
    arrange(desc(fold_increase), q_value, p_value)
  
  fold_top <- fold_df %>% slice_head(n = top_n)
  
  utils::write.csv(
    fold_top,
    file.path(table_dir, paste0(set_slug, "_top", top_n, "_motif_fold_increase.csv")),
    row.names = FALSE,
    fileEncoding = "GBK"
  )
  
  saved <- character(0)
  
  saved <- c(
    saved,
    .plot_one_r_motif_bar(
      plot_df = count_df,
      value_col = "n_foreground_matches",
      xlab = "Number of matched phosphosites",
      title = paste0(set_name, " | Motif matched phosphosite count"),
      out_prefix = file.path(plot_dir, paste0(set_slug, "_top", top_n, "_motif_site_count")),
      top_n = top_n
    )
  )
  
  saved <- c(
    saved,
    .plot_one_r_motif_bar(
      plot_df = fold_df,
      value_col = "fold_increase",
      xlab = "Fold increase",
      title = paste0(set_name, " | Motif enrichment"),
      out_prefix = file.path(plot_dir, paste0(set_slug, "_top", top_n, "_motif_fold_increase")),
      top_n = top_n
    )
  )
  
  unique(saved)
}

# ------------------------------------------------------------
# 3. Main pipeline
# ------------------------------------------------------------

run_r_motif_enrichment_pipeline <- function(motif_table,
                                            trend_motif_table = NULL,
                                            out_dir = "Results/phosphosite_motif_enrichment_Ronly",
                                            flank = 6,
                                            run_all_sites = TRUE,
                                            run_all_sites_by_residue = TRUE,
                                            run_comparison_trends = TRUE,
                                            n_shuffle = 5,
                                            max_order = 2,
                                            top_seed = 30,
                                            min_count = 5,
                                            min_fold = 1.2,
                                            seed_p_cutoff = 0.05,
                                            top_n = 20) {
  .safe_dir_create(out_dir)
  
  motif_tbl <- .get_valid_motif_table(
    motif_table,
    flank = flank
  )
  
  result_list <- list()
  
  # ----------------------------------------------------------
  # A. 所有位点：用 shuffled background
  #    这个回答：整体 phosphosite motif 是否有位置偏好？
  # ----------------------------------------------------------
  if (isTRUE(run_all_sites)) {
    result_list[["all_sites"]] <- run_r_motif_enrichment_one_set(
      foreground_tbl = motif_tbl,
      background_tbl = NULL,
      set_name = "all_sites",
      out_dir = out_dir,
      flank = flank,
      background_mode = "shuffle",
      n_shuffle = n_shuffle,
      max_order = max_order,
      top_seed = top_seed,
      min_count = min_count,
      min_fold = min_fold,
      seed_p_cutoff = seed_p_cutoff,
      top_n = top_n,
      residues = c("S", "T", "Y")
    )
  }
  
  # ----------------------------------------------------------
  # B. 所有位点按 pS / pT / pY 分开
  # ----------------------------------------------------------
  if (isTRUE(run_all_sites_by_residue)) {
    for (res in c("S", "T", "Y")) {
      fg_res <- motif_tbl %>% filter(abs_residue == res)
      
      if (nrow(fg_res) < min_count) {
        next
      }
      
      set_name <- paste0("all_p", res, "_sites")
      
      result_list[[set_name]] <- run_r_motif_enrichment_one_set(
        foreground_tbl = fg_res,
        background_tbl = NULL,
        set_name = set_name,
        out_dir = out_dir,
        flank = flank,
        background_mode = "shuffle",
        n_shuffle = n_shuffle,
        max_order = max_order,
        top_seed = top_seed,
        min_count = min_count,
        min_fold = min_fold,
        seed_p_cutoff = seed_p_cutoff,
        top_n = top_n,
        residues = res
      )
    }
  }
  
  # ----------------------------------------------------------
  # C. 每个 comparison 的 up/down trend
  #
  # foreground:
  #   某 comparison 的 up_trend / down_trend sites
  #
  # background:
  #   本项目所有检测到的 phosphosites
  #
  # 在 run_r_motif_enrichment_one_set 内部会按 pS/pT/pY 分开算，
  # 所以 pS foreground 会和 pS background 比，pT 和 pT 比。
  # ----------------------------------------------------------
  if (isTRUE(run_comparison_trends) && !is.null(trend_motif_table)) {
    trend_tbl <- .get_valid_motif_table(
      trend_motif_table,
      flank = flank,
      require_valid_motif = FALSE
    )
    
    required_cols <- c("comparison", "trend_direction")
    missing_cols <- setdiff(required_cols, colnames(trend_tbl))
    
    if (length(missing_cols) > 0) {
      stop("trend_motif_table 缺少这些列: ", paste(missing_cols, collapse = ", "))
    }
    
    set_info <- trend_tbl %>%
      distinct(comparison, trend_direction, site_id, .keep_all = TRUE) %>%
      count(comparison, trend_direction, name = "n_sites")
    
    utils::write.csv(
      set_info,
      file.path(out_dir, "00_comparison_trend_set_size.csv"),
      row.names = FALSE,
      fileEncoding = "GBK"
    )
    
    for (i in seq_len(nrow(set_info))) {
      cn <- set_info$comparison[i]
      dir_i <- set_info$trend_direction[i]
      
      fg_i <- trend_tbl %>%
        filter(comparison == cn, trend_direction == dir_i) %>%
        distinct(site_id, .keep_all = TRUE)
      
      if (nrow(fg_i) < min_count) {
        next
      }
      
      set_name <- paste(cn, dir_i, sep = "__")
      set_slug <- .safe_slug(set_name)
      
      result_list[[set_slug]] <- run_r_motif_enrichment_one_set(
        foreground_tbl = fg_i,
        background_tbl = motif_tbl,
        set_name = set_name,
        out_dir = out_dir,
        flank = flank,
        background_mode = "provided",
        n_shuffle = n_shuffle,
        max_order = max_order,
        top_seed = top_seed,
        min_count = min_count,
        min_fold = min_fold,
        seed_p_cutoff = seed_p_cutoff,
        top_n = top_n,
        residues = c("S", "T", "Y")
      )
    }
  }
  
  # ----------------------------------------------------------
  # D. 汇总
  # ----------------------------------------------------------
  run_summary <- bind_rows(lapply(result_list, function(x) {
    tibble(
      set_name = x$set_name,
      status = x$status,
      enrichment_file = if (!is.null(x$enrichment_file)) x$enrichment_file else NA_character_,
      n_motifs = if (!is.null(x$enrichment)) nrow(x$enrichment) else 0
    )
  }))
  
  utils::write.csv(
    run_summary,
    file.path(out_dir, "00_r_motif_enrichment_run_summary.csv"),
    row.names = FALSE,
    fileEncoding = "GBK"
  )
  
  enrichment_all <- bind_rows(lapply(result_list, function(x) {
    if (!is.null(x$enrichment) && nrow(x$enrichment) > 0) {
      x$enrichment
    } else {
      data.frame()
    }
  }))
  
  utils::write.csv(
    enrichment_all,
    file.path(out_dir, "01_all_r_motif_enrichment_results.csv"),
    row.names = FALSE,
    fileEncoding = "GBK"
  )
  
  plot_files <- unique(unlist(lapply(result_list, `[[`, "plot_files"), use.names = FALSE))
  
  utils::write.csv(
    data.frame(plot_file = plot_files),
    file.path(out_dir, "02_all_r_motif_enrichment_plot_files.csv"),
    row.names = FALSE,
    fileEncoding = "GBK"
  )
  
  invisible(list(
    run_summary = run_summary,
    enrichment_all = enrichment_all,
    result_list = result_list,
    plot_files = plot_files
  ))
}