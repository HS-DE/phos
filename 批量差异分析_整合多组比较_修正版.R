# ============================================================
# Phosphosite-friendly batch DE / FC-only analysis
#
# Function kept for compatibility:
#   run_batch_de_mixed()
#
# Expected input:
#   expr        : numeric matrix, rownames = feature_id, colnames = sample_id
#   sample_info : data.frame, contains sample_id and group information
#
# Important:
#   For phosphoproteomics, rownames are feature_id, not simple Gene symbols.
#   Feature IDs are preserved exactly.
#
# Single-sample / low-replicate support:
#   If either group has fewer than min_reps_for_limma samples, the function uses
#   group_FC_only. In that mode, logFC is calculated but P.Value and adj.P.Val
#   are intentionally NA.
# ============================================================

run_batch_de_mixed <- function(
  expr,
  sample_info,
  compare_list = NULL,
  multi_compare = NULL,
  group_col = "group",
  exp_col = NULL,
  ctrl_col = NULL,
  multi_compare_group_cols = NULL,
  multi_group_min_n = 3,
  min_reps_for_multi_limma = 2,
  multi_top_n = 50,
  min_reps_for_limma = 3,
  low_rep_strategy = c("group_fc_only", "limma"),
  base_out_dir = "../Results/差异分析",
  create_dir = TRUE,
  p_type = c("adjp", "pvalue"),
  p_thresh = 0.05,
  logfc_thresh = 1,
  adjust_method = "BH",
  save_formats = c("pdf", "png", "tiff"),
  tiff_compression = "lzw",
  make_volcano = TRUE,
  label_extreme_n_each_side = 10,
  plot_width = 10,
  plot_height = 8,
  plot_dpi = 300,
  heat_show_rownames = FALSE,
  heat_show_colnames = FALSE,
  cluster = FALSE,
  assume_log2 = NULL,
  add1 = 1,
  sample_id_col = NULL,
  export_excel = TRUE,
  export_deg_summary = TRUE
) {
  suppressPackageStartupMessages({
    requireNamespace("dplyr", quietly = TRUE)
    requireNamespace("tidyr", quietly = TRUE)
  })

  # -------------------------
  # 0. Basic checks
  # -------------------------
  if (!is.matrix(expr)) expr <- as.matrix(expr)
  if (is.null(rownames(expr))) stop("expr 需要行名；在磷酸化流程中行名应为 feature_id。")
  if (is.null(colnames(expr))) stop("expr 需要列名；列名应为 sample_id。")

  if (!is.numeric(expr)) {
    expr_num <- apply(expr, 2, function(x) suppressWarnings(as.numeric(x)))
    rownames(expr_num) <- rownames(expr)
    colnames(expr_num) <- colnames(expr)
    expr <- as.matrix(expr_num)
  }

  sample_info <- as.data.frame(sample_info, check.names = FALSE, stringsAsFactors = FALSE)

  if (!is.null(sample_id_col)) {
    if (!sample_id_col %in% colnames(sample_info)) stop("sample_id_col 不在 sample_info 中：", sample_id_col)
    rownames(sample_info) <- as.character(sample_info[[sample_id_col]])
  } else if (
    is.null(rownames(sample_info)) ||
      any(is.na(rownames(sample_info))) ||
      any(rownames(sample_info) == "") ||
      identical(rownames(sample_info), as.character(seq_len(nrow(sample_info))))
  ) {
    candidate_cols <- c("sample_id", "sample_name", "SampleID", "Sample_name", "样品编号", "样品信息")
    hit <- intersect(candidate_cols, colnames(sample_info))
    if (length(hit) == 0) {
      stop("sample_info 没有可用行名，也找不到 sample_id/sample_name/样品编号/样品信息 列。")
    }
    rownames(sample_info) <- as.character(sample_info[[hit[1]]])
  }

  if (anyDuplicated(rownames(sample_info)) > 0) {
    stop("sample_info 行名/样本 ID 有重复，请先处理重复样本名。")
  }

  if (!group_col %in% colnames(sample_info)) {
    stop("sample_info 中找不到 group_col='", group_col, "'。当前列：", paste(colnames(sample_info), collapse = ", "))
  }

  low_rep_strategy <- match.arg(low_rep_strategy)
  p_type <- tolower(p_type[1])
  if (p_type %in% c("p", "pval")) p_type <- "pvalue"
  if (!p_type %in% c("adjp", "pvalue")) stop("p_type 只能是 adjp 或 pvalue。")

  if (isTRUE(create_dir)) dir.create(base_out_dir, recursive = TRUE, showWarnings = FALSE)

  if (!exists("limma_de", mode = "function")) {
    local_limma_file <- "limma差异.R"
    if (file.exists(local_limma_file)) source(local_limma_file)
  }

  # -------------------------
  # 1. Helpers
  # -------------------------
  safe_slug <- function(x) {
    x <- as.character(x)
    x <- gsub("[[:space:]]+", "_", x)
    x <- gsub("[/\\\\:;*?\"<>|]", "_", x)
    x <- gsub("\\+", "plus", x)
    x <- gsub("[^[:alnum:]_.-]", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
  }

  get_group_samples <- function(g) {
    s0 <- rownames(sample_info)[as.character(sample_info[[group_col]]) == as.character(g)]
    unique(intersect(colnames(expr), s0))
  }

  guess_is_log2 <- function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) return(TRUE)
    !(max(v) > 1000 || stats::median(v) > 100)
  }

  get_compare_cols <- function(df) {
    if (!is.null(exp_col) && !is.null(ctrl_col)) {
      if (!exp_col %in% names(df) || !ctrl_col %in% names(df)) {
        stop("手动指定的 exp_col/ctrl_col 不在 compare_list 列名中。")
      }
      return(list(exp_col = exp_col, ctrl_col = ctrl_col))
    }

    cn_exp <- intersect(names(df), c("实验", "exp", "group1", "treat", "case", "B"))
    cn_ctrl <- intersect(names(df), c("对照", "ctrl", "group2", "control", "A"))
    if (length(cn_exp) == 0 || length(cn_ctrl) == 0) {
      stop("compare_list 必须包含实验/对照或 exp/ctrl 列。当前列名：", paste(names(df), collapse = ", "))
    }
    list(exp_col = cn_exp[1], ctrl_col = cn_ctrl[1])
  }

  scale_rows <- function(mat) {
    row_mean <- rowMeans(mat, na.rm = TRUE)
    row_sd <- apply(mat, 1, stats::sd, na.rm = TRUE)
    row_sd[is.na(row_sd) | row_sd == 0] <- 1
    out <- sweep(mat, 1, row_mean, FUN = "-")
    out <- sweep(out, 1, row_sd, FUN = "/")
    out[!is.finite(out)] <- 0
    out
  }

  # -------------------------
  # 2. FC-only for low replicate comparisons
  # -------------------------
  group_fc_only_de <- function(group1, group2, save_prefix = NULL, summary_fun = c("mean", "median")) {
    summary_fun <- match.arg(summary_fun)
    s1 <- get_group_samples(group1)
    s2 <- get_group_samples(group2)

    if (length(s1) == 0 || length(s2) == 0) {
      stop("比较组匹配不到样本：", group1, " n=", length(s1), "; ", group2, " n=", length(s2))
    }

    X1 <- expr[, s1, drop = FALSE]
    X2 <- expr[, s2, drop = FALSE]

    is_log2 <- if (is.null(assume_log2)) guess_is_log2(c(as.numeric(X1), as.numeric(X2))) else isTRUE(assume_log2)
    if (!is_log2) {
      X1 <- log2(X1 + add1)
      X2 <- log2(X2 + add1)
    }

    if (summary_fun == "mean") {
      g1 <- rowMeans(X1, na.rm = TRUE)
      g2 <- rowMeans(X2, na.rm = TRUE)
    } else {
      g1 <- apply(X1, 1, function(z) if (all(is.na(z))) NA_real_ else stats::median(z, na.rm = TRUE))
      g2 <- apply(X2, 1, function(z) if (all(is.na(z))) NA_real_ else stats::median(z, na.rm = TRUE))
    }
    g1[is.nan(g1)] <- NA_real_
    g2[is.nan(g2)] <- NA_real_

    # logFC = group1 - group2, same direction as limma_de(group1, group2)
    logFC <- g1 - g2
    AveExpr <- (g1 + g2) / 2

    direction <- ifelse(logFC >= logfc_thresh, "Up", ifelse(logFC <= -logfc_thresh, "Down", "Not Significant"))
    direction <- factor(direction, levels = c("Up", "Not Significant", "Down"))

    tt <- data.frame(
      feature_id = rownames(expr),
      logFC = as.numeric(logFC),
      AveExpr = as.numeric(AveExpr),
      P.Value = NA_real_,
      adj.P.Val = NA_real_,
      direction = direction,
      comparison = paste0(group1, "_vs_", group2),
      group1 = group1,
      group2 = group2,
      n_group1 = length(s1),
      n_group2 = length(s2),
      method = "group_FC_only",
      stringsAsFactors = FALSE
    )
    rownames(tt) <- tt$feature_id

    combined <- cbind(tt, group1_stat = g1, group2_stat = g2, expr[, c(s1, s2), drop = FALSE])
    deg <- tt$feature_id[tt$direction %in% c("Up", "Down")]
    saved_files <- character(0)

    if (!is.null(save_prefix) && nzchar(save_prefix)) {
      if (isTRUE(create_dir)) dir.create(dirname(save_prefix), recursive = TRUE, showWarnings = FALSE)
      f_csv <- paste0(save_prefix, "_FC_only_results.csv")
      f_comb <- paste0(save_prefix, "_combined.csv")
      utils::write.csv(tt, f_csv, row.names = FALSE, fileEncoding = "GBK")
      utils::write.csv(combined, f_comb, row.names = FALSE, fileEncoding = "GBK")
      saved_files <- c(saved_files, f_csv, f_comb)
    }

    volcano_plot <- NULL
    if (isTRUE(make_volcano) && requireNamespace("ggplot2", quietly = TRUE)) {
      volcano_df <- tt
      volcano_df$yval <- abs(volcano_df$logFC)
      volcano_plot <- ggplot2::ggplot(volcano_df, ggplot2::aes(x = logFC, y = yval, colour = direction)) +
        ggplot2::geom_point(alpha = 0.45, size = 2.4) +
        ggplot2::scale_color_manual(values = c("Down" = "#67a9cf", "Not Significant" = "#cccccc", "Up" = "#fd8d3c"), drop = TRUE) +
        ggplot2::geom_vline(xintercept = c(-logfc_thresh, logfc_thresh), linetype = 4, color = "black") +
        ggplot2::geom_hline(yintercept = logfc_thresh, linetype = 4, color = "black") +
        ggplot2::labs(x = paste0("log2FC (", group1, " - ", group2, ")"), y = "|log2FC| (FC-only)", color = "") +
        ggplot2::theme_bw() +
        ggplot2::theme(aspect.ratio = 1, legend.position = "right")

      if (!is.null(save_prefix) && nzchar(save_prefix)) {
        for (fmt in unique(tolower(save_formats))) {
          f <- paste0(save_prefix, "_volcano.", fmt)
          if (fmt %in% c("tif", "tiff")) {
            ggplot2::ggsave(f, volcano_plot, width = plot_width, height = plot_height,
                            dpi = plot_dpi, device = "tiff", compression = tiff_compression)
          } else {
            ggplot2::ggsave(f, volcano_plot, width = plot_width, height = plot_height, dpi = plot_dpi)
          }
          saved_files <- c(saved_files, f)
        }
      }
    }

    heatmap_plot <- NULL
    heatmap_matrix <- NULL
    heatmap_matrix_row_scaled <- NULL
    if (isTRUE(make_volcano) && requireNamespace("pheatmap", quietly = TRUE) && requireNamespace("grid", quietly = TRUE)) {
      sig_features <- tt$feature_id[tt$direction %in% c("Up", "Down")]
      sig_features <- intersect(sig_features, rownames(expr))
      if (length(sig_features) > 1) {
        heat_cols <- c(s2, s1)
        heatmap_matrix <- expr[sig_features, heat_cols, drop = FALSE]
        heatmap_matrix_row_scaled <- scale_rows(heatmap_matrix)
        ann_col <- data.frame(
          Group = factor(c(rep(group2, length(s2)), rep(group1, length(s1))), levels = c(group2, group1)),
          row.names = heat_cols
        )
        my_palette <- grDevices::colorRampPalette(c("#045a8d", "white", "#d7301f"))(1000)
        heatmap_plot <- pheatmap::pheatmap(
          heatmap_matrix_row_scaled,
          scale = "none",
          show_rownames = isTRUE(heat_show_rownames),
          show_colnames = isTRUE(heat_show_colnames),
          cluster_rows = FALSE,
          cluster_cols = isTRUE(cluster),
          color = my_palette,
          annotation_col = ann_col,
          silent = TRUE
        )

        if (!is.null(save_prefix) && nzchar(save_prefix)) {
          f_pdf <- paste0(save_prefix, "_heatmap.pdf")
          grDevices::pdf(f_pdf, width = 5, height = 6)
          grid::grid.newpage(); grid::grid.draw(heatmap_plot$gtable)
          grDevices::dev.off()
          saved_files <- c(saved_files, f_pdf)
        }
      }
    }

    list(
      results = tt,
      summary = table(tt$direction, useNA = "ifany"),
      combined = combined,
      contrast = paste0(group1, " - ", group2),
      samples_used = c(s1, s2),
      deg = deg,
      volcano_plot = volcano_plot,
      heatmap_plot = heatmap_plot,
      heatmap_matrix = heatmap_matrix,
      heatmap_matrix_row_scaled = heatmap_matrix_row_scaled,
      saved_files = unique(saved_files),
      method = "group_FC_only"
    )
  }

  # -------------------------
  # 3. Build pairwise compare list
  # -------------------------
  if (!is.null(compare_list)) {
    compare_list <- as.data.frame(compare_list, check.names = FALSE)
    cc <- get_compare_cols(compare_list)
    pairwise_df <- compare_list[, c(cc$exp_col, cc$ctrl_col), drop = FALSE]
    colnames(pairwise_df) <- c("exp", "ctrl")
    pairwise_df$exp <- trimws(as.character(pairwise_df$exp))
    pairwise_df$ctrl <- trimws(as.character(pairwise_df$ctrl))
    pairwise_df <- pairwise_df[!is.na(pairwise_df$exp) & !is.na(pairwise_df$ctrl) & pairwise_df$exp != "" & pairwise_df$ctrl != "", ]
    pairwise_df <- unique(pairwise_df)
  } else {
    groups <- unique(as.character(sample_info[[group_col]]))
    groups <- groups[!is.na(groups) & groups != ""]
    if (length(groups) < 2) stop("少于两个分组，无法自动生成两两比较。")
    pair_mat <- combn(groups, 2)
    pairwise_df <- data.frame(exp = pair_mat[2, ], ctrl = pair_mat[1, ], stringsAsFactors = FALSE)
  }

  utils::write.csv(pairwise_df, file.path(base_out_dir, "compare_list_used.csv"), row.names = FALSE)
  utils::write.csv(as.data.frame(table(sample_info[[group_col]])), file.path(base_out_dir, "group_sample_count.csv"), row.names = FALSE)

  # -------------------------
  # 4. Run pairwise comparisons
  # -------------------------
  pairwise_results <- list()

  for (i in seq_len(nrow(pairwise_df))) {
    g1 <- pairwise_df$exp[i]
    g2 <- pairwise_df$ctrl[i]
    comp_name <- paste0(safe_slug(g1), "_vs_", safe_slug(g2))
    comp_dir <- file.path(base_out_dir, comp_name)
    if (isTRUE(create_dir)) dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
    save_prefix <- file.path(comp_dir, comp_name)

    s1 <- get_group_samples(g1)
    s2 <- get_group_samples(g2)
    enough_reps <- length(s1) >= min_reps_for_limma && length(s2) >= min_reps_for_limma

    if (enough_reps || low_rep_strategy == "limma") {
      if (!exists("limma_de", mode = "function")) {
        stop("当前环境找不到 limma_de()，请先 source('limma差异.R')。")
      }

      res <- limma_de(
        expr = expr,
        sample_info = sample_info,
        group_col = group_col,
        group1 = g1,
        group2 = g2,
        p_type = p_type,
        p_thresh = p_thresh,
        logfc_thresh = logfc_thresh,
        adjust_method = adjust_method,
        make_volcano = make_volcano,
        make_heatmap = TRUE,
        label_extreme_n_each_side = label_extreme_n_each_side,
        plot_width = plot_width,
        plot_height = plot_height,
        plot_dpi = plot_dpi,
        heat_show_rownames = heat_show_rownames,
        heat_show_colnames = heat_show_colnames,
        cluster_cols = cluster,
        save_prefix = save_prefix,
        save_formats = save_formats,
        tiff_compression = tiff_compression,
        fileEncoding = "GBK",
        create_dir = create_dir
      )
    } else {
      message("低重复比较，使用 group_FC_only：", g1, " vs ", g2,
              "；", g1, " n=", length(s1), ", ", g2, " n=", length(s2))
      res <- group_fc_only_de(g1, g2, save_prefix = save_prefix)
    }

    pairwise_results[[comp_name]] <- res
  }

  # -------------------------
  # 5. Export merged pairwise summaries
  # -------------------------
  all_results <- dplyr::bind_rows(lapply(pairwise_results, `[[`, "results"))
  utils::write.csv(all_results, file.path(base_out_dir, "all_pairwise_DE_results.csv"), row.names = FALSE, fileEncoding = "GBK")

  logfc_wide <- all_results %>%
    dplyr::select(feature_id, comparison, logFC) %>%
    dplyr::distinct(feature_id, comparison, .keep_all = TRUE) %>%
    tidyr::pivot_wider(names_from = comparison, values_from = logFC)

  utils::write.csv(logfc_wide, file.path(base_out_dir, "all_pairwise_log2FC_wide.csv"), row.names = FALSE, fileEncoding = "GBK")

  deg_summary <- all_results %>%
    dplyr::group_by(comparison, method, direction) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  utils::write.csv(deg_summary, file.path(base_out_dir, "DE_direction_summary.csv"), row.names = FALSE, fileEncoding = "GBK")

  invisible(list(
    pairwise = pairwise_results,
    all_results = all_results,
    logfc_wide = logfc_wide,
    deg_summary = deg_summary,
    compare_list = pairwise_df
  ))
}
