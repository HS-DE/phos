# ============================================================
# Phosphosite-friendly limma differential analysis function
#
# This file replaces the older protein/gene-centric limma_de() implementation.
# It is adapted for phosphoproteomics features, where one gene/protein can have
# multiple phosphosites or modified peptides.
#
# Important changes:
#   1) Row names are treated as feature_id, not necessarily Gene.
#   2) Feature IDs are preserved exactly; no sub("_.*$", "", x) truncation.
#   3) The function accepts sample_info with either rownames or a sample_id column.
#   4) It still supports volcano and heatmap outputs.
# ============================================================

limma_de <- function(
  expr,
  sample_info,
  group_col = "group",
  group1,
  group2,
  sample_id_col = NULL,
  p_type = c("adjp", "pvalue"),
  p_thresh = 0.05,
  logfc_thresh = 1,
  adjust_method = "BH",
  make_volcano = TRUE,
  make_heatmap = TRUE,
  label_extreme_n_each_side = 10,
  plot_width = 10,
  plot_height = 8,
  plot_dpi = 300,
  heat_show_rownames = FALSE,
  heat_show_colnames = FALSE,
  cluster_cols = FALSE,
  save_prefix = NULL,
  save_formats = c("pdf", "png"),
  tiff_compression = "lzw",
  fileEncoding = "GBK",
  create_dir = TRUE
) {
  # -------------------------
  # 0. Package checks
  # -------------------------
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("需要 limma 包：BiocManager::install('limma')")
  }
  if (!requireNamespace("tools", quietly = TRUE)) {
    stop("需要 tools 包（R 自带）。")
  }

  p_type <- match.arg(tolower(p_type[1]), c("adjp", "pvalue"))
  p_col <- if (p_type == "adjp") "adj.P.Val" else "P.Value"

  # -------------------------
  # 1. Input checks
  # -------------------------
  expr <- as.matrix(expr)
  if (is.null(rownames(expr))) {
    stop("expr 需要有行名；在磷酸化流程中行名应为 feature_id。")
  }
  if (is.null(colnames(expr))) {
    stop("expr 需要有列名；列名应为 sample_id。")
  }

  if (!is.numeric(expr)) {
    expr_num <- apply(expr, 2, function(x) suppressWarnings(as.numeric(x)))
    rownames(expr_num) <- rownames(expr)
    colnames(expr_num) <- colnames(expr)
    expr <- as.matrix(expr_num)
  }

  sample_info <- as.data.frame(sample_info, check.names = FALSE, stringsAsFactors = FALSE)

  if (!is.null(sample_id_col)) {
    if (!sample_id_col %in% colnames(sample_info)) {
      stop("sample_id_col 不在 sample_info 中：", sample_id_col)
    }
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
    stop("在 sample_info 中找不到分组列：", group_col)
  }

  all_groups <- unique(as.character(sample_info[[group_col]]))
  if (!all(c(group1, group2) %in% all_groups)) {
    stop("group1/group2 不在 sample_info 的分组列中：", group1, " / ", group2)
  }

  # -------------------------
  # 2. Subset samples
  # -------------------------
  in_contrast <- rownames(sample_info)[as.character(sample_info[[group_col]]) %in% c(group1, group2)]
  samples_used <- intersect(colnames(expr), in_contrast)

  if (length(samples_used) == 0) {
    stop("没有找到与 expr 匹配且属于这两个组的样本。")
  }

  X <- expr[, samples_used, drop = FALSE]
  meta <- sample_info[samples_used, , drop = FALSE]

  if (length(unique(as.character(meta[[group_col]]))) < 2) {
    stop("当前比较只匹配到一个分组，请检查 sample_info 和 expr 列名。")
  }

  # -------------------------
  # 3. limma model
  # -------------------------
  # logFC is group1 - group2.
  grp_raw <- as.character(meta[[group_col]])
  grp_fac <- factor(grp_raw, levels = c(group2, group1))
  lvl_safe <- make.names(levels(grp_fac))
  grp_safe <- factor(make.names(as.character(grp_fac)), levels = lvl_safe)

  design <- stats::model.matrix(~ 0 + grp_safe)
  colnames(design) <- lvl_safe

  contrast_str <- paste0(make.names(group1), " - ", make.names(group2))
  contrast_mat <- limma::makeContrasts(contrasts = contrast_str, levels = design)

  fit <- limma::lmFit(X, design)
  fit2 <- limma::contrasts.fit(fit, contrast_mat)
  fit2 <- limma::eBayes(fit2)

  tt <- limma::topTable(fit2, coef = 1, number = Inf, adjust.method = adjust_method, sort.by = "P")
  if ("B" %in% colnames(tt)) {
    colnames(tt)[colnames(tt) == "B"] <- "B_stat"
  }

  tt$feature_id <- rownames(tt)
  tt <- tt[, c("feature_id", setdiff(colnames(tt), "feature_id")), drop = FALSE]
  rownames(tt) <- tt$feature_id

  if (!p_col %in% colnames(tt)) {
    stop("结果中不存在显著性列：", p_col)
  }

  tt$direction <- ifelse(
    !is.na(tt[[p_col]]) & tt[[p_col]] < p_thresh & tt$logFC >= logfc_thresh,
    "Up",
    ifelse(
      !is.na(tt[[p_col]]) & tt[[p_col]] < p_thresh & tt$logFC <= -logfc_thresh,
      "Down",
      "Not Significant"
    )
  )
  tt$direction <- factor(tt$direction, levels = c("Up", "Not Significant", "Down"))

  tt$comparison <- paste0(group1, "_vs_", group2)
  tt$group1 <- group1
  tt$group2 <- group2
  tt$n_group1 <- sum(grp_raw == group1)
  tt$n_group2 <- sum(grp_raw == group2)
  tt$method <- "limma"

  summary_tbl <- table(tt$direction, useNA = "ifany")

  common_features <- intersect(rownames(tt), rownames(X))
  combined <- cbind(
    tt[common_features, , drop = FALSE],
    X[common_features, , drop = FALSE]
  )

  deg <- tt$feature_id[tt$direction %in% c("Up", "Down")]

  # -------------------------
  # 4. Save tables
  # -------------------------
  saved_files <- character(0)
  do_save <- !is.null(save_prefix) && nzchar(save_prefix)

  if (do_save) {
    if (isTRUE(create_dir)) {
      dir.create(dirname(save_prefix), recursive = TRUE, showWarnings = FALSE)
    }

    f_res <- paste0(save_prefix, "_limma_results.csv")
    f_comb <- paste0(save_prefix, "_combined.csv")
    utils::write.csv(tt, f_res, row.names = FALSE, fileEncoding = fileEncoding)
    utils::write.csv(combined, f_comb, row.names = FALSE, fileEncoding = fileEncoding)
    saved_files <- c(saved_files, f_res, f_comb)
  }

  # -------------------------
  # 5. Volcano plot
  # -------------------------
  volcano_plot <- NULL
  volcano_data <- NULL

  if (isTRUE(make_volcano)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("未安装 ggplot2，跳过火山图。")
    } else {
      volcano_df <- tt
      pvals_safe <- pmax(volcano_df[[p_col]], .Machine$double.xmin)
      volcano_df$neg_log10_p <- -log10(pvals_safe)

      sig_df <- volcano_df[
        !is.na(volcano_df[[p_col]]) &
          volcano_df[[p_col]] < p_thresh &
          abs(volcano_df$logFC) >= logfc_thresh,
        ,
        drop = FALSE
      ]
      sig_df <- sig_df[order(sig_df$logFC, decreasing = FALSE), , drop = FALSE]

      most_negative <- head(sig_df, label_extreme_n_each_side)
      most_positive <- tail(sig_df, label_extreme_n_each_side)
      volcano_df$labelFeature <- NA_character_
      volcano_df$labelFeature[volcano_df$feature_id %in% most_negative$feature_id] <- most_negative$feature_id
      volcano_df$labelFeature[volcano_df$feature_id %in% most_positive$feature_id] <- most_positive$feature_id

      volcano_plot <- ggplot2::ggplot(volcano_df, ggplot2::aes(x = logFC, y = neg_log10_p, colour = direction)) +
        ggplot2::geom_point(alpha = 0.45, size = 2.4) +
        ggplot2::scale_color_manual(
          values = c("Down" = "#67a9cf", "Not Significant" = "#cccccc", "Up" = "#fd8d3c"),
          limits = c("Up", "Not Significant", "Down"),
          drop = TRUE
        ) +
        ggplot2::geom_vline(xintercept = c(-logfc_thresh, logfc_thresh), linetype = 4, color = "black") +
        ggplot2::geom_hline(yintercept = -log10(p_thresh), linetype = 4, color = "black") +
        ggplot2::labs(
          x = paste0("log2FC (", group1, " - ", group2, ")"),
          y = paste0("-log10(", p_col, ")"),
          color = ""
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(aspect.ratio = 1, legend.position = "right")

      if (requireNamespace("ggrepel", quietly = TRUE)) {
        volcano_plot <- volcano_plot +
          ggrepel::geom_text_repel(
            data = subset(volcano_df, !is.na(labelFeature)),
            ggplot2::aes(label = labelFeature),
            size = 3,
            show.legend = FALSE,
            max.overlaps = Inf
          )
      }

      volcano_data <- volcano_df

      if (do_save) {
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
  }

  # -------------------------
  # 6. Heatmap for Up/Down features
  # -------------------------
  heatmap_plot <- NULL
  heatmap_matrix <- NULL
  heatmap_matrix_row_scaled <- NULL

  if (isTRUE(make_heatmap)) {
    if (!requireNamespace("pheatmap", quietly = TRUE) || !requireNamespace("grid", quietly = TRUE)) {
      warning("未安装 pheatmap/grid，跳过热图。")
    } else {
      sig_features <- tt$feature_id[tt$direction %in% c("Up", "Down")]
      sig_features <- intersect(sig_features, rownames(X))

      if (length(sig_features) > 1) {
        samples_group2 <- rownames(meta)[grp_raw == group2]
        samples_group1 <- rownames(meta)[grp_raw == group1]
        heat_cols <- c(samples_group2, samples_group1)

        heatmap_matrix <- X[sig_features, heat_cols, drop = FALSE]
        row_sd <- apply(heatmap_matrix, 1, stats::sd, na.rm = TRUE)
        keep <- is.finite(row_sd) & row_sd > 0
        heatmap_matrix <- heatmap_matrix[keep, , drop = FALSE]

        if (nrow(heatmap_matrix) > 1) {
          heatmap_matrix_row_scaled <- t(scale(t(heatmap_matrix)))
          heatmap_matrix_row_scaled[!is.finite(heatmap_matrix_row_scaled)] <- 0

          ann_col <- data.frame(
            Group = factor(c(rep(group2, length(samples_group2)), rep(group1, length(samples_group1))),
                           levels = c(group2, group1)),
            row.names = heat_cols
          )

          my_palette <- grDevices::colorRampPalette(c("#045a8d", "white", "#d7301f"))(1000)
          heatmap_plot <- pheatmap::pheatmap(
            heatmap_matrix_row_scaled,
            scale = "none",
            show_rownames = isTRUE(heat_show_rownames),
            show_colnames = isTRUE(heat_show_colnames),
            cluster_rows = FALSE,
            cluster_cols = isTRUE(cluster_cols),
            color = my_palette,
            annotation_col = ann_col,
            silent = TRUE
          )

          if (do_save) {
            f_raw <- paste0(save_prefix, "_heatmap_matrix_raw.csv")
            f_scaled <- paste0(save_prefix, "_heatmap_matrix_row_scaled.csv")
            utils::write.csv(data.frame(feature_id = rownames(heatmap_matrix), heatmap_matrix, check.names = FALSE),
                             f_raw, row.names = FALSE, fileEncoding = fileEncoding)
            utils::write.csv(data.frame(feature_id = rownames(heatmap_matrix_row_scaled), heatmap_matrix_row_scaled, check.names = FALSE),
                             f_scaled, row.names = FALSE, fileEncoding = fileEncoding)

            f_pdf <- paste0(save_prefix, "_heatmap.pdf")
            grDevices::pdf(f_pdf, width = 5, height = 6)
            grid::grid.newpage(); grid::grid.draw(heatmap_plot$gtable)
            grDevices::dev.off()

            f_png <- paste0(save_prefix, "_heatmap.png")
            grDevices::png(f_png, width = 5, height = 6, units = "in", res = 600)
            grid::grid.newpage(); grid::grid.draw(heatmap_plot$gtable)
            grDevices::dev.off()

            saved_files <- c(saved_files, f_raw, f_scaled, f_pdf, f_png)
          }
        }
      }
    }
  }

  list(
    results = tt,
    summary = summary_tbl,
    combined = combined,
    design = design,
    contrast = contrast_str,
    samples_used = samples_used,
    deg = deg,
    volcano_plot = volcano_plot,
    volcano_data = volcano_data,
    heatmap_plot = heatmap_plot,
    heatmap_matrix = heatmap_matrix,
    heatmap_matrix_row_scaled = heatmap_matrix_row_scaled,
    saved_files = unique(saved_files),
    method = "limma"
  )
}
