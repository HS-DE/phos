limma_de <- function(
  expr,
  sample_info,
  group_col = "type",
  group1,
  group2,
  p_type = c("adjp", "pvalue"),
  p_thresh = 0.05,
  logfc_thresh = 1,
  adjust_method = "BH",
  # ==== Volcano ====
  make_volcano = TRUE,
  label_extreme_n_each_side = 10,
  plot_width = 10,
  plot_height = 8,
  plot_dpi = 300,
  # ==== Heatmap ====
  make_heatmap = TRUE,
  heat_show_rownames = FALSE,
  heat_show_colnames = FALSE,
  cluster_cols = FALSE,
  # ==== Save ====
  save_prefix = NULL, # 目录 or 前缀
  save_formats = c("pdf", "png"), # volcano 输出格式（pdf/png/tiff）
  tiff_compression = "lzw",
  fileEncoding = "GBK",
  create_dir = TRUE
) {
  # -------------------------
  # helper
  # -------------------------
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  normalize_save_prefix <- function(save_prefix, group1, group2) {
    if (is.null(save_prefix) || !nzchar(save_prefix)) {
      return(NULL)
    }

    # 如果是目录（以/或\结尾，或本身是目录），自动拼一个对比名
    is_dir_style <- grepl("[/\\\\]$", save_prefix) || dir.exists(save_prefix)
    if (is_dir_style) {
      dir_path <- sub("[/\\\\]$", "", save_prefix)
      base <- file.path(dir_path, paste0(group1, "_vs_", group2))
      return(base)
    }
    return(save_prefix)
  }

  safe_dir_create <- function(path) {
    if (isTRUE(create_dir)) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }
  }

  safe_ggsave <- function(
    plot,
    filename,
    width,
    height,
    dpi,
    tiff_compression
  ) {
    fmt <- tolower(tools::file_ext(filename))
    if (fmt %in% c("tif", "tiff")) {
      ggplot2::ggsave(
        filename = filename,
        plot = plot,
        width = width,
        height = height,
        dpi = dpi,
        device = "tiff",
        compression = tiff_compression
      )
    } else if (fmt == "png") {
      ggplot2::ggsave(
        filename = filename,
        plot = plot,
        width = width,
        height = height,
        dpi = dpi,
        device = "png"
      )
    } else {
      ggplot2::ggsave(
        filename = filename,
        plot = plot,
        width = width,
        height = height,
        dpi = dpi
      )
    }
  }

  # -------------------------
  # packages
  # -------------------------
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("需要 limma 包：BiocManager::install('limma')")
  }
  if (!requireNamespace("tools", quietly = TRUE)) {
    stop("需要 tools 包（R 自带），请检查环境。")
  }

  p_type <- match.arg(tolower(p_type), c("adjp", "pvalue"))

  # -------------------------
  # input check & coerce
  # -------------------------
  expr <- as.matrix(expr)
  if (is.null(rownames(expr))) {
    stop("expr 需要有行名（基因/蛋白ID）。")
  }
  if (is.null(colnames(expr))) {
    stop("expr 需要有列名（样本ID）。")
  }

  # 强制 expr 数值化（防字符矩阵）
  if (!is.numeric(expr)) {
    expr_num <- apply(expr, 2, function(x) suppressWarnings(as.numeric(x)))
    rownames(expr_num) <- rownames(expr)
    colnames(expr_num) <- colnames(expr)
    expr <- as.matrix(expr_num)
  }

  # sample_info 兼容 tibble
  sample_info <- as.data.frame(
    sample_info,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # 如果没有 rownames，尝试自动设置
  if (
    is.null(rownames(sample_info)) ||
      any(is.na(rownames(sample_info))) ||
      any(rownames(sample_info) == "")
  ) {
    candidate_cols <- c(
      "sample_name",
      "Sample_name",
      "sample_id",
      "SampleID",
      "样品信息",
      "样品编号"
    )
    hit <- intersect(candidate_cols, colnames(sample_info))
    if (length(hit) > 0) {
      rownames(sample_info) <- as.character(sample_info[[hit[1]]])
    } else {
      stop(
        "sample_info 没有行名；也没找到可用列（sample_name/sample_id/样品信息/样品编号）用来设行名。"
      )
    }
  }
  if (anyDuplicated(rownames(sample_info)) > 0) {
    rownames(sample_info) <- make.unique(rownames(sample_info))
    warning("sample_info 行名有重复，已 make.unique()。请检查样本ID是否唯一。")
  }

  if (!group_col %in% colnames(sample_info)) {
    stop(sprintf("在 sample_info 中找不到分组列 '%s'。", group_col))
  }

  # 检查 group1/group2 是否存在
  grps_all <- unique(sample_info[[group_col]])
  if (!all(c(group1, group2) %in% grps_all)) {
    stop(
      "group1 / group2 不在 sample_info 的分组列中：请检查拼写或分组列内容。"
    )
  }

  # -------------------------
  # subset samples
  # -------------------------
  in_contrast <- rownames(sample_info)[
    sample_info[[group_col]] %in% c(group1, group2)
  ]
  samples_used <- intersect(colnames(expr), in_contrast)
  if (length(samples_used) == 0) {
    stop("没有找到与 expr 匹配且属于这两个组的样本。")
  }

  X <- expr[, samples_used, drop = FALSE]
  meta <- sample_info[samples_used, , drop = FALSE]

  # -------------------------
  # design & fit
  # -------------------------
  grp_raw <- meta[[group_col]]
  grp_fac <- factor(grp_raw, levels = c(group1, group2)) # group1 - group2
  lvl_safe <- make.names(levels(grp_fac))
  grp_safe <- factor(make.names(as.character(grp_fac)), levels = lvl_safe)

  design <- stats::model.matrix(~ 0 + grp_safe)
  colnames(design) <- lvl_safe

  fit <- limma::lmFit(X, design)
  contrast_str <- paste0(lvl_safe[1], " - ", lvl_safe[2])
  contrast.mat <- limma::makeContrasts(
    contrasts = contrast_str,
    levels = colnames(design)
  )
  fit2 <- limma::contrasts.fit(fit, contrast.mat)
  fit2 <- limma::eBayes(fit2)

  # -------------------------
  # results table
  # -------------------------
  tt <- limma::topTable(fit2, coef = 1, n = Inf, adjust.method = adjust_method)
  if ("B" %in% colnames(tt)) {
    colnames(tt)[colnames(tt) == "B"] <- "B_stat"
  }

  tt$Gene <- rownames(tt)

  pcol <- if (p_type == "adjp") "adj.P.Val" else "P.Value"
  if (!pcol %in% colnames(tt)) {
    stop(sprintf("结果中不存在列 %s。", pcol))
  }

  tt$direction <- ifelse(
    is.na(tt[[pcol]]),
    "NA",
    ifelse(
      tt[[pcol]] < p_thresh & tt$logFC >= logfc_thresh,
      "Up",
      ifelse(
        tt[[pcol]] < p_thresh & tt$logFC <= -logfc_thresh,
        "Down",
        "Not Significant"
      )
    )
  )

  # 不强依赖 dplyr，避免 select 丢 rownames
  tt <- tt[, c("Gene", setdiff(colnames(tt), "Gene")), drop = FALSE]
  rownames(tt) <- tt$Gene

  tt$direction <- factor(
    tt$direction,
    levels = c("Up", "Not Significant", "Down")
  )
  summ <- table(tt$direction, useNA = "ifany")

  common_genes <- intersect(rownames(tt), rownames(X))
  combined <- cbind(
    tt[common_genes, , drop = FALSE],
    X[common_genes, , drop = FALSE]
  )

  deg <- tt$Gene[tt$direction %in% c("Up", "Down")]

  # -------------------------
  # save prefix handling
  # -------------------------
  file_base <- normalize_save_prefix(save_prefix, group1, group2)
  do_save <- !is.null(file_base) && nzchar(file_base)
  if (do_save) {
    safe_dir_create(dirname(file_base))
  }

  # -------------------------
  # volcano
  # -------------------------
  volcano_plot <- NULL
  volcano_data <- NULL
  saved_files <- character(0)

  if (isTRUE(make_volcano)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("未安装 ggplot2，跳过火山图。")
    } else {
      if (!requireNamespace("grid", quietly = TRUE)) {
        stop("需要 grid 包（R 自带）。")
      }

      pvals_safe <- pmax(tt[[pcol]], .Machine$double.xmin)
      volcano_df <- tt
      volcano_df$neg_log10_p <- -log10(pvals_safe)
      volcano_df$geneName <- rownames(tt)

      sig_df <- volcano_df[
        !is.na(volcano_df[[pcol]]) &
          volcano_df[[pcol]] < p_thresh &
          abs(volcano_df$logFC) > logfc_thresh,
        ,
        drop = FALSE
      ]
      sig_df <- sig_df[order(sig_df$logFC, decreasing = FALSE), , drop = FALSE]
      most_negative <- head(sig_df, label_extreme_n_each_side)
      most_positive <- tail(sig_df, label_extreme_n_each_side)

      volcano_df$labelGene <- NA_character_
      if (nrow(most_negative) > 0) {
        volcano_df$labelGene[
          volcano_df$geneName %in% most_negative$geneName
        ] <- most_negative$geneName
      }
      if (nrow(most_positive) > 0) {
        volcano_df$labelGene[
          volcano_df$geneName %in% most_positive$geneName
        ] <- most_positive$geneName
      }

      col_map <- c(
        "Down" = "#67a9cf",
        "Not Significant" = "#cccccc",
        "Up" = "#fd8d3c"
      )
      y_hline <- -log10(p_thresh)
      x_vlines <- c(-logfc_thresh, logfc_thresh)

      gp <- ggplot2::ggplot(
        volcano_df,
        ggplot2::aes(x = logFC, y = neg_log10_p, colour = direction)
      ) +
        ggplot2::geom_point(alpha = 0.4, size = 2.8) +
        ggplot2::scale_color_manual(
          values = col_map,
          limits = names(col_map),
          drop = TRUE
        ) +
        ggplot2::geom_vline(
          xintercept = x_vlines,
          linetype = 4,
          color = "black",
          linewidth = 0.8
        ) +
        ggplot2::geom_hline(
          yintercept = y_hline,
          linetype = 4,
          color = "black",
          linewidth = 0.8
        ) +
        ggplot2::labs(
          x = paste0("Log2(Fold change, ", contrast_str, ")"),
          y = paste0("-Log10 (", pcol, ")"),
          color = ""
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          aspect.ratio = 1,
          axis.title.y = ggplot2::element_text(
            size = 16,
            face = "plain",
            color = "black"
          ),
          axis.title.x = ggplot2::element_text(
            size = 16,
            face = "plain",
            color = "black"
          ),
          axis.text.x = ggplot2::element_text(
            size = 14,
            face = "plain",
            color = "black"
          ),
          axis.text.y = ggplot2::element_text(
            size = 14,
            face = "plain",
            color = "black"
          ),
          plot.title = ggplot2::element_text(hjust = 0.5),
          legend.position = "right",
          legend.text = ggplot2::element_text(size = 12),
          legend.key.height = grid::unit(0.4, "cm"),
          legend.key.width = grid::unit(0.3, "cm")
        )

      if (requireNamespace("ggrepel", quietly = TRUE)) {
        gp <- gp +
          ggrepel::geom_text_repel(
            data = subset(volcano_df, !is.na(labelGene)),
            ggplot2::aes(label = labelGene),
            size = 4,
            box.padding = grid::unit(0.2, "lines"),
            point.padding = grid::unit(0.4, "lines"),
            segment.color = "black",
            show.legend = FALSE,
            max.overlaps = Inf
          )
      }

      volcano_plot <- gp
      volcano_data <- volcano_df

      if (do_save) {
        # 保存 volcano
        fmts <- unique(tolower(save_formats))
        for (fmt in fmts) {
          f <- paste0(file_base, "_volcano.", fmt)
          safe_ggsave(
            gp,
            f,
            plot_width,
            plot_height,
            plot_dpi,
            tiff_compression
          )
          saved_files <- c(saved_files, f)
        }
      }
    }
  }

  # -------------------------
  # save tables
  # -------------------------
  if (do_save) {
    f_res <- paste0(file_base, "_limma_results.csv")
    f_comb <- paste0(file_base, "_combined.csv")
    utils::write.csv(tt, f_res, row.names = FALSE, fileEncoding = fileEncoding)
    utils::write.csv(
      combined,
      f_comb,
      row.names = TRUE,
      fileEncoding = fileEncoding
    )
    saved_files <- c(saved_files, f_res, f_comb)
  }

  # -------------------------
  # heatmap
  # -------------------------
  ph <- NULL
  heatmap_files <- character(0)

  # 新增：用于保存/返回热图矩阵
  heatmap_matrix <- NULL
  heatmap_matrix_row_scaled <- NULL
  heatmap_matrix_saved_files <- character(0)

  if (isTRUE(make_heatmap)) {
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
      warning("未安装 pheatmap，跳过热图。install.packages('pheatmap') 可用。")
    } else if (!requireNamespace("grid", quietly = TRUE)) {
      warning("缺少 grid（R 自带一般不会缺），跳过热图。")
    } else {
      my_palette <- grDevices::colorRampPalette(c(
        "#045a8d",
        "white",
        "#d7301f"
      ))(1000)

      filtered_degs <- combined[
        combined$direction %in% c("Down", "Up"),
        ,
        drop = FALSE
      ]
      if (nrow(filtered_degs) == 0) {
        warning("没有找到差异基因用于绘制热图（Up/Down 为空）。跳过热图。")
      } else {
        # 排序
        sorted_results <- filtered_degs
        if ("logFC" %in% colnames(sorted_results)) {
          sorted_results <- sorted_results[
            order(sorted_results$logFC, decreasing = FALSE),
            ,
            drop = FALSE
          ]
        }

        # 这次对比用到的样本列
        expr_cols <- intersect(samples_used, colnames(sorted_results))
        if (length(expr_cols) == 0) {
          warning("combined 里找不到样本表达列，跳过热图。")
        } else {
          # 固定列顺序：group2 -> group1（你之前希望红左蓝右/或至少固定）
          grp_vec <- as.character(sample_info[expr_cols, group_col])
          ord <- order(factor(grp_vec, levels = c(group2, group1)))
          expr_cols <- expr_cols[ord]

          # 表达矩阵
          heat_mat <- data.matrix(sorted_results[, expr_cols, drop = FALSE])
          rownames(heat_mat) <- rownames(sorted_results)
          colnames(heat_mat) <- expr_cols

          # 过滤：至少有一个 finite
          keep_non_na <- rowSums(is.finite(heat_mat)) > 0
          heat_mat <- heat_mat[keep_non_na, , drop = FALSE]

          if (nrow(heat_mat) == 0) {
            warning("热图矩阵全是 NA（过滤后为空）。跳过热图。")
          } else {
            # 过滤：去 0 方差（row-scale 会 NaN）
            row_sd <- apply(heat_mat, 1, sd, na.rm = TRUE)
            keep_var <- is.finite(row_sd) & row_sd > 0
            heat_mat <- heat_mat[keep_var, , drop = FALSE]

            if (nrow(heat_mat) == 0) {
              warning(
                "热图矩阵过滤后全是 0 方差行（row-scale 会产生 NaN）。跳过热图。"
              )
            } else {
              # 安全 row-scale
              heat_scaled <- t(scale(t(heat_mat)))
              heat_scaled[!is.finite(heat_scaled)] <- 0

              # =========================
              # 导出热图矩阵
              # =========================
              heatmap_matrix <- heat_mat
              heatmap_matrix_row_scaled <- heat_scaled

              if (do_save) {
                f_heat_raw <- paste0(file_base, "_heatmap_matrix_raw.csv")
                f_heat_scaled <- paste0(
                  file_base,
                  "_heatmap_matrix_row_scaled.csv"
                )

                utils::write.csv(
                  data.frame(
                    Gene = rownames(heatmap_matrix),
                    heatmap_matrix,
                    check.names = FALSE
                  ),
                  f_heat_raw,
                  row.names = FALSE,
                  fileEncoding = fileEncoding
                )

                utils::write.csv(
                  data.frame(
                    Gene = rownames(heatmap_matrix_row_scaled),
                    heatmap_matrix_row_scaled,
                    check.names = FALSE
                  ),
                  f_heat_scaled,
                  row.names = FALSE,
                  fileEncoding = fileEncoding
                )

                heatmap_matrix_saved_files <- c(
                  heatmap_matrix_saved_files,
                  f_heat_raw,
                  f_heat_scaled
                )
              }

              # 注释列
              ann_df <- data.frame(
                Group = factor(
                  sample_info[expr_cols, group_col],
                  levels = c(group2, group1)
                ),
                row.names = expr_cols
              )
              ann_df$Group <- droplevels(ann_df$Group)
              annotation_colors <- list(
                Group = setNames(c("#c74732", "#045a8d"), c(group2, group1))
              )

              ph <- pheatmap::pheatmap(
                heat_scaled,
                scale = "none",
                show_rownames = isTRUE(heat_show_rownames),
                show_colnames = isTRUE(heat_show_colnames),
                cluster_rows = FALSE,
                cluster_cols = isTRUE(cluster_cols),
                color = my_palette,
                annotation_col = ann_df,
                annotation_colors = annotation_colors,
                annotation_legend = TRUE,
                silent = TRUE
              )

              # 会话里画出来
              grid::grid.newpage()
              grid::grid.draw(ph$gtable)

              # 保存（只有 ph 不为空才保存）
              if (do_save && !is.null(ph)) {
                f_pdf <- paste0(file_base, "_heatmap.pdf")
                grDevices::pdf(f_pdf, width = 5, height = 6)
                grid::grid.newpage()
                grid::grid.draw(ph$gtable)
                grDevices::dev.off()

                f_png <- paste0(file_base, "_heatmap.png")
                grDevices::png(
                  f_png,
                  width = 5,
                  height = 6,
                  units = "in",
                  res = 600
                )
                grid::grid.newpage()
                grid::grid.draw(ph$gtable)
                grDevices::dev.off()

                f_tiff <- paste0(file_base, "_heatmap.tiff")
                grDevices::tiff(
                  f_tiff,
                  width = 5,
                  height = 6,
                  units = "in",
                  res = 600,
                  compression = tiff_compression
                )
                grid::grid.newpage()
                grid::grid.draw(ph$gtable)
                grDevices::dev.off()

                heatmap_files <- c(f_pdf, f_png, f_tiff)
              }
            }
          }
        }
      }
    }
  }

  saved_files <- unique(c(
    saved_files,
    heatmap_files,
    heatmap_matrix_saved_files
  ))

  tt$Gene <- sub("_.*$", "", tt$Gene)
  combined$Gene <- sub("_.*$", "", combined$Gene)
  deg <- unique(sub("_.*$", "", deg))
  # -------------------------
  # return
  # -------------------------
  list(
    results = tt,
    summary = summ,
    combined = combined,
    design = design,
    contrast = contrast_str,
    samples_used = samples_used,
    deg = deg,
    volcano_plot = volcano_plot,
    volcano_data = volcano_data,

    # 新增：热图对象和矩阵
    heatmap_plot = ph,
    heatmap_matrix = heatmap_matrix,
    heatmap_matrix_row_scaled = heatmap_matrix_row_scaled,
    heatmap_matrix_saved_files = heatmap_matrix_saved_files,

    saved_files = saved_files
  )
}

if (F) {
  limma_voom_de <- function(
    expr,
    sample_info,
    group_col = "type",
    group1,
    group2,
    p_type = c("adjp", "pvalue"),
    p_thresh = 0.05,
    logfc_thresh = 1,
    adjust_method = "BH",
    # ==== 火山图相关（不默认保存） ====
    make_volcano = TRUE,
    label_extreme_n_each_side = 10,
    plot_width = 10,
    plot_height = 8,
    plot_dpi = 300,
    # ==== 保存相关：调用时指定前缀；默认不保存 ====
    save_prefix = NULL, # 例如 "../Result/Volcano_PE_vs_HC"
    save_formats = c("pdf", "tiff"), # 可选：c("pdf"), c("tiff"), c("pdf","png") 等
    tiff_compression = "lzw",
    cluster = F,
    create_dir = TRUE
  ) {
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop(
        "需要 limma 包：install.packages('limma') 或 BiocManager::install('limma')"
      )
    }
    p_type <- match.arg(tolower(p_type), c("adjp", "pvalue"))

    # —— 基础检查 ----
    expr <- as.matrix(expr)
    if (is.null(rownames(expr))) {
      stop("expr 需要有行名（基因ID）。")
    }
    if (is.null(colnames(expr))) {
      stop("expr 需要有列名（样本ID）。")
    }
    if (is.null(rownames(sample_info))) {
      stop("sample_info 需要以样本ID为行名（必须能与 expr 的列名匹配）。")
    }
    if (!group_col %in% colnames(sample_info)) {
      stop(sprintf("在 sample_info 中找不到分组列 '%s'。", group_col))
    }

    # —— 选样本 ----
    if (!all(c(group1, group2) %in% unique(sample_info[[group_col]]))) {
      stop("group1 / group2 不在 sample_info 的分组列中。")
    }
    in_contrast <- rownames(sample_info)[
      sample_info[[group_col]] %in% c(group1, group2)
    ]
    samples_used <- intersect(colnames(expr), in_contrast)
    if (length(samples_used) == 0) {
      stop("没有找到与 expr 匹配且属于这两个组的样本。")
    }

    X <- expr[, samples_used, drop = FALSE]
    # X 已经是 expr[, samples_used]
    # 检测是否像原始计数（非常简单的启发式判断）
    is_integer_like <- all(round(X) == X)
    max_val <- max(X, na.rm = TRUE)
    likely_counts <- is_integer_like && max_val > 50

    if (!likely_counts) {
      warning(
        "输入矩阵看起来可能不是原始计数（非整数或值较小）。voom 通常用于原始 counts 数据。\n",
        "仍将尝试使用 voom；如果你的数据已经是 log-CPM/TPM 等，请不要使用 voom，改用原始 lmFit 路径。"
      )
    }

    meta <- sample_info[samples_used, , drop = FALSE]
    # —— 设计矩阵 ----
    # 建立 edgeR 对象并归一化
    y <- edgeR::DGEList(counts = X)
    y <- edgeR::calcNormFactors(y)

    grp_raw <- meta[[group_col]]
    grp_fac <- factor(grp_raw, levels = c(group1, group2))
    lvl_safe <- make.names(levels(grp_fac))
    grp_safe <- factor(make.names(as.character(grp_fac)), levels = lvl_safe)
    design <- stats::model.matrix(~ 0 + grp_safe)
    colnames(design) <- lvl_safe

    # voom 计算 log-CPM 和权重
    v <- limma::voom(y, design = design, plot = FALSE)

    # —— 拟合与对比 ----
    fit <- limma::lmFit(v, design)

    contrast_str <- paste0(lvl_safe[1], " - ", lvl_safe[2])
    contrast.mat <- limma::makeContrasts(
      contrasts = contrast_str,
      levels = colnames(design)
    )
    fit2 <- limma::contrasts.fit(fit, contrast.mat)
    fit2 <- limma::eBayes(fit2)

    # —— 结果 ----
    tt <- limma::topTable(
      fit2,
      coef = 1,
      n = Inf,
      adjust.method = adjust_method
    )
    tt$Gene <- rownames(tt)

    pcol <- if (p_type == "adjp") "adj.P.Val" else "P.Value"
    if (!pcol %in% colnames(tt)) {
      stop(sprintf("结果中不存在列 %s。", pcol))
    }

    tt$direction <- ifelse(
      is.na(tt[[pcol]]),
      "NA",
      ifelse(
        tt[[pcol]] < p_thresh & tt$logFC >= logfc_thresh,
        "Up",
        ifelse(
          tt[[pcol]] < p_thresh & tt$logFC <= -logfc_thresh,
          "Down",
          "Not Significant"
        )
      )
    )
    summ <- table(tt$direction, useNA = "ifany")

    common_genes <- intersect(rownames(tt), rownames(X))
    combined <- cbind(
      tt[common_genes, , drop = FALSE],
      X[common_genes, , drop = FALSE]
    )

    deg <- rownames(tt[tt$direction %in% c("Up", "Down"), ])

    # =========================
    #        火山图
    # =========================
    volcano_plot <- NULL
    volcano_data <- NULL
    volcano_saved_files <- character(0)

    if (make_volcano && requireNamespace("ggplot2", quietly = TRUE)) {
      if (!requireNamespace("ggrepel", quietly = TRUE)) {
        warning(
          "未安装 ggrepel，火山图将不进行标签标注。install.packages('ggrepel')"
        )
      }
      if (!requireNamespace("grid", quietly = TRUE)) {
        stop("需要 grid 包（R 自带），请检查环境。")
      }

      pvals_safe <- pmax(tt[[pcol]], .Machine$double.xmin)
      volcano_df <- tt
      volcano_df$neg_log10_p <- -log10(pvals_safe)
      volcano_df$geneName <- rownames(tt)
      volcano_df$p_col_used <- pcol

      sig_df <- volcano_df[
        volcano_df[[pcol]] < p_thresh & !is.na(volcano_df[[pcol]]),
        ,
        drop = FALSE
      ]
      sig_df <- sig_df[order(sig_df$logFC, decreasing = FALSE), , drop = FALSE]

      most_negative <- head(sig_df, label_extreme_n_each_side)
      most_positive <- tail(sig_df, label_extreme_n_each_side)

      volcano_df$labelGene <- NA_character_
      if (nrow(most_negative) > 0) {
        volcano_df$labelGene[
          volcano_df$geneName %in% most_negative$geneName
        ] <- most_negative$geneName
      }
      if (nrow(most_positive) > 0) {
        volcano_df$labelGene[
          volcano_df$geneName %in% most_positive$geneName
        ] <- most_positive$geneName
      }

      col_map <- c(
        "Down" = "#67a9cf",
        "Not Significant" = "#cccccc",
        "Up" = "#fd8d3c"
      )
      y_hline <- -log10(p_thresh)
      x_vlines <- c(-logfc_thresh, logfc_thresh)

      gp <- ggplot2::ggplot(
        volcano_df,
        ggplot2::aes(x = logFC, y = neg_log10_p, colour = direction)
      ) +
        ggplot2::geom_point(alpha = 0.4, size = 2.8) +
        ggplot2::scale_color_manual(
          values = col_map,
          limits = names(col_map),
          drop = TRUE
        ) +
        ggplot2::geom_vline(
          xintercept = x_vlines,
          linetype = 4,
          color = "black",
          linewidth = 0.8
        ) +
        ggplot2::geom_hline(
          yintercept = y_hline,
          linetype = 4,
          color = "black",
          linewidth = 0.8
        ) +
        ggplot2::labs(
          x = paste0("Log2(Fold change, ", contrast_str, ")"),
          y = paste0("-Log10 (", pcol, ")"),
          color = ""
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          aspect.ratio = 1,
          axis.title.y = ggplot2::element_text(
            size = 16,
            face = "plain",
            color = "black"
          ),
          axis.title.x = ggplot2::element_text(
            size = 16,
            face = "plain",
            color = "black"
          ),
          axis.text.x = ggplot2::element_text(
            size = 14,
            face = "plain",
            color = "black"
          ),
          axis.text.y = ggplot2::element_text(
            size = 14,
            face = "plain",
            color = "black"
          ),
          plot.title = ggplot2::element_text(hjust = 0.5),
          legend.position = "right",
          legend.text = ggplot2::element_text(size = 12),
          legend.key.height = grid::unit(0.4, "cm"),
          legend.key.width = grid::unit(0.3, "cm")
        )

      if (requireNamespace("ggrepel", quietly = TRUE)) {
        gp <- gp +
          ggrepel::geom_text_repel(
            data = subset(volcano_df, !is.na(labelGene)),
            ggplot2::aes(x = logFC, y = neg_log10_p, label = labelGene),
            size = 4,
            box.padding = grid::unit(0.2, "lines"),
            point.padding = grid::unit(0.4, "lines"),
            segment.color = "black",
            show.legend = FALSE,
            max.overlaps = Inf
          )
      }

      volcano_plot <- gp
      volcano_data <- volcano_df

      # —— 仅当调用时显式提供 save_prefix 才保存 ----
      if (!is.null(save_prefix) && nzchar(save_prefix)) {
        if (create_dir) {
          dir.create(
            dirname(save_prefix),
            recursive = TRUE,
            showWarnings = FALSE
          )
        }

        # 若 save_prefix 自带扩展名，则只保存该格式；否则按 save_formats 追加扩展名
        if (requireNamespace("tools", quietly = TRUE)) {
          ext <- tools::file_ext(save_prefix)
        } else {
          ext <- sub(".*\\.(.+)$", "\\1", save_prefix)
        }

        if (nzchar(ext)) {
          # 有扩展名：只保存这一个文件
          fmt <- tolower(ext)
          path <- save_prefix
          if (fmt %in% c("tif", "tiff")) {
            ggplot2::ggsave(
              filename = path,
              plot = gp,
              width = plot_width,
              height = plot_height,
              dpi = plot_dpi,
              device = "tiff",
              compression = tiff_compression
            )
          } else {
            ggplot2::ggsave(
              filename = path,
              plot = gp,
              width = plot_width,
              height = plot_height,
              dpi = plot_dpi
            )
          }
          volcano_saved_files <- c(volcano_saved_files, path)
        } else {
          # 无扩展名：按 save_formats 一次性导出
          for (fmt in tolower(save_formats)) {
            fpath <- paste0(save_prefix, "_volcano", ".", fmt)
            if (fmt %in% c("tif", "tiff")) {
              ggplot2::ggsave(
                filename = fpath,
                plot = gp,
                width = plot_width,
                height = plot_height,
                dpi = plot_dpi,
                device = "tiff",
                compression = tiff_compression
              )
            } else {
              ggplot2::ggsave(
                filename = fpath,
                plot = gp,
                width = plot_width,
                height = plot_height,
                dpi = plot_dpi
              )
            }
            volcano_saved_files <- c(volcano_saved_files, fpath)
          }
        }
      }
    } else if (make_volcano) {
      warning(
        "未安装 ggplot2，无法生成火山图。请先 install.packages('ggplot2')"
      )
    }
    write.csv(combined, paste0(save_prefix, ".csv"), fileEncoding = "GBK")
    # =========================
    #        Heatmap
    # =========================
    heatmap_saved_files <- character(0)

    if (requireNamespace("pheatmap", quietly = TRUE)) {
      # 颜色面板：行 z-score 用调色板
      my_palette <- grDevices::colorRampPalette(c(
        "#045a8d",
        "white",
        "#d7301f"
      ))(1000)

      # 从 combined 中筛选差异基因（Up / Down），并按 logFC 排序（从小到大）
      filtered_degs <- combined[
        combined$direction %in% c("Down", "Up"),
        ,
        drop = FALSE
      ]
      if (nrow(filtered_degs) == 0) {
        warning(
          "没有找到差异基因用于绘制热图（filtered_degs 为空）。跳过热图。"
        )
      } else {
        if ("logFC" %in% colnames(filtered_degs)) {
          sorted_results <- filtered_degs[
            order(filtered_degs$logFC, decreasing = FALSE),
            ,
            drop = FALSE
          ]
        } else {
          sorted_results <- filtered_degs
        }

        # 确定表达量矩阵的列（样本）——优先使用 samples_used（函数中已定义）
        if (exists("samples_used") && length(samples_used) > 0) {
          expr_cols <- intersect(colnames(sorted_results), samples_used)
        } else {
          numeric_mask <- vapply(sorted_results, is.numeric, logical(1))
          expr_cols <- colnames(sorted_results)[numeric_mask]
        }

        if (length(expr_cols) == 0) {
          warning("无法在 combined 中识别出表达列（samples）。热图跳过。")
        } else {
          heat_mat <- as.matrix(sorted_results[, expr_cols, drop = FALSE])
          rownames(heat_mat) <- rownames(sorted_results)

          # ======= 构建 annotation_col（优先使用 sample_info） =======
          annotation_col <- NULL
          groups_present <- NULL

          if (
            !is.null(sample_info) &&
              !is.null(rownames(sample_info)) &&
              all(expr_cols %in% rownames(sample_info)) &&
              group_col %in% colnames(sample_info)
          ) {
            # 直接用 sample_info 的分组列（这是最可靠的方式）
            ann_df <- data.frame(
              Group = factor(sample_info[expr_cols, group_col])
            )
            rownames(ann_df) <- expr_cols
            annotation_col <- ann_df
            groups_present <- levels(ann_df$Group)
          } else {
            # 回退方式：按列名中是否包含 group1/group2 文本来判断（不区分大小写）
            # 先尝试匹配 group1 -> group1, group2 -> group2；若都不匹配则全部标为 "Unknown"
            label_vec <- rep(NA_character_, length(expr_cols))
            g1_match <- grepl(group1, expr_cols, ignore.case = TRUE)
            g2_match <- grepl(group2, expr_cols, ignore.case = TRUE)
            label_vec[g1_match] <- as.character(group1)
            label_vec[g2_match & is.na(label_vec)] <- as.character(group2)
            # 未匹配到的标为 Unknown
            label_vec[is.na(label_vec)] <- "Unknown"
            ann_df <- data.frame(Group = factor(label_vec))
            rownames(ann_df) <- expr_cols
            annotation_col <- ann_df
            groups_present <- levels(ann_df$Group)
          }

          # ======= 构建 annotation_colors（基于 group1/group2 动态） =======
          # 默认颜色池（如果只有两组，按 order: group1, group2），如果有 Unknown 追加灰色
          default_two <- c("#c74732", "#045a8d") # 第一个对应 group1, 第二个对应 group2
          unique_groups <- unique(groups_present)
          # 为 unique_groups 分配颜色：若包含 group1/group2，则按它们顺序赋默认颜色
          ann_colors <- setNames(
            rep("#808080", length(unique_groups)), # 默认灰
            c(unique_groups[2], unique_groups[1])
          )
          if (group1 %in% unique_groups) {
            ann_colors[group2] <- default_two[1]
          }
          if (group2 %in% unique_groups) {
            ann_colors[group1] <- default_two[2]
          }
          # 若还有其他组（不太可能），用 RColorBrewer 补色
          other_groups <- setdiff(unique_groups, c(group1, group2))
          if (length(other_groups) > 0) {
            if (requireNamespace("RColorBrewer", quietly = TRUE)) {
              pal_extra <- RColorBrewer::brewer.pal(
                max(3, length(other_groups) + 1),
                "Set3"
              )
              ann_colors[other_groups] <- pal_extra[seq_len(length(
                other_groups
              ))]
            } else {
              # fallback 随机颜色（但确定性）
              ann_colors[other_groups] <- grDevices::rainbow(length(
                other_groups
              ))
            }
          }

          annotation_colors <- list(Group = ann_colors)
          if (cluster == F) {
            # ======= 绘制热图（在当前会话先绘图） =======
            try(
              {
                pheatmap::pheatmap(
                  heat_mat,
                  scale = "row",
                  show_rownames = FALSE,
                  show_colnames = FALSE,
                  cluster_rows = FALSE,
                  cluster_cols = F,
                  color = my_palette,
                  annotation_col = annotation_col,
                  annotation_colors = annotation_colors,
                  annotation_legend = TRUE
                )
              },
              silent = TRUE
            )

            # ======= 若 save_prefix 存在则保存（按照 save_formats） =======
            fpath_heat_tiff <- paste0(save_prefix, "_heatmap", ".tiff")
            tiff(
              fpath_heat_tiff,
              width = 5,
              height = 6,
              units = "in",
              res = 600,
              compression = tiff_compression
            )
            pheatmap::pheatmap(
              heat_mat,
              scale = "row",
              show_rownames = FALSE,
              show_colnames = FALSE,
              cluster_rows = FALSE,
              cluster_cols = F,
              color = my_palette,
              annotation_col = annotation_col,
              annotation_colors = annotation_colors,
              annotation_legend = TRUE
            )
            dev.off()
            fpath_heat_pdf <- paste0(save_prefix, "_heatmap", ".pdf")
            #pdf(fpath_heat_pdf, width = 5, height = 6)
            pheatmap::pheatmap(
              heat_mat,
              scale = "row",
              show_rownames = FALSE,
              show_colnames = FALSE,
              cluster_rows = FALSE,
              cluster_cols = F,
              color = my_palette,
              filename = fpath_heat_pdf,
              width = 5,
              height = 6,
              annotation_col = annotation_col,
              annotation_colors = annotation_colors,
              annotation_legend = TRUE
            )
            dev.off()

            heatmap_saved_files <- c(fpath_heat_tiff, fpath_heat_pdf)
          } else {
            # ======= 绘制热图（在当前会话先绘图） =======
            try(
              {
                pheatmap::pheatmap(
                  heat_mat,
                  scale = "row",
                  show_rownames = FALSE,
                  show_colnames = FALSE,
                  cluster_rows = FALSE,
                  cluster_cols = TRUE,
                  color = my_palette,
                  annotation_col = annotation_col,
                  annotation_colors = annotation_colors,
                  annotation_legend = TRUE
                )
              },
              silent = TRUE
            )

            # ======= 若 save_prefix 存在则保存（按照 save_formats） =======
            fpath_heat_tiff <- paste0(save_prefix, "_heatmap", ".tiff")
            tiff(
              fpath_heat_tiff,
              width = 5,
              height = 6,
              units = "in",
              res = 600,
              compression = tiff_compression
            )
            pheatmap::pheatmap(
              heat_mat,
              scale = "row",
              show_rownames = FALSE,
              show_colnames = FALSE,
              cluster_rows = FALSE,
              cluster_cols = TRUE,
              color = my_palette,
              annotation_col = annotation_col,
              annotation_colors = annotation_colors,
              annotation_legend = TRUE
            )
            dev.off()
            fpath_heat_pdf <- paste0(save_prefix, "_heatmap", ".pdf")
            #pdf(fpath_heat_pdf, width = 5, height = 6)
            pheatmap::pheatmap(
              heat_mat,
              scale = "row",
              show_rownames = FALSE,
              show_colnames = FALSE,
              cluster_rows = FALSE,
              cluster_cols = TRUE,
              color = my_palette,
              filename = fpath_heat_pdf,
              width = 5,
              height = 6,
              annotation_col = annotation_col,
              annotation_colors = annotation_colors,
              annotation_legend = TRUE
            )
            dev.off()

            heatmap_saved_files <- c(fpath_heat_tiff, fpath_heat_pdf)
          }
        } # end have expr_cols
      } # end have degs
    } else {
      warning(
        "未安装 pheatmap 包，热图功能被跳过。install.packages('pheatmap') 可用。"
      )
    }

    # 将热图保存路径并入火山图保存列表（如果有）
    volcano_saved_files <- unique(c(volcano_saved_files, heatmap_saved_files))
    dev.off()
    list(
      results = tt,
      summary = summ,
      combined = combined,
      design = design,
      contrast = contrast_str,
      samples_used = samples_used,
      deg = deg,
      volcano_plot = volcano_plot,
      volcano_data = volcano_data,
      volcano_saved_files = volcano_saved_files
    )
  }

  limma_de_sele_label <- function(
    expr,
    sample_info,
    group_col = "type",
    group1,
    group2,
    p_type = c("adjp", "pvalue"),
    p_thresh = 0.05,
    logfc_thresh = 1,
    adjust_method = "BH",
    # ==== 火山图相关（不默认保存） ====
    make_volcano = TRUE,
    label_extreme_n_each_side = 10,
    plot_width = 10,
    plot_height = 8,
    plot_dpi = 300,
    # ==== 保存相关：调用时指定前缀；默认不保存 ====
    save_prefix = NULL, # 例如 "../Result/Volcano_PE_vs_HC"
    save_formats = c("pdf", "tiff"), # 可选：c("pdf"), c("tiff"), c("pdf","png") 等
    tiff_compression = "lzw",
    cluster = F,
    create_dir = TRUE
  ) {
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop(
        "需要 limma 包：install.packages('limma') 或 BiocManager::install('limma')"
      )
    }
    p_type <- match.arg(tolower(p_type), c("adjp", "pvalue"))

    # —— 基础检查 ----
    expr <- as.matrix(expr)
    if (is.null(rownames(expr))) {
      stop("expr 需要有行名（基因ID）。")
    }
    if (is.null(colnames(expr))) {
      stop("expr 需要有列名（样本ID）。")
    }
    if (is.null(rownames(sample_info))) {
      stop("sample_info 需要以样本ID为行名（必须能与 expr 的列名匹配）。")
    }
    if (!group_col %in% colnames(sample_info)) {
      stop(sprintf("在 sample_info 中找不到分组列 '%s'。", group_col))
    }

    # —— 选样本 ----
    if (!all(c(group1, group2) %in% unique(sample_info[[group_col]]))) {
      stop("group1 / group2 不在 sample_info 的分组列中。")
    }
    in_contrast <- rownames(sample_info)[
      sample_info[[group_col]] %in% c(group1, group2)
    ]
    samples_used <- intersect(colnames(expr), in_contrast)
    if (length(samples_used) == 0) {
      stop("没有找到与 expr 匹配且属于这两个组的样本。")
    }
    X <- expr[, samples_used, drop = FALSE]

    meta <- sample_info[samples_used, , drop = FALSE]
    # —— 设计矩阵 ----
    grp_raw <- meta[[group_col]]
    grp_fac <- factor(grp_raw, levels = c(group1, group2))
    lvl_safe <- make.names(levels(grp_fac))
    grp_safe <- factor(make.names(as.character(grp_fac)), levels = lvl_safe)
    design <- stats::model.matrix(~ 0 + grp_safe)
    colnames(design) <- lvl_safe

    # —— 拟合与对比 ----
    fit <- limma::lmFit(X, design)
    contrast_str <- paste0(lvl_safe[1], " - ", lvl_safe[2])
    contrast.mat <- limma::makeContrasts(
      contrasts = contrast_str,
      levels = colnames(design)
    )
    fit2 <- limma::contrasts.fit(fit, contrast.mat)
    fit2 <- limma::eBayes(fit2)

    # —— 结果 ----
    tt <- limma::topTable(
      fit2,
      coef = 1,
      n = Inf,
      adjust.method = adjust_method
    )
    tt$Gene <- rownames(tt)

    pcol <- if (p_type == "adjp") "adj.P.Val" else "P.Value"
    if (!pcol %in% colnames(tt)) {
      stop(sprintf("结果中不存在列 %s。", pcol))
    }

    tt$direction <- ifelse(
      is.na(tt[[pcol]]),
      "NA",
      ifelse(
        tt[[pcol]] < p_thresh & tt$logFC >= logfc_thresh,
        "Up",
        ifelse(
          tt[[pcol]] < p_thresh & tt$logFC <= -logfc_thresh,
          "Down",
          "Not Significant"
        )
      )
    )
    summ <- table(tt$direction, useNA = "ifany")

    common_genes <- intersect(rownames(tt), rownames(X))
    combined <- cbind(
      tt[common_genes, , drop = FALSE],
      X[common_genes, , drop = FALSE]
    )

    deg <- rownames(tt[tt$direction %in% c("Up", "Down"), ])

    # =========================
    #        火山图
    # =========================
    volcano_plot <- NULL
    volcano_data <- NULL
    volcano_saved_files <- character(0)

    if (make_volcano && requireNamespace("ggplot2", quietly = TRUE)) {
      if (!requireNamespace("ggrepel", quietly = TRUE)) {
        warning(
          "未安装 ggrepel，火山图将不进行标签标注。install.packages('ggrepel')"
        )
      }
      if (!requireNamespace("grid", quietly = TRUE)) {
        stop("需要 grid 包（R 自带），请检查环境。")
      }

      pvals_safe <- pmax(tt[[pcol]], .Machine$double.xmin)
      volcano_df <- tt
      volcano_df$neg_log10_p <- -log10(pvals_safe)
      volcano_df$geneName <- rownames(tt)
      volcano_df$p_col_used <- pcol

      sig_df <- volcano_df[
        volcano_df[[pcol]] < p_thresh & !is.na(volcano_df[[pcol]]),
        ,
        drop = FALSE
      ]
      sig_df <- sig_df[order(sig_df$logFC, decreasing = FALSE), , drop = FALSE]

      most_negative <- head(sig_df, label_extreme_n_each_side)
      most_positive <- tail(sig_df, label_extreme_n_each_side)

      volcano_df$labelGene <- NA_character_
      volcano_df$labelGene <- ifelse(
        grepl("RAB7A", volcano_df$geneName),
        "RAB7A",
        ifelse(grepl("MFN2", volcano_df$geneName), "MFN2", "")
      )
      col_map <- c(
        "Down" = "#67a9cf",
        "Not Significant" = "#cccccc",
        "Up" = "#fd8d3c"
      )
      y_hline <- -log10(p_thresh)
      x_vlines <- c(-logfc_thresh, logfc_thresh)

      gp <- ggplot2::ggplot(
        volcano_df,
        ggplot2::aes(x = logFC, y = neg_log10_p, colour = direction)
      ) +
        ggplot2::geom_point(alpha = 0.4, size = 2.8) +
        ggplot2::scale_color_manual(
          values = col_map,
          limits = names(col_map),
          drop = TRUE
        ) +
        ggplot2::geom_vline(
          xintercept = x_vlines,
          linetype = 4,
          color = "black",
          linewidth = 0.8
        ) +
        ggplot2::geom_hline(
          yintercept = y_hline,
          linetype = 4,
          color = "black",
          linewidth = 0.8
        ) +
        ggplot2::labs(
          x = paste0("Log2(Fold change, ", contrast_str, ")"),
          y = paste0("-Log10 (", pcol, ")"),
          color = ""
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          aspect.ratio = 1,
          axis.title.y = ggplot2::element_text(
            size = 16,
            face = "plain",
            color = "black"
          ),
          axis.title.x = ggplot2::element_text(
            size = 16,
            face = "plain",
            color = "black"
          ),
          axis.text.x = ggplot2::element_text(
            size = 14,
            face = "plain",
            color = "black"
          ),
          axis.text.y = ggplot2::element_text(
            size = 14,
            face = "plain",
            color = "black"
          ),
          plot.title = ggplot2::element_text(hjust = 0.5),
          legend.position = "right",
          legend.text = ggplot2::element_text(size = 12),
          legend.key.height = grid::unit(0.4, "cm"),
          legend.key.width = grid::unit(0.3, "cm")
        )

      if (requireNamespace("ggrepel", quietly = TRUE)) {
        gp <- gp +
          ggrepel::geom_text_repel(
            data = subset(volcano_df, !is.na(labelGene)),
            ggplot2::aes(x = logFC, y = neg_log10_p, label = labelGene),
            size = 5,
            box.padding = grid::unit(0.2, "lines"),
            point.padding = grid::unit(0.4, "lines"),
            segment.color = "black",
            show.legend = FALSE,
            color = "black",
            max.overlaps = Inf
          )
      }

      volcano_plot <- gp
      volcano_data <- volcano_df

      # —— 仅当调用时显式提供 save_prefix 才保存 ----
      if (!is.null(save_prefix) && nzchar(save_prefix)) {
        if (create_dir) {
          dir.create(
            dirname(save_prefix),
            recursive = TRUE,
            showWarnings = FALSE
          )
        }

        # 若 save_prefix 自带扩展名，则只保存该格式；否则按 save_formats 追加扩展名
        if (requireNamespace("tools", quietly = TRUE)) {
          ext <- tools::file_ext(save_prefix)
        } else {
          ext <- sub(".*\\.(.+)$", "\\1", save_prefix)
        }

        if (nzchar(ext)) {
          # 有扩展名：只保存这一个文件
          fmt <- tolower(ext)
          path <- save_prefix
          if (fmt %in% c("tif", "tiff")) {
            ggplot2::ggsave(
              filename = path,
              plot = gp,
              width = plot_width,
              height = plot_height,
              dpi = plot_dpi,
              device = "tiff",
              compression = tiff_compression
            )
          } else {
            ggplot2::ggsave(
              filename = path,
              plot = gp,
              width = plot_width,
              height = plot_height,
              dpi = plot_dpi
            )
          }
          volcano_saved_files <- c(volcano_saved_files, path)
        } else {
          # 无扩展名：按 save_formats 一次性导出
          for (fmt in tolower(save_formats)) {
            fpath <- paste0(save_prefix, "_volcano", ".", fmt)
            if (fmt %in% c("tif", "tiff")) {
              ggplot2::ggsave(
                filename = fpath,
                plot = gp,
                width = plot_width,
                height = plot_height,
                dpi = plot_dpi,
                device = "tiff",
                compression = tiff_compression
              )
            } else {
              ggplot2::ggsave(
                filename = fpath,
                plot = gp,
                width = plot_width,
                height = plot_height,
                dpi = plot_dpi
              )
            }
            volcano_saved_files <- c(volcano_saved_files, fpath)
          }
        }
      }
    } else if (make_volcano) {
      warning(
        "未安装 ggplot2，无法生成火山图。请先 install.packages('ggplot2')"
      )
    }
    write.csv(combined, paste0(save_prefix, ".csv"))
    # =========================
    #        Heatmap
    # =========================
    heatmap_saved_files <- character(0)

    if (requireNamespace("pheatmap", quietly = TRUE)) {
      # 颜色面板：行 z-score 用调色板
      my_palette <- grDevices::colorRampPalette(c(
        "#045a8d",
        "white",
        "#d7301f"
      ))(1000)

      # 从 combined 中筛选差异基因（Up / Down），并按 logFC 排序（从小到大）
      filtered_degs <- combined[
        combined$direction %in% c("Down", "Up"),
        ,
        drop = FALSE
      ]
      if (nrow(filtered_degs) == 0) {
        warning(
          "没有找到差异基因用于绘制热图（filtered_degs 为空）。跳过热图。"
        )
      } else {
        if ("logFC" %in% colnames(filtered_degs)) {
          sorted_results <- filtered_degs[
            order(filtered_degs$logFC, decreasing = FALSE),
            ,
            drop = FALSE
          ]
        } else {
          sorted_results <- filtered_degs
        }

        # 确定表达量矩阵的列（样本）——优先使用 samples_used（函数中已定义）
        if (exists("samples_used") && length(samples_used) > 0) {
          expr_cols <- intersect(colnames(sorted_results), samples_used)
        } else {
          numeric_mask <- vapply(sorted_results, is.numeric, logical(1))
          expr_cols <- colnames(sorted_results)[numeric_mask]
        }

        if (length(expr_cols) == 0) {
          warning("无法在 combined 中识别出表达列（samples）。热图跳过。")
        } else {
          heat_mat <- as.matrix(sorted_results[, expr_cols, drop = FALSE])
          rownames(heat_mat) <- rownames(sorted_results)

          # ======= 构建 annotation_col（优先使用 sample_info） =======
          annotation_col <- NULL
          groups_present <- NULL

          if (
            !is.null(sample_info) &&
              !is.null(rownames(sample_info)) &&
              all(expr_cols %in% rownames(sample_info)) &&
              group_col %in% colnames(sample_info)
          ) {
            # 直接用 sample_info 的分组列（这是最可靠的方式）
            ann_df <- data.frame(
              Group = factor(sample_info[expr_cols, group_col])
            )
            rownames(ann_df) <- expr_cols
            annotation_col <- ann_df
            groups_present <- levels(ann_df$Group)
          } else {
            # 回退方式：按列名中是否包含 group1/group2 文本来判断（不区分大小写）
            # 先尝试匹配 group1 -> group1, group2 -> group2；若都不匹配则全部标为 "Unknown"
            label_vec <- rep(NA_character_, length(expr_cols))
            g1_match <- grepl(group1, expr_cols, ignore.case = TRUE)
            g2_match <- grepl(group2, expr_cols, ignore.case = TRUE)
            label_vec[g1_match] <- as.character(group1)
            label_vec[g2_match & is.na(label_vec)] <- as.character(group2)
            # 未匹配到的标为 Unknown
            label_vec[is.na(label_vec)] <- "Unknown"
            ann_df <- data.frame(Group = factor(label_vec))
            rownames(ann_df) <- expr_cols
            annotation_col <- ann_df
            groups_present <- levels(ann_df$Group)
          }

          # ======= 构建 annotation_colors（基于 group1/group2 动态） =======
          # 默认颜色池（如果只有两组，按 order: group1, group2），如果有 Unknown 追加灰色
          default_two <- c("#c74732", "#045a8d") # 第一个对应 group1, 第二个对应 group2
          unique_groups <- unique(groups_present)
          # 为 unique_groups 分配颜色：若包含 group1/group2，则按它们顺序赋默认颜色
          ann_colors <- setNames(
            rep("#808080", length(unique_groups)), # 默认灰
            c(unique_groups[2], unique_groups[1])
          )
          if (group1 %in% unique_groups) {
            ann_colors[group2] <- default_two[1]
          }
          if (group2 %in% unique_groups) {
            ann_colors[group1] <- default_two[2]
          }
          # 若还有其他组（不太可能），用 RColorBrewer 补色
          other_groups <- setdiff(unique_groups, c(group1, group2))
          if (length(other_groups) > 0) {
            if (requireNamespace("RColorBrewer", quietly = TRUE)) {
              pal_extra <- RColorBrewer::brewer.pal(
                max(3, length(other_groups) + 1),
                "Set3"
              )
              ann_colors[other_groups] <- pal_extra[seq_len(length(
                other_groups
              ))]
            } else {
              # fallback 随机颜色（但确定性）
              ann_colors[other_groups] <- grDevices::rainbow(length(
                other_groups
              ))
            }
          }

          annotation_colors <- list(Group = ann_colors)

          if (cluster == F) {
            # ======= 绘制热图（在当前会话先绘图） =======
            try(
              {
                pheatmap::pheatmap(
                  heat_mat,
                  scale = "row",
                  show_rownames = FALSE,
                  show_colnames = FALSE,
                  cluster_rows = FALSE,
                  cluster_cols = F,
                  color = my_palette,
                  annotation_col = annotation_col,
                  annotation_colors = annotation_colors,
                  annotation_legend = TRUE
                )
              },
              silent = TRUE
            )

            # ======= 若 save_prefix 存在则保存（按照 save_formats） =======
            fpath_heat_tiff <- paste0(save_prefix, "_heatmap", ".tiff")
            tiff(
              fpath_heat_tiff,
              width = 5,
              height = 6,
              units = "in",
              res = 600,
              compression = tiff_compression
            )
            pheatmap::pheatmap(
              heat_mat,
              scale = "row",
              show_rownames = FALSE,
              show_colnames = FALSE,
              cluster_rows = FALSE,
              cluster_cols = F,
              color = my_palette,
              annotation_col = annotation_col,
              annotation_colors = annotation_colors,
              annotation_legend = TRUE
            )
            dev.off()
            fpath_heat_pdf <- paste0(save_prefix, "_heatmap", ".pdf")
            #pdf(fpath_heat_pdf, width = 5, height = 6)
            pheatmap::pheatmap(
              heat_mat,
              scale = "row",
              show_rownames = FALSE,
              show_colnames = FALSE,
              cluster_rows = FALSE,
              cluster_cols = F,
              color = my_palette,
              filename = fpath_heat_pdf,
              width = 5,
              height = 6,
              annotation_col = annotation_col,
              annotation_colors = annotation_colors,
              annotation_legend = TRUE
            )
            dev.off()

            heatmap_saved_files <- c(fpath_heat_tiff, fpath_heat_pdf)
          } else {
            # ======= 绘制热图（在当前会话先绘图） =======
            try(
              {
                pheatmap::pheatmap(
                  heat_mat,
                  scale = "row",
                  show_rownames = FALSE,
                  show_colnames = FALSE,
                  cluster_rows = FALSE,
                  cluster_cols = TRUE,
                  color = my_palette,
                  annotation_col = annotation_col,
                  annotation_colors = annotation_colors,
                  annotation_legend = TRUE
                )
              },
              silent = TRUE
            )

            # ======= 若 save_prefix 存在则保存（按照 save_formats） =======
            fpath_heat_tiff <- paste0(save_prefix, "_heatmap", ".tiff")
            tiff(
              fpath_heat_tiff,
              width = 5,
              height = 6,
              units = "in",
              res = 600,
              compression = tiff_compression
            )
            pheatmap::pheatmap(
              heat_mat,
              scale = "row",
              show_rownames = FALSE,
              show_colnames = FALSE,
              cluster_rows = FALSE,
              cluster_cols = TRUE,
              color = my_palette,
              annotation_col = annotation_col,
              annotation_colors = annotation_colors,
              annotation_legend = TRUE
            )
            dev.off()
            fpath_heat_pdf <- paste0(save_prefix, "_heatmap", ".pdf")
            #pdf(fpath_heat_pdf, width = 5, height = 6)
            pheatmap::pheatmap(
              heat_mat,
              scale = "row",
              show_rownames = FALSE,
              show_colnames = FALSE,
              cluster_rows = FALSE,
              cluster_cols = TRUE,
              color = my_palette,
              filename = fpath_heat_pdf,
              width = 5,
              height = 6,
              annotation_col = annotation_col,
              annotation_colors = annotation_colors,
              annotation_legend = TRUE
            )
            dev.off()

            heatmap_saved_files <- c(fpath_heat_tiff, fpath_heat_pdf)
          }
        } # end have expr_cols
      } # end have degs
    } else {
      warning(
        "未安装 pheatmap 包，热图功能被跳过。install.packages('pheatmap') 可用。"
      )
    }

    # 将热图保存路径并入火山图保存列表（如果有）
    volcano_saved_files <- unique(c(volcano_saved_files, heatmap_saved_files))
    dev.off()
    list(
      results = tt,
      summary = summ,
      combined = combined,
      design = design,
      contrast = contrast_str,
      samples_used = samples_used,
      deg = deg,
      volcano_plot = volcano_plot,
      volcano_data = volcano_data,
      volcano_saved_files = volcano_saved_files
    )
  }
}


# 使用方法
if (F) {
  out <- limma_de_pair(
    expr = log_data,
    sample_info = data1_group,
    group_col = "type",
    group1 = "PE",
    group2 = "HC",
    p_type = "adjp",
    p_thresh = 0.05,
    logfc_thresh = 1,
    # 火山图可选参数
    make_volcano = TRUE,
    label_extreme_n_each_side = 10,
    save_prefix = "../Result/bulkRNA差异分析/Volcano_PE_vs_HC", # 只在这里指定
    save_formats = c("pdf", "tiff") # 想要的格式
  )

  # 查看火山图
  out$volcano_plot

  # 已保存的文件路径
  out$volcano_saved_files

  head(out$results)
  out$summary
  head(out$combined)
}
