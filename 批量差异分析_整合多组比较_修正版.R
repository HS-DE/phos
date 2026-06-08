run_batch_de_mixed <- function(
  expr,
  sample_info,
  compare_list = NULL,
  multi_compare = NULL,
  group_col = "group",

  # compare_list 列名（不填就自动识别中文/英文）
  exp_col = NULL,
  ctrl_col = NULL,

  # multi_compare：每行一个多组整体比较，至少 3 组
  # 例如：
  #   data.frame(group1 = "sg2", group2 = "sg1", group3 = "WT")
  multi_compare_group_cols = NULL,
  multi_group_min_n = 3,
  min_reps_for_multi_limma = 2,
  multi_top_n = 50,

  # 关键：limma 最小重复数阈值（默认 3；想允许 2v2 就设 2）
  min_reps_for_limma = 3,

  # 低重复组间时（< min_reps_for_limma）采用的策略：默认 group_FC_only
  low_rep_strategy = c("group_fc_only", "limma"),

  # 输出
  base_out_dir = "../Results/差异分析",
  create_dir = TRUE,

  # limma 参数（组间）
  p_type = c("adjp", "pvalue"),
  p_thresh = 0.05,
  logfc_thresh = 1,
  adjust_method = "BH",

  # 图/保存
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

  # FC-only（单样本 & 低重复组间）参数
  assume_log2 = NULL,
  add1 = 1,

  # sample_info 没有 rownames 时，指定哪一列作为样本名（必须能匹配 expr 列名）
  sample_id_col = NULL,

  # 汇总导出
  export_excel = TRUE,
  export_deg_summary = TRUE
) {
  if (exists("limma_de", mode = "function")) {
    rm(limma_de)
  }
  source("C:/Work/SH/code/source/limma差异.R")
  # =========================
  # 0) 基础检查
  # =========================
  if (!is.matrix(expr)) {
    expr <- as.matrix(expr)
  }
  if (is.null(rownames(expr))) {
    stop("expr 需要行名（feature/gene id）。")
  }
  if (is.null(colnames(expr))) {
    stop("expr 需要列名（样本名）。")
  }

  if (!exists("limma_de", mode = "function")) {
    stop("当前环境找不到 limma_de()。请先 source 你的 limma差异.R。")
  }

  if (!is.null(sample_id_col)) {
    if (!sample_id_col %in% colnames(sample_info)) {
      stop("sample_id_col 不在 sample_info 列名中：", sample_id_col)
    }
    rownames(sample_info) <- as.character(sample_info[[sample_id_col]])
  }

  if (!group_col %in% colnames(sample_info)) {
    stop(
      "sample_info 中找不到 group_col='",
      group_col,
      "'。当前列：",
      paste(colnames(sample_info), collapse = ", ")
    )
  }

  low_rep_strategy <- match.arg(low_rep_strategy)

  # p_type 兼容
  p_type <- tolower(p_type[1])
  if (p_type %in% c("pval", "p")) {
    p_type <- "pvalue"
  }
  if (!p_type %in% c("adjp", "pvalue")) {
    stop("p_type 只能是 'adjp' 或 'pvalue'。")
  }

  # =========================
  # 1) 小工具
  # =========================
  safe_slug <- function(x) {
    x <- as.character(x)
    x <- gsub("[[:space:]]+", "_", x)
    x <- gsub("[/\\\\:;\\*\\?\"<>\\|]", "_", x)
    x <- gsub("[\\+]", "plus", x)
    x <- gsub("[^[:alnum:]_\\-\\.]", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
  }

  # compare_list 列名自动识别
  get_compare_cols <- function(df) {
    if (!is.null(exp_col) && !is.null(ctrl_col)) {
      if (!exp_col %in% names(df) || !ctrl_col %in% names(df)) {
        stop("手动指定的 exp_col/ctrl_col 不在 compare_list 列名中。")
      }
      return(list(exp_col = exp_col, ctl_col = ctrl_col))
    }
    cn_exp <- intersect(names(df), c("实验", "exp", "group1", "treat"))
    cn_ctl <- intersect(names(df), c("对照", "ctrl", "group2", "control"))
    if (length(cn_exp) == 0 || length(cn_ctl) == 0) {
      stop(
        "compare_list 必须包含两列：实验/对照（或 exp/ctrl）。当前列名：",
        paste(names(df), collapse = ", ")
      )
    }
    list(exp_col = cn_exp[1], ctl_col = cn_ctl[1])
  }

  # 取某组在 expr 里实际能用的样本列（交集）
  get_group_samples <- function(g) {
    s0 <- rownames(sample_info)[
      as.character(sample_info[[group_col]]) == as.character(g)
    ]
    unique(intersect(colnames(expr), s0))
  }

  # 启发式判断是否已 log2
  guess_is_log2 <- function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) {
      return(TRUE)
    }
    vmax <- max(v)
    vmed <- stats::median(v)
    !(vmax > 1000 || vmed > 100)
  }

  make_unique_id <- function(candidate, existing_names) {
    candidate <- as.character(candidate)
    if (!candidate %in% existing_names) {
      return(candidate)
    }
    idx <- 2L
    out <- paste0(candidate, "__", idx)
    while (out %in% existing_names) {
      idx <- idx + 1L
      out <- paste0(candidate, "__", idx)
    }
    out
  }

  # 解析 multi_compare：每一行取所有非空组名，形成一个多组整体比较
  normalize_multi_compare <- function(df) {
    if (is.null(df)) {
      return(list())
    }
    if (!is.data.frame(df)) {
      stop("multi_compare 必须是 data.frame。")
    }
    if (nrow(df) == 0) {
      return(list())
    }

    cols_use <- if (!is.null(multi_compare_group_cols)) {
      intersect(as.character(multi_compare_group_cols), colnames(df))
    } else {
      colnames(df)
    }
    if (length(cols_use) == 0) {
      return(list())
    }

    out <- vector("list", nrow(df))
    for (i in seq_len(nrow(df))) {
      vals <- unlist(df[i, cols_use, drop = FALSE], use.names = FALSE)
      vals <- trimws(as.character(vals))
      vals <- vals[!is.na(vals) & nzchar(vals)]
      vals <- vals[!duplicated(vals)]
      out[[i]] <- vals
    }
    out
  }

  plot_heatmap_multi_group_overall <- function(
    res_all,
    expr_impute,
    sample_info,
    groups,
    group_col = "group",
    pcol = "adj.P.Val",
    p_thresh = 0.05,
    top_n = 50,
    save_prefix = NULL,
    create_dir = TRUE,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    heat_show_rownames = FALSE,
    heat_show_colnames = TRUE,
    tiff_compression = "lzw",
    fileEncoding = "GBK"
  ) {
    # =========================
    # 0) 初始化：一定要放在前面
    # =========================
    saved_files <- character(0)
    heatmap_matrix_saved_files <- character(0)
    heat_mat <- NULL
    heat_mat_scaled <- NULL
    do_save <- !is.null(save_prefix) && nzchar(save_prefix)

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
      warning("未安装 pheatmap，跳过多组整体比较热图。")
      return(invisible(NULL))
    }
    if (!requireNamespace("grid", quietly = TRUE)) {
      warning("缺少 grid，跳过多组整体比较热图。")
      return(invisible(NULL))
    }

    expr2 <- as.matrix(expr_impute)
    if (is.null(rownames(expr2)) || is.null(colnames(expr2))) {
      warning("expr_impute 缺少行名或列名，跳过多组整体比较热图。")
      return(invisible(NULL))
    }

    sample_info <- as.data.frame(sample_info, check.names = FALSE, stringsAsFactors = FALSE)
    if (is.null(rownames(sample_info)) || !group_col %in% colnames(sample_info)) {
      warning("sample_info 缺少行名或 group_col，跳过多组整体比较热图。")
      return(invisible(NULL))
    }
    if (!"feature" %in% colnames(res_all)) {
      warning("res_all 缺少 feature 列，跳过多组整体比较热图。")
      return(invisible(NULL))
    }
    if (!pcol %in% colnames(res_all)) {
      warning("res_all 缺少显著性列：", pcol, "，跳过多组整体比较热图。")
      return(invisible(NULL))
    }

    # =========================
    # 1) 筛选显著特征
    # =========================
    res2 <- res_all
    rownames(res2) <- as.character(res2$feature)

    sig_tbl <- res2[
      !is.na(res2[[pcol]]) & res2[[pcol]] < p_thresh,
      ,
      drop = FALSE
    ]

    if (nrow(sig_tbl) == 0) {
      warning("多组整体比较没有显著基因/蛋白，跳过热图。")
      return(invisible(NULL))
    }

    # top_n：先限制展示的显著特征数量
    # 优先按 F 从大到小；如果没有 F，则按 p 值从小到大。
    if (!is.null(top_n) && is.finite(top_n) && nrow(sig_tbl) > top_n) {
      if ("F" %in% colnames(sig_tbl)) {
        sig_tbl <- sig_tbl[order(sig_tbl$F, decreasing = TRUE), , drop = FALSE]
      } else {
        sig_tbl <- sig_tbl[order(sig_tbl[[pcol]], decreasing = FALSE), , drop = FALSE]
      }
      sig_tbl <- sig_tbl[seq_len(top_n), , drop = FALSE]
    }

    # =========================
    # 2) 确定样本列顺序：按 groups 的顺序排列
    # =========================
    samples_used <- rownames(sample_info)[
      as.character(sample_info[[group_col]]) %in% groups
    ]
    samples_used <- intersect(colnames(expr2), samples_used)

    if (length(samples_used) < length(groups)) {
      warning("多组整体比较可匹配样本过少，跳过热图。")
      return(invisible(NULL))
    }

    grp_vec <- as.character(sample_info[samples_used, group_col])
    grp_fac <- factor(grp_vec, levels = groups)
    ord <- order(grp_fac)
    samples_used <- samples_used[ord]
    grp_fac <- grp_fac[ord]

    # =========================
    # 3) 构建热图矩阵并确定行顺序
    # =========================
    genes_use <- intersect(rownames(sig_tbl), rownames(expr2))
    if (length(genes_use) == 0) {
      warning("显著基因/蛋白在表达矩阵中找不到，跳过热图。")
      return(invisible(NULL))
    }

    heat_mat <- expr2[genes_use, samples_used, drop = FALSE]

    # 行排序：
    # 多组整体比较没有单一 logFC，因此默认按 mean_range 从大到小；
    # 如果没有 mean_range，则按 F 从大到小。
    if ("mean_range" %in% colnames(sig_tbl)) {
      heat_mat <- heat_mat[
        order(sig_tbl[rownames(heat_mat), "mean_range"], decreasing = TRUE),
        ,
        drop = FALSE
      ]
    } else if ("F" %in% colnames(sig_tbl)) {
      heat_mat <- heat_mat[
        order(sig_tbl[rownames(heat_mat), "F"], decreasing = TRUE),
        ,
        drop = FALSE
      ]
    }

    # 至少保留有一个 finite 值的行
    keep_non_na <- rowSums(is.finite(heat_mat)) > 0
    heat_mat <- heat_mat[keep_non_na, , drop = FALSE]

    if (nrow(heat_mat) == 0) {
      warning("热图矩阵过滤后为空，跳过热图。")
      return(invisible(NULL))
    }

    # =========================
    # 4) 手动 row-scale
    #    后续 pheatmap 使用 scale = "none"，避免重复标准化。
    # =========================
    scale_rows <- function(mat) {
      row_mean <- rowMeans(mat, na.rm = TRUE)
      row_sd <- apply(mat, 1, stats::sd, na.rm = TRUE)

      # 如果某行所有值相同，row_sd 为 0；这里设为 1，标准化后该行全为 0。
      row_sd[is.na(row_sd) | row_sd == 0] <- 1

      mat_scaled <- sweep(mat, 1, row_mean, FUN = "-")
      mat_scaled <- sweep(mat_scaled, 1, row_sd, FUN = "/")
      mat_scaled[!is.finite(mat_scaled)] <- 0

      mat_scaled
    }

    heat_mat_scaled <- scale_rows(heat_mat)

    # =========================
    # 5) 导出热图矩阵
    # =========================
    if (do_save) {
      if (isTRUE(create_dir)) {
        dir.create(dirname(save_prefix), recursive = TRUE, showWarnings = FALSE)
      }

      f_heat_raw <- paste0(save_prefix, "_heatmap_matrix_raw.csv")
      f_heat_scaled <- paste0(save_prefix, "_heatmap_matrix_row_scaled.csv")

      utils::write.csv(
        data.frame(
          Gene = rownames(heat_mat),
          heat_mat,
          check.names = FALSE
        ),
        f_heat_raw,
        row.names = FALSE,
        fileEncoding = fileEncoding
      )

      utils::write.csv(
        data.frame(
          Gene = rownames(heat_mat_scaled),
          heat_mat_scaled,
          check.names = FALSE
        ),
        f_heat_scaled,
        row.names = FALSE,
        fileEncoding = fileEncoding
      )

      heatmap_matrix_saved_files <- c(f_heat_raw, f_heat_scaled)
    }

    # =========================
    # 6) 列注释
    # =========================
    ann_df <- data.frame(
      Group = factor(grp_fac, levels = groups),
      row.names = samples_used
    )

    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      pal_n <- max(3, min(8, length(groups)))
      base_pal <- RColorBrewer::brewer.pal(pal_n, "Set2")
      if (length(groups) > length(base_pal)) {
        ann_colors <- grDevices::colorRampPalette(base_pal)(length(groups))
      } else {
        ann_colors <- base_pal[seq_along(groups)]
      }
    } else {
      ann_colors <- grDevices::rainbow(length(groups))
    }

    annotation_colors <- list(Group = stats::setNames(ann_colors, groups))
    my_palette <- grDevices::colorRampPalette(c("#045a8d", "white", "#d7301f"))(1000)

    # =========================
    # 7) 画热图
    #    注意：如果 cluster_rows = TRUE，图上行顺序会被聚类重排；
    #    当前批量函数调用时传 cluster_rows = FALSE，所以图行顺序与导出矩阵一致。
    # =========================
    ph <- pheatmap::pheatmap(
      heat_mat_scaled,
      scale = "none",
      show_rownames = isTRUE(heat_show_rownames),
      show_colnames = isTRUE(heat_show_colnames),
      cluster_rows = isTRUE(cluster_rows),
      cluster_cols = isTRUE(cluster_cols),
      color = my_palette,
      annotation_col = ann_df,
      annotation_colors = annotation_colors,
      annotation_legend = TRUE,
      silent = TRUE
    )

    grid::grid.newpage()
    grid::grid.draw(ph$gtable)

    # =========================
    # 8) 保存图片
    # =========================
    if (do_save) {
      f_pdf <- paste0(save_prefix, "_heatmap.pdf")
      grDevices::pdf(f_pdf, width = 5, height = 6)
      grid::grid.newpage()
      grid::grid.draw(ph$gtable)
      grDevices::dev.off()

      f_png <- paste0(save_prefix, "_heatmap.png")
      grDevices::png(f_png, width = 5, height = 6, units = "in", res = 600)
      grid::grid.newpage()
      grid::grid.draw(ph$gtable)
      grDevices::dev.off()

      f_tiff <- paste0(save_prefix, "_heatmap.tiff")
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

      saved_files <- c(f_pdf, f_png, f_tiff)
    }

    invisible(list(
      pheatmap_obj = ph,
      heatmap_matrix = heat_mat,
      heatmap_matrix_row_scaled = heat_mat_scaled,
      genes_used = rownames(heat_mat),
      samples_used = colnames(heat_mat),
      heatmap_matrix_saved_files = heatmap_matrix_saved_files,
      saved_files = unique(c(saved_files, heatmap_matrix_saved_files))
    ))
  }

  multi_group_limma_de <- function(
    expr,
    sample_info,
    group_col = "group",
    groups,
    p_type = c("adjp", "pvalue"),
    p_thresh = 0.05,
    adjust_method = "BH",
    top_n = 50,
    save_prefix = NULL,
    save_formats = c("pdf", "png", "tiff"),
    tiff_compression = "lzw",
    create_dir = TRUE,
    heat_show_rownames = FALSE,
    heat_show_colnames = FALSE,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fileEncoding = "GBK"
  ) {
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("需要 limma 包：BiocManager::install('limma')")
    }

    p_type <- match.arg(tolower(p_type[1]), c("adjp", "pvalue"))
    pcol <- if (p_type == "adjp") "adj.P.Val" else "P.Value"

    expr2 <- as.matrix(expr)
    if (is.null(rownames(expr2))) {
      stop("expr 需要行名（feature/gene id）。")
    }
    if (is.null(colnames(expr2))) {
      stop("expr 需要列名（样本名）。")
    }

    sample_info2 <- as.data.frame(
      sample_info,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    if (is.null(rownames(sample_info2))) {
      stop("sample_info 必须有行名（样本名），且要能和 expr 列名对应。")
    }

    if (!group_col %in% colnames(sample_info2)) {
      stop("sample_info 中找不到 group_col='", group_col, "'。")
    }

    groups <- trimws(as.character(groups))
    groups <- groups[!is.na(groups) & nzchar(groups)]
    groups <- groups[!duplicated(groups)]

    if (length(groups) < multi_group_min_n) {
      stop("多组整体比较至少需要 ", multi_group_min_n, " 个组。")
    }

    missing_groups <- setdiff(
      groups,
      unique(as.character(sample_info2[[group_col]]))
    )

    if (length(missing_groups) > 0) {
      stop(
        "以下组不在 sample_info 的分组列中：",
        paste(missing_groups, collapse = ", ")
      )
    }

    samples_used <- rownames(sample_info2)[
      as.character(sample_info2[[group_col]]) %in% groups
    ]
    samples_used <- intersect(colnames(expr2), samples_used)

    if (length(samples_used) < length(groups)) {
      stop("可匹配样本过少，无法完成多组整体比较。")
    }

    X <- expr2[, samples_used, drop = FALSE]
    meta <- sample_info2[samples_used, , drop = FALSE]

    grp_raw <- as.character(meta[[group_col]])
    grp_fac <- factor(grp_raw, levels = groups)
    grp_safe <- factor(
      make.names(as.character(grp_fac)),
      levels = make.names(groups)
    )

    # =========================
    # 1) 多组整体 limma F-test
    # =========================
    design <- stats::model.matrix(~ grp_safe)
    fit <- limma::lmFit(X, design)
    fit <- limma::eBayes(fit)

    coef_idx <- seq_len(ncol(design))
    coef_idx <- coef_idx[coef_idx != 1]

    ttF <- limma::topTable(
      fit,
      coef = coef_idx,
      number = Inf,
      sort.by = "F",
      adjust.method = adjust_method
    )

    ttF$Gene <- rownames(ttF)
    rownames(ttF) <- ttF$Gene

    # =========================
    # 2) 各组均值和组间均值范围
    # =========================
    mean_by_group <- sapply(groups, function(g) {
      rowMeans(X[, grp_raw == g, drop = FALSE], na.rm = TRUE)
    })

    mean_by_group <- as.data.frame(mean_by_group, check.names = FALSE)
    colnames(mean_by_group) <- paste0("mean_", groups)

    mean_range <- apply(mean_by_group, 1, function(x) {
      max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    })

    which_max <- apply(mean_by_group, 1, function(x) {
      paste0("mean_", groups[which.max(x)])
    })

    which_min <- apply(mean_by_group, 1, function(x) {
      paste0("mean_", groups[which.min(x)])
    })

    # =========================
    # 3) 额外输出：其他组 vs baseline 的系数
    # =========================
    baseline <- groups[1]
    coef_tbl <- NULL

    for (g in groups[-1]) {
      coef_name <- paste0("grp_safe", make.names(g))
      out_name <- paste0("logFC_", g, "_vs_", baseline)

      if (coef_name %in% colnames(fit$coefficients)) {
        coef_tbl <- cbind(coef_tbl, fit$coefficients[, coef_name])
        colnames(coef_tbl)[ncol(coef_tbl)] <- out_name
      }
    }

    if (is.null(coef_tbl)) {
      coef_tbl <- data.frame(row.names = rownames(ttF))
    } else {
      coef_tbl <- as.data.frame(coef_tbl, check.names = FALSE)
    }

    # =========================
    # 4) 汇总结果表
    # =========================
    res_all <- cbind(
      feature = rownames(ttF),
      ttF,
      coef_tbl[rownames(ttF), , drop = FALSE],
      mean_by_group[rownames(ttF), , drop = FALSE],
      mean_range = mean_range[rownames(ttF)],
      highest_group_mean = which_max[rownames(ttF)],
      lowest_group_mean = which_min[rownames(ttF)]
    )

    rownames(res_all) <- as.character(res_all$feature)

    res_all$direction <- ifelse(
      !is.na(res_all[[pcol]]) & res_all[[pcol]] < p_thresh,
      "Sig",
      "Not Significant"
    )
    res_all$direction <- factor(
      res_all$direction,
      levels = c("Sig", "Not Significant")
    )

    sig_mask <- !is.na(res_all[[pcol]]) & res_all[[pcol]] < p_thresh
    res_sig <- res_all[sig_mask, , drop = FALSE]
    deg <- rownames(res_sig)

    summ <- table(res_all$direction, useNA = "ifany")

    common_features <- intersect(rownames(res_all), rownames(X))
    combined <- cbind(
      res_all[common_features, , drop = FALSE],
      X[common_features, , drop = FALSE]
    )

    # =========================
    # 5) 保存结果表 + 热图
    # =========================
    saved_files <- character(0)
    do_save <- !is.null(save_prefix) && nzchar(save_prefix)
    hm_out <- NULL

    if (do_save) {
      if (isTRUE(create_dir)) {
        dir.create(dirname(save_prefix), recursive = TRUE, showWarnings = FALSE)
      }

      f_all <- paste0(save_prefix, "_overall_F_test.csv")
      f_sig <- paste0(save_prefix, "_significant.csv")
      f_comb <- paste0(save_prefix, "_combined.csv")

      utils::write.csv(
        res_all,
        f_all,
        row.names = FALSE,
        fileEncoding = fileEncoding
      )

      utils::write.csv(
        res_sig,
        f_sig,
        row.names = FALSE,
        fileEncoding = fileEncoding
      )

      utils::write.csv(
        combined,
        f_comb,
        row.names = TRUE,
        fileEncoding = fileEncoding
      )

      saved_files <- c(saved_files, f_all, f_sig, f_comb)
    }

    # 不管是否保存表格，只要 save_prefix 有效就画图并保存热图；
    # 当前批量主函数会传入 save_prefix，所以正常会执行。
    if (do_save) {
      hm_out <- tryCatch(
        {
          plot_heatmap_multi_group_overall(
            res_all = res_all,
            expr_impute = X,
            sample_info = meta,
            groups = groups,
            group_col = group_col,
            pcol = pcol,
            p_thresh = p_thresh,
            top_n = top_n,
            save_prefix = paste0(save_prefix, "_overall"),
            create_dir = create_dir,
            cluster_rows = cluster_rows,
            cluster_cols = cluster_cols,
            heat_show_rownames = heat_show_rownames,
            heat_show_colnames = heat_show_colnames,
            tiff_compression = tiff_compression,
            fileEncoding = fileEncoding
          )
        },
        error = function(e) {
          warning("多组整体比较热图绘制失败：", conditionMessage(e))
          NULL
        }
      )

      if (!is.null(hm_out) && !is.null(hm_out$saved_files)) {
        saved_files <- unique(c(saved_files, hm_out$saved_files))
      }
    }

    group_sizes <- stats::setNames(
      as.integer(table(factor(grp_raw, levels = groups))),
      groups
    )

    list(
      results = res_all,
      summary = summ,
      combined = combined,
      design = design,
      contrast = paste(groups, collapse = " - "),
      samples_used = samples_used,
      deg = deg,
      volcano_plot = NULL,
      volcano_data = NULL,
      saved_files = unique(saved_files),
      method = "limma_F_test",

      heatmap_plot = if (!is.null(hm_out) && !is.null(hm_out$pheatmap_obj)) {
        hm_out$pheatmap_obj
      } else {
        NULL
      },

      heatmap_matrix = if (!is.null(hm_out) && !is.null(hm_out$heatmap_matrix)) {
        hm_out$heatmap_matrix
      } else {
        NULL
      },

      heatmap_matrix_row_scaled = if (
        !is.null(hm_out) && !is.null(hm_out$heatmap_matrix_row_scaled)
      ) {
        hm_out$heatmap_matrix_row_scaled
      } else {
        NULL
      },

      heatmap_matrix_saved_files = if (
        !is.null(hm_out) && !is.null(hm_out$heatmap_matrix_saved_files)
      ) {
        hm_out$heatmap_matrix_saved_files
      } else {
        character(0)
      },

      comparison_type = "multi_group",
      group_names = groups,
      n_groups = length(groups),
      group_sizes = group_sizes
    )
  }

  # =========================
  # 2) FC-only：单样本比较（两列）
  # =========================
  fc_only_de <- function(
    expr,
    sample1,
    sample2,
    logfc_thresh = 1,
    make_volcano = TRUE,
    label_extreme_n_each_side = 10,
    plot_width = 10,
    plot_height = 8,
    plot_dpi = 300,
    save_prefix = NULL,
    save_formats = c("pdf", "tiff", "png"),
    tiff_compression = "lzw",
    create_dir = TRUE,
    assume_log2 = NULL,
    add1 = 1
  ) {
    expr <- as.matrix(expr)
    if (!sample1 %in% colnames(expr)) {
      stop("sample1 不在 expr 列名中：", sample1)
    }
    if (!sample2 %in% colnames(expr)) {
      stop("sample2 不在 expr 列名中：", sample2)
    }

    x1 <- expr[, sample1]
    x2 <- expr[, sample2]

    if (is.null(assume_log2)) {
      assume_log2 <- guess_is_log2(c(x1, x2))
    }

    if (!assume_log2) {
      x1_log <- log2(x1 + add1)
      x2_log <- log2(x2 + add1)
    } else {
      x1_log <- x1
      x2_log <- x2
    }

    logFC <- x1_log - x2_log
    AveExpr <- (x1_log + x2_log) / 2

    direction <- ifelse(
      logFC >= logfc_thresh,
      "Up",
      ifelse(logFC <= -logfc_thresh, "Down", "Not Significant")
    )
    direction <- factor(direction, levels = c("Up", "Not Significant", "Down"))

    tt <- data.frame(
      Gene = rownames(expr),
      logFC = as.numeric(logFC),
      AveExpr = as.numeric(AveExpr),
      P.Value = NA_real_,
      adj.P.Val = NA_real_,
      direction = direction,
      stringsAsFactors = FALSE
    )
    rownames(tt) <- tt$Gene

    summ <- table(tt$direction, useNA = "ifany")
    combined <- cbind(tt, expr[, c(sample1, sample2), drop = FALSE])
    deg <- rownames(tt)[tt$direction %in% c("Up", "Down")]

    if (!is.null(save_prefix) && nzchar(save_prefix)) {
      if (create_dir) {
        dir.create(dirname(save_prefix), recursive = TRUE, showWarnings = FALSE)
      }
      utils::write.csv(
        combined,
        paste0(save_prefix, ".csv"),
        fileEncoding = "GBK"
      )
    }

    volcano_plot <- NULL
    volcano_saved_files <- character(0)

    if (isTRUE(make_volcano) && requireNamespace("ggplot2", quietly = TRUE)) {
      volcano_df <- tt
      volcano_df$yval <- abs(volcano_df$logFC) # FC-only：y 轴用 |log2FC|
      volcano_df$geneName <- volcano_df$Gene

      sig_df <- volcano_df[abs(volcano_df$logFC) > logfc_thresh, , drop = FALSE]
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

      gp <- ggplot2::ggplot(
        volcano_df,
        ggplot2::aes(x = logFC, y = yval, colour = direction)
      ) +
        ggplot2::geom_point(alpha = 0.4, size = 2.8) +
        ggplot2::scale_color_manual(
          values = col_map,
          limits = names(col_map),
          drop = TRUE
        ) +
        ggplot2::geom_vline(
          xintercept = c(-logfc_thresh, logfc_thresh),
          linetype = 4,
          color = "black",
          linewidth = 0.8
        ) +
        ggplot2::geom_hline(
          yintercept = logfc_thresh,
          linetype = 4,
          color = "black",
          linewidth = 0.8
        ) +
        ggplot2::labs(
          x = paste0("Log2FC (", sample1, " - ", sample2, ")"),
          y = "|Log2FC| (FC-only)",
          color = ""
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(aspect.ratio = 1, legend.position = "right")

      if (requireNamespace("ggrepel", quietly = TRUE)) {
        gp <- gp +
          ggrepel::geom_text_repel(
            data = subset(volcano_df, !is.na(labelGene)),
            ggplot2::aes(label = labelGene),
            size = 4,
            show.legend = FALSE,
            max.overlaps = Inf
          )
      }

      volcano_plot <- gp

      if (!is.null(save_prefix) && nzchar(save_prefix)) {
        if (create_dir) {
          dir.create(
            dirname(save_prefix),
            recursive = TRUE,
            showWarnings = FALSE
          )
        }
        for (fmt in tolower(save_formats)) {
          fpath <- paste0(save_prefix, "_volcano.", fmt)
          if (fmt %in% c("tif", "tiff")) {
            ggplot2::ggsave(
              fpath,
              gp,
              width = plot_width,
              height = plot_height,
              dpi = plot_dpi,
              device = "tiff",
              compression = tiff_compression
            )
          } else {
            ggplot2::ggsave(
              fpath,
              gp,
              width = plot_width,
              height = plot_height,
              dpi = plot_dpi
            )
          }
          volcano_saved_files <- c(volcano_saved_files, fpath)
        }
      }
    }

    list(
      results = tt,
      summary = summ,
      combined = combined,
      contrast = paste0(sample1, " - ", sample2),
      samples_used = c(sample1, sample2),
      deg = deg,
      volcano_plot = volcano_plot,
      volcano_saved_files = volcano_saved_files,
      method = "FC_only_single"
    )
  }

  # =========================
  # 3) 低重复组间：group_FC_only（用组均值/中位数做 log2FC，不算 p）
  # =========================
  group_fc_only_de <- function(
    expr,
    group1,
    group2,
    summary_fun = c("mean", "median"),
    logfc_thresh = 1,
    make_volcano = TRUE,
    label_extreme_n_each_side = 10,
    plot_width = 10,
    plot_height = 8,
    plot_dpi = 300,
    save_prefix = NULL,
    save_formats = c("pdf", "tiff", "png"),
    tiff_compression = "lzw",
    create_dir = TRUE,
    assume_log2 = NULL,
    add1 = 1,
    cluster = FALSE,
    heat_show_rownames = FALSE,
    heat_show_colnames = FALSE
  ) {
    summary_fun <- match.arg(summary_fun)

    s1 <- get_group_samples(group1)
    s2 <- get_group_samples(group2)
    if (length(s1) == 0 || length(s2) == 0) {
      stop(
        "group_fc_only_de：组在 expr 中匹配不到样本：",
        group1,
        " n=",
        length(s1),
        " ; ",
        group2,
        " n=",
        length(s2)
      )
    }

    X1 <- expr[, s1, drop = FALSE]
    X2 <- expr[, s2, drop = FALSE]

    if (is.null(assume_log2)) {
      assume_log2 <- guess_is_log2(c(as.numeric(X1), as.numeric(X2)))
    }

    if (!assume_log2) {
      X1 <- log2(X1 + add1)
      X2 <- log2(X2 + add1)
    }

    if (summary_fun == "mean") {
      g1 <- rowMeans(X1, na.rm = TRUE)
      g1[is.nan(g1)] <- NA
      g2 <- rowMeans(X2, na.rm = TRUE)
      g2[is.nan(g2)] <- NA
    } else {
      g1 <- apply(X1, 1, function(z) {
        if (all(is.na(z))) NA_real_ else stats::median(z, na.rm = TRUE)
      })
      g2 <- apply(X2, 1, function(z) {
        if (all(is.na(z))) NA_real_ else stats::median(z, na.rm = TRUE)
      })
    }

    logFC <- g1 - g2
    AveExpr <- (g1 + g2) / 2

    direction <- ifelse(
      logFC >= logfc_thresh,
      "Up",
      ifelse(logFC <= -logfc_thresh, "Down", "Not Significant")
    )
    direction <- factor(direction, levels = c("Up", "Not Significant", "Down"))

    tt <- data.frame(
      Gene = rownames(expr),
      logFC = as.numeric(logFC),
      AveExpr = as.numeric(AveExpr),
      P.Value = NA_real_,
      adj.P.Val = NA_real_,
      direction = direction,
      stringsAsFactors = FALSE
    )
    rownames(tt) <- tt$Gene

    summ <- table(tt$direction, useNA = "ifany")
    combined <- cbind(
      tt,
      group1_stat = g1,
      group2_stat = g2,
      expr[, c(s1, s2), drop = FALSE]
    )
    deg <- rownames(tt)[tt$direction %in% c("Up", "Down")]

    if (!is.null(save_prefix) && nzchar(save_prefix)) {
      if (create_dir) {
        dir.create(dirname(save_prefix), recursive = TRUE, showWarnings = FALSE)
      }
      utils::write.csv(
        combined,
        paste0(save_prefix, ".csv"),
        fileEncoding = "GBK"
      )
    }

    volcano_plot <- NULL
    volcano_saved_files <- character(0)

    # 火山图：y 用 |log2FC|
    if (isTRUE(make_volcano) && requireNamespace("ggplot2", quietly = TRUE)) {
      volcano_df <- tt
      volcano_df$yval <- abs(volcano_df$logFC)
      volcano_df$geneName <- volcano_df$Gene

      sig_df <- volcano_df[abs(volcano_df$logFC) > logfc_thresh, , drop = FALSE]
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

      gp <- ggplot2::ggplot(
        volcano_df,
        ggplot2::aes(x = logFC, y = yval, colour = direction)
      ) +
        ggplot2::geom_point(alpha = 0.4, size = 2.8) +
        ggplot2::scale_color_manual(
          values = col_map,
          limits = names(col_map),
          drop = TRUE
        ) +
        ggplot2::geom_vline(
          xintercept = c(-logfc_thresh, logfc_thresh),
          linetype = 4,
          color = "black",
          linewidth = 0.8
        ) +
        ggplot2::geom_hline(
          yintercept = logfc_thresh,
          linetype = 4,
          color = "black",
          linewidth = 0.8
        ) +
        ggplot2::labs(
          x = paste0(
            "Log2FC (",
            group1,
            " - ",
            group2,
            " ; ",
            summary_fun,
            ")"
          ),
          y = "|Log2FC| (Group FC-only)",
          color = ""
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(aspect.ratio = 1, legend.position = "right")

      if (requireNamespace("ggrepel", quietly = TRUE)) {
        gp <- gp +
          ggrepel::geom_text_repel(
            data = subset(volcano_df, !is.na(labelGene)),
            ggplot2::aes(label = labelGene),
            size = 4,
            show.legend = FALSE,
            max.overlaps = Inf
          )
      }

      volcano_plot <- gp

      if (!is.null(save_prefix) && nzchar(save_prefix)) {
        if (create_dir) {
          dir.create(
            dirname(save_prefix),
            recursive = TRUE,
            showWarnings = FALSE
          )
        }
        for (fmt in tolower(save_formats)) {
          fpath <- paste0(save_prefix, "_volcano.", fmt)
          if (fmt %in% c("tif", "tiff")) {
            ggplot2::ggsave(
              fpath,
              gp,
              width = plot_width,
              height = plot_height,
              dpi = plot_dpi,
              device = "tiff",
              compression = tiff_compression
            )
          } else {
            ggplot2::ggsave(
              fpath,
              gp,
              width = plot_width,
              height = plot_height,
              dpi = plot_dpi
            )
          }
          volcano_saved_files <- c(volcano_saved_files, fpath)
        }
      }
    }

    # 可选：热图（按样本展示，更直观；依然不涉及 p 值）
    heatmap_saved_files <- character(0)
    heatmap_matrix_saved_files <- character(0)
    heat_mat <- NULL
    heat_mat_scaled <- NULL

    scale_rows <- function(mat) {
      row_mean <- rowMeans(mat, na.rm = TRUE)
      row_sd <- apply(mat, 1, stats::sd, na.rm = TRUE)

      row_sd[is.na(row_sd) | row_sd == 0] <- 1

      mat_scaled <- sweep(mat, 1, row_mean, FUN = "-")
      mat_scaled <- sweep(mat_scaled, 1, row_sd, FUN = "/")
      mat_scaled[!is.finite(mat_scaled)] <- 0

      mat_scaled
    }

    if (
      requireNamespace("pheatmap", quietly = TRUE) &&
        requireNamespace("grid", quietly = TRUE)
    ) {
      do_save <- !is.null(save_prefix) && nzchar(save_prefix)
      if (do_save && create_dir) {
        dir.create(dirname(save_prefix), recursive = TRUE, showWarnings = FALSE)
      }

      filtered <- combined[
        combined$direction %in% c("Down", "Up"),
        ,
        drop = FALSE
      ]
      if (nrow(filtered) > 0) {
        # 按 logFC 排序
        filtered <- filtered[
          order(filtered$logFC, decreasing = FALSE),
          ,
          drop = FALSE
        ]

        # 样本列：group2 在左，group1 在右
        expr_cols <- c(s2, s1)
        heat_mat <- as.matrix(filtered[, expr_cols, drop = FALSE])
        rownames(heat_mat) <- filtered$Gene

        heat_mat_scaled <- scale_rows(heat_mat)

        if (do_save) {
          f_heat_raw <- paste0(save_prefix, "_heatmap_matrix_raw.csv")
          f_heat_scaled <- paste0(save_prefix, "_heatmap_matrix_row_scaled.csv")

          utils::write.csv(
            data.frame(
              Gene = rownames(heat_mat),
              heat_mat,
              check.names = FALSE
            ),
            f_heat_raw,
            row.names = FALSE,
            fileEncoding = "GBK"
          )

          utils::write.csv(
            data.frame(
              Gene = rownames(heat_mat_scaled),
              heat_mat_scaled,
              check.names = FALSE
            ),
            f_heat_scaled,
            row.names = FALSE,
            fileEncoding = "GBK"
          )

          heatmap_matrix_saved_files <- c(
            heatmap_matrix_saved_files,
            f_heat_raw,
            f_heat_scaled
          )
        }

        ann_df <- data.frame(
          Group = factor(
            c(rep(group2, length(s2)), rep(group1, length(s1))),
            levels = c(group2, group1)
          ),
          row.names = expr_cols
        )
        annotation_colors <- list(
          Group = setNames(c("#c74732", "#045a8d"), c(group2, group1))
        )

        my_palette <- grDevices::colorRampPalette(c(
          "#045a8d",
          "white",
          "#d7301f"
        ))(1000)

        ph <- pheatmap::pheatmap(
          heat_mat_scaled,
          scale = "none",
          show_rownames = isTRUE(heat_show_rownames),
          show_colnames = isTRUE(heat_show_colnames),
          cluster_rows = FALSE,
          cluster_cols = isTRUE(cluster),
          color = my_palette,
          annotation_col = ann_df,
          annotation_colors = annotation_colors,
          silent = TRUE
        )

        # 当前会话画一下
        grid::grid.newpage()
        grid::grid.draw(ph$gtable)

        if (do_save) {
          f_pdf <- paste0(save_prefix, "_heatmap.pdf")
          grDevices::pdf(f_pdf, width = 5, height = 6)
          grid::grid.newpage()
          grid::grid.draw(ph$gtable)
          grDevices::dev.off()

          f_png <- paste0(save_prefix, "_heatmap.png")
          grDevices::png(f_png, width = 5, height = 6, units = "in", res = 600)
          grid::grid.newpage()
          grid::grid.draw(ph$gtable)
          grDevices::dev.off()

          f_tiff <- paste0(save_prefix, "_heatmap.tiff")
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

          heatmap_saved_files <- c(f_pdf, f_png, f_tiff)
        }
      }
    }

    volcano_saved_files <- unique(c(
      volcano_saved_files,
      heatmap_saved_files,
      heatmap_matrix_saved_files
    ))

    list(
      results = tt,
      summary = summ,
      combined = combined,
      contrast = paste0(group1, " - ", group2),
      samples_used = c(s1, s2),
      deg = deg,
      volcano_plot = volcano_plot,
      volcano_saved_files = volcano_saved_files,

      heatmap_matrix = heat_mat,
      heatmap_matrix_row_scaled = heat_mat_scaled,
      heatmap_matrix_saved_files = heatmap_matrix_saved_files,

      method = paste0("FC_only_group_", summary_fun)
    )
  }

  # sheet 名安全
  safe_sheet_name <- function(x, max_len = 31) {
    x <- as.character(x)
    x <- gsub("[:\\\\/\\?\\*\\[\\]]", "_", x)
    x <- gsub("[[:space:]]+", " ", x)
    x <- gsub("^\\s+|\\s+$", "", x)
    if (nchar(x) > max_len) {
      x <- substr(x, 1, max_len)
    }
    if (nchar(x) == 0) {
      x <- "sheet"
    }
    x
  }

  make_unique <- function(x, max_len = 31) {
    out <- character(length(x))
    seen <- list()
    for (i in seq_along(x)) {
      nm <- x[i]
      if (is.null(seen[[nm]])) {
        seen[[nm]] <- 1
        out[i] <- nm
      } else {
        seen[[nm]] <- seen[[nm]] + 1
        suffix <- paste0("_", seen[[nm]])
        base <- nm
        if (nchar(base) + nchar(suffix) > max_len) {
          base <- substr(base, 1, max_len - nchar(suffix))
        }
        out[i] <- paste0(base, suffix)
      }
    }
    out
  }

  # =========================
  # 4) 自动判定逻辑（你要求的新版）
  # =========================
  group_levels <- unique(as.character(sample_info[[group_col]]))
  sample_names <- colnames(expr)

  detect_type <- function(a, b) {
    a <- as.character(a)
    b <- as.character(b)

    # 先判断是不是“组名对组名”
    is_group_pair <- (a %in% group_levels) && (b %in% group_levels)
    if (is_group_pair) {
      s1 <- get_group_samples(a)
      s2 <- get_group_samples(b)
      n1 <- length(s1)
      n2 <- length(s2)

      if (n1 == 0 || n2 == 0) {
        return(list(type = "unknown", n1 = n1, n2 = n2))
      }

      # ★你的规则：两组样本数都 >= min_reps_for_limma 才跑 limma
      if (min(n1, n2) >= min_reps_for_limma) {
        return(list(type = "group_limma", n1 = n1, n2 = n2))
      } else {
        # 低重复：默认走 group_fc_only（也可改 low_rep_strategy="limma"）
        if (low_rep_strategy == "limma" && min(n1, n2) >= 2) {
          return(list(type = "group_limma", n1 = n1, n2 = n2))
        } else {
          return(list(type = "group_fc_only", n1 = n1, n2 = n2))
        }
      }
    }

    # 再判断是不是“样本名对样本名”
    is_sample_pair <- (a %in% sample_names) && (b %in% sample_names)
    if (is_sample_pair) {
      return(list(type = "single_fc_only", n1 = 1, n2 = 1))
    }

    list(type = "unknown", n1 = NA_integer_, n2 = NA_integer_)
  }

  # =========================
  # 5) 主循环
  # =========================
  if (isTRUE(create_dir)) {
    dir.create(base_out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  results_list <- list()
  log_df <- data.frame(
    contrast_id = character(0),
    exp = character(0),
    ctrl = character(0),
    groups = character(0),
    type = character(0),
    n_exp = integer(0),
    n_ctrl = integer(0),
    n_groups = integer(0),
    n_samples = integer(0),
    status = character(0),
    message = character(0),
    stringsAsFactors = FALSE
  )

  # ---------- 5A) 两组 / 单样本比较 ----------
  if (
    !is.null(compare_list) &&
      is.data.frame(compare_list) &&
      nrow(compare_list) > 0
  ) {
    cols <- get_compare_cols(compare_list)
    exp_col2 <- cols$exp_col
    ctl_col2 <- cols$ctl_col

    for (i in seq_len(nrow(compare_list))) {
      a <- as.character(compare_list[[exp_col2]][i])
      b <- as.character(compare_list[[ctl_col2]][i])
      base_contrast_id <- paste0(a, "_vs_", b)
      contrast_id <- make_unique_id(base_contrast_id, names(results_list))

      dt <- detect_type(a, b)
      ctype <- dt$type

      subdir <- if (ctype == "group_limma") {
        "01_group_limma"
      } else if (ctype == "group_fc_only") {
        "01b_group_FC_only"
      } else if (ctype == "single_fc_only") {
        "02_single_FC_only"
      } else {
        "99_unknown"
      }

      contrast_dir <- file.path(base_out_dir, subdir, safe_slug(contrast_id))
      if (isTRUE(create_dir)) {
        dir.create(contrast_dir, recursive = TRUE, showWarnings = FALSE)
      }
      save_prefix <- file.path(contrast_dir, "DE")

      message(
        "[",
        i,
        "/",
        nrow(compare_list),
        "] ",
        ctype,
        ": ",
        contrast_id,
        " (n_exp=",
        dt$n1,
        ", n_ctrl=",
        dt$n2,
        ")"
      )

      if (ctype == "group_limma") {
        out_i <- tryCatch(
          {
            limma_de(
              expr = expr,
              sample_info = sample_info,
              group_col = group_col,
              group1 = a,
              group2 = b,
              p_type = p_type,
              p_thresh = p_thresh,
              logfc_thresh = logfc_thresh,
              adjust_method = adjust_method,
              make_volcano = make_volcano,
              label_extreme_n_each_side = label_extreme_n_each_side,
              plot_width = plot_width,
              plot_height = plot_height,
              plot_dpi = plot_dpi,
              heat_show_rownames = heat_show_rownames,
              heat_show_colnames = heat_show_colnames,
              save_prefix = save_prefix,
              save_formats = save_formats,
              tiff_compression = tiff_compression,
              cluster_cols = cluster,
              create_dir = create_dir
            )
          },
          error = function(e) e
        )

        if (inherits(out_i, "error")) {
          log_df <- rbind(
            log_df,
            data.frame(
              contrast_id = contrast_id,
              exp = a,
              ctrl = b,
              groups = paste(c(a, b), collapse = ";"),
              type = ctype,
              n_exp = dt$n1,
              n_ctrl = dt$n2,
              n_groups = 2L,
              n_samples = dt$n1 + dt$n2,
              status = "ERROR",
              message = conditionMessage(out_i),
              stringsAsFactors = FALSE
            )
          )
          next
        }

        out_i$method <- "limma"
        results_list[[contrast_id]] <- out_i

        log_df <- rbind(
          log_df,
          data.frame(
            contrast_id = contrast_id,
            exp = a,
            ctrl = b,
            groups = paste(c(a, b), collapse = ";"),
            type = ctype,
            n_exp = dt$n1,
            n_ctrl = dt$n2,
            n_groups = 2L,
            n_samples = dt$n1 + dt$n2,
            status = "OK",
            message = save_prefix,
            stringsAsFactors = FALSE
          )
        )
      } else if (ctype == "group_fc_only") {
        out_i <- tryCatch(
          {
            group_fc_only_de(
              expr = expr,
              group1 = a,
              group2 = b,
              summary_fun = "mean",
              logfc_thresh = logfc_thresh,
              make_volcano = make_volcano,
              label_extreme_n_each_side = label_extreme_n_each_side,
              plot_width = plot_width,
              plot_height = plot_height,
              plot_dpi = plot_dpi,
              save_prefix = save_prefix,
              save_formats = save_formats,
              tiff_compression = tiff_compression,
              create_dir = create_dir,
              assume_log2 = assume_log2,
              add1 = add1,
              cluster = cluster,
              heat_show_rownames = heat_show_rownames,
              heat_show_colnames = heat_show_colnames
            )
          },
          error = function(e) e
        )

        if (inherits(out_i, "error")) {
          log_df <- rbind(
            log_df,
            data.frame(
              contrast_id = contrast_id,
              exp = a,
              ctrl = b,
              groups = paste(c(a, b), collapse = ";"),
              type = ctype,
              n_exp = dt$n1,
              n_ctrl = dt$n2,
              n_groups = 2L,
              n_samples = dt$n1 + dt$n2,
              status = "ERROR",
              message = conditionMessage(out_i),
              stringsAsFactors = FALSE
            )
          )
          next
        }

        results_list[[contrast_id]] <- out_i
        log_df <- rbind(
          log_df,
          data.frame(
            contrast_id = contrast_id,
            exp = a,
            ctrl = b,
            groups = paste(c(a, b), collapse = ";"),
            type = ctype,
            n_exp = dt$n1,
            n_ctrl = dt$n2,
            n_groups = 2L,
            n_samples = dt$n1 + dt$n2,
            status = "OK",
            message = save_prefix,
            stringsAsFactors = FALSE
          )
        )
      } else if (ctype == "single_fc_only") {
        out_i <- tryCatch(
          {
            fc_only_de(
              expr = expr,
              sample1 = a,
              sample2 = b,
              logfc_thresh = logfc_thresh,
              make_volcano = make_volcano,
              label_extreme_n_each_side = label_extreme_n_each_side,
              plot_width = plot_width,
              plot_height = plot_height,
              plot_dpi = plot_dpi,
              save_prefix = save_prefix,
              save_formats = save_formats,
              tiff_compression = tiff_compression,
              create_dir = create_dir,
              assume_log2 = assume_log2,
              add1 = add1
            )
          },
          error = function(e) e
        )

        if (inherits(out_i, "error")) {
          log_df <- rbind(
            log_df,
            data.frame(
              contrast_id = contrast_id,
              exp = a,
              ctrl = b,
              groups = paste(c(a, b), collapse = ";"),
              type = ctype,
              n_exp = dt$n1,
              n_ctrl = dt$n2,
              n_groups = 2L,
              n_samples = dt$n1 + dt$n2,
              status = "ERROR",
              message = conditionMessage(out_i),
              stringsAsFactors = FALSE
            )
          )
          next
        }

        results_list[[contrast_id]] <- out_i
        log_df <- rbind(
          log_df,
          data.frame(
            contrast_id = contrast_id,
            exp = a,
            ctrl = b,
            groups = paste(c(a, b), collapse = ";"),
            type = ctype,
            n_exp = dt$n1,
            n_ctrl = dt$n2,
            n_groups = 2L,
            n_samples = dt$n1 + dt$n2,
            status = "OK",
            message = save_prefix,
            stringsAsFactors = FALSE
          )
        )
      } else {
        log_df <- rbind(
          log_df,
          data.frame(
            contrast_id = contrast_id,
            exp = a,
            ctrl = b,
            groups = paste(c(a, b), collapse = ";"),
            type = ctype,
            n_exp = dt$n1,
            n_ctrl = dt$n2,
            n_groups = 2L,
            n_samples = ifelse(
              all(is.finite(c(dt$n1, dt$n2))),
              dt$n1 + dt$n2,
              NA_integer_
            ),
            status = "SKIP",
            message = "无法判定（既不是组名对，也不是样本名对，或组在 expr 中无样本）",
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }

  # ---------- 5B) 多组整体比较（3组及以上） ----------
  multi_groups_list <- normalize_multi_compare(multi_compare)
  if (length(multi_groups_list) > 0) {
    for (i in seq_along(multi_groups_list)) {
      groups_i <- multi_groups_list[[i]]
      groups_i <- groups_i[!is.na(groups_i) & nzchar(groups_i)]
      groups_i <- groups_i[!duplicated(groups_i)]

      if (length(groups_i) < multi_group_min_n) {
        log_df <- rbind(
          log_df,
          data.frame(
            contrast_id = paste(groups_i, collapse = "_vs_"),
            exp = if (length(groups_i) >= 1) groups_i[1] else NA_character_,
            ctrl = if (length(groups_i) >= 2) {
              paste(groups_i[-1], collapse = "|")
            } else {
              NA_character_
            },
            groups = paste(groups_i, collapse = ";"),
            type = "multi_group_limma",
            n_exp = NA_integer_,
            n_ctrl = NA_integer_,
            n_groups = length(groups_i),
            n_samples = NA_integer_,
            status = "SKIP",
            message = paste0(
              "该行有效组数少于 ",
              multi_group_min_n,
              "，跳过多组整体比较"
            ),
            stringsAsFactors = FALSE
          )
        )
        next
      }

      missing_groups <- setdiff(groups_i, group_levels)
      sample_counts <- sapply(groups_i, function(g) {
        length(get_group_samples(g))
      })
      base_contrast_id <- paste(groups_i, collapse = "_vs_")
      contrast_id <- make_unique_id(base_contrast_id, names(results_list))
      contrast_dir <- file.path(
        base_out_dir,
        "03_multi_group_limma",
        safe_slug(contrast_id)
      )
      if (isTRUE(create_dir)) {
        dir.create(contrast_dir, recursive = TRUE, showWarnings = FALSE)
      }
      save_prefix <- file.path(contrast_dir, "DE_multi")

      message(
        "[multi ",
        i,
        "/",
        length(multi_groups_list),
        "] multi_group_limma: ",
        contrast_id,
        " (groups=",
        length(groups_i),
        ", n_samples=",
        sum(sample_counts),
        ")"
      )

      if (length(missing_groups) > 0) {
        log_df <- rbind(
          log_df,
          data.frame(
            contrast_id = contrast_id,
            exp = groups_i[1],
            ctrl = paste(groups_i[-1], collapse = "|"),
            groups = paste(groups_i, collapse = ";"),
            type = "multi_group_limma",
            n_exp = sample_counts[1],
            n_ctrl = sum(sample_counts[-1]),
            n_groups = length(groups_i),
            n_samples = sum(sample_counts),
            status = "ERROR",
            message = paste0(
              "以下组不在 sample_info 的分组列中：",
              paste(missing_groups, collapse = ", ")
            ),
            stringsAsFactors = FALSE
          )
        )
        next
      }

      if (any(sample_counts < min_reps_for_multi_limma)) {
        log_df <- rbind(
          log_df,
          data.frame(
            contrast_id = contrast_id,
            exp = groups_i[1],
            ctrl = paste(groups_i[-1], collapse = "|"),
            groups = paste(groups_i, collapse = ";"),
            type = "multi_group_limma",
            n_exp = sample_counts[1],
            n_ctrl = sum(sample_counts[-1]),
            n_groups = length(groups_i),
            n_samples = sum(sample_counts),
            status = "SKIP",
            message = paste0(
              "至少有一组样本数 < min_reps_for_multi_limma(",
              min_reps_for_multi_limma,
              ")：",
              paste(paste0(groups_i, "=", sample_counts), collapse = ", ")
            ),
            stringsAsFactors = FALSE
          )
        )
        next
      }

      out_i <- tryCatch(
        {
          multi_group_limma_de(
            expr = expr,
            sample_info = sample_info,
            group_col = group_col,
            groups = groups_i,
            p_type = p_type,
            p_thresh = p_thresh,
            adjust_method = adjust_method,
            top_n = multi_top_n,
            save_prefix = save_prefix,
            save_formats = save_formats,
            tiff_compression = tiff_compression,
            create_dir = create_dir,
            heat_show_rownames = heat_show_rownames,
            heat_show_colnames = heat_show_colnames,
            cluster_cols = cluster
          )
        },
        error = function(e) e
      )

      if (inherits(out_i, "error")) {
        log_df <- rbind(
          log_df,
          data.frame(
            contrast_id = contrast_id,
            exp = groups_i[1],
            ctrl = paste(groups_i[-1], collapse = "|"),
            groups = paste(groups_i, collapse = ";"),
            type = "multi_group_limma",
            n_exp = sample_counts[1],
            n_ctrl = sum(sample_counts[-1]),
            n_groups = length(groups_i),
            n_samples = sum(sample_counts),
            status = "ERROR",
            message = conditionMessage(out_i),
            stringsAsFactors = FALSE
          )
        )
        next
      }

      results_list[[contrast_id]] <- out_i
      log_df <- rbind(
        log_df,
        data.frame(
          contrast_id = contrast_id,
          exp = groups_i[1],
          ctrl = paste(groups_i[-1], collapse = "|"),
          groups = paste(groups_i, collapse = ";"),
          type = "multi_group_limma",
          n_exp = sample_counts[1],
          n_ctrl = sum(sample_counts[-1]),
          n_groups = length(groups_i),
          n_samples = sum(sample_counts),
          status = "OK",
          message = save_prefix,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  if (nrow(log_df) == 0) {
    warning(
      "未执行任何比较：compare_list 为空，且 multi_compare 也为空或无有效行。"
    )
  }

  status_table <- table(log_df$status)
  colnames(log_df) <- c(
    "比较对",
    "实验组",
    "对照组",
    "组集合",
    "比较类型(组间/单样本/多组整体)",
    "实验组样本量",
    "对照组样本量",
    "组数",
    "总样本量",
    "分析状态",
    "结果路径前缀"
  )
  # =========================
  # 6) 汇总导出 Excel（每对比一个 sheet + log）
  # =========================
  out_xlsx <- NULL
  if (isTRUE(export_excel)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      warning(
        "缺少 openxlsx：install.packages('openxlsx') 后可输出 Excel 汇总。"
      )
    } else {
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "log")
      openxlsx::writeDataTable(wb, "log", log_df)

      contrast_names <- names(results_list)
      sheet_names <- make_unique(vapply(
        contrast_names,
        safe_sheet_name,
        FUN.VALUE = character(1)
      ))

      for (k in seq_along(contrast_names)) {
        cn <- contrast_names[k]
        sn <- sheet_names[k]
        out_i <- results_list[[cn]]
        df_i <- out_i$combined

        openxlsx::addWorksheet(wb, sn)
        if (is.null(df_i) || !is.data.frame(df_i) || nrow(df_i) == 0) {
          openxlsx::writeData(
            wb,
            sn,
            data.frame(note = paste0("No combined table for: ", cn))
          )
        } else {
          if (!is.null(rownames(df_i)) && any(nzchar(rownames(df_i)))) {
            df_i <- cbind(feature_id = rownames(df_i), df_i)
            rownames(df_i) <- NULL
          }
          openxlsx::writeDataTable(wb, sn, df_i)
        }
      }

      out_xlsx <- file.path(base_out_dir, "差异分析结果汇总表.xlsx")
      openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
    }
  }

  # =========================
  # 7) 汇总 Up/Down 数量表
  # =========================
  diff_summary_df <- NULL
  out_deg_xlsx <- NULL

  if (isTRUE(export_deg_summary) && length(results_list) > 0) {
    diff_summary_df <- do.call(
      rbind,
      lapply(names(results_list), function(contrast_id) {
        out_i <- results_list[[contrast_id]]
        summ <- out_i$summary
        up_n <- if ("Up" %in% names(summ)) as.integer(summ[["Up"]]) else 0L
        down_n <- if ("Down" %in% names(summ)) {
          as.integer(summ[["Down"]])
        } else {
          0L
        }
        overall_sig_n <- if ("Sig" %in% names(summ)) {
          as.integer(summ[["Sig"]])
        } else {
          0L
        }
        deg_n <- if (!is.null(out_i$deg)) {
          length(unique(na.omit(as.character(out_i$deg))))
        } else {
          0L
        }
        sig_total <- max(up_n + down_n, overall_sig_n, deg_n, na.rm = TRUE)

        data.frame(
          contrast = contrast_id,
          up = up_n,
          down = down_n,
          overall_sig = overall_sig_n,
          sig_total = sig_total,
          method = if (!is.null(out_i$method)) out_i$method else NA_character_,
          stringsAsFactors = FALSE
        )
      })
    )

    diff_summary_df <- diff_summary_df[
      order(diff_summary_df$sig_total, decreasing = TRUE),
    ]
    diff_summary_df$contrast <- gsub("_vs_", " vs ", diff_summary_df$contrast)
    colnames(diff_summary_df) <- c(
      "比较对",
      "上调差异基因/蛋白数量",
      "下调差异基因/蛋白数量",
      "整体显著差异基因/蛋白数量",
      "差异基因/蛋白总数",
      "比较方法"
    )
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      wb2 <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb2, "DEG_count_summary")
      openxlsx::writeDataTable(wb2, "DEG_count_summary", diff_summary_df)
      out_deg_xlsx <- file.path(base_out_dir, "DEG_count_summary.xlsx")
      openxlsx::saveWorkbook(wb2, out_deg_xlsx, overwrite = TRUE)
    }
  }

  list(
    results_list = results_list,
    log_df = log_df,
    status_table = status_table,
    diff_summary_df = diff_summary_df,
    out_xlsx = out_xlsx,
    out_deg_xlsx = out_deg_xlsx,
    base_out_dir = base_out_dir
  )
}


if (F) {
  out <- run_batch_de_mixed(
    expr = expr_impute,
    sample_info = dep_id_map,
    group_col = "group",
    compare_list = compare_list2,
    multi_compare = data.frame(group1 = "sg2", group2 = "sg1", group3 = "WT"),
    base_out_dir = "../Results/差异分析",
    logfc_thresh = 1,

    # 按你说的：≥3 才 limma；<3 走 group_fc_only
    min_reps_for_limma = 3,
    low_rep_strategy = "group_fc_only"
  )

  out$status_table
  head(out$log_df)
  out$diff_summary_df
}
