# ============================================================
# 04_plot_ksea_kinase_substrate_chord.R
#
# 目的：
#   基于 03_run_phosphosite_ksea.R 的中文输出表，
#   绘制 kinase-substrate phosphosite 关系弦图。
#
# 适配当前 03 脚本输出文件名：
#   *_KSEA 磷酸化位点 logFC 向量表.csv
#   *_KSEA 上调方向富集结果.csv
#   *_KSEA 下调方向富集结果.csv
#   *_KSEA 激酶活性趋势评分表.csv
#
# 输出：
#   ./demo/Results/phosphosite_KSEA/kinase_substrate_chord/
#
# 注意：
#   这里画的是 “激酶 - 底物磷酸化位点” 关系图。
#   KSEA 使用的是 site_label，例如 RPS6;S240;，
#   不是直接使用肽段序列。
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(circlize)
})

# ------------------------------------------------------------
# 0. 小工具函数
# ------------------------------------------------------------

.safe_dir_create <- function(path) {
  if (!dir.exists(path)) {
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

.read_csv_base <- function(file) {
  read.csv(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# PhosR 的 substrates 字符串通常类似：
#   RPS6;S240;;RPS6;S236;;EIF4B;S422;
#
# 不能直接按 ; 切分，因为每个 site_label 内部本身就有分号。
# 所以用正则提取完整 site_label：
#   GENE;S123;
#   GENE;T456;
#   GENE;Y789;
.parse_phosr_substrates <- function(x) {
  if (is.na(x) || x == "") {
    return(character(0))
  }
  
  sites <- stringr::str_extract_all(
    x,
    pattern = "[^;]+;[STY][0-9]+;"
  )[[1]]
  
  unique(sites)
}

# ------------------------------------------------------------
# 1. 在某个 comparison 文件夹中自动找中文 KSEA 文件
# ------------------------------------------------------------

.find_one_file <- function(dir, pattern, required = TRUE) {
  files <- list.files(
    dir,
    pattern = pattern,
    full.names = TRUE
  )
  
  if (length(files) == 0) {
    if (required) {
      stop("在目录中找不到文件: ", dir, "\npattern: ", pattern)
    } else {
      return(NA_character_)
    }
  }
  
  files[1]
}

.find_ksea_files_in_comp_dir <- function(comp_dir) {
  list(
    input_logfc = .find_one_file(
      comp_dir,
      pattern = "KSEA 磷酸化位点 logFC 向量表\\.csv$"
    ),
    up_result = .find_one_file(
      comp_dir,
      pattern = "KSEA 上调方向富集结果\\.csv$"
    ),
    down_result = .find_one_file(
      comp_dir,
      pattern = "KSEA 下调方向富集结果\\.csv$"
    ),
    signed_score = .find_one_file(
      comp_dir,
      pattern = "KSEA 激酶活性趋势评分表\\.csv$"
    )
  )
}

# ------------------------------------------------------------
# 2. 从一个 comparison 的 4 个 KSEA 文件生成边表
# ------------------------------------------------------------

make_ksea_chord_edge_table_from_chinese_outputs <- function(
    input_logfc_file,
    up_result_file,
    down_result_file,
    signed_score_file,
    top_n_kinases = 8,
    max_sites_per_kinase = 10,
    min_abs_signed_score = 0
) {
  input_df <- .read_csv_base(input_logfc_file)
  up_df <- .read_csv_base(up_result_file)
  down_df <- .read_csv_base(down_result_file)
  signed_df <- .read_csv_base(signed_score_file)
  
  # ----------------------------------------------------------
  # 检查必要列
  # ----------------------------------------------------------
  
  required_input_cols <- c("site_label", "logFC")
  missing_input_cols <- setdiff(required_input_cols, colnames(input_df))
  if (length(missing_input_cols) > 0) {
    stop(
      "KSEA 输入 logFC 表缺少列: ",
      paste(missing_input_cols, collapse = ", ")
    )
  }
  
  required_signed_cols <- c("Kinase", "signed_score", "dominant_direction")
  missing_signed_cols <- setdiff(required_signed_cols, colnames(signed_df))
  if (length(missing_signed_cols) > 0) {
    stop(
      "KSEA 激酶活性趋势评分表缺少列: ",
      paste(missing_signed_cols, collapse = ", ")
    )
  }
  
  required_enrich_cols <- c("Kinase", "substrates")
  missing_up_cols <- setdiff(required_enrich_cols, colnames(up_df))
  missing_down_cols <- setdiff(required_enrich_cols, colnames(down_df))
  
  if (length(missing_up_cols) > 0) {
    stop(
      "KSEA 上调方向富集结果表缺少列: ",
      paste(missing_up_cols, collapse = ", ")
    )
  }
  
  if (length(missing_down_cols) > 0) {
    stop(
      "KSEA 下调方向富集结果表缺少列: ",
      paste(missing_down_cols, collapse = ", ")
    )
  }
  
  # comparison 名称优先从 signed_df 里取。
  comparison_name <- if ("comparison" %in% colnames(signed_df)) {
    unique(as.character(signed_df$comparison))[1]
  } else {
    basename(dirname(signed_score_file))
  }
  
  # ----------------------------------------------------------
  # 选择要展示的 kinase/signature
  # ----------------------------------------------------------
  
  top_kinases <- signed_df %>%
    filter(
      !is.na(Kinase),
      Kinase != "",
      !is.na(signed_score),
      abs(signed_score) >= min_abs_signed_score
    ) %>%
    arrange(desc(abs(signed_score))) %>%
    slice_head(n = top_n_kinases)
  
  if (nrow(top_kinases) == 0) {
    warning("没有可用于弦图的 kinase/signature: ", comparison_name)
    return(data.frame())
  }
  
  # ----------------------------------------------------------
  # 根据 dominant_direction 选择对应方向的 substrates
  #
  # 如果某个 kinase 是 UP_or_greater，
  #   就从 上调方向富集结果表 中取 substrates。
  #
  # 如果某个 kinase 是 DOWN_or_less，
  #   就从 下调方向富集结果表 中取 substrates。
  #
  # 这样画图时，连线展示的是支撑当前主方向的底物位点。
  # ----------------------------------------------------------
  
  up_sub <- up_df %>%
    select(Kinase, substrates_up = substrates)
  
  down_sub <- down_df %>%
    select(Kinase, substrates_down = substrates)
  
  top_kinases2 <- top_kinases %>%
    left_join(up_sub, by = "Kinase") %>%
    left_join(down_sub, by = "Kinase") %>%
    mutate(
      substrates_for_plot = dplyr::case_when(
        dominant_direction == "UP_or_greater" ~ substrates_up,
        dominant_direction == "DOWN_or_less" ~ substrates_down,
        TRUE ~ substrates_up
      )
    )
  
  # ----------------------------------------------------------
  # 把 substrates 拆成长表：Kinase - substrate_site_label
  # ----------------------------------------------------------
  
  edge_df <- bind_rows(lapply(seq_len(nrow(top_kinases2)), function(i) {
    sites <- .parse_phosr_substrates(top_kinases2$substrates_for_plot[i])
    
    if (length(sites) == 0) {
      return(NULL)
    }
    
    data.frame(
      comparison = comparison_name,
      Kinase = top_kinases2$Kinase[i],
      kinase_signed_score = top_kinases2$signed_score[i],
      kinase_direction = top_kinases2$dominant_direction[i],
      substrate_site_label = sites,
      stringsAsFactors = FALSE
    )
  }))
  
  if (nrow(edge_df) == 0) {
    warning("substrates 解析后为空: ", comparison_name)
    return(data.frame())
  }
  
  # ----------------------------------------------------------
  # 接回每个 substrate phosphosite 的 logFC
  # ----------------------------------------------------------
  
  input_logfc <- input_df %>%
    select(site_label, logFC) %>%
    distinct(site_label, .keep_all = TRUE)
  
  edge_df <- edge_df %>%
    left_join(
      input_logfc,
      by = c("substrate_site_label" = "site_label")
    ) %>%
    mutate(
      substrate_gene = stringr::str_replace(substrate_site_label, ";.*$", ""),
      substrate_residue_site = stringr::str_match(
        substrate_site_label,
        ";([STY][0-9]+);"
      )[, 2],
      substrate_display = paste0(substrate_gene, " ", substrate_residue_site)
    )
  
  # 每个 kinase 最多展示 max_sites_per_kinase 个底物位点。
  # 如果底物太多，弦图会完全糊掉。
  # 优先展示 abs(logFC) 大的位点。
  edge_df <- edge_df %>%
    group_by(Kinase) %>%
    arrange(desc(!is.na(logFC)), desc(abs(logFC)), .by_group = TRUE) %>%
    slice_head(n = max_sites_per_kinase) %>%
    ungroup()
  
  edge_df
}

# ------------------------------------------------------------
# 3. 绘制单个 comparison 的弦图
# ------------------------------------------------------------

plot_ksea_kinase_substrate_chord <- function(
    edge_df,
    out_prefix,
    plot_title = NULL,
    fc_limit = NULL,
    width = 8,
    height = 8
) {
  if (nrow(edge_df) == 0) {
    warning("edge_df 为空，跳过绘图。")
    return(invisible(NULL))
  }
  
  .safe_dir_create(dirname(out_prefix))
  
  # 为避免 kinase 和 substrate 同名冲突，内部节点名加前缀。
  edge_df <- edge_df %>%
    mutate(
      kinase_node = paste0("KINASE|", Kinase),
      substrate_node = paste0(
        "SUBSTRATE|",
        substrate_gene,
        "_",
        substrate_residue_site
      )
    )
  
  kinase_nodes <- unique(edge_df$kinase_node)
  substrate_nodes <- unique(edge_df$substrate_node)
  sector_order <- c(kinase_nodes, substrate_nodes)
  
  kinases <- unique(edge_df$Kinase)
  kinase_cols <- setNames(
    grDevices::hcl.colors(length(kinases), palette = "Dark 3"),
    kinases
  )
  
  substrate_genes <- unique(edge_df$substrate_gene)
  substrate_gene_cols <- setNames(
    grDevices::hcl.colors(length(substrate_genes), palette = "Set 3"),
    substrate_genes
  )
  
  # 第一圈：红色 kinase，绿色 substrate
  role_col <- setNames(
    ifelse(
      startsWith(sector_order, "KINASE|"),
      "#D7191C",
      "#1A9641"
    ),
    sector_order
  )
  
  # 第三圈：kinase 或 substrate gene 色块
  name_block_col <- setNames(rep("grey90", length(sector_order)), sector_order)
  
  for (k in kinases) {
    name_block_col[paste0("KINASE|", k)] <- kinase_cols[k]
  }
  
  substrate_node_to_gene <- edge_df %>%
    select(substrate_node, substrate_gene) %>%
    distinct()
  
  for (i in seq_len(nrow(substrate_node_to_gene))) {
    node <- substrate_node_to_gene$substrate_node[i]
    gene <- substrate_node_to_gene$substrate_gene[i]
    name_block_col[node] <- substrate_gene_cols[gene]
  }
  
  # 第二圈：底物磷酸化位点 logFC 颜色
  substrate_fc <- edge_df %>%
    group_by(substrate_node) %>%
    summarise(
      logFC = if (all(is.na(logFC))) NA_real_ else median(logFC, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (is.null(fc_limit)) {
    fc_limit <- suppressWarnings(max(abs(substrate_fc$logFC), na.rm = TRUE))
    if (!is.finite(fc_limit) || fc_limit == 0) {
      fc_limit <- 1
    }
  }
  
  fc_fun <- circlize::colorRamp2(
    c(-fc_limit, 0, fc_limit),
    c("#2C7BB6", "white", "#D7191C")
  )
  
  fc_col <- setNames(rep(NA_character_, length(sector_order)), sector_order)
  
  for (i in seq_len(nrow(substrate_fc))) {
    node <- substrate_fc$substrate_node[i]
    val <- substrate_fc$logFC[i]
    
    if (!is.na(val)) {
      fc_col[node] <- fc_fun(val)
    }
  }
  
  # 标签
  label_map <- setNames(rep("", length(sector_order)), sector_order)
  
  for (k in kinases) {
    label_map[paste0("KINASE|", k)] <- k
  }
  
  substrate_label_df <- edge_df %>%
    select(substrate_node, substrate_display) %>%
    distinct()
  
  for (i in seq_len(nrow(substrate_label_df))) {
    label_map[substrate_label_df$substrate_node[i]] <-
      substrate_label_df$substrate_display[i]
  }
  
  # 弦图连接关系
  link_df <- edge_df %>%
    select(kinase_node, substrate_node, Kinase) %>%
    distinct() %>%
    mutate(value = 1)
  
  link_col <- kinase_cols[link_df$Kinase]
  
  # kinase 区和 substrate 区之间留大 gap
  gap_after <- rep(1, length(sector_order))
  gap_after[length(kinase_nodes)] <- 8
  gap_after[length(sector_order)] <- 8
  
  draw_chord <- function() {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    
    # 减少图形四周留白，让标题和图例离圆环更近
    graphics::par(
      mar = c(0, 0, 0.2, 0),
      oma = c(0, 0, 0, 0),
      xaxs = "i",
      yaxs = "i",
      xpd = NA
    )
    
    circlize::circos.clear()
    
    circlize::circos.par(
      start.degree = 90,
      gap.after = gap_after,
      track.margin = c(0.002, 0.002),
      canvas.xlim = c(-1.12, 1.12),
      canvas.ylim = c(-1.12, 1.12)
    )
    
    circlize::chordDiagram(
      x = link_df[, c("kinase_node", "substrate_node", "value")],
      order = sector_order,
      col = link_col,
      transparency = 0.35,
      annotationTrack = NULL,
      preAllocateTracks = list(
        list(track.height = 0.025),  # 第一圈：Kinase/Substrate 类型
        list(track.height = 0.035),  # 第二圈：logFC
        list(track.height = 0.080)   # 第三圈：名称/基因色块
      )
    )
    
    # 第一圈：节点类型
    circlize::circos.trackPlotRegion(
      track.index = 1,
      bg.border = NA,
      panel.fun = function(x, y) {
        sector <- circlize::get.cell.meta.data("sector.index")
        xlim <- circlize::get.cell.meta.data("xlim")
        ylim <- circlize::get.cell.meta.data("ylim")
        
        circlize::circos.rect(
          xlim[1], ylim[1], xlim[2], ylim[2],
          col = role_col[sector],
          border = NA
        )
      }
    )
    
    # 第二圈：logFC
    # kinase 没有 phosphosite logFC，所以保持白色。
    # 底物位点如果 logFC 缺失，也保持白色。
    circlize::circos.trackPlotRegion(
      track.index = 2,
      bg.border = NA,
      panel.fun = function(x, y) {
        sector <- circlize::get.cell.meta.data("sector.index")
        xlim <- circlize::get.cell.meta.data("xlim")
        ylim <- circlize::get.cell.meta.data("ylim")
        
        fill_col <- fc_col[sector]
        if (is.na(fill_col)) {
          fill_col <- "white"
        }
        
        circlize::circos.rect(
          xlim[1], ylim[1], xlim[2], ylim[2],
          col = fill_col,
          border = "grey85",
          lwd = 0.3
        )
      }
    )
    
    # 第三圈：名称和色块
    circlize::circos.trackPlotRegion(
      track.index = 3,
      bg.border = NA,
      panel.fun = function(x, y) {
        sector <- circlize::get.cell.meta.data("sector.index")
        xlim <- circlize::get.cell.meta.data("xlim")
        ylim <- circlize::get.cell.meta.data("ylim")
        
        circlize::circos.rect(
          xlim[1], ylim[1], xlim[2], ylim[2],
          col = name_block_col[sector],
          border = "white",
          lwd = 0.4
        )
        
        circlize::circos.text(
          x = mean(xlim),
          y = mean(ylim),
          labels = label_map[sector],
          facing = "clockwise",
          niceFacing = TRUE,
          adj = c(0.5, 0.5),
          cex = 0.45
        )
      }
    )
    
    if (!is.null(plot_title)) {
      graphics::mtext(
        text = plot_title,
        side = 3,
        line = -1.5,
        cex = 0.9,
        font = 2
      )
    }
    
    legend(
      x = -1.00,
      y = 1.00,
      legend = c("Kinase", "Substrate", "logFC high", "logFC low"),
      fill = c("#D7191C", "#1A9641", "#D7191C", "#2C7BB6"),
      border = NA,
      bty = "n",
      cex = 0.7,
      x.intersp = 0.6,
      y.intersp = 0.8
    )
    
    circlize::circos.clear()
  }
  
  pdf_file <- paste0(out_prefix, "_kinase_substrate_chord.pdf")
  png_file <- paste0(out_prefix, "_kinase_substrate_chord.png")
  edge_file <- paste0(out_prefix, "_kinase_substrate_edge_table.csv")
  
  grDevices::pdf(pdf_file, width = width, height = height)
  draw_chord()
  grDevices::dev.off()
  
  grDevices::png(png_file, width = width, height = height, units = "in", res = 600)
  draw_chord()
  grDevices::dev.off()
  
  utils::write.csv(edge_df, edge_file, row.names = FALSE)
  
  invisible(list(
    pdf = pdf_file,
    png = png_file,
    edge_table = edge_file
  ))
}

# ------------------------------------------------------------
# 4. 批量绘制所有 comparison 的弦图
# ------------------------------------------------------------

plot_all_ksea_kinase_substrate_chord <- function(
    ksea_dir = "./demo/Results/phosphosite_KSEA",
    out_dir = file.path(ksea_dir, "kinase_substrate_chord"),
    top_n_kinases = 8,
    max_sites_per_kinase = 10,
    min_abs_signed_score = 0,
    comparisons = NULL
) {
  .safe_dir_create(out_dir)
  
  # 当前中文 03 脚本每个 comparison 有一个子目录。
  # 子目录里有：
  #   *_KSEA 激酶活性趋势评分表.csv
  comp_dirs <- list.dirs(ksea_dir, recursive = FALSE, full.names = TRUE)
  
  comp_dirs <- comp_dirs[
    vapply(comp_dirs, function(d) {
      length(list.files(d, pattern = "KSEA 激酶活性趋势评分表\\.csv$")) > 0
    }, logical(1))
  ]
  
  if (length(comp_dirs) == 0) {
    stop(
      "没有找到 comparison 子目录下的中文 KSEA 结果表。\n",
      "请确认已经运行 03 脚本，并且目录为: ", ksea_dir
    )
  }
  
  all_edges <- list()
  all_outputs <- list()
  
  for (comp_dir in comp_dirs) {
    files <- .find_ksea_files_in_comp_dir(comp_dir)
    
    signed_tmp <- .read_csv_base(files$signed_score)
    comparison_name <- if ("comparison" %in% colnames(signed_tmp)) {
      unique(as.character(signed_tmp$comparison))[1]
    } else {
      basename(comp_dir)
    }
    
    if (!is.null(comparisons) && !(comparison_name %in% comparisons)) {
      next
    }
    
    message("绘制 KSEA kinase-substrate chord: ", comparison_name)
    
    edge_df <- make_ksea_chord_edge_table_from_chinese_outputs(
      input_logfc_file = files$input_logfc,
      up_result_file = files$up_result,
      down_result_file = files$down_result,
      signed_score_file = files$signed_score,
      top_n_kinases = top_n_kinases,
      max_sites_per_kinase = max_sites_per_kinase,
      min_abs_signed_score = min_abs_signed_score
    )
    
    if (nrow(edge_df) == 0) {
      next
    }
    
    comp_slug <- .safe_slug(comparison_name)
    comp_out_dir <- file.path(out_dir, comp_slug)
    .safe_dir_create(comp_out_dir)
    
    out_prefix <- file.path(
      comp_out_dir,
      paste0(comp_slug, "_top", top_n_kinases, "_kinases")
    )
    
    outputs <- plot_ksea_kinase_substrate_chord(
      edge_df = edge_df,
      out_prefix = out_prefix,
      plot_title = comparison_name
    )
    
    all_edges[[comparison_name]] <- edge_df
    all_outputs[[comparison_name]] <- outputs
  }
  
  if (length(all_edges) > 0) {
    all_edge_df <- bind_rows(all_edges)
    
    utils::write.csv(
      all_edge_df,
      file.path(out_dir, "KSEA_激酶_底物磷酸化位点关系边表.csv"),
      row.names = FALSE
    )
  } else {
    all_edge_df <- data.frame()
  }
  
  invisible(list(
    edges = all_edges,
    outputs = all_outputs,
    all_edge_table = all_edge_df
  ))
}