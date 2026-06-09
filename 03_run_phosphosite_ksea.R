# ============================================================
# 03. Rank-based KSEA using verified PhosR native output structure
#
# Input:
#   de_results : de$all_results from 02_prepare_metadata_and_phosphosite_de.R
#   site_mat2  : prep$site_mat2 from 01_build_phosphosite_matrix_from_diann.R
#
# External biological annotation:
#   PhosR::PhosphoSitePlus
#     species = "mouse" -> PhosphoSite.mouse
#     species = "human" -> PhosphoSite.human
#
# Verified local PhosR::pathwayRankBasedEnrichment() return structure:
#   class      : matrix / array
#   rownames   : kinase / kinase-signature names
#   colnames   : pvalue, # of substrates, substrates
#
# Important:
#   site_mat2 is a site-level table. Multi-phosphorylated peptides can appear
#   more than once because one modified peptide maps to several phosphosites.
#   For KSEA, keep_only_single_site = TRUE filters these multi-site rows before
#   feature_id uniqueness is checked.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

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

.read_table_if_path <- function(x) {
  if (is.character(x) && length(x) == 1 && file.exists(x)) {
    utils::read.csv(x, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    x
  }
}

# ------------------------------------------------------------
# 1. Load PhosR kinase-substrate annotation
# ------------------------------------------------------------
load_phosr_kinase_annotation <- function(species = c("mouse", "human")) {
  species <- match.arg(species)

  if (!requireNamespace("PhosR", quietly = TRUE)) {
    stop("需要 PhosR 包。请先安装/加载 PhosR。")
  }

  env <- new.env(parent = emptyenv())
  utils::data("PhosphoSitePlus", package = "PhosR", envir = env)

  obj_name <- if (species == "mouse") "PhosphoSite.mouse" else "PhosphoSite.human"
  if (!exists(obj_name, envir = env, inherits = FALSE)) {
    stop("PhosR::PhosphoSitePlus 中找不到对象: ", obj_name)
  }

  get(obj_name, envir = env, inherits = FALSE)
}

normalize_kinase_annotation <- function(annotation,
                                        annotation_kinase_col = NULL,
                                        annotation_site_col = NULL) {
  if (is.list(annotation) && !is.data.frame(annotation)) {
    kinase_names <- names(annotation)
    if (is.null(kinase_names) || any(is.na(kinase_names)) || any(kinase_names == "")) {
      stop("annotation 是 list，但 names(annotation) 不完整；无法明确 kinase/signature 名称。")
    }

    anno_tbl <- tibble::tibble(
      Kinase = rep(kinase_names, lengths(annotation)),
      site_label = as.character(unlist(annotation, use.names = FALSE))
    )
  } else {
    anno_df <- as.data.frame(annotation, check.names = FALSE, stringsAsFactors = FALSE)

    if (is.null(annotation_kinase_col) || is.null(annotation_site_col)) {
      stop(
        "自定义 annotation 是 data.frame 时，必须显式指定 annotation_kinase_col 和 annotation_site_col。\n",
        "当前 annotation 列名：", paste(colnames(anno_df), collapse = ", ")
      )
    }
    if (!annotation_kinase_col %in% colnames(anno_df)) stop("annotation_kinase_col 不在 annotation 中：", annotation_kinase_col)
    if (!annotation_site_col %in% colnames(anno_df)) stop("annotation_site_col 不在 annotation 中：", annotation_site_col)

    anno_tbl <- anno_df %>%
      transmute(
        Kinase = as.character(.data[[annotation_kinase_col]]),
        site_label = as.character(.data[[annotation_site_col]])
      )
  }

  anno_tbl %>%
    mutate(
      Kinase = trimws(as.character(Kinase)),
      site_label = trimws(as.character(site_label))
    ) %>%
    filter(!is.na(Kinase), Kinase != "", !is.na(site_label), site_label != "") %>%
    distinct(Kinase, site_label)
}

# ------------------------------------------------------------
# 2. Build site annotation matching step 02 feature_id
# ------------------------------------------------------------
build_ksea_site_annotation <- function(site_mat2,
                                       feature_id_cols = c("Protein.Id", "Genes", "Residue.Both", "Modified.Sequence"),
                                       keep_only_single_site = TRUE,
                                       uppercase_gene = TRUE,
                                       stop_if_ambiguous_feature_id = TRUE) {
  site_mat2 <- .read_table_if_path(site_mat2)
  site_mat2 <- as.data.frame(site_mat2, check.names = FALSE, stringsAsFactors = FALSE)

  required_cols <- c(feature_id_cols, "abs_pos", "abs_residue")
  missing_cols <- setdiff(required_cols, colnames(site_mat2))
  if (length(missing_cols) > 0) {
    stop("site_mat2 缺少这些列: ", paste(missing_cols, collapse = ", "))
  }

  site_anno <- site_mat2 %>%
    mutate(
      feature_id = do.call(paste, c(across(all_of(feature_id_cols)), sep = "|")),
      Residue = as.character(abs_residue),
      Site = suppressWarnings(as.numeric(abs_pos)),
      GeneSymbol = as.character(Genes),
      n_site_in_row = ifelse(
        "Residue.Both" %in% colnames(.),
        stringr::str_count(as.character(Residue.Both), fixed(";")) + 1L,
        1L
      )
    )

  # Important: multi-site peptides are expected to generate repeated feature_id
  # rows in site_mat2. They must be filtered before feature_id ambiguity checks.
  if (isTRUE(keep_only_single_site)) {
    site_anno <- site_anno %>% filter(n_site_in_row == 1L)
  }

  if (isTRUE(uppercase_gene)) {
    site_anno <- site_anno %>% mutate(GeneSymbol = toupper(GeneSymbol))
  }

  site_anno <- site_anno %>%
    filter(!is.na(Site), !is.na(Residue), Residue %in% c("S", "T", "Y")) %>%
    mutate(site_label = paste0(GeneSymbol, ";", Residue, Site, ";"))

  # Repeated identical mappings are harmless and are collapsed. Only one feature_id
  # mapping to multiple different site_label values is truly ambiguous.
  ambiguity_tbl <- site_anno %>%
    group_by(feature_id) %>%
    summarise(
      n_rows = n(),
      n_site_label = n_distinct(site_label),
      site_labels = paste(unique(site_label), collapse = ";"),
      .groups = "drop"
    ) %>%
    filter(n_site_label > 1)

  if (isTRUE(stop_if_ambiguous_feature_id) && nrow(ambiguity_tbl) > 0) {
    stop(
      "KSEA 位点注释中存在同一个 feature_id 对应多个不同 site_label 的情况。请检查上游映射。示例: ",
      paste(utils::head(ambiguity_tbl$feature_id, 10), collapse = "; ")
    )
  }

  site_anno %>%
    arrange(feature_id, site_label) %>%
    distinct(feature_id, site_label, .keep_all = TRUE) %>%
    distinct(feature_id, .keep_all = TRUE)
}

# ------------------------------------------------------------
# 3. Parse verified PhosR pathwayRankBasedEnrichment() output
# ------------------------------------------------------------
parse_phosr_rank_enrichment <- function(enr,
                                        comparison,
                                        direction,
                                        drop_source_prefixed_signatures = TRUE,
                                        source_signature_prefixes = c("Yang", "Humphrey")) {
  enr_df <- as.data.frame(enr, check.names = FALSE, stringsAsFactors = FALSE)

  required_cols <- c("pvalue", "# of substrates", "substrates")
  missing_cols <- setdiff(required_cols, colnames(enr_df))
  if (length(missing_cols) > 0) {
    stop(
      "PhosR::pathwayRankBasedEnrichment() 返回结构与已验证结构不一致，缺少列: ",
      paste(missing_cols, collapse = ", "),
      "。当前列: ", paste(colnames(enr_df), collapse = ", ")
    )
  }

  if (is.null(rownames(enr_df)) || any(rownames(enr_df) == "")) {
    stop("PhosR::pathwayRankBasedEnrichment() 返回结果没有完整 rownames，无法确定 kinase/signature 名称。")
  }

  out <- enr_df %>%
    tibble::rownames_to_column("Kinase") %>%
    transmute(
      comparison = comparison,
      direction = direction,
      Kinase = as.character(Kinase),
      pvalue = suppressWarnings(as.numeric(pvalue)),
      n_substrates = suppressWarnings(as.integer(.data[["# of substrates"]])),
      substrates = as.character(substrates)
    )

  out$source_prefixed_signature <- FALSE
  if (length(source_signature_prefixes) > 0) {
    prefix_pattern <- paste0("^(?:", paste(source_signature_prefixes, collapse = "|"), ")\\.")
    out$source_prefixed_signature <- grepl(prefix_pattern, out$Kinase)
  }

  if (isTRUE(drop_source_prefixed_signatures)) {
    out <- out %>% filter(!source_prefixed_signature)
  }

  out %>% arrange(pvalue)
}

make_signed_score_from_phosr <- function(enr_up_tbl, enr_down_tbl, comparison) {
  up2 <- enr_up_tbl %>%
    select(Kinase, p_up = pvalue, n_substrates_up = n_substrates, substrates_up = substrates, source_prefixed_signature)

  down2 <- enr_down_tbl %>%
    select(Kinase, p_down = pvalue, n_substrates_down = n_substrates, substrates_down = substrates, source_prefixed_signature)

  full_join(up2, down2, by = c("Kinase", "source_prefixed_signature")) %>%
    mutate(
      comparison = comparison,
      p_best = pmin(p_up, p_down, na.rm = TRUE),
      p_best = ifelse(is.infinite(p_best), NA_real_, p_best),
      dominant_direction = case_when(
        is.na(p_up) & is.na(p_down) ~ NA_character_,
        is.na(p_down) ~ "UP_or_greater",
        is.na(p_up) ~ "DOWN_or_less",
        p_up <= p_down ~ "UP_or_greater",
        TRUE ~ "DOWN_or_less"
      ),
      signed_score = case_when(
        is.na(p_best) ~ NA_real_,
        dominant_direction == "UP_or_greater" ~ -log10(pmax(p_best, .Machine$double.xmin)),
        dominant_direction == "DOWN_or_less" ~ log10(pmax(p_best, .Machine$double.xmin)),
        TRUE ~ NA_real_
      ),
      n_substrates = ifelse(dominant_direction == "UP_or_greater", n_substrates_up, n_substrates_down),
      substrates = ifelse(dominant_direction == "UP_or_greater", substrates_up, substrates_down)
    ) %>%
    select(comparison, Kinase, signed_score, dominant_direction, p_best, p_up, p_down,
           n_substrates, substrates, source_prefixed_signature, everything()) %>%
    arrange(desc(abs(signed_score)))
}

# ------------------------------------------------------------
# 4. Plot per-comparison KSEA results
# ------------------------------------------------------------
plot_one_ksea_signed_score <- function(signed_score,
                                       comparison,
                                       comp_dir,
                                       top_n = 30,
                                       make_heatmap = TRUE,
                                       make_barplot = TRUE) {
  saved_files <- character(0)
  if (is.null(signed_score) || nrow(signed_score) == 0) return(saved_files)

  plot_df <- signed_score %>%
    filter(!is.na(signed_score), !is.na(Kinase), Kinase != "") %>%
    arrange(desc(abs(signed_score))) %>%
    slice_head(n = top_n) %>%
    mutate(
      Kinase = as.character(Kinase),
      Kinase_plot = factor(Kinase, levels = rev(Kinase)),
      activity_direction = ifelse(signed_score >= 0, "Higher", "Lower")
    )

  if (nrow(plot_df) == 0) return(saved_files)
  comp_slug <- .safe_slug(comparison)

  if (isTRUE(make_barplot) && requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Kinase_plot, y = signed_score, fill = activity_direction)) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::coord_flip() +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
      ggplot2::labs(
        x = NULL,
        y = "PhosR rank-based KSEA signed score",
        fill = "Direction",
        title = comparison
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        axis.text.y = ggplot2::element_text(size = 7)
      )

    f_pdf <- file.path(comp_dir, paste0(comp_slug, "_signed_score_barplot.pdf"))
    f_png <- file.path(comp_dir, paste0(comp_slug, "_signed_score_barplot.png"))
    ggplot2::ggsave(f_pdf, p, width = 7, height = max(5, 0.18 * nrow(plot_df) + 2), dpi = 300)
    ggplot2::ggsave(f_png, p, width = 7, height = max(5, 0.18 * nrow(plot_df) + 2), dpi = 300)
    saved_files <- c(saved_files, f_pdf, f_png)
  }

  if (isTRUE(make_heatmap) && requireNamespace("pheatmap", quietly = TRUE)) {
    mat <- as.matrix(plot_df$signed_score)
    rownames(mat) <- plot_df$Kinase
    colnames(mat) <- comparison

    f_pdf <- file.path(comp_dir, paste0(comp_slug, "_signed_score_heatmap.pdf"))
    f_png <- file.path(comp_dir, paste0(comp_slug, "_signed_score_heatmap.png"))

    grDevices::pdf(f_pdf, width = 4.8, height = max(5, 0.16 * nrow(mat) + 2))
    pheatmap::pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_row = 7, fontsize_col = 8, main = comparison)
    grDevices::dev.off()

    grDevices::png(f_png, width = 4.8, height = max(5, 0.16 * nrow(mat) + 2), units = "in", res = 600)
    pheatmap::pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_row = 7, fontsize_col = 8, main = comparison)
    grDevices::dev.off()

    saved_files <- c(saved_files, f_pdf, f_png)
  }

  saved_files
}

# ------------------------------------------------------------
# 5. Run one comparison with PhosR native function
# ------------------------------------------------------------
run_one_phosr_native_ksea <- function(stats_tbl,
                                      comparison,
                                      annotation,
                                      out_dir,
                                      min_sites = 10,
                                      drop_source_prefixed_signatures = TRUE,
                                      source_signature_prefixes = c("Yang", "Humphrey"),
                                      make_per_comparison_heatmap = TRUE,
                                      make_per_comparison_barplot = TRUE,
                                      top_n_per_comparison = 30) {
  comp_slug <- .safe_slug(comparison)
  comp_dir <- file.path(out_dir, comp_slug)
  .safe_dir_create(comp_dir)

  stats_tbl <- stats_tbl %>%
    filter(!is.na(logFC), !is.na(site_label), site_label != "") %>%
    group_by(site_label) %>%
    summarise(logFC = median(logFC, na.rm = TRUE), n_features_collapsed = n(), .groups = "drop") %>%
    arrange(desc(logFC))

  utils::write.csv(stats_tbl, file.path(comp_dir, paste0(comp_slug, "_KSEA_input_logFC_vector.csv")), row.names = FALSE)

  if (nrow(stats_tbl) < min_sites) {
    warning("comparison ", comparison, " 可用于 KSEA 的位点数少于 min_sites=", min_sites, "，跳过。")
    return(list(comparison = comparison, input_stats = stats_tbl, enrichment_up = data.frame(), enrichment_down = data.frame(), signed_score = data.frame(), plot_files = character(0)))
  }

  geneStats <- stats_tbl$logFC
  names(geneStats) <- stats_tbl$site_label
  geneStats <- geneStats[!is.na(geneStats)]

  enr_up <- PhosR::pathwayRankBasedEnrichment(geneStats, annotation = annotation, alter = "greater")
  enr_down <- PhosR::pathwayRankBasedEnrichment(geneStats, annotation = annotation, alter = "less")

  enr_up_tbl <- parse_phosr_rank_enrichment(
    enr = enr_up,
    comparison = comparison,
    direction = "UP_or_greater",
    drop_source_prefixed_signatures = drop_source_prefixed_signatures,
    source_signature_prefixes = source_signature_prefixes
  )

  enr_down_tbl <- parse_phosr_rank_enrichment(
    enr = enr_down,
    comparison = comparison,
    direction = "DOWN_or_less",
    drop_source_prefixed_signatures = drop_source_prefixed_signatures,
    source_signature_prefixes = source_signature_prefixes
  )

  signed_score <- make_signed_score_from_phosr(enr_up_tbl, enr_down_tbl, comparison = comparison)

  utils::write.csv(enr_up_tbl, file.path(comp_dir, paste0(comp_slug, "_rank_KSEA_UP_greater.csv")), row.names = FALSE)
  utils::write.csv(enr_down_tbl, file.path(comp_dir, paste0(comp_slug, "_rank_KSEA_DOWN_less.csv")), row.names = FALSE)
  utils::write.csv(signed_score, file.path(comp_dir, paste0(comp_slug, "_rank_KSEA_signed_score.csv")), row.names = FALSE)

  plot_files <- plot_one_ksea_signed_score(
    signed_score = signed_score,
    comparison = comparison,
    comp_dir = comp_dir,
    top_n = top_n_per_comparison,
    make_heatmap = make_per_comparison_heatmap,
    make_barplot = make_per_comparison_barplot
  )

  list(
    comparison = comparison,
    input_stats = stats_tbl,
    enrichment_up = enr_up_tbl,
    enrichment_down = enr_down_tbl,
    signed_score = signed_score,
    plot_files = plot_files
  )
}

# ------------------------------------------------------------
# 6. Main wrapper
# ------------------------------------------------------------
run_phosphosite_ksea <- function(de_results,
                                 site_mat2,
                                 out_dir = "Results/phosphosite_KSEA",
                                 species = c("mouse", "human"),
                                 annotation = NULL,
                                 annotation_kinase_col = NULL,
                                 annotation_site_col = NULL,
                                 comparison_col = "comparison",
                                 logfc_col = "logFC",
                                 keep_only_single_site = TRUE,
                                 uppercase_gene = TRUE,
                                 min_sites = 10,
                                 make_heatmap = TRUE,
                                 top_n_kinases_heatmap = 80,
                                 make_per_comparison_heatmap = TRUE,
                                 make_per_comparison_barplot = TRUE,
                                 top_n_per_comparison = 30,
                                 drop_source_prefixed_signatures = TRUE,
                                 source_signature_prefixes = c("Yang", "Humphrey")) {
  species <- match.arg(species)
  .safe_dir_create(out_dir)

  if (!requireNamespace("PhosR", quietly = TRUE)) stop("需要 PhosR 包。")

  de_results <- .read_table_if_path(de_results)
  de_results <- as.data.frame(de_results, check.names = FALSE, stringsAsFactors = FALSE)

  if (!comparison_col %in% colnames(de_results)) stop("de_results 缺少 comparison 列: ", comparison_col)
  if (!logfc_col %in% colnames(de_results)) stop("de_results 缺少 logFC 列: ", logfc_col)
  if (!"feature_id" %in% colnames(de_results)) stop("de_results 缺少 feature_id 列。请使用第二步脚本输出的 all_results。")

  site_anno <- build_ksea_site_annotation(
    site_mat2 = site_mat2,
    keep_only_single_site = keep_only_single_site,
    uppercase_gene = uppercase_gene
  )
  utils::write.csv(site_anno, file.path(out_dir, "01_ksea_site_annotation.csv"), row.names = FALSE)

  if (is.null(annotation)) annotation <- load_phosr_kinase_annotation(species = species)

  annotation_tbl <- normalize_kinase_annotation(
    annotation = annotation,
    annotation_kinase_col = annotation_kinase_col,
    annotation_site_col = annotation_site_col
  )
  utils::write.csv(annotation_tbl, file.path(out_dir, "00_ksea_kinase_substrate_annotation_used.csv"), row.names = FALSE)

  kinase_name_check <- annotation_tbl %>%
    count(Kinase, name = "n_substrate_sites") %>%
    mutate(
      source_prefixed_signature = ifelse(
        length(source_signature_prefixes) > 0,
        grepl(paste0("^(?:", paste(source_signature_prefixes, collapse = "|"), ")\\."), Kinase),
        FALSE
      )
    ) %>%
    arrange(Kinase)
  utils::write.csv(kinase_name_check, file.path(out_dir, "00_ksea_kinase_names_used.csv"), row.names = FALSE)

  de2 <- de_results %>%
    mutate(logFC = suppressWarnings(as.numeric(.data[[logfc_col]]))) %>%
    left_join(site_anno %>% select(feature_id, site_label, GeneSymbol, Residue, Site), by = "feature_id")
  utils::write.csv(de2, file.path(out_dir, "02_ksea_de_results_with_site_label.csv"), row.names = FALSE)

  match_summary <- data.frame(
    total_de_rows = nrow(de2),
    rows_with_site_label = sum(!is.na(de2$site_label) & de2$site_label != ""),
    rows_without_site_label = sum(is.na(de2$site_label) | de2$site_label == ""),
    unique_site_labels = length(unique(de2$site_label[!is.na(de2$site_label) & de2$site_label != ""])),
    stringsAsFactors = FALSE
  )
  utils::write.csv(match_summary, file.path(out_dir, "02_ksea_site_label_match_summary.csv"), row.names = FALSE)

  comparisons <- unique(as.character(de2[[comparison_col]]))
  comparisons <- comparisons[!is.na(comparisons) & comparisons != ""]
  if (length(comparisons) == 0) stop("de_results 中没有可用 comparison。")

  result_list <- list()
  for (cn in comparisons) {
    stats_tbl <- de2 %>%
      filter(.data[[comparison_col]] == cn) %>%
      select(feature_id, site_label, logFC, GeneSymbol, Residue, Site)

    result_list[[cn]] <- run_one_phosr_native_ksea(
      stats_tbl = stats_tbl,
      comparison = cn,
      annotation = annotation,
      out_dir = out_dir,
      min_sites = min_sites,
      drop_source_prefixed_signatures = drop_source_prefixed_signatures,
      source_signature_prefixes = source_signature_prefixes,
      make_per_comparison_heatmap = make_per_comparison_heatmap,
      make_per_comparison_barplot = make_per_comparison_barplot,
      top_n_per_comparison = top_n_per_comparison
    )
  }

  signed_long <- bind_rows(lapply(result_list, `[[`, "signed_score"))
  utils::write.csv(signed_long, file.path(out_dir, "all_rank_based_KSEA_signed_score_long.csv"), row.names = FALSE)

  plot_files <- unlist(lapply(result_list, `[[`, "plot_files"), use.names = FALSE)
  utils::write.csv(data.frame(plot_file = unique(plot_files), stringsAsFactors = FALSE), file.path(out_dir, "per_comparison_plot_files.csv"), row.names = FALSE)

  if (nrow(signed_long) > 0) {
    signed_wide <- signed_long %>%
      select(Kinase, comparison, signed_score) %>%
      distinct(Kinase, comparison, .keep_all = TRUE) %>%
      pivot_wider(names_from = comparison, values_from = signed_score)

    utils::write.csv(signed_wide, file.path(out_dir, "all_rank_based_KSEA_signed_score_wide.csv"), row.names = FALSE)

    if (isTRUE(make_heatmap) && requireNamespace("pheatmap", quietly = TRUE)) {
      mat <- signed_wide %>% tibble::column_to_rownames("Kinase") %>% as.matrix()
      mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]

      if (nrow(mat) > 1 && ncol(mat) > 0) {
        max_abs <- apply(abs(mat), 1, max, na.rm = TRUE)
        max_abs[!is.finite(max_abs)] <- NA_real_
        keep <- order(max_abs, decreasing = TRUE, na.last = NA)
        keep <- keep[seq_len(min(length(keep), top_n_kinases_heatmap))]
        mat_plot <- mat[keep, , drop = FALSE]

        pdf(file.path(out_dir, "all_rank_based_KSEA_signed_score_heatmap.pdf"), width = 10, height = 10)
        pheatmap::pheatmap(mat_plot, fontsize_row = 6, fontsize_col = 8, main = "PhosR rank-based KSEA signed score")
        grDevices::dev.off()

        png(file.path(out_dir, "all_rank_based_KSEA_signed_score_heatmap.png"), width = 10, height = 10, units = "in", res = 600)
        pheatmap::pheatmap(mat_plot, fontsize_row = 6, fontsize_col = 8, main = "PhosR rank-based KSEA signed score")
        grDevices::dev.off()
      }
    }
  } else {
    signed_wide <- data.frame()
  }

  saveRDS(result_list, file.path(out_dir, "rank_based_KSEA_result_list.rds"))

  invisible(list(
    site_annotation = site_anno,
    annotation_tbl = annotation_tbl,
    kinase_name_check = kinase_name_check,
    de_results_with_site_label = de2,
    match_summary = match_summary,
    result_list = result_list,
    signed_long = signed_long,
    signed_wide = signed_wide,
    per_comparison_plot_files = unique(plot_files),
    annotation = annotation
  ))
}

# ============================================================
# Example usage
# ============================================================
# source("03_run_phosphosite_ksea.R")
#
# ksea <- run_phosphosite_ksea(
#   de_results = de$all_results,
#   site_mat2 = prep$site_mat2,
#   out_dir = "./demo/Results/phosphosite_KSEA",
#   species = "mouse",
#   keep_only_single_site = TRUE,
#   min_sites = 10,
#   drop_source_prefixed_signatures = TRUE,
#   source_signature_prefixes = c("Yang", "Humphrey")
# )
