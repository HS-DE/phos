# ============================================================
# Debug script: inspect PhosR::pathwayRankBasedEnrichment() return structure
#
# Purpose:
#   This script does NOT perform the final KSEA workflow.
#   It only helps inspect the real output structure of
#   PhosR::pathwayRankBasedEnrichment() in your local R/PhosR version.
#
# Why this matters:
#   Different package versions or annotation objects may return objects with
#   different column names / row names. We should inspect the actual returned
#   structure before writing a wrapper that extracts Kinase / P value / score.
#
# Required objects before running:
#   prep$site_mat2
#   de$all_results
#
# Typical usage:
#   source("03_debug_phosr_pathwayRankBasedEnrichment_structure.R")
#
#   dbg <- debug_phosr_pathwayRankBasedEnrichment_structure(
#     de_results = de$all_results,
#     site_mat2 = prep$site_mat2,
#     out_dir = "./demo/Results/debug_PhosR_pathwayRankBasedEnrichment",
#     species = "mouse"
#   )
#
# Then send these files / console outputs for wrapper development:
#   debug_summary.txt
#   enrichment_up_head.csv
#   enrichment_down_head.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
})

.debug_read_table_if_path <- function(x) {
  if (is.character(x) && length(x) == 1 && file.exists(x)) {
    utils::read.csv(x, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    x
  }
}

.debug_build_site_annotation <- function(site_mat2,
                                         feature_id_cols = c("Protein.Id", "Genes", "Residue.Both", "Modified.Sequence"),
                                         keep_only_single_site = TRUE,
                                         uppercase_gene = TRUE) {
  site_mat2 <- .debug_read_table_if_path(site_mat2)
  site_mat2 <- as.data.frame(site_mat2, check.names = FALSE, stringsAsFactors = FALSE)

  required_cols <- c(feature_id_cols, "abs_pos", "abs_residue")
  missing_cols <- setdiff(required_cols, colnames(site_mat2))
  if (length(missing_cols) > 0) {
    stop("site_mat2 缺少这些列: ", paste(missing_cols, collapse = ", "))
  }

  site_anno <- site_mat2 %>%
    mutate(
      feature_id_raw = do.call(paste, c(across(all_of(feature_id_cols)), sep = "|")),
      feature_id = make.unique(feature_id_raw),
      Residue = as.character(abs_residue),
      Site = suppressWarnings(as.numeric(abs_pos)),
      GeneSymbol = as.character(Genes),
      n_site_in_row = ifelse(
        "Residue.Both" %in% colnames(.),
        stringr::str_count(as.character(Residue.Both), fixed(";")) + 1L,
        1L
      )
    )

  if (uppercase_gene) {
    site_anno <- site_anno %>% mutate(GeneSymbol = toupper(GeneSymbol))
  }

  if (keep_only_single_site) {
    site_anno <- site_anno %>% filter(n_site_in_row == 1L)
  }

  site_anno %>%
    filter(!is.na(Site), !is.na(Residue), Residue %in% c("S", "T", "Y")) %>%
    mutate(site_label = paste0(GeneSymbol, ";", Residue, Site, ";")) %>%
    distinct(feature_id, .keep_all = TRUE)
}

debug_phosr_pathwayRankBasedEnrichment_structure <- function(de_results,
                                                             site_mat2,
                                                             out_dir = "Results/debug_PhosR_pathwayRankBasedEnrichment",
                                                             species = c("mouse", "human"),
                                                             comparison = NULL,
                                                             comparison_col = "comparison",
                                                             logfc_col = "logFC",
                                                             keep_only_single_site = TRUE,
                                                             max_sites = Inf) {
  species <- match.arg(species)

  if (!requireNamespace("PhosR", quietly = TRUE)) {
    stop("需要 PhosR 包。")
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  de_results <- .debug_read_table_if_path(de_results)
  de_results <- as.data.frame(de_results, check.names = FALSE, stringsAsFactors = FALSE)

  if (!comparison_col %in% colnames(de_results)) stop("de_results 缺少 comparison 列: ", comparison_col)
  if (!logfc_col %in% colnames(de_results)) stop("de_results 缺少 logFC 列: ", logfc_col)
  if (!"feature_id" %in% colnames(de_results)) stop("de_results 缺少 feature_id 列。")

  site_anno <- .debug_build_site_annotation(
    site_mat2 = site_mat2,
    keep_only_single_site = keep_only_single_site
  )

  de2 <- de_results %>%
    mutate(logFC = suppressWarnings(as.numeric(.data[[logfc_col]]))) %>%
    left_join(site_anno %>% select(feature_id, site_label), by = "feature_id")

  comparisons <- unique(as.character(de2[[comparison_col]]))
  comparisons <- comparisons[!is.na(comparisons) & comparisons != ""]

  if (length(comparisons) == 0) stop("没有可用 comparison。")

  if (is.null(comparison)) {
    comparison <- comparisons[1]
  }

  if (!comparison %in% comparisons) {
    stop("指定 comparison 不存在: ", comparison, "。可用 comparison: ", paste(comparisons, collapse = ", "))
  }

  stats_tbl <- de2 %>%
    filter(.data[[comparison_col]] == comparison) %>%
    filter(!is.na(logFC), !is.na(site_label), site_label != "") %>%
    group_by(site_label) %>%
    summarise(logFC = median(logFC, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(logFC))

  if (is.finite(max_sites) && nrow(stats_tbl) > max_sites) {
    stats_tbl <- stats_tbl %>% slice_head(n = max_sites)
  }

  geneStats <- stats_tbl$logFC
  names(geneStats) <- stats_tbl$site_label
  geneStats <- geneStats[!is.na(geneStats)]

  utils::write.csv(stats_tbl, file.path(out_dir, "geneStats_input_vector.csv"), row.names = FALSE)

  env <- new.env(parent = emptyenv())
  utils::data("PhosphoSitePlus", package = "PhosR", envir = env)
  annotation_name <- if (species == "mouse") "PhosphoSite.mouse" else "PhosphoSite.human"
  annotation <- get(annotation_name, envir = env, inherits = FALSE)

  enr_up <- PhosR::pathwayRankBasedEnrichment(
    geneStats,
    annotation = annotation,
    alter = "greater"
  )

  enr_down <- PhosR::pathwayRankBasedEnrichment(
    geneStats,
    annotation = annotation,
    alter = "less"
  )

  enr_up_df <- as.data.frame(enr_up, check.names = FALSE, stringsAsFactors = FALSE)
  enr_down_df <- as.data.frame(enr_down, check.names = FALSE, stringsAsFactors = FALSE)

  utils::write.csv(utils::head(enr_up_df, 50), file.path(out_dir, "enrichment_up_head.csv"), row.names = TRUE)
  utils::write.csv(utils::head(enr_down_df, 50), file.path(out_dir, "enrichment_down_head.csv"), row.names = TRUE)

  summary_file <- file.path(out_dir, "debug_summary.txt")
  zz <- file(summary_file, open = "wt", encoding = "UTF-8")
  on.exit(close(zz), add = TRUE)

  cat("# PhosR::pathwayRankBasedEnrichment() debug summary\n", file = zz)
  cat("\n## sessionInfo()\n", file = zz)
  capture.output(utils::sessionInfo(), file = zz, append = TRUE)

  cat("\n## selected comparison\n", file = zz)
  cat(comparison, "\n", file = zz)

  cat("\n## geneStats\n", file = zz)
  cat("length: ", length(geneStats), "\n", file = zz)
  cat("head names:\n", file = zz)
  capture.output(utils::head(names(geneStats), 20), file = zz, append = TRUE)
  cat("head values:\n", file = zz)
  capture.output(utils::head(geneStats, 20), file = zz, append = TRUE)

  cat("\n## annotation object\n", file = zz)
  cat("class: ", paste(class(annotation), collapse = ", "), "\n", file = zz)
  cat("length: ", length(annotation), "\n", file = zz)
  cat("head(names(annotation)):\n", file = zz)
  capture.output(utils::head(names(annotation), 30), file = zz, append = TRUE)

  cat("\n## enrichment_up object\n", file = zz)
  cat("class: ", paste(class(enr_up), collapse = ", "), "\n", file = zz)
  cat("dim as data.frame: ", paste(dim(enr_up_df), collapse = " x "), "\n", file = zz)
  cat("colnames:\n", file = zz)
  capture.output(colnames(enr_up_df), file = zz, append = TRUE)
  cat("rownames head:\n", file = zz)
  capture.output(utils::head(rownames(enr_up_df), 30), file = zz, append = TRUE)
  cat("str(enr_up):\n", file = zz)
  capture.output(str(enr_up), file = zz, append = TRUE)

  cat("\n## enrichment_down object\n", file = zz)
  cat("class: ", paste(class(enr_down), collapse = ", "), "\n", file = zz)
  cat("dim as data.frame: ", paste(dim(enr_down_df), collapse = " x "), "\n", file = zz)
  cat("colnames:\n", file = zz)
  capture.output(colnames(enr_down_df), file = zz, append = TRUE)
  cat("rownames head:\n", file = zz)
  capture.output(utils::head(rownames(enr_down_df), 30), file = zz, append = TRUE)
  cat("str(enr_down):\n", file = zz)
  capture.output(str(enr_down), file = zz, append = TRUE)

  message("Debug files written to: ", out_dir)

  invisible(list(
    comparison = comparison,
    geneStats = geneStats,
    annotation = annotation,
    enrichment_up = enr_up,
    enrichment_down = enr_down,
    enrichment_up_df = enr_up_df,
    enrichment_down_df = enr_down_df,
    debug_summary_file = summary_file,
    out_dir = out_dir
  ))
}
