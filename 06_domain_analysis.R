# ============================================================
# 06. Protein domain analysis and domain enrichment
#
# Purpose:
#   Reproduce company-style domain analysis:
#
#   1) Domain Analysis:
#      Differential / trend modified proteins -> domain count Top20 bar plot
#
#   2) Domain Enrichment:
#      Fisher exact test using all identified modified proteins as background
#
# Main input:
#   de$all_results
#   domain_annotation_from_uniprot.csv
#
# Important:
#   This is protein-level domain analysis.
#   It counts proteins carrying each domain, not phosphosite counts.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(tibble)
})

# ------------------------------------------------------------
# Helper functions
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

.write_csv_gbk <- function(x, file) {
  .safe_dir_create(dirname(file))
  write.csv(x, file, row.names = FALSE, fileEncoding = "GBK")
}

read_domain_annotation <- function(domain_annotation_file) {
  if (!file.exists(domain_annotation_file)) {
    stop("Domain annotation file not found: ", domain_annotation_file)
  }
  
  ext <- tolower(tools::file_ext(domain_annotation_file))
  
  if (ext == "csv") {
    anno <- read.csv(domain_annotation_file, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    anno <- read.delim(domain_annotation_file, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    stop("Only csv / tsv / txt domain annotation files are supported.")
  }
  
  anno
}

standardize_domain_annotation <- function(domain_annotation) {
  domain_annotation <- as.data.frame(domain_annotation, check.names = FALSE)
  
  if (!"Protein.Id" %in% colnames(domain_annotation)) {
    stop("domain annotation must contain column: Protein.Id")
  }
  
  if (!"Domain.Name" %in% colnames(domain_annotation)) {
    stop("domain annotation must contain column: Domain.Name")
  }
  
  if (!"Genes" %in% colnames(domain_annotation)) {
    domain_annotation$Genes <- NA_character_
  }
  
  if (!"Domain.ID" %in% colnames(domain_annotation)) {
    domain_annotation$Domain.ID <- NA_character_
  }
  
  if (!"Domain.Source" %in% colnames(domain_annotation)) {
    domain_annotation$Domain.Source <- NA_character_
  }
  
  if (!"Domain.Start" %in% colnames(domain_annotation)) {
    domain_annotation$Domain.Start <- NA_integer_
  }
  
  if (!"Domain.End" %in% colnames(domain_annotation)) {
    domain_annotation$Domain.End <- NA_integer_
  }
  
  domain_annotation %>%
    transmute(
      Protein.Id = trimws(as.character(Protein.Id)),
      Genes = trimws(as.character(Genes)),
      Domain.ID = as.character(Domain.ID),
      Domain.Name = trimws(as.character(Domain.Name)),
      Domain.Source = as.character(Domain.Source),
      Domain.Start = suppressWarnings(as.integer(Domain.Start)),
      Domain.End = suppressWarnings(as.integer(Domain.End))
    ) %>%
    filter(!is.na(Protein.Id), Protein.Id != "") %>%
    filter(!is.na(Domain.Name), Domain.Name != "") %>%
    distinct(
      Protein.Id,
      Domain.Name,
      Domain.ID,
      Domain.Source,
      Domain.Start,
      Domain.End,
      .keep_all = TRUE
    )
}

# ------------------------------------------------------------
# 1. Define background proteins from de$all_results
# ------------------------------------------------------------

make_background_modified_protein_table <- function(de_results) {
  required_cols <- c(
    "feature_id",
    "Protein.Id",
    "Genes",
    "Modified.Sequence",
    "Residue.Both",
    "logFC",
    "comparison",
    "direction"
  )
  
  missing_cols <- setdiff(required_cols, colnames(de_results))
  if (length(missing_cols) > 0) {
    stop("de_results missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  de_results %>%
    as.data.frame(check.names = FALSE) %>%
    as_tibble() %>%
    mutate(
      Protein.Id = as.character(Protein.Id),
      Genes = as.character(Genes),
      comparison = as.character(comparison),
      direction = as.character(direction),
      logFC = suppressWarnings(as.numeric(logFC))
    ) %>%
    group_by(comparison, Protein.Id, Genes) %>%
    summarise(
      n_all_modified_features = n_distinct(feature_id),
      n_all_modified_sequences = n_distinct(Modified.Sequence),
      all_residues = paste(sort(unique(Residue.Both)), collapse = ";"),
      .groups = "drop"
    )
}

# ------------------------------------------------------------
# 2. Select differential / trend modified features
# ------------------------------------------------------------

make_diff_modified_feature_table <- function(de_results,
                                             diff_mode = c("direction", "abs_logFC"),
                                             keep_directions = c("Up", "Down"),
                                             logfc_cutoff = 1) {
  diff_mode <- match.arg(diff_mode)
  
  required_cols <- c(
    "feature_id",
    "Protein.Id",
    "Genes",
    "Modified.Sequence",
    "Residue.Both",
    "logFC",
    "comparison",
    "direction"
  )
  
  missing_cols <- setdiff(required_cols, colnames(de_results))
  if (length(missing_cols) > 0) {
    stop("de_results missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  x <- de_results %>%
    as.data.frame(check.names = FALSE) %>%
    as_tibble() %>%
    mutate(
      Protein.Id = as.character(Protein.Id),
      Genes = as.character(Genes),
      comparison = as.character(comparison),
      direction = as.character(direction),
      logFC = suppressWarnings(as.numeric(logFC))
    )
  
  if (diff_mode == "direction") {
    x <- x %>%
      filter(direction %in% keep_directions)
  }
  
  if (diff_mode == "abs_logFC") {
    x <- x %>%
      filter(!is.na(logFC), abs(logFC) >= logfc_cutoff) %>%
      mutate(
        direction = case_when(
          logFC >= logfc_cutoff ~ "Up",
          logFC <= -logfc_cutoff ~ "Down",
          TRUE ~ "Not Significant"
        )
      ) %>%
      filter(direction %in% keep_directions)
  }
  
  x
}

make_diff_modified_protein_table <- function(diff_feature_table) {
  if (nrow(diff_feature_table) == 0) {
    return(tibble())
  }
  
  diff_feature_table %>%
    group_by(comparison, Protein.Id, Genes) %>%
    summarise(
      n_modified_features = n_distinct(feature_id),
      n_modified_sequences = n_distinct(Modified.Sequence),
      residues = paste(sort(unique(Residue.Both)), collapse = ";"),
      directions = paste(sort(unique(direction)), collapse = ";"),
      max_abs_logFC = suppressWarnings(max(abs(logFC), na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      max_abs_logFC = ifelse(is.infinite(max_abs_logFC), NA_real_, max_abs_logFC)
    )
}

# ------------------------------------------------------------
# 3. Attach domains
# ------------------------------------------------------------

attach_domain_annotation <- function(protein_table,
                                     domain_annotation_std) {
  if (nrow(protein_table) == 0) {
    return(tibble())
  }
  
  ## domain_annotation_std 里也有 Genes，
  ## 为了避免 left_join 后变成 Genes.x / Genes.y，
  ## 这里不使用 annotation 里的 Genes。
  domain_anno_join <- domain_annotation_std %>%
    select(
      Protein.Id,
      Domain.ID,
      Domain.Name,
      Domain.Source,
      Domain.Start,
      Domain.End
    ) %>%
    distinct()
  
  protein_table %>%
    left_join(
      domain_anno_join,
      by = "Protein.Id",
      relationship = "many-to-many"
    )
}

make_domain_count_table <- function(diff_protein_domain) {
  if (nrow(diff_protein_domain) == 0) {
    return(tibble())
  }
  
  diff_protein_domain %>%
    filter(!is.na(Domain.Name), Domain.Name != "") %>%
    distinct(comparison, Protein.Id, Genes, Domain.Name, .keep_all = TRUE) %>%
    group_by(comparison, Domain.Name) %>%
    summarise(
      protein_count = n_distinct(Protein.Id),
      protein_ids = paste(sort(unique(Protein.Id)), collapse = ";"),
      genes = paste(sort(unique(Genes[!is.na(Genes) & Genes != ""])), collapse = ";"),
      domain_ids = paste(sort(unique(Domain.ID[!is.na(Domain.ID) & Domain.ID != ""])), collapse = ";"),
      domain_sources = paste(sort(unique(Domain.Source[!is.na(Domain.Source) & Domain.Source != ""])), collapse = ";"),
      .groups = "drop"
    ) %>%
    group_by(comparison) %>%
    mutate(
      total_domain_protein_count = sum(protein_count),
      protein_percent = protein_count / total_domain_protein_count * 100
    ) %>%
    ungroup() %>%
    arrange(comparison, desc(protein_count), Domain.Name)
}

# ------------------------------------------------------------
# 4. Domain enrichment by Fisher exact test
# ------------------------------------------------------------

run_domain_enrichment_one_comparison <- function(comparison_name,
                                                 bg_protein_table,
                                                 diff_protein_table,
                                                 domain_annotation_std) {
  bg_all <- bg_protein_table %>%
    filter(comparison == comparison_name) %>%
    distinct(Protein.Id, Genes)
  
  diff_all <- diff_protein_table %>%
    filter(comparison == comparison_name) %>%
    distinct(Protein.Id, Genes)
  
  bg_total <- n_distinct(bg_all$Protein.Id)
  diff_total <- n_distinct(diff_all$Protein.Id)
  
  if (bg_total == 0 || diff_total == 0) {
    return(tibble())
  }
  
  ## 关键修复：
  ## domain_annotation_std 里也有 Genes，
  ## 这里先去掉，避免 join 后 Genes 变成 Genes.x / Genes.y。
  domain_anno_join <- domain_annotation_std %>%
    select(
      Protein.Id,
      Domain.ID,
      Domain.Name,
      Domain.Source,
      Domain.Start,
      Domain.End
    ) %>%
    filter(!is.na(Domain.Name), Domain.Name != "") %>%
    distinct()
  
  bg_domain <- bg_all %>%
    left_join(
      domain_anno_join,
      by = "Protein.Id",
      relationship = "many-to-many"
    ) %>%
    filter(!is.na(Domain.Name), Domain.Name != "") %>%
    distinct(Protein.Id, Genes, Domain.Name, .keep_all = TRUE)
  
  diff_domain <- diff_all %>%
    left_join(
      domain_anno_join,
      by = "Protein.Id",
      relationship = "many-to-many"
    ) %>%
    filter(!is.na(Domain.Name), Domain.Name != "") %>%
    distinct(Protein.Id, Genes, Domain.Name, .keep_all = TRUE)
  
  domain_names <- sort(unique(bg_domain$Domain.Name))
  
  if (length(domain_names) == 0) {
    return(tibble())
  }
  
  out <- lapply(domain_names, function(domain_i) {
    bg_domain_proteins <- bg_domain %>%
      filter(Domain.Name == domain_i) %>%
      pull(Protein.Id) %>%
      unique()
    
    diff_domain_proteins <- diff_domain %>%
      filter(Domain.Name == domain_i) %>%
      pull(Protein.Id) %>%
      unique()
    
    a <- length(intersect(diff_all$Protein.Id, diff_domain_proteins))
    bg_domain_count <- length(bg_domain_proteins)
    
    b <- diff_total - a
    c <- bg_domain_count - a
    d <- bg_total - diff_total - c
    
    if (any(c(a, b, c, d) < 0)) {
      p_value <- NA_real_
      odds_ratio <- NA_real_
    } else {
      ft <- fisher.test(
        matrix(c(a, b, c, d), nrow = 2, byrow = TRUE),
        alternative = "greater"
      )
      
      p_value <- ft$p.value
      odds_ratio <- unname(ft$estimate)
    }
    
    domain_rows <- bg_domain %>%
      filter(Domain.Name == domain_i)
    
    diff_rows <- diff_domain %>%
      filter(Domain.Name == domain_i)
    
    tibble(
      comparison = comparison_name,
      Domain.Name = domain_i,
      Domain.ID = paste(
        sort(unique(domain_rows$Domain.ID[!is.na(domain_rows$Domain.ID) & domain_rows$Domain.ID != ""])),
        collapse = ";"
      ),
      Domain.Source = paste(
        sort(unique(domain_rows$Domain.Source[!is.na(domain_rows$Domain.Source) & domain_rows$Domain.Source != ""])),
        collapse = ";"
      ),
      diff_domain_protein_count = a,
      bg_domain_protein_count = bg_domain_count,
      diff_total_protein_count = diff_total,
      bg_total_protein_count = bg_total,
      rich_factor = ifelse(bg_domain_count > 0, a / bg_domain_count, NA_real_),
      p_value = p_value,
      odds_ratio = odds_ratio,
      protein_ids = paste(sort(unique(diff_rows$Protein.Id)), collapse = ";"),
      genes = paste(sort(unique(diff_rows$Genes[!is.na(diff_rows$Genes) & diff_rows$Genes != ""])), collapse = ";")
    )
  })
  
  bind_rows(out) %>%
    filter(diff_domain_protein_count > 0) %>%
    mutate(
      minus_log10_p = -log10(pmax(p_value, .Machine$double.xmin))
    )
}

run_domain_enrichment <- function(bg_protein_table,
                                  diff_protein_table,
                                  domain_annotation_std) {
  comps <- unique(bg_protein_table$comparison)
  
  res <- bind_rows(lapply(
    comps,
    run_domain_enrichment_one_comparison,
    bg_protein_table = bg_protein_table,
    diff_protein_table = diff_protein_table,
    domain_annotation_std = domain_annotation_std
  ))
  
  if (nrow(res) == 0) {
    return(tibble())
  }
  
  res %>%
    group_by(comparison) %>%
    mutate(
      adj_p_value = p.adjust(p_value, method = "BH")
    ) %>%
    ungroup() %>%
    arrange(comparison, p_value, desc(diff_domain_protein_count))
}

# ------------------------------------------------------------
# 5. Plot Domain Analysis Top20 barplot
# ------------------------------------------------------------

plot_domain_analysis_bar <- function(domain_count,
                                     comparison,
                                     out_prefix,
                                     top_n = 20,
                                     save_formats = c("pdf", "png"),
                                     width = 8,
                                     height = 6,
                                     dpi = 300) {
  df <- domain_count %>%
    filter(comparison == !!comparison) %>%
    arrange(desc(protein_count), Domain.Name) %>%
    slice_head(n = top_n)
  
  if (nrow(df) == 0) {
    warning("No domain count data to plot for comparison: ", comparison)
    return(invisible(NULL))
  }
  
  df <- df %>%
    mutate(
      Domain.Name = factor(Domain.Name, levels = rev(Domain.Name))
    )
  
  p <- ggplot(df, aes(x = protein_count, y = Domain.Name)) +
    geom_col(width = 0.75, fill = "#4B3B91") +
    geom_text(
      aes(label = protein_count),
      hjust = -0.2,
      size = 3.5
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = "Domain Analysis",
      x = "The number of Proteins",
      y = paste0("Domain Name(Top ", top_n, ")")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  for (fmt in tolower(save_formats)) {
    outfile <- paste0(out_prefix, ".", fmt)
    ggsave(outfile, p, width = width, height = height, dpi = dpi)
  }
  
  invisible(p)
}

# ------------------------------------------------------------
# 6. Plot Domain Enrichment bubble plot
# ------------------------------------------------------------

plot_domain_enrichment_bubble <- function(enrich_df,
                                          comparison,
                                          out_prefix,
                                          top_n = 20,
                                          save_formats = c("pdf", "png"),
                                          width = 8,
                                          height = 8,
                                          dpi = 600,
                                          shorten_max_len = 45,
                                          description_mode = c("wrap", "truncate")) {
  description_mode <- match.arg(description_mode)
  
  wrap_text_by_width <- function(x, width) {
    x <- as.character(x)
    width_num <- suppressWarnings(as.numeric(width)[1])
    
    if (is.na(width_num) || !is.finite(width_num) || width_num <= 0) {
      return(x)
    }
    
    width_num <- as.integer(width_num)
    
    vapply(
      x,
      function(s) {
        if (is.na(s) || s == "") return(s)
        paste(base::strwrap(s, width = width_num), collapse = "\n")
      },
      FUN.VALUE = character(1),
      USE.NAMES = FALSE
    )
  }
  
  truncate_text_by_width <- function(x, width) {
    x <- as.character(x)
    width_num <- suppressWarnings(as.numeric(width)[1])
    
    if (is.na(width_num) || !is.finite(width_num) || width_num <= 0) {
      return(x)
    }
    
    width_num <- as.integer(width_num)
    
    vapply(
      x,
      function(s) {
        if (is.na(s) || s == "") return(s)
        if (nchar(s) <= width_num) return(s)
        paste0(substr(s, 1, width_num), "...")
      },
      FUN.VALUE = character(1),
      USE.NAMES = FALSE
    )
  }
  
  format_domain_description <- function(x) {
    if (identical(description_mode, "wrap")) {
      return(wrap_text_by_width(x, shorten_max_len))
    }
    
    truncate_text_by_width(x, shorten_max_len)
  }
  
  df <- enrich_df %>%
    dplyr::filter(comparison == !!comparison) %>%
    dplyr::filter(!is.na(rich_factor)) %>%
    dplyr::filter(!is.na(p_value)) %>%
    dplyr::filter(!is.na(Domain.Name), Domain.Name != "") %>%
    dplyr::filter(diff_domain_protein_count > 0)
  
  if (nrow(df) == 0) {
    warning("No domain enrichment data to plot for comparison: ", comparison)
    return(invisible(NULL))
  }
  
  plot_df <- df %>%
    dplyr::arrange(p_value) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      pvalue_safe = pmax(p_value, .Machine$double.xmin),
      minus_log10_p = -log10(pvalue_safe),
      Count = diff_domain_protein_count,
      RichFactor = rich_factor,
      Description = format_domain_description(Domain.Name)
    )
  
  if (nrow(plot_df) == 0) {
    warning("No top domain enrichment data to plot for comparison: ", comparison)
    return(invisible(NULL))
  }
  
  ## 这里和你的 KEGG 脚本保持一致：
  ## x 轴 RichFactor，y 轴按 RichFactor 降序
  plot_df$Description <- factor(
    plot_df$Description,
    levels = unique(plot_df$Description)
  )
  
  min_val <- min(plot_df$minus_log10_p, na.rm = TRUE)
  max_val <- max(plot_df$minus_log10_p, na.rm = TRUE)
  
  if (!is.finite(min_val) || !is.finite(max_val)) {
    min_val <- 0
    max_val <- 1
  }
  
  if (identical(min_val, max_val)) {
    min_val <- min_val - 0.1
    max_val <- max_val + 0.1
  }
  
  min_count <- min(plot_df$Count, na.rm = TRUE)
  med_count <- round(stats::median(plot_df$Count, na.rm = TRUE))
  max_count <- max(plot_df$Count, na.rm = TRUE)
  
  max_size_val <- max(min(max_count, 11), 6)
  
  vals <- sort(unique(c(min_count, med_count, max_count)))
  
  if (length(vals) == 3) {
    breaks_vals <- c(min_count, med_count, max_count)
    labels_vals <- c(min_count, med_count, max_count)
  } else if (length(vals) == 2) {
    breaks_vals <- c(min_count, max_count)
    labels_vals <- c(min_count, max_count)
  } else {
    breaks_vals <- c(max_count)
    labels_vals <- c(max_count)
  }
  
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = RichFactor,
      y = reorder(Description, -RichFactor)
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        fill = minus_log10_p,
        size = Count
      ),
      shape = 21,
      color = "grey40"
    ) +
    ggplot2::scale_fill_gradientn(
      colors = c("#d7301f", "#045a8d"),
      limits = c(min_val, max_val),
      breaks = pretty(c(min_val, max_val), n = 3),
      labels = scales::number_format(accuracy = 0.1),
      name = "-Log10(pvalue)"
    ) +
    ggplot2::scale_size_area(
      max_size = max_size_val,
      breaks = breaks_vals,
      labels = labels_vals,
      name = "Counts"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(linewidth = 0),
      panel.grid.minor.x = ggplot2::element_line(linewidth = 0),
      axis.text.x = ggplot2::element_text(size = 12, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
      legend.title = ggplot2::element_text(
        size = 12,
        colour = "black",
        margin = ggplot2::margin(b = 8)
      ),
      legend.text = ggplot2::element_text(size = 12),
      legend.key.size = grid::unit(0.4, "cm"),
      legend.spacing = grid::unit(0.4, "cm")
    ) +
    ggplot2::labs(
      x = "RichFactor",
      y = "",
      title = ""
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colourbar(order = 1),
      size = ggplot2::guide_legend(order = 2)
    )
  
  for (fmt in tolower(save_formats)) {
    outfile <- paste0(out_prefix, ".", fmt)
    
    if (fmt %in% c("tiff", "tif")) {
      ggplot2::ggsave(
        outfile,
        plot = p,
        width = width,
        height = height,
        dpi = dpi,
        units = "in",
        compression = "lzw"
      )
    } else {
      ggplot2::ggsave(
        outfile,
        plot = p,
        width = width,
        height = height,
        dpi = dpi,
        units = "in"
      )
    }
  }
  
  invisible(p)
}

# ------------------------------------------------------------
# 7. Main wrapper
# ------------------------------------------------------------

run_phosphosite_domain_analysis <- function(
    de_results,
    domain_annotation_file,
    out_dir = "./demo/Results/domain_analysis_uniprot",
    diff_mode = c("direction", "abs_logFC"),
    keep_directions = c("Up", "Down"),
    logfc_cutoff = 1,
    top_n = 20,
    save_formats = c("pdf", "png")
) {
  diff_mode <- match.arg(diff_mode)
  
  .safe_dir_create(out_dir)
  
  domain_raw <- read_domain_annotation(domain_annotation_file)
  domain_annotation_std <- standardize_domain_annotation(domain_raw)
  
  .write_csv_gbk(
    domain_annotation_std,
    file.path(out_dir, "domain_annotation_used.csv")
  )
  
  bg_protein_table <- make_background_modified_protein_table(de_results)
  
  .write_csv_gbk(
    bg_protein_table,
    file.path(out_dir, "background_all_identified_modified_proteins.csv")
  )
  
  diff_feature_table <- make_diff_modified_feature_table(
    de_results = de_results,
    diff_mode = diff_mode,
    keep_directions = keep_directions,
    logfc_cutoff = logfc_cutoff
  )
  
  .write_csv_gbk(
    diff_feature_table,
    file.path(out_dir, "differential_modified_features_used.csv")
  )
  
  if (nrow(diff_feature_table) == 0) {
    warning(
      "No differential / trend modified features found. ",
      "Check table(de_results$direction), or use diff_mode = 'abs_logFC'."
    )
    
    return(invisible(list(
      domain_annotation = domain_annotation_std,
      bg_protein_table = bg_protein_table,
      diff_feature_table = diff_feature_table,
      diff_protein_table = tibble(),
      diff_protein_domain = tibble(),
      domain_count = tibble(),
      domain_enrichment = tibble()
    )))
  }
  
  diff_protein_table <- make_diff_modified_protein_table(diff_feature_table)
  
  .write_csv_gbk(
    diff_protein_table,
    file.path(out_dir, "differential_modified_proteins.csv")
  )
  
  diff_protein_domain <- attach_domain_annotation(
    protein_table = diff_protein_table,
    domain_annotation_std = domain_annotation_std
  )
  
  .write_csv_gbk(
    diff_protein_domain,
    file.path(out_dir, "differential_modified_protein_domain.csv")
  )
  
  unmatched <- diff_protein_domain %>%
    filter(is.na(Domain.Name) | Domain.Name == "") %>%
    select(comparison, Protein.Id, Genes) %>%
    distinct()
  
  .write_csv_gbk(
    unmatched,
    file.path(out_dir, "check_diff_proteins_without_domain_annotation.csv")
  )
  
  domain_count <- make_domain_count_table(diff_protein_domain)
  
  .write_csv_gbk(
    domain_count,
    file.path(out_dir, "domain_count_by_comparison.csv")
  )
  
  domain_enrichment <- run_domain_enrichment(
    bg_protein_table = bg_protein_table,
    diff_protein_table = diff_protein_table,
    domain_annotation_std = domain_annotation_std
  )
  
  .write_csv_gbk(
    domain_enrichment,
    file.path(out_dir, "domain_enrichment_by_comparison.csv")
  )
  
  plot_dir <- file.path(out_dir, "plots")
  .safe_dir_create(plot_dir)
  
  comps <- unique(bg_protein_table$comparison)
  
  for (comp in comps) {
    comp_slug <- .safe_slug(comp)
    
    plot_domain_analysis_bar(
      domain_count = domain_count,
      comparison = comp,
      out_prefix = file.path(
        plot_dir,
        paste0("domain_analysis_top", top_n, "_", comp_slug)
      ),
      top_n = top_n,
      save_formats = save_formats
    )
    
    plot_domain_enrichment_bubble(
      domain_enrich = domain_enrichment,
      comparison = comp,
      out_prefix = file.path(
        plot_dir,
        paste0("domain_enrichment_top", top_n, "_", comp_slug)
      ),
      top_n = top_n,
      save_formats = save_formats
    )
  }
  
  summary_table <- data.frame(
    item = c(
      "n_background_comparison_protein_rows",
      "n_diff_modified_features",
      "n_diff_comparison_protein_rows",
      "n_unique_diff_proteins_across_all_comparisons",
      "n_domain_annotation_rows",
      "n_unique_domain_names",
      "n_diff_proteins_without_domain_annotation"
    ),
    value = c(
      nrow(bg_protein_table),
      nrow(diff_feature_table),
      nrow(diff_protein_table),
      dplyr::n_distinct(diff_protein_table$Protein.Id),
      nrow(domain_annotation_std),
      dplyr::n_distinct(domain_annotation_std$Domain.Name),
      nrow(unmatched)
    )
  )
  
  .write_csv_gbk(
    summary_table,
    file.path(out_dir, "domain_analysis_summary.csv")
  )
  
  message("Domain analysis finished.")
  message("Output dir: ", out_dir)
  
  invisible(list(
    domain_annotation = domain_annotation_std,
    bg_protein_table = bg_protein_table,
    diff_feature_table = diff_feature_table,
    diff_protein_table = diff_protein_table,
    diff_protein_domain = diff_protein_domain,
    domain_count = domain_count,
    domain_enrichment = domain_enrichment,
    unmatched = unmatched,
    summary = summary_table
  ))
}

# ============================================================
# Example usage
# ============================================================
#
# source("./demo/R/06_domain_analysis.R")
#
# domain_res <- run_phosphosite_domain_analysis(
#   de_results = de$all_results,
#   domain_annotation_file = "./demo/data/domain_annotation_from_uniprot.csv",
#   out_dir = "./demo/Results/domain_analysis_uniprot",
#   diff_mode = "direction",
#   keep_directions = c("Up", "Down"),
#   top_n = 20
# )
#
# If direction has too few Up / Down rows:
#
# domain_res <- run_phosphosite_domain_analysis(
#   de_results = de$all_results,
#   domain_annotation_file = "./demo/data/domain_annotation_from_uniprot.csv",
#   out_dir = "./demo/Results/domain_analysis_uniprot_logFC0.5",
#   diff_mode = "abs_logFC",
#   logfc_cutoff = 0.5,
#   keep_directions = c("Up", "Down"),
#   top_n = 20
# )