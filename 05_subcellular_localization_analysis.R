# ============================================================
# 05. Subcellular localization analysis
#
# Input:
#   1) de$all_results
#   2) subcellular_annotation_from_uniprot.csv
#
# Output:
#   - differential_modified_features_used.csv
#   - differential_modified_proteins.csv
#   - differential_modified_protein_subcellular_location.csv
#   - subcellular_location_count_by_comparison.csv
#   - pie charts for each comparison
#
# Meaning:
#   This module counts subcellular localization of proteins
#   that contain differential / trend phosphosites.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
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

read_subcellular_annotation <- function(annotation_file) {
  if (!file.exists(annotation_file)) {
    stop("Annotation file not found: ", annotation_file)
  }
  
  ext <- tolower(tools::file_ext(annotation_file))
  
  if (ext == "csv") {
    anno <- read.csv(annotation_file, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    anno <- read.delim(annotation_file, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    stop("Only csv / tsv / txt annotation files are supported here.")
  }
  
  anno
}

# ------------------------------------------------------------
# Standardize annotation table
# ------------------------------------------------------------

standardize_subcellular_annotation <- function(annotation) {
  annotation <- as.data.frame(annotation, check.names = FALSE)
  
  if (!"Protein.Id" %in% colnames(annotation)) {
    stop("subcellular annotation must contain column: Protein.Id")
  }
  
  if (!"Location" %in% colnames(annotation)) {
    stop("subcellular annotation must contain column: Location")
  }
  
  if (!"Genes" %in% colnames(annotation)) {
    annotation$Genes <- NA_character_
  }
  
  if (!"Location.Raw" %in% colnames(annotation)) {
    annotation$Location.Raw <- annotation$Location
  }
  
  annotation %>%
    transmute(
      Protein.Id = as.character(Protein.Id),
      Genes = as.character(Genes),
      Location.Raw = as.character(Location.Raw),
      Location = as.character(Location)
    ) %>%
    mutate(
      Protein.Id = trimws(Protein.Id),
      Genes = trimws(Genes),
      Location.Raw = trimws(Location.Raw),
      Location = trimws(Location),
      Location = ifelse(is.na(Location) | Location == "", "Others", Location)
    ) %>%
    distinct(Protein.Id, .keep_all = TRUE)
}

# ------------------------------------------------------------
# Select differential / trend modified features
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

# ------------------------------------------------------------
# Convert modified features to protein-level table
# ------------------------------------------------------------

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
      max_abs_logFC = max(abs(logFC), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      max_abs_logFC = ifelse(is.infinite(max_abs_logFC), NA_real_, max_abs_logFC)
    )
}

# ------------------------------------------------------------
# Attach localization annotation
# ------------------------------------------------------------

attach_subcellular_annotation <- function(diff_protein_table,
                                          annotation_std) {
  if (nrow(diff_protein_table) == 0) {
    return(tibble())
  }
  
  diff_protein_table %>%
    left_join(
      annotation_std %>%
        select(Protein.Id, Location.Raw, Location),
      by = "Protein.Id"
    ) %>%
    mutate(
      Location.Raw = ifelse(is.na(Location.Raw) | Location.Raw == "", "Others", Location.Raw),
      Location = ifelse(is.na(Location) | Location == "", "Others", Location)
    )
}

# ------------------------------------------------------------
# Count localization
# ------------------------------------------------------------

count_subcellular_location <- function(diff_protein_loc) {
  if (nrow(diff_protein_loc) == 0) {
    return(tibble())
  }
  
  diff_protein_loc %>%
    group_by(comparison, Location) %>%
    summarise(
      protein_count = n_distinct(Protein.Id),
      .groups = "drop"
    ) %>%
    group_by(comparison) %>%
    mutate(
      total_protein_count = sum(protein_count),
      percent = protein_count / total_protein_count * 100
    ) %>%
    ungroup() %>%
    arrange(comparison, desc(protein_count))
}

count_subcellular_location_by_direction <- function(diff_feature_table,
                                                    annotation_std) {
  if (nrow(diff_feature_table) == 0) {
    return(tibble())
  }
  
  protein_direction <- diff_feature_table %>%
    group_by(comparison, direction, Protein.Id, Genes) %>%
    summarise(
      n_modified_features = n_distinct(feature_id),
      n_modified_sequences = n_distinct(Modified.Sequence),
      residues = paste(sort(unique(Residue.Both)), collapse = ";"),
      max_abs_logFC = max(abs(logFC), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      max_abs_logFC = ifelse(is.infinite(max_abs_logFC), NA_real_, max_abs_logFC)
    )
  
  protein_direction_loc <- protein_direction %>%
    left_join(
      annotation_std %>%
        select(Protein.Id, Location.Raw, Location),
      by = "Protein.Id"
    ) %>%
    mutate(
      Location.Raw = ifelse(is.na(Location.Raw) | Location.Raw == "", "Others", Location.Raw),
      Location = ifelse(is.na(Location) | Location == "", "Others", Location)
    )
  
  protein_direction_loc %>%
    group_by(comparison, direction, Location) %>%
    summarise(
      protein_count = n_distinct(Protein.Id),
      .groups = "drop"
    ) %>%
    group_by(comparison, direction) %>%
    mutate(
      total_protein_count = sum(protein_count),
      percent = protein_count / total_protein_count * 100
    ) %>%
    ungroup() %>%
    arrange(comparison, direction, desc(protein_count))
}

# ------------------------------------------------------------
# Plot pie chart
# ------------------------------------------------------------

.location_colors <- c(
  Nuclear = "#8DD3C7",
  Cytoplasmic = "#FFFFB3",
  PlasmaMembrane = "#BEBADA",
  Mitochondrial = "#FB8072",
  Extracellular = "#FDB462",
  Others = "#B3B3B3"
)
plot_subcellular_pie <- function(count_df,
                                 comparison,
                                 out_prefix,
                                 save_formats = c("pdf", "png"),
                                 width = 7,
                                 height = 5,
                                 dpi = 300) {
  df <- count_df %>%
    filter(comparison == !!comparison) %>%
    arrange(desc(protein_count))
  
  if (nrow(df) == 0) {
    warning("No data to plot for comparison: ", comparison)
    return(invisible(NULL))
  }
  
  df <- df %>%
    mutate(
      Location = as.character(Location),
      percent_label = sprintf("%.1f%%", percent),
      label = paste0(Location, ", ", protein_count),
      legend_label = paste0(Location, " (", protein_count, ", ", percent_label, ")")
    )
  
  ## 为了图例也显示数量和比例
  df$legend_label <- factor(df$legend_label, levels = df$legend_label)
  
  color_map <- .location_colors[df$Location]
  names(color_map) <- df$legend_label
  color_map[is.na(color_map)] <- "#B3B3B3"
  
  p <- ggplot(df, aes(x = "", y = protein_count, fill = legend_label)) +
    geom_col(width = 1, color = "white", linewidth = 0.4) +
    coord_polar(theta = "y") +
    labs(
      title = paste0(comparison, "\nSubcellular localization"),
      fill = "Location"
    ) +
    scale_fill_manual(values = color_map, drop = FALSE) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "right",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )
  
  for (fmt in tolower(save_formats)) {
    outfile <- paste0(out_prefix, ".", fmt)
    ggsave(outfile, p, width = width, height = height, dpi = dpi)
  }
  
  invisible(p)
}

# ------------------------------------------------------------
# Plot stacked bar
# ------------------------------------------------------------

plot_subcellular_stacked_bar <- function(count_df,
                                         out_prefix,
                                         save_formats = c("pdf", "png"),
                                         width = 8,
                                         height = 5,
                                         dpi = 300) {
  if (nrow(count_df) == 0) {
    warning("No data to plot stacked bar.")
    return(invisible(NULL))
  }
  
  p <- count_df %>%
    mutate(
      comparison = factor(comparison, levels = unique(comparison)),
      Location = factor(Location, levels = names(.location_colors))
    ) %>%
    ggplot(aes(x = comparison, y = protein_count, fill = Location)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = .location_colors, drop = FALSE) +
    labs(
      x = NULL,
      y = "Number of differential modified proteins",
      fill = "Location"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank()
    )
  
  for (fmt in tolower(save_formats)) {
    outfile <- paste0(out_prefix, ".", fmt)
    ggsave(outfile, p, width = width, height = height, dpi = dpi)
  }
  
  invisible(p)
}

# ------------------------------------------------------------
# Main function
# ------------------------------------------------------------

run_phosphosite_subcellular_localization <- function(
    de_results,
    subcellular_annotation_file,
    out_dir = "./demo/Results/subcellular_localization_uniprot",
    diff_mode = c("direction", "abs_logFC"),
    keep_directions = c("Up", "Down"),
    logfc_cutoff = 1,
    save_formats = c("pdf", "png")
) {
  diff_mode <- match.arg(diff_mode)
  
  .safe_dir_create(out_dir)
  
  annotation_raw <- read_subcellular_annotation(subcellular_annotation_file)
  
  annotation_std <- standardize_subcellular_annotation(annotation_raw)
  
  .write_csv_gbk(
    annotation_std,
    file.path(out_dir, "subcellular_annotation_used.csv")
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
      "No differential modified features found. ",
      "Check table(de_results$direction), or use diff_mode = 'abs_logFC'."
    )
    
    return(invisible(list(
      annotation = annotation_std,
      diff_feature_table = diff_feature_table,
      diff_protein_table = tibble(),
      diff_protein_loc = tibble(),
      subcell_count = tibble(),
      subcell_count_by_direction = tibble()
    )))
  }
  
  diff_protein_table <- make_diff_modified_protein_table(diff_feature_table)
  
  .write_csv_gbk(
    diff_protein_table,
    file.path(out_dir, "differential_modified_proteins.csv")
  )
  
  diff_protein_loc <- attach_subcellular_annotation(
    diff_protein_table = diff_protein_table,
    annotation_std = annotation_std
  )
  
  .write_csv_gbk(
    diff_protein_loc,
    file.path(out_dir, "differential_modified_protein_subcellular_location.csv")
  )
  
  subcell_count <- count_subcellular_location(diff_protein_loc)
  
  .write_csv_gbk(
    subcell_count,
    file.path(out_dir, "subcellular_location_count_by_comparison.csv")
  )
  
  subcell_count_by_direction <- count_subcellular_location_by_direction(
    diff_feature_table = diff_feature_table,
    annotation_std = annotation_std
  )
  
  .write_csv_gbk(
    subcell_count_by_direction,
    file.path(out_dir, "subcellular_location_count_by_comparison_direction.csv")
  )
  
  unmatched <- diff_protein_loc %>%
    filter(Location == "Others", Location.Raw == "Others") %>%
    select(comparison, Protein.Id, Genes) %>%
    distinct()
  
  .write_csv_gbk(
    unmatched,
    file.path(out_dir, "check_unmatched_or_others_proteins.csv")
  )
  
  pie_dir <- file.path(out_dir, "pie_charts")
  .safe_dir_create(pie_dir)
  
  for (comp in unique(subcell_count$comparison)) {
    comp_slug <- .safe_slug(comp)
    
    plot_subcellular_pie(
      count_df = subcell_count,
      comparison = comp,
      out_prefix = file.path(
        pie_dir,
        paste0("subcellular_location_pie_", comp_slug)
      ),
      save_formats = save_formats
    )
  }
  
  plot_subcellular_stacked_bar(
    count_df = subcell_count,
    out_prefix = file.path(
      out_dir,
      "subcellular_location_stacked_bar_all_comparisons"
    ),
    save_formats = save_formats
  )
  
  summary_table <- data.frame(
    item = c(
      "n_differential_modified_features",
      "n_comparison_protein_rows",
      "n_unique_proteins_across_all_comparisons",
      "n_unmatched_or_others_rows"
    ),
    value = c(
      nrow(diff_feature_table),
      nrow(diff_protein_table),
      dplyr::n_distinct(diff_protein_table$Protein.Id),
      nrow(unmatched)
    )
  )
  
  .write_csv_gbk(
    summary_table,
    file.path(out_dir, "subcellular_localization_summary.csv")
  )
  
  message("Subcellular localization analysis finished.")
  message("Output dir: ", out_dir)
  
  invisible(list(
    annotation = annotation_std,
    diff_feature_table = diff_feature_table,
    diff_protein_table = diff_protein_table,
    diff_protein_loc = diff_protein_loc,
    subcell_count = subcell_count,
    subcell_count_by_direction = subcell_count_by_direction,
    unmatched = unmatched,
    summary = summary_table
  ))
}

# ============================================================
# Example usage
# ============================================================
#
# source("./demo/R/05_subcellular_localization_analysis.R")
#
# subcell <- run_phosphosite_subcellular_localization(
#   de_results = de$all_results,
#   subcellular_annotation_file = "./demo/data/subcellular_annotation_from_uniprot.csv",
#   out_dir = "./demo/Results/subcellular_localization_uniprot",
#   diff_mode = "direction",
#   keep_directions = c("Up", "Down")
# )
#
# If direction has too few Up / Down rows:
#
# subcell <- run_phosphosite_subcellular_localization(
#   de_results = de$all_results,
#   subcellular_annotation_file = "./demo/data/subcellular_annotation_from_uniprot.csv",
#   out_dir = "./demo/Results/subcellular_localization_uniprot_logFC0.5",
#   diff_mode = "abs_logFC",
#   logfc_cutoff = 0.5,
#   keep_directions = c("Up", "Down")
# )