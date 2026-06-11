# ============================================================
# 05a. Get subcellular localization annotation from UniProt
#
# Purpose:
#   Download UniProtKB subcellular location annotation
#   and convert it to:
#     Protein.Id, Genes, Location.Raw, GO.Cellular.Component, Location
#
# This file generates:
#   ./demo/data/subcellular_annotation_from_uniprot.csv
#
# Then this file can be used by:
#   run_phosphosite_subcellular_localization()
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(tidyr)
})

.safe_dir_create <- function(path) {
  if (!is.null(path) && nzchar(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

.detect_col <- function(df, patterns, required = TRUE, what = "column") {
  nms <- names(df)
  
  for (pat in patterns) {
    hit <- grep(pat, nms, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  
  if (required) {
    stop(
      "Cannot detect ", what, ". Available columns are:\n",
      paste(nms, collapse = "\n")
    )
  }
  
  NA_character_
}

# ------------------------------------------------------------
# 1. Download UniProtKB annotation table
# ------------------------------------------------------------

download_uniprot_subcellular_table <- function(
    organism_id = 10090,
    reviewed_only = TRUE,
    out_file = "./demo/Results/subcellular_localization/uniprot_annotation/uniprot_mouse_subcellular_raw.tsv"
) {
  .safe_dir_create(dirname(out_file))
  
  query <- paste0("(organism_id:", organism_id, ")")
  if (isTRUE(reviewed_only)) {
    query <- paste0(query, " AND (reviewed:true)")
  }
  
  query_encoded <- utils::URLencode(query, reserved = TRUE)
  
  # fields:
  # accession                 -> UniProt accession
  # id                        -> UniProt entry name
  # gene_primary              -> primary gene name
  # protein_name              -> protein name
  # cc_subcellular_location   -> Subcellular location [CC]
  # go_c                      -> Gene Ontology cellular component
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/stream?",
    "compressed=false",
    "&format=tsv",
    "&fields=accession,id,gene_primary,protein_name,cc_subcellular_location,go_c",
    "&query=", query_encoded
  )
  
  message("Downloading UniProt annotation...")
  message(url)
  
  x <- readr::read_tsv(url, show_col_types = FALSE, progress = TRUE)
  
  readr::write_tsv(x, out_file)
  
  message("Saved raw UniProt table: ", out_file)
  message("Rows: ", nrow(x))
  message("Columns: ", paste(names(x), collapse = " | "))
  
  x
}

# ------------------------------------------------------------
# 2. Convert UniProt raw location text to report-style categories
# ------------------------------------------------------------

normalize_uniprot_location_to_report_category <- function(location_cc,
                                                          go_cc = NA_character_) {
  text <- paste(location_cc, go_cc, sep = " ; ")
  text <- tolower(text)
  text[is.na(text)] <- ""
  
  dplyr::case_when(
    str_detect(text, "nucleus|nuclear|nucleoplasm|nucleolus|chromosome") ~ "Nuclear",
    
    str_detect(text, "cytoplasm|cytosol|cytoskeleton|cell cortex") ~ "Cytoplasmic",
    
    str_detect(text, "plasma membrane|cell membrane|cell surface") ~ "PlasmaMembrane",
    
    str_detect(text, "mitochondrion|mitochondrial") ~ "Mitochondrial",
    
    str_detect(text, "secreted|extracellular|cell exterior|extracellular space|extracellular region") ~ "Extracellular",
    
    TRUE ~ "Others"
  )
}

# ------------------------------------------------------------
# 3. Standardize UniProt table
# ------------------------------------------------------------

standardize_uniprot_subcellular_table <- function(uniprot_raw,
                                                  keep_only_accessions = NULL,
                                                  out_file = "./demo/data/subcellular_annotation_from_uniprot.csv",
                                                  file_encoding = "GBK") {
  .safe_dir_create(dirname(out_file))
  
  uniprot_raw <- as.data.frame(uniprot_raw, check.names = FALSE)
  
  entry_col <- .detect_col(
    uniprot_raw,
    patterns = c("^Entry$", "^Accession$", "accession"),
    what = "UniProt accession column"
  )
  
  gene_col <- .detect_col(
    uniprot_raw,
    patterns = c("Gene Names \\(primary\\)", "Gene.*primary", "^Gene$", "Gene Names"),
    required = FALSE,
    what = "gene column"
  )
  
  loc_col <- .detect_col(
    uniprot_raw,
    patterns = c("Subcellular location", "cc_subcellular", "Location"),
    required = FALSE,
    what = "subcellular location column"
  )
  
  go_col <- .detect_col(
    uniprot_raw,
    patterns = c("Gene Ontology \\(cellular component\\)", "GO.*cellular", "go_c"),
    required = FALSE,
    what = "GO cellular component column"
  )
  
  anno <- uniprot_raw %>%
    transmute(
      Protein.Id = as.character(.data[[entry_col]]),
      Genes = if (!is.na(gene_col)) as.character(.data[[gene_col]]) else NA_character_,
      Location.Raw = if (!is.na(loc_col)) as.character(.data[[loc_col]]) else NA_character_,
      GO.Cellular.Component = if (!is.na(go_col)) as.character(.data[[go_col]]) else NA_character_
    ) %>%
    mutate(
      Protein.Id = trimws(Protein.Id),
      Protein.Id.base = sub("-.*$", "", Protein.Id),
      Genes = trimws(Genes),
      Location.Raw = trimws(Location.Raw),
      GO.Cellular.Component = trimws(GO.Cellular.Component),
      Location = normalize_uniprot_location_to_report_category(
        Location.Raw,
        GO.Cellular.Component
      )
    )
  
  if (!is.null(keep_only_accessions)) {
    keep_df <- tibble(
      Protein.Id.input = unique(as.character(keep_only_accessions)),
      Protein.Id.base = sub("-.*$", "", Protein.Id.input)
    )
    
    anno <- keep_df %>%
      left_join(anno, by = "Protein.Id.base") %>%
      mutate(
        Protein.Id = ifelse(is.na(Protein.Id.input), Protein.Id, Protein.Id.input),
        Genes = ifelse(is.na(Genes), NA_character_, Genes),
        Location.Raw = ifelse(is.na(Location.Raw), NA_character_, Location.Raw),
        GO.Cellular.Component = ifelse(is.na(GO.Cellular.Component), NA_character_, GO.Cellular.Component),
        Location = ifelse(is.na(Location), "Others", Location)
      ) %>%
      select(
        Protein.Id,
        Genes,
        Location.Raw,
        GO.Cellular.Component,
        Location
      )
  } else {
    anno <- anno %>%
      select(
        Protein.Id,
        Genes,
        Location.Raw,
        GO.Cellular.Component,
        Location
      )
  }
  
  anno <- anno %>%
    distinct(Protein.Id, .keep_all = TRUE)
  
  write.csv(
    anno,
    out_file,
    row.names = FALSE,
    fileEncoding = file_encoding
  )
  
  message("Saved standardized annotation: ", out_file)
  message("Rows: ", nrow(anno))
  print(table(anno$Location, useNA = "ifany"))
  
  anno
}

# ------------------------------------------------------------
# 4. One-step wrapper for your de$all_results
# ------------------------------------------------------------

make_uniprot_subcellular_annotation_for_phos <- function(
    de_results,
    organism_id = 10090,
    reviewed_only = TRUE,
    use_only_proteins_in_de_results = TRUE,
    out_dir = "./demo/Results/subcellular_localization/uniprot_annotation",
    final_out_file = "./demo/data/subcellular_annotation_from_uniprot.csv"
) {
  .safe_dir_create(out_dir)
  .safe_dir_create(dirname(final_out_file))
  
  if (!"Protein.Id" %in% names(de_results)) {
    stop("de_results must contain Protein.Id")
  }
  
  protein_ids <- unique(as.character(de_results$Protein.Id))
  protein_ids <- protein_ids[!is.na(protein_ids) & protein_ids != ""]
  
  raw_file <- file.path(out_dir, "uniprot_subcellular_raw.tsv")
  
  uniprot_raw <- download_uniprot_subcellular_table(
    organism_id = organism_id,
    reviewed_only = reviewed_only,
    out_file = raw_file
  )
  
  anno <- standardize_uniprot_subcellular_table(
    uniprot_raw = uniprot_raw,
    keep_only_accessions = if (isTRUE(use_only_proteins_in_de_results)) protein_ids else NULL,
    out_file = final_out_file
  )
  
  check_file <- file.path(out_dir, "check_uniprot_subcellular_annotation_summary.csv")
  
  check_summary <- data.frame(
    item = c(
      "n_protein_ids_in_de_results",
      "n_annotation_rows",
      "n_others",
      "n_non_others"
    ),
    value = c(
      length(protein_ids),
      nrow(anno),
      sum(anno$Location == "Others", na.rm = TRUE),
      sum(anno$Location != "Others", na.rm = TRUE)
    )
  )
  
  write.csv(check_summary, check_file, row.names = FALSE, fileEncoding = "GBK")
  
  message("Saved check summary: ", check_file)
  
  invisible(list(
    raw = uniprot_raw,
    annotation = anno,
    check_summary = check_summary
  ))
}

# ============================================================
# Example usage
# ============================================================
#
# source("./demo/R/05a_get_uniprot_subcellular_annotation.R")
#
# uniprot_loc <- make_uniprot_subcellular_annotation_for_phos(
#   de_results = de$all_results,
#   organism_id = 10090, # 人：organism_id = 9606；大鼠：organism_id = 10116
#   reviewed_only = TRUE,
#   use_only_proteins_in_de_results = TRUE,
#   out_dir = "./demo/Results/subcellular_localization/uniprot_annotation",
#   final_out_file = "./demo/data/subcellular_annotation_from_uniprot.csv"
# )
#
# head(uniprot_loc$annotation)
#