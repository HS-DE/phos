# ============================================================
# 06a. Get protein domain annotation from UniProt
#
# Purpose:
#   Download UniProt domain / region / repeat / InterPro / Pfam annotation
#   and standardize it to:
#
#   Protein.Id
#   Genes
#   Domain.ID
#   Domain.Name
#   Domain.Source
#   Domain.Start
#   Domain.End
#
# Output:
#   ./demo/data/domain_annotation_from_uniprot.csv
#
# Note:
#   This is a UniProt-annotation based domain workflow,
#   not InterProScan prediction.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(tibble)
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

.clean_domain_name <- function(x) {
  x <- as.character(x)
  x <- gsub("^\"|\"$", "", x)
  x <- gsub("\\.$", "", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  
  x <- ifelse(
    is.na(x) | x == "" | x %in% c("NA", "NULL"),
    NA_character_,
    x
  )
  
  x
}

# ------------------------------------------------------------
# 1. Download UniProt raw domain table
# ------------------------------------------------------------

.download_uniprot_tsv <- function(query,
                                  fields,
                                  out_file) {
  query_encoded <- utils::URLencode(query, reserved = TRUE)
  
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/stream?",
    "compressed=false",
    "&format=tsv",
    "&fields=", paste(fields, collapse = ","),
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

download_uniprot_domain_table <- function(
    organism_id = 10090,
    reviewed_only = TRUE,
    out_file = "./demo/Results/domain_analysis/uniprot_annotation/uniprot_domain_raw.tsv"
) {
  .safe_dir_create(dirname(out_file))
  
  query <- paste0("(organism_id:", organism_id, ")")
  if (isTRUE(reviewed_only)) {
    query <- paste0(query, " AND (reviewed:true)")
  }
  
  fields_full <- c(
    "accession",
    "id",
    "gene_primary",
    "protein_name",
    "ft_domain",
    "ft_region",
    "ft_repeat",
    "ft_motif",
    "xref_interpro",
    "xref_pfam"
  )
  
  fields_fallback <- c(
    "accession",
    "id",
    "gene_primary",
    "protein_name",
    "ft_domain",
    "ft_region",
    "ft_repeat",
    "ft_motif"
  )
  
  x <- tryCatch(
    {
      .download_uniprot_tsv(
        query = query,
        fields = fields_full,
        out_file = out_file
      )
    },
    error = function(e) {
      warning(
        "Downloading with InterPro/Pfam xref fields failed. ",
        "Retrying with feature fields only.\nOriginal error: ",
        conditionMessage(e)
      )
      
      .download_uniprot_tsv(
        query = query,
        fields = fields_fallback,
        out_file = out_file
      )
    }
  )
  
  x
}

# ------------------------------------------------------------
# 2. Parse UniProt feature field
# ------------------------------------------------------------

.parse_feature_cell <- function(cell,
                                feature_type = "DOMAIN") {
  if (is.na(cell) || trimws(cell) == "") {
    return(tibble())
  }
  
  s <- as.character(cell)
  s <- gsub("\\r?\\n", " ", s)
  s <- gsub("\\s+", " ", s)
  s <- trimws(s)
  
  # Split by repeated feature labels, such as DOMAIN / REGION / REPEAT / MOTIF.
  pattern <- paste0("(?i)(?=", feature_type, "\\s+)")
  parts <- unlist(strsplit(s, pattern, perl = TRUE))
  parts <- trimws(parts)
  parts <- parts[parts != ""]
  parts <- parts[grepl(paste0("(?i)^", feature_type, "\\s+"), parts, perl = TRUE)]
  
  if (length(parts) == 0) {
    return(tibble())
  }
  
  out <- lapply(parts, function(p) {
    pos <- stringr::str_match(
      p,
      paste0("(?i)^", feature_type, "\\s+<?([0-9]+)\\.\\.>?([0-9]+)")
    )
    
    start <- suppressWarnings(as.integer(pos[, 2]))
    end <- suppressWarnings(as.integer(pos[, 3]))
    
    note <- stringr::str_match(p, "/note=\"([^\"]+)\"")[, 2]
    
    if (is.na(note)) {
      note <- stringr::str_match(p, "Note=([^;]+)")[, 2]
    }
    
    if (is.na(note)) {
      note <- p
      note <- gsub(paste0("(?i)^", feature_type, "\\s+<?[0-9]+\\.\\.>?[0-9]+;?"), "", note, perl = TRUE)
      note <- gsub("/evidence=.*$", "", note)
      note <- gsub(";.*$", "", note)
      note <- trimws(note)
    }
    
    note <- .clean_domain_name(note)
    
    tibble(
      Domain.ID = NA_character_,
      Domain.Name = note,
      Domain.Source = paste0("UniProtFT_", toupper(feature_type)),
      Domain.Start = start,
      Domain.End = end
    )
  })
  
  bind_rows(out) %>%
    filter(!is.na(Domain.Name), Domain.Name != "")
}

# ------------------------------------------------------------
# 3. Parse InterPro / Pfam xref fields
# ------------------------------------------------------------

.parse_xref_cell <- function(cell,
                             db = c("InterPro", "Pfam")) {
  db <- match.arg(db)
  
  if (is.na(cell) || trimws(cell) == "") {
    return(tibble())
  }
  
  s <- as.character(cell)
  s <- gsub("\\r?\\n", ";", s)
  s <- gsub("\\s+", " ", s)
  s <- trimws(s)
  
  tokens <- unlist(strsplit(s, ";"))
  tokens <- trimws(tokens)
  tokens <- tokens[tokens != ""]
  
  if (length(tokens) == 0) {
    return(tibble())
  }
  
  id_pattern <- if (db == "InterPro") "^IPR[0-9]+$" else "^PF[0-9]+$"
  
  id_idx <- grep(id_pattern, tokens)
  
  if (length(id_idx) == 0) {
    return(tibble())
  }
  
  out <- lapply(id_idx, function(i) {
    id <- tokens[i]
    
    next_token <- if (i < length(tokens)) tokens[i + 1] else NA_character_
    
    name <- if (!is.na(next_token) && !grepl(id_pattern, next_token)) {
      next_token
    } else {
      id
    }
    
    name <- .clean_domain_name(name)
    
    tibble(
      Domain.ID = id,
      Domain.Name = name,
      Domain.Source = db,
      Domain.Start = NA_integer_,
      Domain.End = NA_integer_
    )
  })
  
  bind_rows(out) %>%
    filter(!is.na(Domain.Name), Domain.Name != "")
}

# ------------------------------------------------------------
# 4. Standardize UniProt raw domain table
# ------------------------------------------------------------

standardize_uniprot_domain_table <- function(uniprot_raw,
                                             keep_only_accessions = NULL,
                                             include_feature_types = c("DOMAIN", "REPEAT", "REGION", "MOTIF"),
                                             use_xref_interpro = TRUE,
                                             use_xref_pfam = TRUE,
                                             out_file = "./demo/data/domain_annotation_from_uniprot.csv",
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
  
  domain_col <- .detect_col(
    uniprot_raw,
    patterns = c("^Domain \\[FT\\]$", "Domain.*\\[FT\\]", "^ft_domain$"),
    required = FALSE,
    what = "domain feature column"
  )
  
  region_col <- .detect_col(
    uniprot_raw,
    patterns = c("^Region \\[FT\\]$", "Region.*\\[FT\\]", "^ft_region$"),
    required = FALSE,
    what = "region feature column"
  )
  
  repeat_col <- .detect_col(
    uniprot_raw,
    patterns = c("^Repeat \\[FT\\]$", "Repeat.*\\[FT\\]", "^ft_repeat$"),
    required = FALSE,
    what = "repeat feature column"
  )
  
  motif_col <- .detect_col(
    uniprot_raw,
    patterns = c("^Motif \\[FT\\]$", "Motif.*\\[FT\\]", "^ft_motif$"),
    required = FALSE,
    what = "motif feature column"
  )
  
  interpro_col <- .detect_col(
    uniprot_raw,
    patterns = c("^InterPro$", "InterPro"),
    required = FALSE,
    what = "InterPro column"
  )
  
  pfam_col <- .detect_col(
    uniprot_raw,
    patterns = c("^Pfam$", "Pfam"),
    required = FALSE,
    what = "Pfam column"
  )
  
  rows <- vector("list", nrow(uniprot_raw))
  
  for (i in seq_len(nrow(uniprot_raw))) {
    protein_id <- as.character(uniprot_raw[[entry_col]][i])
    gene <- if (!is.na(gene_col)) as.character(uniprot_raw[[gene_col]][i]) else NA_character_
    
    pieces <- list()
    
    if ("DOMAIN" %in% include_feature_types && !is.na(domain_col)) {
      pieces[["DOMAIN"]] <- .parse_feature_cell(uniprot_raw[[domain_col]][i], "DOMAIN")
    }
    
    if ("REPEAT" %in% include_feature_types && !is.na(repeat_col)) {
      pieces[["REPEAT"]] <- .parse_feature_cell(uniprot_raw[[repeat_col]][i], "REPEAT")
    }
    
    if ("REGION" %in% include_feature_types && !is.na(region_col)) {
      pieces[["REGION"]] <- .parse_feature_cell(uniprot_raw[[region_col]][i], "REGION")
    }
    
    if ("MOTIF" %in% include_feature_types && !is.na(motif_col)) {
      pieces[["MOTIF"]] <- .parse_feature_cell(uniprot_raw[[motif_col]][i], "MOTIF")
    }
    
    if (isTRUE(use_xref_interpro) && !is.na(interpro_col)) {
      pieces[["InterPro"]] <- .parse_xref_cell(uniprot_raw[[interpro_col]][i], "InterPro")
    }
    
    if (isTRUE(use_xref_pfam) && !is.na(pfam_col)) {
      pieces[["Pfam"]] <- .parse_xref_cell(uniprot_raw[[pfam_col]][i], "Pfam")
    }
    
    one <- bind_rows(pieces)
    
    if (nrow(one) > 0) {
      rows[[i]] <- one %>%
        mutate(
          Protein.Id = protein_id,
          Genes = gene,
          .before = 1
        )
    }
  }
  
  domain_anno <- bind_rows(rows)
  
  if (nrow(domain_anno) == 0) {
    warning("No domain annotation parsed from UniProt raw table.")
    domain_anno <- tibble(
      Protein.Id = character(),
      Genes = character(),
      Domain.ID = character(),
      Domain.Name = character(),
      Domain.Source = character(),
      Domain.Start = integer(),
      Domain.End = integer()
    )
  }
  
  domain_anno <- domain_anno %>%
    mutate(
      Protein.Id = trimws(as.character(Protein.Id)),
      Protein.Id.base = sub("-.*$", "", Protein.Id),
      Genes = trimws(as.character(Genes)),
      Domain.ID = as.character(Domain.ID),
      Domain.Name = .clean_domain_name(Domain.Name),
      Domain.Source = as.character(Domain.Source)
    ) %>%
    filter(!is.na(Protein.Id), Protein.Id != "") %>%
    filter(!is.na(Domain.Name), Domain.Name != "")
  
  # If input contains isoform IDs, map canonical UniProt annotation back to input IDs.
  if (!is.null(keep_only_accessions)) {
    keep_df <- tibble(
      Protein.Id.input = unique(as.character(keep_only_accessions)),
      Protein.Id.base = sub("-.*$", "", Protein.Id.input)
    ) %>%
      filter(!is.na(Protein.Id.input), Protein.Id.input != "")
    
    domain_anno <- keep_df %>%
      left_join(domain_anno, by = "Protein.Id.base") %>%
      filter(!is.na(Domain.Name), Domain.Name != "") %>%
      mutate(
        Protein.Id = Protein.Id.input
      ) %>%
      select(
        Protein.Id,
        Genes,
        Domain.ID,
        Domain.Name,
        Domain.Source,
        Domain.Start,
        Domain.End
      )
  } else {
    domain_anno <- domain_anno %>%
      select(
        Protein.Id,
        Genes,
        Domain.ID,
        Domain.Name,
        Domain.Source,
        Domain.Start,
        Domain.End
      )
  }
  
  # Keep one protein-domain-source row.
  domain_anno <- domain_anno %>%
    distinct(
      Protein.Id,
      Domain.Name,
      Domain.Source,
      Domain.ID,
      Domain.Start,
      Domain.End,
      .keep_all = TRUE
    )
  
  write.csv(
    domain_anno,
    out_file,
    row.names = FALSE,
    fileEncoding = file_encoding
  )
  
  message("Saved standardized domain annotation: ", out_file)
  message("Rows: ", nrow(domain_anno))
  message("Proteins with at least one domain: ", dplyr::n_distinct(domain_anno$Protein.Id))
  message("Unique domain names: ", dplyr::n_distinct(domain_anno$Domain.Name))
  
  invisible(domain_anno)
}

# ------------------------------------------------------------
# 5. One-step wrapper for this phosphoproteomics workflow
# ------------------------------------------------------------

make_uniprot_domain_annotation_for_phos <- function(
    de_results,
    organism_id = 10090,
    reviewed_only = TRUE,
    use_only_proteins_in_de_results = TRUE,
    include_feature_types = c("DOMAIN", "REPEAT", "REGION", "MOTIF"),
    use_xref_interpro = TRUE,
    use_xref_pfam = TRUE,
    out_dir = "./demo/Results/domain_analysis/uniprot_annotation",
    final_out_file = "./demo/data/domain_annotation_from_uniprot.csv"
) {
  .safe_dir_create(out_dir)
  .safe_dir_create(dirname(final_out_file))
  
  if (!"Protein.Id" %in% names(de_results)) {
    stop("de_results must contain Protein.Id")
  }
  
  protein_ids <- unique(as.character(de_results$Protein.Id))
  protein_ids <- protein_ids[!is.na(protein_ids) & protein_ids != ""]
  
  raw_file <- file.path(out_dir, "uniprot_domain_raw.tsv")
  
  uniprot_raw <- download_uniprot_domain_table(
    organism_id = organism_id,
    reviewed_only = reviewed_only,
    out_file = raw_file
  )
  
  domain_anno <- standardize_uniprot_domain_table(
    uniprot_raw = uniprot_raw,
    keep_only_accessions = if (isTRUE(use_only_proteins_in_de_results)) protein_ids else NULL,
    include_feature_types = include_feature_types,
    use_xref_interpro = use_xref_interpro,
    use_xref_pfam = use_xref_pfam,
    out_file = final_out_file
  )
  
  check_summary <- data.frame(
    item = c(
      "n_protein_ids_in_de_results",
      "n_domain_annotation_rows",
      "n_proteins_with_domain",
      "n_unique_domain_names"
    ),
    value = c(
      length(protein_ids),
      nrow(domain_anno),
      dplyr::n_distinct(domain_anno$Protein.Id),
      dplyr::n_distinct(domain_anno$Domain.Name)
    )
  )
  
  check_file <- file.path(out_dir, "check_uniprot_domain_annotation_summary.csv")
  
  write.csv(
    check_summary,
    check_file,
    row.names = FALSE,
    fileEncoding = "GBK"
  )
  
  message("Saved check summary: ", check_file)
  
  invisible(list(
    raw = uniprot_raw,
    annotation = domain_anno,
    check_summary = check_summary
  ))
}

# ============================================================
# Example usage
# ============================================================
#
# source("./demo/R/06a_get_uniprot_domain_annotation.R")
#
# domain_anno <- make_uniprot_domain_annotation_for_phos(
#   de_results = de$all_results,
#   organism_id = 10090,
#   reviewed_only = TRUE,
#   use_only_proteins_in_de_results = TRUE,
#   out_dir = "./demo/Results/domain_analysis/uniprot_annotation",
#   final_out_file = "./demo/data/domain_annotation_from_uniprot.csv"
# )
#
# head(domain_anno$annotation)