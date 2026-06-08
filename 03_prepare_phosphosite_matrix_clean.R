# ============================================================
# Clean upstream phosphoproteomics preparation workflow
#
# Purpose:
#   Convert DIA-NN pr_matrix.tsv + matching FASTA into clean phosphosite-level
#   tables for downstream FC analysis / KSEA / PhosR kinase activity analysis.
#
# Main outputs:
#   1) phospho_peptide_table.csv
#   2) phosphosite_mapping_all_candidates.csv
#   3) phosphosite_mapping_selected.csv
#   4) phosphosite_intensity_long.csv
#   5) phosphosite_intensity_matrix.csv
#   6) site_mat2.csv                 # compatible with old step1.R / KSEA流程.R
#   7) final_mat.csv                 # subset for FC analysis
#
# Design principles:
#   - Do not depend on variables typed manually in the R console.
#   - Do not silently truncate site mapping rows.
#   - Keep old object names such as site_mat2/final_mat available for compatibility.
#   - Use max precursor intensity by default to collapse precursor-level values to
#     modified-peptide-level values, matching the old exploratory script.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(readr)
})

# ------------------------------------------------------------
# 1. Read UniProt FASTA as named sequence vector
# ------------------------------------------------------------
read_fasta_as_map <- function(fasta_file) {
  x <- readLines(fasta_file, warn = FALSE)

  header_idx <- which(startsWith(x, ">"))
  if (length(header_idx) == 0) {
    stop("No FASTA headers found in: ", fasta_file)
  }

  end_idx <- c(header_idx[-1] - 1, length(x))
  headers <- x[header_idx]

  seqs <- purrr::map2_chr(header_idx, end_idx, function(st, ed) {
    paste0(x[(st + 1):ed], collapse = "")
  })

  # UniProt headers usually look like:
  # >sp|Q3UH06|RREB1_MOUSE Ras-responsive element-binding protein 1 OS=Mus musculus ...
  acc <- vapply(headers, function(h) {
    h <- sub("^>", "", h)
    parts <- strsplit(h, "\\|")[[1]]
    if (length(parts) >= 2) {
      parts[2]
    } else {
      strsplit(h, " ")[[1]][1]
    }
  }, character(1))

  stats::setNames(seqs, acc)
}

# ------------------------------------------------------------
# 2. Detect and optionally clean DIA-NN sample column names
# ------------------------------------------------------------
detect_sample_cols <- function(expr,
                               meta_cols = c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
                                             "First.Protein.Description", "Proteotypic",
                                             "Stripped.Sequence", "Modified.Sequence",
                                             "Precursor.Charge", "Precursor.Id")) {
  if ("Precursor.Id" %in% names(expr)) {
    start_col <- match("Precursor.Id", names(expr)) + 1
    return(names(expr)[start_col:ncol(expr)])
  }

  # Fallback: all columns outside known annotation columns are treated as samples.
  setdiff(names(expr), meta_cols)
}

clean_sample_colnames <- function(expr,
                                  sample_cols,
                                  sample_prefix_pattern = NULL,
                                  remove_raw_ext = TRUE) {
  old_names <- names(expr)
  idx <- match(sample_cols, old_names)

  if (anyNA(idx)) {
    stop("Some sample_cols are not found in expr: ", paste(sample_cols[is.na(idx)], collapse = ", "))
  }

  new_sample_names <- basename(old_names[idx])

  if (remove_raw_ext) {
    new_sample_names <- sub("\\.raw$", "", new_sample_names, ignore.case = TRUE)
  }

  # Example: if raw name is YAS202601070072-3-P-Lvs-26-1-1-1.raw
  # and sample_prefix_pattern = "YAS202601070072-3-",
  # the sample name becomes P-Lvs-26-1-1-1.
  if (!is.null(sample_prefix_pattern) && nzchar(sample_prefix_pattern)) {
    new_sample_names <- sub(paste0("^.*", sample_prefix_pattern), "", new_sample_names)
  }

  if (any(duplicated(new_sample_names))) {
    dup <- unique(new_sample_names[duplicated(new_sample_names)])
    stop("Sample names are duplicated after cleaning: ", paste(dup, collapse = ", "))
  }

  old_names[idx] <- new_sample_names
  names(expr) <- old_names

  list(expr = expr, sample_cols = new_sample_names)
}

# ------------------------------------------------------------
# 3. Parse phosphorylation sites inside Modified.Sequence
# ------------------------------------------------------------
get_phos_sites_in_peptide <- function(modseq, phos_tag = "(UniMod:21)") {
  n <- nchar(modseq)
  i <- 1
  aa_idx <- 0
  last_aa <- NA_character_
  last_pos <- NA_integer_

  phos_index <- integer(0)
  pep_pos <- integer(0)
  residue <- character(0)

  while (i <= n) {
    ch <- substr(modseq, i, i)

    # Amino acids are represented by uppercase letters in DIA-NN Modified.Sequence.
    if (grepl("^[A-Z]$", ch)) {
      aa_idx <- aa_idx + 1
      last_aa <- ch
      last_pos <- aa_idx
      i <- i + 1
      next
    }

    # Modification tag, e.g. (UniMod:21) or (UniMod:4)
    if (ch == "(") {
      j <- regexpr("\\)", substr(modseq, i, n))[1]
      if (j == -1) break

      tag <- substr(modseq, i, i + j - 1)
      if (identical(tag, phos_tag)) {
        phos_index <- c(phos_index, length(phos_index) + 1)
        pep_pos <- c(pep_pos, last_pos)
        residue <- c(residue, last_aa)
      }

      i <- i + j
      next
    }

    i <- i + 1
  }

  tibble(
    phos_index = phos_index,
    pep_pos = pep_pos,
    residue = residue
  )
}

find_all_matches <- function(prot_seq, peptide, il_equiv = FALSE) {
  if (il_equiv) {
    prot_seq <- chartr("I", "L", prot_seq)
    peptide <- chartr("I", "L", peptide)
  }

  hits <- gregexpr(peptide, prot_seq, fixed = TRUE)[[1]]
  if (length(hits) == 1 && hits[1] == -1) integer(0) else as.integer(hits)
}

# ------------------------------------------------------------
# 4. Map one modified peptide to candidate phosphosite rows
# ------------------------------------------------------------
map_one_modified_peptide <- function(protein_group,
                                     protein_ids,
                                     gene,
                                     stripped_sequence,
                                     modified_sequence,
                                     fasta_map,
                                     phos_tag = "(UniMod:21)",
                                     il_equiv = FALSE) {
  sites <- get_phos_sites_in_peptide(modified_sequence, phos_tag = phos_tag)
  if (nrow(sites) == 0) return(tibble())

  ids <- strsplit(as.character(protein_ids), ";")[[1]] |> trimws()
  lead_id <- strsplit(as.character(protein_group), ";")[[1]][1] |> trimws()

  out <- list()

  for (id_i in seq_along(ids)) {
    pid <- ids[id_i]
    pid_used <- pid

    # If isoform is absent in FASTA, try canonical accession.
    if (!pid_used %in% names(fasta_map) && grepl("-", pid_used)) {
      pid_base <- sub("-.*$", "", pid_used)
      if (pid_base %in% names(fasta_map)) pid_used <- pid_base
    }

    if (!pid_used %in% names(fasta_map)) next

    prot_seq <- as.character(fasta_map[[pid_used]])
    starts <- find_all_matches(prot_seq, stripped_sequence, il_equiv = il_equiv)
    if (length(starts) == 0) next

    for (st in starts) {
      abs_pos <- st + sites$pep_pos - 1
      abs_residue_from_fasta <- substr(prot_seq, abs_pos, abs_pos)

      out[[length(out) + 1]] <- tibble(
        Protein.Group = protein_group,
        Protein.Ids = protein_ids,
        Genes = gene,
        leading_protein = lead_id,
        Protein.Id.original = pid,
        Protein.Id = pid_used,
        protein_id_order = id_i,
        Stripped.Sequence = stripped_sequence,
        Modified.Sequence = modified_sequence,
        pep_start = st,
        pep_end = st + nchar(stripped_sequence) - 1,
        phos_index = sites$phos_index,
        pep_pos = sites$pep_pos,
        residue = sites$residue,
        abs_pos = abs_pos,
        abs_residue = sites$residue,
        abs_residue_from_fasta = abs_residue_from_fasta,
        residue_match = sites$residue == abs_residue_from_fasta,
        site_id = paste0(pid_used, "_", sites$residue, abs_pos)
      )
    }
  }

  if (length(out) == 0) tibble() else bind_rows(out)
}

# ------------------------------------------------------------
# 5. Select representative mapping without silent truncation
# ------------------------------------------------------------
select_representative_site_mapping <- function(site_map_all, out_dir = NULL) {
  if (nrow(site_map_all) == 0) {
    stop("No phosphosite mapping rows were generated. Please check FASTA and Protein.Ids.")
  }

  # Report mapping ambiguity instead of silently truncating rows.
  mapping_check <- site_map_all %>%
    group_by(Protein.Group, Protein.Ids, Genes, Stripped.Sequence, Modified.Sequence, phos_index) %>%
    summarise(
      n_candidate_rows = n(),
      n_candidate_proteins = n_distinct(Protein.Id),
      n_candidate_starts = n_distinct(pep_start),
      residue_match_all = all(residue_match),
      candidate_sites = paste(unique(site_id), collapse = ";"),
      .groups = "drop"
    ) %>%
    mutate(is_ambiguous = n_candidate_rows > 1 | !residue_match_all)

  if (!is.null(out_dir)) {
    write.csv(mapping_check, file.path(out_dir, "check_phosphosite_mapping_candidates.csv"), row.names = FALSE)
    write.csv(mapping_check %>% filter(is_ambiguous),
              file.path(out_dir, "check_phosphosite_mapping_ambiguous.csv"), row.names = FALSE)
  }

  if (any(mapping_check$is_ambiguous)) {
    warning(
      sum(mapping_check$is_ambiguous),
      " phosphosite mapping entries are ambiguous or have residue mismatch. ",
      "See check_phosphosite_mapping_ambiguous.csv. Representative rows will still be selected."
    )
  }

  site_map_selected <- site_map_all %>%
    mutate(
      # Prefer the leading protein / canonical accession / original order.
      is_leading = Protein.Id == leading_protein | Protein.Id.original == leading_protein,
      is_canonical = !str_detect(Protein.Id, "-"),
      mapping_rank = case_when(
        is_leading & is_canonical ~ 1L,
        is_leading ~ 2L,
        is_canonical ~ 3L,
        TRUE ~ 4L
      )
    ) %>%
    group_by(Protein.Group, Protein.Ids, Genes, Stripped.Sequence, Modified.Sequence, phos_index) %>%
    arrange(mapping_rank, protein_id_order, pep_start, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    select(-is_leading, -is_canonical, -mapping_rank)

  site_map_selected
}

# ------------------------------------------------------------
# 6. Collapse precursor-level intensity to modified-peptide-level intensity
# ------------------------------------------------------------
collapse_precursor_intensity <- function(expr,
                                         sample_cols,
                                         collapse_method = c("max", "sum", "mean")) {
  collapse_method <- match.arg(collapse_method)

  long_int <- expr %>%
    select(Protein.Group, Protein.Ids, Genes, Stripped.Sequence, Modified.Sequence,
           Precursor.Id, Precursor.Charge, all_of(sample_cols)) %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "intensity"
    ) %>%
    mutate(intensity = as.numeric(intensity))

  pep_int <- long_int %>%
    group_by(sample, Protein.Group, Protein.Ids, Genes, Stripped.Sequence, Modified.Sequence) %>%
    summarise(
      intensity = case_when(
        all(is.na(intensity)) ~ NA_real_,
        collapse_method == "max"  ~ max(intensity, na.rm = TRUE),
        collapse_method == "sum"  ~ sum(intensity, na.rm = TRUE),
        collapse_method == "mean" ~ mean(intensity, na.rm = TRUE),
        TRUE ~ max(intensity, na.rm = TRUE)
      ),
      n_precursors = sum(!is.na(intensity)),
      .groups = "drop"
    )

  pep_int
}

# ------------------------------------------------------------
# 7. Main workflow
# ------------------------------------------------------------
prepare_phosphosite_matrix <- function(pr_matrix_file,
                                       fasta_file,
                                       out_dir = "Results/prepare_phosphosite_matrix",
                                       sample_cols = NULL,
                                       sample_prefix_pattern = NULL,
                                       collapse_method = c("max", "sum", "mean"),
                                       phos_tag = "(UniMod:21)",
                                       il_equiv = FALSE,
                                       keep_only_detected_phos = TRUE) {
  collapse_method <- match.arg(collapse_method)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ----------------------------
  # Step 1. Read input files
  # ----------------------------
  expr <- read.delim(pr_matrix_file, header = TRUE, check.names = FALSE,
                     quote = "", comment.char = "", sep = "\t")

  if (is.null(sample_cols)) {
    sample_cols_raw <- detect_sample_cols(expr)
  } else {
    sample_cols_raw <- sample_cols
  }

  cleaned <- clean_sample_colnames(
    expr = expr,
    sample_cols = sample_cols_raw,
    sample_prefix_pattern = sample_prefix_pattern,
    remove_raw_ext = TRUE
  )
  expr <- cleaned$expr
  sample_cols <- cleaned$sample_cols

  expr <- expr %>% mutate(across(all_of(sample_cols), as.numeric))

  fasta_map <- read_fasta_as_map(fasta_file)

  write.csv(
    data.frame(sample_col = sample_cols),
    file = file.path(out_dir, "sample_columns_used.csv"),
    row.names = FALSE
  )

  # ----------------------------
  # Step 2. Keep phosphopeptides
  # ----------------------------
  phos_regex <- stringr::fixed(phos_tag)

  df_phos <- expr %>%
    filter(str_detect(Modified.Sequence, phos_regex))

  if (keep_only_detected_phos) {
    df_phos <- df_phos %>%
      filter(if_any(all_of(sample_cols), ~ !is.na(.) & . > 0))
  }

  # One row per modified peptide sequence / protein group for mapping.
  phos_peptide_table <- df_phos %>%
    distinct(Protein.Group, Protein.Ids, Protein.Names, Genes,
             Stripped.Sequence, Modified.Sequence, .keep_all = TRUE) %>%
    select(Protein.Group, Protein.Ids, Protein.Names, Genes,
           Stripped.Sequence, Modified.Sequence)

  write.csv(phos_peptide_table, file.path(out_dir, "phospho_peptide_table.csv"), row.names = FALSE)

  # ----------------------------
  # Step 3. Map modified peptides to absolute protein sites
  # ----------------------------
  site_map_all <- pmap_dfr(
    list(
      phos_peptide_table$Protein.Group,
      phos_peptide_table$Protein.Ids,
      phos_peptide_table$Genes,
      phos_peptide_table$Stripped.Sequence,
      phos_peptide_table$Modified.Sequence
    ),
    ~ map_one_modified_peptide(
      protein_group = ..1,
      protein_ids = ..2,
      gene = ..3,
      stripped_sequence = ..4,
      modified_sequence = ..5,
      fasta_map = fasta_map,
      phos_tag = phos_tag,
      il_equiv = il_equiv
    )
  )

  write.csv(site_map_all, file.path(out_dir, "phosphosite_mapping_all_candidates.csv"), row.names = FALSE)

  site_map_selected <- select_representative_site_mapping(site_map_all, out_dir = out_dir)
  write.csv(site_map_selected, file.path(out_dir, "phosphosite_mapping_selected.csv"), row.names = FALSE)

  # Check whether each modified sequence has the expected number of selected site rows.
  mapping_count_check <- site_map_selected %>%
    group_by(Modified.Sequence) %>%
    summarise(
      n_phos_in_sequence = first(str_count(Modified.Sequence, fixed(phos_tag))),
      n_selected_sites = n(),
      ok = n_phos_in_sequence == n_selected_sites,
      sites = paste(site_id, collapse = ";"),
      .groups = "drop"
    )

  write.csv(mapping_count_check, file.path(out_dir, "check_modified_sequence_site_counts.csv"), row.names = FALSE)

  if (any(!mapping_count_check$ok)) {
    warning(
      sum(!mapping_count_check$ok),
      " Modified.Sequence entries have inconsistent selected site counts. ",
      "See check_modified_sequence_site_counts.csv."
    )
  }

  # ----------------------------
  # Step 4. Collapse precursor intensities to modified-peptide intensities
  # ----------------------------
  pep_int <- collapse_precursor_intensity(
    expr = expr,
    sample_cols = sample_cols,
    collapse_method = collapse_method
  ) %>%
    semi_join(phos_peptide_table,
              by = c("Protein.Group", "Protein.Ids", "Genes", "Stripped.Sequence", "Modified.Sequence"))

  write.csv(pep_int, file.path(out_dir, "phospho_peptide_intensity_long.csv"), row.names = FALSE)

  pep_mat <- pep_int %>%
    select(Protein.Group, Protein.Ids, Genes, Stripped.Sequence, Modified.Sequence, sample, intensity) %>%
    pivot_wider(names_from = sample, values_from = intensity)

  write.csv(pep_mat, file.path(out_dir, "phospho_peptide_intensity_matrix.csv"), row.names = FALSE)

  # ----------------------------
  # Step 5. Build phosphosite-level intensity table
  # ----------------------------
  site_table_with_int <- site_map_selected %>%
    left_join(
      pep_int,
      by = c("Protein.Group", "Protein.Ids", "Genes", "Stripped.Sequence", "Modified.Sequence")
    )

  phosphosite_intensity_long <- site_table_with_int %>%
    select(Protein.Group, Protein.Id, Genes, site_id, abs_pos, abs_residue,
           Stripped.Sequence, Modified.Sequence, phos_index,
           sample, intensity, n_precursors)

  write.csv(phosphosite_intensity_long,
            file.path(out_dir, "phosphosite_intensity_long.csv"), row.names = FALSE)

  site_mat <- phosphosite_intensity_long %>%
    select(Protein.Group, Protein.Id, Genes, site_id, abs_pos, abs_residue,
           Stripped.Sequence, Modified.Sequence, phos_index, sample, intensity) %>%
    distinct(Protein.Group, Protein.Id, Genes, site_id, abs_pos, abs_residue,
             Stripped.Sequence, Modified.Sequence, phos_index, sample, .keep_all = TRUE) %>%
    pivot_wider(names_from = sample, values_from = intensity)

  # Residue.Both is retained for compatibility with old KSEAapp/step scripts.
  residue_both_df <- site_map_selected %>%
    group_by(Modified.Sequence) %>%
    summarise(
      Residue.Both = paste0(unique(paste0(abs_residue, abs_pos)) %>% sort(), collapse = ";"),
      .groups = "drop"
    )

  site_mat2 <- site_mat %>%
    left_join(residue_both_df, by = "Modified.Sequence") %>%
    relocate(Residue.Both, .after = Modified.Sequence)

  final_mat <- site_mat2 %>%
    select(Protein.Id, Genes, Stripped.Sequence, Modified.Sequence, Residue.Both, all_of(sample_cols))

  write.csv(site_mat2, file.path(out_dir, "site_mat2.csv"), row.names = FALSE)
  write.csv(final_mat, file.path(out_dir, "final_mat.csv"), row.names = FALSE)
  write.csv(site_mat2, file.path(out_dir, "phosphosite_intensity_matrix.csv"), row.names = FALSE)

  # ----------------------------
  # Step 6. Return objects for interactive use
  # ----------------------------
  invisible(list(
    expr = expr,
    fasta_map = fasta_map,
    sample_cols = sample_cols,
    phos_peptide_table = phos_peptide_table,
    site_map_all = site_map_all,
    site_map_selected = site_map_selected,
    pep_int = pep_int,
    pep_mat = pep_mat,
    phosphosite_intensity_long = phosphosite_intensity_long,
    site_mat = site_mat,
    site_mat2 = site_mat2,
    final_mat = final_mat,
    mapping_count_check = mapping_count_check
  ))
}

# ============================================================
# Example usage
# ============================================================
# prep <- prepare_phosphosite_matrix(
#   pr_matrix_file = "4_Phos_Mouse_report.pr_matrix.tsv",
#   fasta_file = "UP000000589_Mouse_Reviewed_17246_include_Isoform_2026_03_13.fasta",
#   out_dir = "Results/prepare_phosphosite_matrix",
#   sample_prefix_pattern = "YAS202601070072-3-",
#   collapse_method = "max"
# )
#
# # Keep old downstream object names available:
# expr <- prep$expr
# fasta_map <- prep$fasta_map
# sample_cols <- prep$sample_cols
# site_mat2 <- prep$site_mat2
# final_mat <- prep$final_mat
