# ============================================================
# 01. Build phosphosite matrix from DIA-NN pr_matrix
#
# This wrapper loads the core implementation and overrides the precursor
# collapsing function with the safe version. Therefore users only need:
#
#   source("01_build_phosphosite_matrix_from_diann.R")
#
# No extra local patch script is required.
# ============================================================

.this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
.this_dir <- if (!is.na(.this_file)) dirname(.this_file) else getwd()

source(file.path(.this_dir, "01_build_phosphosite_matrix_from_diann_core.R"), encoding = "UTF-8")

# ------------------------------------------------------------
# Safe replacement for collapse_precursor_intensity()
# ------------------------------------------------------------
# The previous core implementation used dplyr::case_when() inside summarise().
# case_when() evaluates all right-hand-side expressions, so max(..., na.rm=TRUE)
# was still evaluated for all-NA groups and produced many warnings:
#   max里所有的参数都不存在；返回-Inf
# This function avoids that by explicitly checking each group first.
collapse_precursor_intensity <- function(expr,
                                         sample_cols,
                                         collapse_method = c("max", "sum", "mean")) {
  collapse_method <- match.arg(collapse_method)

  collapse_one <- function(x) {
    x <- as.numeric(x)
    x <- x[!is.na(x)]

    if (length(x) == 0) {
      return(NA_real_)
    }

    if (collapse_method == "max") {
      return(max(x))
    }

    if (collapse_method == "sum") {
      return(sum(x))
    }

    if (collapse_method == "mean") {
      return(mean(x))
    }

    max(x)
  }

  long_int <- expr %>%
    dplyr::select(
      Protein.Group, Protein.Ids, Genes,
      Stripped.Sequence, Modified.Sequence,
      Precursor.Id, Precursor.Charge,
      dplyr::all_of(sample_cols)
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(sample_cols),
      names_to = "sample",
      values_to = "intensity"
    ) %>%
    dplyr::mutate(intensity = as.numeric(intensity))

  pep_int <- long_int %>%
    dplyr::group_by(sample, Protein.Group, Protein.Ids, Genes, Stripped.Sequence, Modified.Sequence) %>%
    dplyr::summarise(
      intensity = collapse_one(intensity),
      n_precursors = sum(!is.na(intensity)),
      .groups = "drop"
    )

  pep_int
}
