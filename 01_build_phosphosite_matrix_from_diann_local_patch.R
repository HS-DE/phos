# ============================================================
# Local patch for 01_build_phosphosite_matrix_from_diann.R
#
# Usage:
#   source("01_build_phosphosite_matrix_from_diann.R")
#   source("01_build_phosphosite_matrix_from_diann_local_patch.R")
#   prep <- prepare_phosphosite_matrix(...)
#
# Why this patch exists:
#   The first version used dplyr::case_when() inside summarise().
#   case_when() evaluates all right-hand-side expressions, so max(..., na.rm=TRUE)
#   is still evaluated even when all intensities are NA, producing many warnings:
#     max里所有的参数都不存在；返回-Inf
#   This replacement avoids that warning by using an explicit helper function.
# ============================================================

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
