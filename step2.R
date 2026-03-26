# ============================================================
# Phospho-site/peptide matrix: log2 + lenient missing + log2FC
# 4 groups, 1 sample each
# ============================================================

library(dplyr)
library(tidyr)

# -----------------------------
# 0) Inputs you already have:
#   final_mat (data.frame/tibble)
#   sample_metadata (data.frame)
#   sample_cols (character vector)
# -----------------------------

# -----------------------------
# 1) Basic checks
# -----------------------------
stopifnot(all(sample_cols %in% colnames(final_mat)))
stopifnot(all(sample_metadata$sample_id %in% sample_cols))

# 你这个矩阵的注释列（按你当前 final_mat 来）
stopifnot(all(meta_cols %in% colnames(final_mat)))

# -----------------------------
# 2) Ensure numeric + log2 transform
#    (treat <=0 as NA)
# -----------------------------
mat0 <- final_mat %>%
  mutate(across(all_of(sample_cols), ~ as.numeric(.))) %>%
  mutate(across(all_of(sample_cols), ~ ifelse(is.na(.) | . <= 0, NA_real_, .)))

mat_log2 <- mat0 %>%
  mutate(across(all_of(sample_cols), ~ log2(.)))

# -----------------------------
# 3) Optional: median scaling in log2 space
#    (recommended; comment out if you don't want)
# -----------------------------
sample_medians <- sapply(mat_log2[, sample_cols, drop = FALSE], median, na.rm = TRUE)

mat_log2_scaled <- mat_log2
mat_log2_scaled[, sample_cols] <- sweep(
  as.matrix(mat_log2_scaled[, sample_cols, drop = FALSE]),
  2, sample_medians, "-"
)

# -----------------------------
# 4) Lenient missing imputation (per-sample 1% quantile in log2 space)
#    NOTE: if a sample has too few values, quantile may be NA.
#          We'll fallback to (min - 1) for that sample.
# -----------------------------
q01 <- sapply(mat_log2_scaled[, sample_cols, drop = FALSE], function(v) {
  v2 <- v[is.finite(v)]
  if (length(v2) < 20) return(NA_real_)
  as.numeric(stats::quantile(v2, probs = 0.01, na.rm = TRUE, type = 7))
})

mins <- sapply(mat_log2_scaled[, sample_cols, drop = FALSE], function(v) {
  v2 <- v[is.finite(v)]
  if (length(v2) == 0) return(NA_real_)
  min(v2, na.rm = TRUE)
})

fill_val <- q01
fill_val[is.na(fill_val)] <- mins[is.na(fill_val)] - 1  # fallback

mat_imp <- mat_log2_scaled %>%
  mutate(across(all_of(sample_cols), ~ ifelse(is.na(.), fill_val[cur_column()], .)))

# -----------------------------
# 5) Helper: compute log2FC = (B - A)
# -----------------------------
compute_log2fc <- function(df_imp, sample_A, sample_B, comparison_name) {
  df_imp %>%
    transmute(
      Protein.Id, Genes, Stripped.Sequence, Residue.Both,
      comparison = comparison_name,
      sample_A = sample_A,
      sample_B = sample_B,
      log2FC = .data[[sample_B]] - .data[[sample_A]]
    )
}

# -----------------------------
# 6) Map group -> sample_id
#    (your "样品信息" column stores group name)
# -----------------------------
group2sample <- sample_metadata %>%
  dplyr::select(sample_id, 样品信息) %>%
  distinct()

get_sample <- function(group_name) {
  sid <- group2sample %>% filter(样品信息 == group_name) %>% pull(sample_id)
  if (length(sid) != 1) stop(paste0("Group '", group_name, "' not uniquely mapped to a single sample_id."))
  sid
}

s_control <- get_sample("Control")
s_lps     <- get_sample("LPS")
s_lpsikk  <- get_sample("LPS+IKK")
s_ml162   <- get_sample("ML162")

# -----------------------------
# 7) Compute all requested comparisons
# -----------------------------
res_list <- list(
  LPS_Control     = compute_log2fc(mat_imp, s_control, s_lps,    "LPS vs Control"),
  LPSIKK_Control  = compute_log2fc(mat_imp, s_control, s_lpsikk, "LPS+IKK vs Control"),
  ML162_Control   = compute_log2fc(mat_imp, s_control, s_ml162,  "ML162 vs Control"),
  LPSIKK_LPS      = compute_log2fc(mat_imp, s_lps,     s_lpsikk, "LPS+IKK vs LPS"),
  ML162_LPS       = compute_log2fc(mat_imp, s_lps,     s_ml162,  "ML162 vs LPS"),
  ML162_LPSIKK    = compute_log2fc(mat_imp, s_lpsikk,  s_ml162,  "ML162 vs LPS+IKK")
)

names(res_list)


res_fc_long <- bind_rows(res_list)

# -----------------------------
# 8) Wide output: one row per site/peptide, multiple comparison columns
# -----------------------------
res_fc_wide <- res_fc_long %>%
  select(-sample_A, -sample_B) %>%
  pivot_wider(names_from = comparison, values_from = log2FC)

# -----------------------------
# 9) Optional: save results
# -----------------------------
# write.csv(res_fc_long, "log2FC_results_long.csv", row.names = FALSE)
# write.csv(res_fc_wide, "log2FC_results_wide.csv", row.names = FALSE)

# -----------------------------
# 10) Quick sanity checks / previews
# -----------------------------
print(head(res_fc_long, 10))
print(head(res_fc_wide, 5))
