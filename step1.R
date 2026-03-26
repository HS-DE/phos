setwd("C:/Work/SH/浙大磷酸化4例_小鼠/3.25售后/R/")
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)

read_fasta_as_map <- function(fasta_file) {
  x <- readLines(fasta_file, warn = FALSE)
  
  hdr_idx <- which(startsWith(x, ">"))
  end_idx <- c(hdr_idx[-1] - 1, length(x))
  
  headers <- x[hdr_idx]
  seqs <- map2(hdr_idx, end_idx, ~ paste0(x[(.x + 1):.y], collapse = ""))
  
  # 你的header形如：>sp|A0A087X1C5|CP2D7_HUMAN ...
  acc <- vapply(headers, function(h) {
    h <- sub("^>", "", h)
    parts <- strsplit(h, "\\|")[[1]]
    if (length(parts) >= 2) parts[2] else strsplit(h, " ")[[1]][1]
  }, character(1))
  
  setNames(seqs, acc)
}

fasta_file <- "C:/Work/SH/Pub_database/phos/uniprot-download_true_format_fasta_20422_includeIsoform_true_query__28homo-2023.03.17-05.33.29.22.fasta"
fasta_map <- read_fasta_as_map(fasta_file)


get_phos_sites_in_peptide <- function(modseq) {
  n <- nchar(modseq)
  i <- 1
  aa_idx <- 0
  last_aa <- NA_character_
  last_pos <- NA_integer_
  
  pos <- integer(0)
  res <- character(0)
  
  while (i <= n) {
    ch <- substr(modseq, i, i)
    
    # 氨基酸（大写字母）
    if (grepl("^[A-Z]$", ch)) {
      aa_idx <- aa_idx + 1
      last_aa <- ch
      last_pos <- aa_idx
      i <- i + 1
      next
    }
    
    # 修饰括号：跳到 )，若是 UniMod:21 记录位点
    if (ch == "(") {
      j <- regexpr("\\)", substr(modseq, i, n))[1]
      if (j == -1) break
      tag <- substr(modseq, i, i + j - 1)  # 包含括号
      
      if (tag == "(UniMod:21)") {
        pos <- c(pos, last_pos)
        res <- c(res, last_aa)
      }
      
      i <- i + j
      next
    }
    
    # 其它字符（一般不会有），跳过
    i <- i + 1
  }
  
  tibble(pep_pos = pos, residue = res)
}

find_all_matches <- function(prot_seq, pep, il_equiv = FALSE) {
  if (il_equiv) {
    prot_seq <- chartr("I", "L", prot_seq)
    pep <- chartr("I", "L", pep)
  }
  hits <- gregexpr(pep, prot_seq, fixed = TRUE)[[1]]
  if (length(hits) == 1 && hits[1] == -1) integer(0) else as.integer(hits)
}

map_one_peptide <- function(protein_ids, stripped, modified, fasta_map, il_equiv = FALSE) {
  # 取肽段内磷酸化位点
  sites <- get_phos_sites_in_peptide(modified)
  if (nrow(sites) == 0) return(tibble())
  
  ids <- strsplit(protein_ids, ";")[[1]] |> trimws()
  out <- list()
  
  for (pid in ids) {
    pid2 <- pid
    if (!pid2 %in% names(fasta_map) && grepl("-", pid2)) {
      # 若FASTA里没isoform，尝试去掉 -2
      pid_base <- sub("-.*$", "", pid2)
      if (pid_base %in% names(fasta_map)) pid2 <- pid_base
    }
    if (!pid2 %in% names(fasta_map)) next
    
    prot_seq <- fasta_map[[pid2]][[1]]
    starts <- find_all_matches(prot_seq, stripped, il_equiv = il_equiv)
    if (length(starts) == 0) next
    
    for (st in starts) {
      abs_pos <- st + sites$pep_pos - 1
      abs_res <- substr(prot_seq, abs_pos, abs_pos)
      
      out[[length(out) + 1]] <- tibble(
        Protein.Id = pid2,
        Stripped.Sequence = stripped,
        Modified.Sequence = modified,
        pep_start = st,
        pep_end = st + nchar(stripped) - 1,
        pep_pos = sites$pep_pos,
        residue = sites$residue,
        abs_pos = abs_pos,
        abs_residue = sites$residue,
        site_id = paste0(pid2, "_", sites$residue, abs_pos)
      )
    }
  }
  
  if (length(out) == 0) tibble() else bind_rows(out)
}


# 自动识别样本列：把前面注释列排除掉
meta_cols <- c("Protein.Group","Protein.Ids","Protein.Names","Genes",
               "First.Protein.Description","Proteotypic",
               "Stripped.Sequence","Modified.Sequence",
               "Precursor.Charge","Precursor.Id")
sample_cols <- setdiff(names(expr), meta_cols)

# 只保留磷酸化肽段 + 至少一个样本检出
df_phos <- expr %>%
  filter(str_detect(Modified.Sequence, "\\(UniMod:21\\)")) %>%
  filter(if_any(all_of(sample_cols), ~ !is.na(.) & . > 0)) %>%
  distinct(Protein.Ids, Stripped.Sequence, Modified.Sequence)  # 同一肽段不同charge只算一次

# 批量映射（生成 site-level 长表）
# df_phos: 你那三列（Protein.Ids, Stripped.Sequence, Modified.Sequence）

# 1) 生成映射结果（site-level 长表）
site_map <- pmap_dfr(
  list(df_phos$Protein.Ids, df_phos$Stripped.Sequence, df_phos$Modified.Sequence),
  ~ map_one_peptide(..1, ..2, ..3, fasta_map, il_equiv = FALSE)
)

site_map2 <- site_map %>%
  mutate(.rowid = seq_len(n())) %>%                         # 记录原始顺序（保证“最上面”）
  group_by(Modified.Sequence) %>%
  mutate(n_phos = str_count(Modified.Sequence, "\\(UniMod:21\\)")) %>%  # 该序列磷酸化数
  arrange(.rowid, .by_group = TRUE) %>%                     # 按原顺序排
  filter(row_number() <= n_phos) %>%                        # 只保留前 n_phos 行
  ungroup() %>%
  select(-.rowid, -n_phos)

# 检查：处理后每个 Modified.Sequence 的行数是否等于 UniMod:21 个数
check2 <- site_map2 %>%
  group_by(Modified.Sequence) %>%
  summarise(
    n_phos = str_count(Modified.Sequence, "\\(UniMod:21\\)"),
    n_rows = n(),
    ok = (n_phos == n_rows),
    .groups = "drop"
  )

table(check2$ok)
check2 %>% filter(!ok)



# 注意：site_map 里蛋白列叫 Protein.Id（单个 accession），不是 Protein.Ids（可能多个）
# 2) 合并回 df_phos（按三列键）
site_table <- df_phos %>%
  left_join(
    site_map2 %>% dplyr::select(Stripped.Sequence, Modified.Sequence, Protein.Id, pep_start, pep_end,
                        pep_pos, residue, abs_pos, abs_residue, site_id),
    by = c("Stripped.Sequence", "Modified.Sequence")
  )

# 1) precursor-level -> long
long_int <- expr %>%
  dplyr::select(Protein.Ids, Stripped.Sequence, Modified.Sequence,
         Precursor.Id, Precursor.Charge, all_of(sample_cols)) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "intensity"
  ) %>%
  mutate(intensity = as.numeric(intensity))

# 2) collapse 到 “peptide-level”（同一肽段同一样本只保留 1 个强度）
pep_int <- long_int %>%
  group_by(sample, Protein.Ids, Stripped.Sequence, Modified.Sequence) %>%
  summarise(
    intensity = if (all(is.na(intensity))) NA_real_ else max(intensity, na.rm = TRUE),
    n_precursors = sum(!is.na(intensity)),     # 这个能让你看到该肽段在该样本由多少条precursor贡献
    .groups = "drop"
  )

pep_mat <- pep_int %>%
  dplyr::select(Protein.Ids, Stripped.Sequence, Modified.Sequence, sample, intensity) %>%
  pivot_wider(names_from = sample, values_from = intensity)

head(pep_mat)

pep_mat2 <- pep_mat %>% filter(pep_mat$Modified.Sequence %in% site_table$Modified.Sequence)


site_table_with_int <- site_table %>%
  left_join(
    pep_int,
    by = c("Protein.Ids", "Stripped.Sequence", "Modified.Sequence")
  )

site_mat <- site_table_with_int %>%
  dplyr::select(Protein.Id, site_id, abs_pos, abs_residue,
                Stripped.Sequence, Modified.Sequence,
                sample, intensity) %>%
  dplyr::distinct(Protein.Id, site_id, abs_pos, abs_residue,
                  Stripped.Sequence, Modified.Sequence,
                  sample, .keep_all = TRUE) %>%
  tidyr::pivot_wider(names_from = sample, values_from = intensity)
site_mat <- expr %>% dplyr::select(Genes,Modified.Sequence) %>% distinct() %>% 
  left_join(site_mat, by = "Modified.Sequence") %>% filter(! is.na(.$Stripped.Sequence))
head(site_mat)

site_mat_one <- site_mat %>%
  group_by(Genes, Modified.Sequence) %>%
  mutate(
    is_canonical = !str_detect(Protein.Id, "-"),
    Protein.Id.keep = if (any(is_canonical)) Protein.Id[which(is_canonical)[1]] else Protein.Id[1]
  ) %>%
  filter(Protein.Id == Protein.Id.keep) %>%
  ungroup() %>%
  dplyr::select(-is_canonical, -Protein.Id.keep)

check_df <- site_mat_one %>%
  mutate(
    n_phos = str_count(Modified.Sequence, "\\(UniMod:21\\)")
  ) %>%
  group_by(Genes, Modified.Sequence) %>%
  summarise(
    n_phos = dplyr::first(n_phos),                  # 每个 Modified.Sequence 固定
    n_sites = n_distinct(site_id),           # 该序列对应多少个位点行
    n_rows = n(),                            # 实际行数（一般应等于 n_sites）
    .groups = "drop"
  ) %>%
  mutate(
    ok = (n_sites == n_phos)
  )

# 1) 总体概览
table(check_df$ok)

# 2) 看看哪些不一致
diff_df <- check_df %>%
  filter(!ok) %>%
  arrange(desc(abs(n_sites - n_phos)))

diff_df

residue_both_df <- site_mat %>%
  group_by(Modified.Sequence) %>%
  summarise(
    Residue.Both = paste0(
      unique(paste0(abs_residue, abs_pos)) %>% sort(),
      collapse = ";"
    ),
    .groups = "drop"
  )

site_mat2 <- site_mat %>%
  left_join(residue_both_df, by = "Modified.Sequence")



############ FC + P ##############
final_mat <- site_mat2 %>% dplyr::select(Protein.Id, Genes, Stripped.Sequence, 
                                         Modified.Sequence, Residue.Both, all_of(sample_cols))

# ============================================================
# Phospho-site/peptide matrix: log2 + lenient missing + log2FC
# 4 groups, 1 sample each
# ============================================================

library(dplyr)
library(tidyr)
sample_metadata <- readxl::read_xlsx("../data/sampleinfo.xlsx")
sample_metadata$sample_id <- sample_metadata$样品编号
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
meta_cols <- c("Protein.Id", "Genes", "Stripped.Sequence", "Residue.Both")
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
  dplyr::select(sample_id, group) %>%
  distinct()

get_sample <- function(group_name) {
  sid <- group2sample %>% filter(group == group_name) %>% pull(sample_id)
  if (length(sid) != 1) stop(paste0("Group '", group_name, "' not uniquely mapped to a single sample_id."))
  sid
}

s_control <- get_sample("TiO2")
s_lps     <- get_sample("MOAC")
s_lpsikk  <- get_sample("LPS+IKK")
s_ml162   <- get_sample("ML162")

# -----------------------------
# 7) Compute all requested comparisons
# -----------------------------
res_list <- list(
  compute_log2fc(mat_imp, s_control, s_lps,    "LPS vs Control"),
  compute_log2fc(mat_imp, s_control, s_lpsikk, "LPS+IKK vs Control"),
  compute_log2fc(mat_imp, s_control, s_ml162,  "ML162 vs Control"),
  compute_log2fc(mat_imp, s_lps,     s_lpsikk, "LPS+IKK vs LPS"),
  compute_log2fc(mat_imp, s_lps,     s_ml162,  "ML162 vs LPS"),
  compute_log2fc(mat_imp, s_lpsikk,  s_ml162,  "ML162 vs LPS+IKK")
)

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


