library(dplyr)
library(purrr)
library(stringr)
library(readr)

# =========================================================
# 1) 从 fasta_map 取蛋白序列（带一个小兜底：找不到 isoform 就去掉 -2/-3）
#    fasta_map: 你之前构建的 map，key = Uniprot ID，value = 序列字符串（或 list）
# =========================================================
get_prot_seq <- function(pid, fasta_map) {
  pid <- as.character(pid)
  
  # 你之前的 fasta_map 可能是 fasta_map[[id]][[1]] 结构
  fetch <- function(x) {
    if (is.null(x)) return(NA_character_)
    if (is.list(x)) return(as.character(x[[1]]))
    as.character(x)
  }
  
  s <- fetch(fasta_map[[pid]])
  if (!is.na(s)) return(s)
  
  # fallback：去掉 -2/-3 之类
  if (str_detect(pid, "-")) {
    pid_base <- sub("-.*$", "", pid)
    s2 <- fetch(fasta_map[[pid_base]])
    if (!is.na(s2)) return(s2)
  }
  NA_character_
}

# =========================================================
# 2) 生成 15aa 窗口（±7），超出边界用 "_" 补齐
# =========================================================
make_15mer <- function(prot_seq, abs_pos, pad_char = "_") {
  L <- nchar(prot_seq)
  left  <- abs_pos - 7
  right <- abs_pos + 7
  
  chars <- vapply(left:right, function(i) {
    if (i < 1 || i > L) pad_char else substr(prot_seq, i, i)
  }, character(1))
  
  paste0(chars, collapse = "")
}

# =========================================================
# 3) 生成 PSP 输入序列
#    format = "asterisk"  -> 中心位点写成 S* / T* / Y*
#    format = "central"   -> 不加星号，但中心必须是位点（长度必须奇数，这里固定15没问题）
# =========================================================
make_psp_seq <- function(window15, center_res, format = c("asterisk","central"),
                         use_lowercase_center = FALSE) {
  format <- match.arg(format)
  center_res <- as.character(center_res)
  if (use_lowercase_center) center_res <- tolower(center_res)
  
  # window15 长度固定 15，中心位置是第 8 个字符
  if (nchar(window15) != 15) stop("window15 必须是 15 aa")
  
  if (format == "central") {
    # 把中心强行改成你的位点氨基酸（避免因 isoform / 定位误差导致中心不是 S/T/Y）
    paste0(substr(window15, 1, 7), toupper(center_res), substr(window15, 9, 15))
  } else {
    # asterisk 格式：中心位点后加 "*"
    paste0(substr(window15, 1, 7), toupper(center_res), "*", substr(window15, 9, 15))
  }
}

# =========================================================
# 4) 批量生成（输入：site_mat_one 或你最终的 site 表）
#    需要列：Protein.Id, abs_pos, abs_residue, site_id（或你想用的id列）
# =========================================================
build_psp_input <- function(site_df, fasta_map,
                            id_col = "site_id",
                            pid_col = "Protein.Id",
                            pos_col = "abs_pos",
                            res_col = "abs_residue",
                            format = c("asterisk","central"),
                            use_lowercase_center = FALSE) {
  format <- match.arg(format)
  
  site_df %>%
    mutate(
      Protein.Id = as.character(.data[[pid_col]]),
      abs_pos    = as.integer(.data[[pos_col]]),
      abs_res    = as.character(.data[[res_col]]),
      prot_seq   = map_chr(Protein.Id, ~ get_prot_seq(.x, fasta_map)),
      window15   = if_else(!is.na(prot_seq), map2_chr(prot_seq, abs_pos, make_15mer), NA_character_),
      site_seq   = if_else(!is.na(window15),
                           map2_chr(window15, abs_res, ~ make_psp_seq(.x, .y, format = format,
                                                                      use_lowercase_center = use_lowercase_center)),
                           NA_character_)
    ) %>%
    # 输出 PSP 需要的两列：identifier(可选) + sequence(必须)
    transmute(
      identifier = .data[[id_col]],
      sequence   = site_seq
    ) %>%
    filter(!is.na(sequence))
}

# =======================
# 用法示例（你按需改对象名）
# =======================
# 假设你的位点表叫 site_mat_one，且已经有 Protein.Id / abs_pos / abs_residue / site_id
psp_df <- build_psp_input(
  site_df = site_mat_one,
  fasta_map = fasta_map,
  id_col = "site_id",
  format = "asterisk"      # 推荐：更直观
)

# 导出给 PhosphoSitePlus（tsv）
write_tsv(psp_df, "PSP_input_15aa.tsv")

# 简单检查：每条序列（去掉*）长度应该是15
table(nchar(gsub("\\*", "", psp_df$sequence)))

