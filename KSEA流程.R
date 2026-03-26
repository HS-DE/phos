# ============================================================
# PhosR phosphoproteomics pipeline (mouse) for Kinase activity (KSEA-like)
# Input: site-level table (推荐 site_mat2) + fasta_map + sample_info
# Output: processed matrices + limma DE + kinase activity matrix + heatmaps
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  
  library(PhosR)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(limma)
  
  # 画热图用（你也可以只用 PhosR 自带 heatmap）
  library(pheatmap)
})

# ----------------------------
# Step 0. 用户需要改的参数
# ----------------------------

out_dir <- "C:/Work/SH/浙大磷酸化4例_小鼠/Results_2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1) 选择你要作为 PhosR 输入的表：
# 推荐：site_mat2（应包含 abs_pos / abs_residue / Protein.Id / Genes + 样本列）
# 这里假设你环境里已经有 site_mat2
phos_df <- site_mat2

# 2) 样本列（你前面已经算过 sample_cols，这里直接复用）
# sample_cols <- setdiff(names(expr), meta_cols)

# 3) 样本分组信息：你必须给出每个 sample 属于哪个 group
#    ——如果你的列名有规律，也可以用正则自动提取
sample_info <- tibble(
  sample = sample_cols,
  # TODO: 改成你真实的组名提取规则
  # 例：如果你列名像 "CTRL_1" "CTRL_2" "TRT_1" ...
  group  = gsub("_[0-9]+$", "", sample_cols)
)

# 4) 指定对照组（用于把 intensity 转成 log2FC/ratio）
#    ——KSEA/激酶活性一般在 ratio 上更合理
control_group <- unique(sample_info$group)[1]  # TODO: 改成你的真实对照组名（如 "CTRL"）

# 5) 是否丢弃多位点肽段（强烈建议：先只做单一位点，结果更干净）
keep_only_single_site <- TRUE


# ----------------------------
# Step 1. 清洗 & 构造“单一位点”行（PhosR 更偏向 site-centric）
# ----------------------------

# (1) 只保留必要列：Protein.Id / Genes / abs_pos / abs_residue / sample_cols
need_cols <- c("Protein.Id", "Genes", "abs_pos", "abs_residue", sample_cols)
miss_cols <- setdiff(need_cols, names(phos_df))
if (length(miss_cols) > 0) {
  stop("phos_df 缺少这些列，请用 site_mat2 作为输入，缺少列：\n", paste(miss_cols, collapse = ", "))
}

phos_df <- phos_df %>%
  mutate(
    GeneSymbol = toupper(Genes),                 # 为了更好匹配 PhosphoSitePlus（常用大写）
    Residue    = as.character(abs_residue),
    Site       = as.numeric(abs_pos)
  )

# (2) 去掉 NA Site / Residue
phos_df <- phos_df %>%
  filter(!is.na(Site), !is.na(Residue), Residue %in% c("S","T","Y"))

# (3) 如果你已经有 Residue.Both / Modified.Sequence，可用它判断多位点
if (keep_only_single_site && "Residue.Both" %in% names(phos_df)) {
  phos_df <- phos_df %>%
    mutate(n_site_in_row = str_count(Residue.Both, ";") + 1) %>%
    filter(n_site_in_row == 1) %>%
    select(-n_site_in_row)
}

# (4) 构造 PhosR 常用的 site label： "GENE;S123;"
#     这个格式在 PhosR 的很多函数/示例里会用到
phos_df <- phos_df %>%
  mutate(site_label = paste0(GeneSymbol, ";", Residue, Site, ";"))


# ----------------------------
# Step 2. 生成 15-mer 序列（-7, p, +7），供 kinaseSubstrateScore 用
# ----------------------------
# kinaseSubstrateScore() 要求每个位点一个 15aa 序列（中心是被磷酸化的位点） :contentReference[oaicite:4]{index=4}

get_flank15 <- function(prot_seq, pos, win = 7, pad = "X") {
  # prot_seq: 蛋白全长序列
  # pos: 1-based 位点位置（abs_pos）
  # 返回长度=15 的序列，中心第 win+1 位为该位点
  if (is.na(prot_seq) || is.na(pos)) return(NA_character_)
  L <- nchar(prot_seq)
  if (pos < 1 || pos > L) return(NA_character_)
  
  st <- pos - win
  ed <- pos + win
  
  left_pad  <- if (st < 1) paste(rep(pad, 1 - st), collapse = "") else ""
  right_pad <- if (ed > L) paste(rep(pad, ed - L), collapse = "") else ""
  
  st2 <- max(1, st)
  ed2 <- min(L, ed)
  
  core <- substr(prot_seq, st2, ed2)
  seq15 <- paste0(left_pad, core, right_pad)
  
  if (nchar(seq15) != (2 * win + 1)) return(NA_character_)
  seq15
}

# 从 fasta_map 按 Protein.Id 取序列（你前面读出来的 fasta_map 是 setNames(seqs, acc)）
phos_df <- phos_df %>%
  mutate(
    prot_seq = map_chr(Protein.Id, ~{
      if (.x %in% names(fasta_map)) fasta_map[[.x]][[1]] else NA_character_
    }),
    Sequence15 = map2_chr(prot_seq, Site, ~ get_flank15(.x, .y, win = 7))
  ) %>%
  select(-prot_seq)

# 丢弃没有 15-mer 的位点（否则 kinaseSubstrateScore 会很痛苦）
phos_df <- phos_df %>% filter(!is.na(Sequence15))


# ----------------------------
# Step 3. 构造 quant matrix，并把 intensity -> log2 ratio（相对 control）
# ----------------------------

quant_mat <- phos_df %>%
  select(all_of(sample_cols)) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

# 行名用 site_label，保证与 PhosR 注释/匹配一致
rownames(quant_mat) <- phos_df$site_label
colnames(quant_mat) <- sample_cols

# (1) log2 transform（加一个很小的伪计数，避免 0）
quant_log2 <- log2(quant_mat + 1)

# 用所有样本的行中位数当参照（更稳健）
row_ref <- apply(quant_log2, 1, median, na.rm = TRUE)
quant_ratio <- sweep(quant_log2, 1, row_ref, "-")

if(F) {# 以第一个样本作为对照
# (2) 计算 ratio：每个位点减去 control 组均值 => log2FC
grp_vec <- sample_info$group[match(colnames(quant_log2), sample_info$sample)]
if (anyNA(grp_vec)) stop("sample_info$sample 与 sample_cols 对不上，请检查 sample_info")

ctl_cols <- grp_vec == control_group
if (!any(ctl_cols)) stop("control_group 在 sample_info$group 中不存在：", control_group)

ctl_mean <- rowMeans(quant_log2[, ctl_cols, drop = FALSE], na.rm = TRUE)
quant_ratio <- sweep(quant_log2, 1, ctl_mean, FUN = "-")
}
# 保存一下原始 ratio 矩阵
write.csv(
  data.frame(site_label = rownames(quant_ratio), quant_ratio, check.names = FALSE),
  file = file.path(out_dir, "01_quant_log2ratio_raw.csv"),
  row.names = FALSE
)


# ----------------------------
# Step 4. 构造 PhosphoExperiment 对象（PhosR 的核心容器）
# ----------------------------
# PhosphoExperiment(...) 可同时放 assay + UniprotID/GeneSymbol/Site/Residue/Sequence 等注释 :contentReference[oaicite:5]{index=5}

coldata <- S4Vectors::DataFrame(
  group = factor(grp_vec),
  row.names = colnames(quant_ratio)
)

ppe <- PhosR::PhosphoExperiment(
  assays    = list(Quantification = quant_ratio),
  colData   = coldata,
  UniprotID = phos_df$Protein.Id,
  GeneSymbol= phos_df$GeneSymbol,
  Site      = phos_df$Site,
  Residue   = phos_df$Residue,
  Sequence  = phos_df$Sequence15
)

saveRDS(ppe, file.path(out_dir, "02_ppe_raw.rds"))


# ----------------------------
# Step 5. 过滤 + 插补 + median scaling（PhosR 常规预处理）
# ----------------------------

# (A) 过滤：至少在某组内 >=70% replicate 被定量，并且在 >=1 个组满足
# 让同一组的重复共享同一个 grps 标签
grps <- sub("-[0-9]+$", "", colnames(ppe))   # 去掉末尾 -数字
table(grps)                                  # 看看每组有几个重复（必须>=2 才适合 selectGrps）

ppe_filt <- PhosR::selectGrps(ppe, grps, percent = 0.7, n = 1, assay = "Quantification")

# selectGrps 的定义：按组内定量率筛位点 :contentReference[oaicite:6]{index=6}

mat_filt <- SummarizedExperiment::assay(ppe_filt, "Quantification")

# (B) scImpute：组内缺失（site- & condition-specific） :contentReference[oaicite:7]{index=7}
set.seed(123)
mat_sc <- PhosR::scImpute(mat_filt, percent = 0.7, grps)

# (C) tImpute：尾部分布插补（Perseus 风格） :contentReference[oaicite:8]{index=8}
set.seed(123)
mat_t <- PhosR::tImpute(mat_sc)

# (D) medianScaling：样本中位数中心化（可选 scale=FALSE） :contentReference[oaicite:9]{index=9}
mat_scaled <- PhosR::medianScaling(mat_t, scale = FALSE)

# 写回 PhosphoExperiment
SummarizedExperiment::assay(ppe_filt, "imputed") <- mat_t
SummarizedExperiment::assay(ppe_filt, "scaled")  <- mat_scaled

saveRDS(ppe_filt, file.path(out_dir, "03_ppe_scaled.rds"))

# 输出处理后的矩阵
write.csv(
  data.frame(site_label = rownames(mat_scaled), mat_scaled, check.names = FALSE),
  file = file.path(out_dir, "03_quant_scaled.csv"),
  row.names = FALSE
)

# ----------------------------
# Step 6. （可选）RUVphospho 去批次/去不想要的变异
# ----------------------------
# RUVphospho 是 RUVIII 的 wrapper，需要 design matrix + stable sites(ctl) :contentReference[oaicite:10]{index=10}

do_ruv <- FALSE  # TODO: 如你有明显批次效应再开 TRUE
mat_norm <- mat_scaled

if (do_ruv) {
  data("SPSs", package = "PhosR")  # 稳定磷酸化位点列表
  ctl_idx <- which(rownames(mat_scaled) %in% SPSs)
  
  design <- model.matrix(~ group - 1, data = as.data.frame(coldata))
  
  mat_norm <- PhosR::RUVphospho(mat_scaled, M = design, ctl = ctl_idx, k = 3)
  
  SummarizedExperiment::assay(ppe_filt, "normalised") <- mat_norm
  saveRDS(ppe_filt, file.path(out_dir, "04_ppe_normalised.rds"))
  write.csv(
    data.frame(site_label = rownames(mat_norm), mat_norm, check.names = FALSE),
    file = file.path(out_dir, "04_quant_normalised.csv"),
    row.names = FALSE
  )
}

# ============================================================
# Step 7-9 (REPLACE ALL) 适配：组名非法 + 每组1个样本（只看FC）
# 依赖：mat_norm, ppe_filt, out_dir 已存在
# ============================================================

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(PhosR)
  library(SummarizedExperiment)
  library(pheatmap)
})

# ----------------------------
# Step 7. 差异磷酸化：有重复用 limma；无重复只输出 logFC
# ----------------------------

# (0) 从 ppe_filt 的 colData 获取 group（比你手写 grp_vec 更稳）
grp_vec_raw <- as.character(SummarizedExperiment::colData(ppe_filt)$group)
names(grp_vec_raw) <- colnames(mat_norm)

if (length(grp_vec_raw) != ncol(mat_norm)) {
  stop("grp_vec_raw 长度与 mat_norm 列数不一致，请检查 ppe_filt 与 mat_norm 是否同一批样本。")
}

# (1) 让组名变成 R 合法变量名（关键：解决 P-Lvs-... 这种含 '-' 的错误）
grp_fac <- factor(grp_vec_raw)
old_lvls <- levels(grp_fac)
new_lvls <- make.names(old_lvls)            # '-' 会变成 '.'
levels(grp_fac) <- new_lvls

# 记录原名->新名，便于你看结果
lvl_map <- data.frame(group_raw = old_lvls, group_safe = new_lvls)
write.csv(lvl_map, file.path(out_dir, "05_group_name_map.csv"), row.names = FALSE)

# (2) 判断是否“有重复”
tab_n <- table(grp_fac)
has_replicate <- any(tab_n >= 2)

# (3) 自动生成所有两两对比（基于 safe 名字）
grp_levels <- levels(grp_fac)

# 组合所有 pair： A_vs_B 与公式 A-B
pair_mat <- combn(grp_levels, 2)  # 每列是一对：c(g1, g2)
contrast_formula <- apply(pair_mat, 2, function(x) paste0(x[2], "-", x[1]))
contrast_names   <- apply(pair_mat, 2, function(x) paste0(x[2], "_vs_", x[1]))

de_list <- setNames(vector("list", length(contrast_names)), contrast_names)

# --------- 情况 A：有重复 -> limma 正常跑 ---------
if (has_replicate) {
  
  design <- model.matrix(~ 0 + grp_fac)
  colnames(design) <- levels(grp_fac)
  
  fit <- limma::lmFit(mat_norm, design)
  cont_mat <- limma::makeContrasts(contrasts = contrast_formula, levels = design)
  
  # 给 contrasts 设置名字（可读）
  colnames(cont_mat) <- contrast_names
  
  fit2 <- limma::contrasts.fit(fit, cont_mat)
  fit2 <- limma::eBayes(fit2)
  
  for (cn in contrast_names) {
    tt <- limma::topTable(fit2, coef = cn, number = Inf, sort.by = "P")
    tt$site_label <- rownames(tt)
    tt$contrast <- cn
    de_list[[cn]] <- tt
    
    write.csv(tt, file.path(out_dir, paste0("DE_", cn, ".csv")),
              row.names = FALSE)
  }
  
  saveRDS(de_list, file.path(out_dir, "05_DE_list.rds"))
  
} else {
  
  # --------- 情况 B：无重复 -> 只算 logFC(=样本差值)，P/FDR 全部 NA ---------
  message("检测到所有组都没有重复（每组=1），将仅输出 logFC（P.Value/adj.P.Val 设为 NA）。")
  
  # 为了避免误解，这里把 “每个组” 当作一个样本（因为组内只有1列）
  for (k in seq_along(contrast_names)) {
    cn <- contrast_names[k]
    
    g1 <- pair_mat[1, k]   # baseline
    g2 <- pair_mat[2, k]   # numerator
    
    # 找到对应列（每组只有1列，取第一列即可）
    col1 <- which(grp_fac == g1)[1]
    col2 <- which(grp_fac == g2)[1]
    
    logFC <- mat_norm[, col2] - mat_norm[, col1]
    
    tt <- data.frame(
      logFC = as.numeric(logFC),
      AveExpr = rowMeans(mat_norm, na.rm = TRUE),
      t = NA_real_,
      P.Value = NA_real_,
      adj.P.Val = NA_real_,
      B = NA_real_,
      site_label = rownames(mat_norm),
      contrast = cn,
      stringsAsFactors = FALSE
    )
    
    de_list[[cn]] <- tt
    write.csv(tt, file.path(out_dir, paste0("DE_", cn, ".csv")),
              row.names = FALSE)
  }
  
  saveRDS(de_list, file.path(out_dir, "05_DE_list.rds"))
}

# ----------------------------
# Step 8. 激酶活性 / KSEA-like（kinaseSubstrateScore）
# 关键修改：不再用 matANOVA/padj 筛位点（因为你无重复）
# 改为：按“跨样本变异/变化幅度”选 topN 位点，再 standardise
# ----------------------------

data("PhosphoSitePlus") 

# (1) 选“变化大的位点”
#     推荐用：行标准差（跨样本）作为“动态强度”
row_sd <- apply(mat_norm, 1, sd, na.rm = TRUE)

# 同时过滤：至少在 min_non_na 个样本中非NA（避免全缺失导致 sd=NA）
min_non_na <- max(2, floor(0.5 * ncol(mat_norm)))
non_na_n <- rowSums(!is.na(mat_norm))
ok_keep <- non_na_n >= min_non_na

row_sd2 <- row_sd
row_sd2[!ok_keep] <- NA_real_

topN <- 5000# 原本用的2000  # 你可调：1000~4000 常用
idx <- order(row_sd2, decreasing = TRUE, na.last = NA)
idx <- idx[seq_len(min(topN, length(idx)))]

mat_reg <- mat_norm[idx, , drop = FALSE]

# (2) standardise（行标准化）
mat_std <- PhosR::standardise(mat_reg)

# (3) 对齐 seqs（保证 rownames 一致）
# 注意：Sequence(ppe_filt) 返回每个位点的 15-mer
seq_all <- PhosR::Sequence(ppe_filt)
# 把 seq_all 命名成和 ppe_filt assay 的行名一致
names(seq_all) <- rownames(SummarizedExperiment::assay(ppe_filt, "Quantification"))

seq_reg <- seq_all[rownames(mat_reg)]

# 若有 NA 的序列（没对齐上），就去掉对应行，避免 kinaseSubstrateScore 报错
keep_seq <- !is.na(seq_reg) & nzchar(seq_reg)
mat_std2 <- mat_std[keep_seq, , drop = FALSE]
seq_reg2 <- seq_reg[keep_seq]

# (4) KSEA-like 打分 + 激酶活性矩阵
KSR <- PhosR::kinaseSubstrateScore(
  substrate.list = PhosphoSite.mouse,
  mat  = mat_std2,
  seqs = seq_reg2,
  numMotif = 5,
  numSub   = 1,
  species  = "mouse",
  verbose  = TRUE
)
ksActivityMatrix <- as.data.frame(KSR[["ksActivityMatrix"]])
saveRDS(KSR, file.path(out_dir, "06_KSR_kinaseSubstrateScore.rds"))

# (5) PhosR 自带热图
pdf(file.path(out_dir, "06_kinaseSubstrateHeatmap.pdf"), width = 10, height = 8)
PhosR::kinaseSubstrateHeatmap(KSR)
dev.off()
png(file.path(out_dir, "06_kinaseSubstrateHeatmap.png"), width = 10, height = 8, units = "in", res = 600)
PhosR::kinaseSubstrateHeatmap(KSR)
dev.off()

# (6) 提取 kinase activity matrix
ks_act <- KSR[["ksActivityMatrix"]]
write.csv(
  data.frame(Kinase = rownames(ks_act), ks_act, check.names = FALSE),
  file = file.path(out_dir, "06_kinase_activity_matrix.csv"),
  row.names = FALSE
)

pdf(file.path(out_dir, "06_kinase_activity_heatmap_pheatmap.pdf"), width = 10, height = 10)
pheatmap::pheatmap(ks_act, fontsize_row = 6, fontsize_col = 8)
dev.off()
png(file.path(out_dir, "06_kinase_activity_heatmap_pheatmap.png"), width = 10, height = 10, units = "in", res = 600)
pheatmap::pheatmap(ks_act, fontsize_row = 6, fontsize_col = 8)
dev.off()

data("PhosphoSitePlus") 

kinase_anno_tbl <- tibble::enframe(PhosphoSite.mouse, name = "Kinase", value = "Site") %>%
  tidyr::unnest_longer(Site)

# 看看长什么样
head(kinase_anno_tbl)

# 导出(这个是注释表)
write.csv(kinase_anno_tbl,"./PhosphoSite_mouse_kinase_substrate_annotation.csv",
          row.names = FALSE)
# Kinase: 只首字母大写，其余小写（保留数字）
kinase_anno_tbl2 <- kinase_anno_tbl %>%
  mutate(
    Kinase = paste0(toupper(substr(Kinase, 1, 1)), tolower(substr(Kinase, 2, nchar(Kinase))))
  )
tpm_data2 <- tpm_data %>% dplyr::filter(Genes %in% unique(kinase_anno_tbl2$Kinase))
write.csv(tpm_data2,"../Results/表格/激酶肽段定量表.csv", row.names = F)

data2 <- data %>% dplyr::filter(Genes %in% unique(kinase_anno_tbl2$Kinase))
write.csv(data2,"../Results/表格/激酶蛋白定量表.csv", row.names = F)


# (7) 可选：预测 kinase-substrate network
predMat <- PhosR::kinaseSubstratePred(KSR, top = 30)
saveRDS(predMat, file.path(out_dir, "06_predMat_top30.rds"))


# ----------------------------
# Step 9. rank-based KSEA（每个 contrast 一张表）
# 这里完全可以用 Step7 的 logFC（无重复也OK）
# ----------------------------

ksea_rank_dir <- file.path(out_dir, "07_rank_based_KSEA")
dir.create(ksea_rank_dir, showWarnings = FALSE, recursive = TRUE)

for (cn in names(de_list)) {
  tt <- de_list[[cn]]
  
  geneStats <- tt$logFC
  names(geneStats) <- tt$site_label
  
  # 去掉 NA，避免 enrichment 内部出错
  geneStats <- geneStats[!is.na(geneStats)]
  
  enr_up <- PhosR::pathwayRankBasedEnrichment(
    geneStats, annotation = PhosphoSite.mouse, alter = "greater"
  )
  enr_dn <- PhosR::pathwayRankBasedEnrichment(
    geneStats, annotation = PhosphoSite.mouse, alter = "less"
  )
  
  write.csv(enr_up, file.path(ksea_rank_dir, paste0("KSEA_rank_UP_", cn, ".csv")), row.names = FALSE)
  write.csv(enr_dn, file.path(ksea_rank_dir, paste0("KSEA_rank_DN_", cn, ".csv")), row.names = FALSE)
}

message("All done! Outputs in: ", out_dir)




# ============================================================
# A1. 由 ks_act 生成：每个 contrast 的 kinase activity 差值矩阵（Δ = B - A）
# ============================================================

# ks_act: 行=Kinase, 列=sample（你已经有了）
stopifnot(is.matrix(ks_act) || is.data.frame(ks_act))

ks_act <- as.matrix(ks_act)

# 生成所有两两对比（B_vs_A）
pairs <- combn(colnames(ks_act), 2, simplify = FALSE)

delta_mat <- do.call(
  cbind,
  lapply(pairs, function(p) {
    # p[1] = A, p[2] = B
    ks_act[, p[2]] - ks_act[, p[1]]
  })
)

colnames(delta_mat) <- vapply(pairs, function(p) paste0(p[2], "_vs_", p[1]), character(1))
rownames(delta_mat) <- rownames(ks_act)

# 导出
write.csv(
  data.frame(Kinase = rownames(delta_mat), delta_mat, check.names = FALSE),
  file = file.path(out_dir, "06_kinase_activity_delta_all_contrasts.csv"),
  row.names = FALSE
)

# 对比热图（哪组 vs 哪组）——一眼看全局
pdf(file.path(out_dir, "06_kinase_activity_delta_heatmap.pdf"), width = 10, height = 10)
pheatmap::pheatmap(delta_mat, fontsize_row = 6, fontsize_col = 8)
dev.off()

png(file.path(out_dir, "06_kinase_activity_delta_heatmap.png"),
    width = 10, height = 10, units = "in", res = 300)
pheatmap::pheatmap(delta_mat, fontsize_row = 6, fontsize_col = 8)
dev.off()


# ============================================================
# A2. 针对单个 contrast 画 Top up/down kinase barplot
# ============================================================

contrast_pick <- colnames(delta_mat)[1]  # TODO: 换成你想看的对比

delta_vec <- delta_mat[, contrast_pick]
delta_vec <- delta_vec[!is.na(delta_vec)]

# 先拿上调：从大到小
top_n <- 20
up_kin <- names(sort(delta_vec, decreasing = TRUE))[seq_len(min(top_n, length(delta_vec)))]

# 下调：先把上调的 kinase 去掉，再从小到大取 top_n
delta_vec_down_pool <- delta_vec[setdiff(names(delta_vec), up_kin)]
dn_kin <- names(sort(delta_vec_down_pool, decreasing = FALSE))[seq_len(min(top_n, length(delta_vec_down_pool)))]

top_up <- delta_vec[up_kin]
top_dn <- delta_vec[dn_kin]

plot_df <- rbind(
  data.frame(Kinase = names(top_up), delta = as.numeric(top_up), direction = "Up"),
  data.frame(Kinase = names(top_dn), delta = as.numeric(top_dn), direction = "Down")
)

# levels 必须唯一：按 delta 从小到大排列
levs <- plot_df$Kinase[order(plot_df$delta)]
levs <- unique(as.character(levs))   # 关键：去重
plot_df$Kinase <- factor(plot_df$Kinase, levels = levs)

# 画图：优先 ggplot2；没有就 base
# ============================================================
# Batch: for each contrast in delta_mat, draw Top up/down barplots
# Outputs: PDF + PNG for every contrast
# ============================================================
delta_mat
plot_kinase_top_ud <- function(delta_mat,
                               out_dir,
                               top_n = 20,
                               abs_min = 0,          # 可设 0.05 之类过滤几乎不变的 kinase
                               use_ggplot = TRUE) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  has_gg <- requireNamespace("ggplot2", quietly = TRUE)
  if (use_ggplot && has_gg) {
    library(ggplot2)
  }
  
  for (contrast_pick in colnames(delta_mat)) {
    
    delta_vec <- delta_mat[, contrast_pick]
    delta_vec <- delta_vec[!is.na(delta_vec)]
    
    # 可选：过滤掉几乎不变的 kinase，避免上下调全是 0
    if (abs_min > 0) {
      delta_vec <- delta_vec[abs(delta_vec) >= abs_min]
    }
    
    # 如果过滤后没有足够 kinase，跳过
    if (length(delta_vec) < 2) next
    
    # Up: 最大的 top_n
    up_kin <- names(sort(delta_vec, decreasing = TRUE))[seq_len(min(top_n, length(delta_vec)))]
    
    # Down: 排除 up_kin 后，最小的 top_n
    down_pool <- delta_vec[setdiff(names(delta_vec), up_kin)]
    if (length(down_pool) == 0) next
    
    dn_kin <- names(sort(down_pool, decreasing = FALSE))[seq_len(min(top_n, length(down_pool)))]
    
    plot_df <- rbind(
      data.frame(Kinase = up_kin, delta = as.numeric(delta_vec[up_kin]), direction = "Up"),
      data.frame(Kinase = dn_kin, delta = as.numeric(delta_vec[dn_kin]), direction = "Down")
    )
    
    # factor levels 唯一 + 按 delta 排序
    levs <- unique(as.character(plot_df$Kinase[order(plot_df$delta)]))
    plot_df$Kinase <- factor(plot_df$Kinase, levels = levs)
    
    # 文件名里避免奇怪字符
    safe_name <- gsub("[^A-Za-z0-9_.-]", "_", contrast_pick)
    
    pdf_file <- file.path(out_dir, paste0("06_kinase_delta_top", top_n, "_", safe_name, ".pdf"))
    png_file <- file.path(out_dir, paste0("06_kinase_delta_top", top_n, "_", safe_name, ".png"))
    
    if (use_ggplot && has_gg) {
      
      p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Kinase, y = delta)) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = paste0("Kinase activity change: ", contrast_pick),
          x = NULL,
          y = "Δ kinase activity (B - A)"
        )
      
      ggplot2::ggsave(pdf_file, p, width = 10, height = 8)
      ggplot2::ggsave(png_file, p, width = 10, height = 8, dpi = 300)
      
    } else {
      
      pdf(pdf_file, width = 10, height = 8)
      par(mar = c(5, 10, 4, 2))
      barplot(plot_df$delta,
              names.arg = as.character(plot_df$Kinase),
              horiz = TRUE, las = 1,
              main = paste0("Kinase activity change: ", contrast_pick),
              xlab = "Δ kinase activity (B - A)")
      dev.off()
      
      png(png_file, width = 10, height = 8, units = "in", res = 300)
      par(mar = c(5, 10, 4, 2))
      barplot(plot_df$delta,
              names.arg = as.character(plot_df$Kinase),
              horiz = TRUE, las = 1,
              main = paste0("Kinase activity change: ", contrast_pick),
              xlab = "Δ kinase activity (B - A)")
      dev.off()
    }
  }
  
  message("Batch plots done. Output dir: ", out_dir)
}

# 运行：每个对比都会出图
plot_dir <- file.path(out_dir, "06_kinase_topUD_each_contrast")
plot_kinase_top_ud(delta_mat, plot_dir, top_n = 20, abs_min = 0)

#------------------------------------
tpm_data2
library(dplyr)
library(stringr)
ksea_tab <- read.csv("../Results/表格/磷酸化肽段表格(可以拿来做KSEA).csv") %>% 
  dplyr::select("Modified.Sequence", "Residue.Both") %>%
  distinct(Modified.Sequence, Residue.Both, .keep_all = TRUE)
tpm_data3 <- tpm_data2 %>% left_join(ksea_tab, by = "Modified.Sequence")
write.csv(tpm_data3, "../Results/表格/激酶肽段定量表(附加磷酸化位点信息).csv",row.names = F,
          na = "0")


