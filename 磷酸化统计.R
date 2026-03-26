setwd("C:/Work/SH/浙大磷酸化4例_小鼠/3.25售后/R/")

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

df <- read.csv("../data/Phos_Test_report-first-pass.pr_matrix.tsv",
               sep = "\t", header = TRUE, check.names = FALSE,
               quote = "", comment.char = "")

# 先获取列名并进行处理
source("C:/Work/SH/code/source/质谱样本名处理.R")
names(df)
df <- extract_sample_names(expr = df, prefix_pattern = "20260119")
df$sample_cols
expr    <- df$expr
samples <- df$sample_cols   # ✅ 不要再用 names(df_dis[,-1:-10]) 这种写法

# （可选但推荐）确保样本列是 numeric，避免读入成字符导致 >0 出问题
expr <- expr %>% mutate(across(all_of(samples), as.numeric))

# ============================================================
# A) 统一构建“长表 + 检出过滤”（第二段逻辑的核心）
#    以后所有计数/分类都从这个 long_detect 出发
# ============================================================
long_detect <- expr %>%
  dplyr::select(Stripped.Sequence, Modified.Sequence, all_of(samples)) %>%
  pivot_longer(
    cols = all_of(samples),
    names_to = "sample",
    values_to = "intensity"
  ) %>%
  filter(!is.na(intensity) & intensity > 0) %>%           # ✅ 检出定义（not NA 且 >0）
  distinct(sample, Modified.Sequence, .keep_all = TRUE)   # ✅ 同一样本同一肽段多行（不同precursor）只算一次

# ============================================================
# B) 每个样本：总肽段数（按 Modified.Sequence 去重）
# ============================================================
total_peptides <- long_detect %>%
  distinct(sample, Modified.Sequence) %>%
  dplyr::count(sample, name = "Total_Peptides") %>%
  as.data.frame()

row.names(total_peptides) <- total_peptides$sample


# ============================================================
# C) 每个样本：磷酸化肽段数（含 UniMod:21）
# ============================================================
phos_peptides <- long_detect %>%
  filter(str_detect(Modified.Sequence, "\\(UniMod:21\\)")) %>%
  distinct(sample, Modified.Sequence) %>%
  dplyr::count(sample, name = "Phos_Peptides") %>%
  as.data.frame()

row.names(phos_peptides) <- phos_peptides$sample


# ============================================================
# D) 做“肽段注释表”（只用于解析序列/修饰数量，不用于计数）
#    这里 distinct(Modified.Sequence) 是合理的：只为得到每条 Modified.Sequence 的注释字段
# ============================================================
seq_anno <- expr %>%
  distinct(Modified.Sequence, .keep_all = TRUE) %>%
  select(Stripped.Sequence, Modified.Sequence)

# 只取磷酸化肽段做注释（⚠️ 精确匹配 (UniMod:21)）
df_phos_anno <- seq_anno %>%
  filter(str_detect(Modified.Sequence, "\\(UniMod:21\\)")) %>%   # ✅ 不要用 "UniMod:21"
  mutate(
    # 1) 精确计算 UniMod:21 的个数（最稳，绝不会把 LC(...)(UniMod:21) 算成 0）
    count = str_count(Modified.Sequence, "\\(UniMod:21\\)"),
    
    # 2) 提取 UniMod:21 之前的“骨架”（你原来那个 sub 会截断到第一个21之后，
    #    这里保留也行，但注意：有多磷酸化时它只保留第一个21之前）
    Modified.Sequence_UniMod_21 = sub("\\(UniMod:21\\).*", "", Modified.Sequence),
    
    # 3) 如果你后面还要用 split 来取位点残基（S/T/Y），保留 list-column（可选）
    Modified.Sequence_UniMod_list = map(
      Modified.Sequence,
      ~ strsplit(.x, "\\(UniMod:21\\)")[[1]]
    )
  )


# ============================================================
# E) 构建“每个样本检出的磷酸肽段清单”，并 join 注释（count 等）
#    之后你所有按样本的 count 分类/位置统计都从这里做
# ============================================================
df_phos_detect <- long_detect %>%
  filter(str_detect(Modified.Sequence, "\\(UniMod:21\\)")) %>%
  distinct(sample, Modified.Sequence) %>%                       # ✅ 每样本每肽段只算一次
  left_join(
    df_phos_anno %>% select(Modified.Sequence, Stripped.Sequence,
                            Modified.Sequence_UniMod_list, count),
    by = "Modified.Sequence"
  )

# ============================================================
# F) “磷酸化位点数（count）分类统计” ——替代你原来的 phos_list 循环
# ============================================================
count_by_sample <- df_phos_detect %>%
  dplyr::count(sample, count, name = "frequency") %>%
  arrange(sample, count)
count_by_sample <- df_phos_detect %>%
  dplyr::count(sample, count, name = "frequency") %>%
  arrange(sample, count)
count_by_sample

# 如果你想要和你原来 sample_count 那种“宽表矩阵”（列=1..max_count）
sample_count <- count_by_sample %>%
  tidyr::pivot_wider(
    names_from = count,
    values_from = frequency,
    values_fill = 0
  ) %>%
  as.data.frame()

# ============================================================
# G) “磷酸化修饰位置/氨基酸类型统计”
#    你原代码 unlist(x[-length(x)]) 只是拿到了片段，不是位点字母
#    正确做法：对每个分割片段（除最后一段）取最后一个字符（S/T/Y）
# ============================================================
get_phos_residues <- function(parts) {
  # parts: 例如 strsplit by "(UniMod:21)" 得到的向量
  if (length(parts) <= 1) return(character(0))
  pre <- parts[-length(parts)]                       # 去掉最后一段
  substr(pre, nchar(pre), nchar(pre))                # 取每段末尾那个氨基酸字母
}

df_phos_detect2 <- df_phos_detect %>%
  mutate(phos_residue = map(Modified.Sequence_UniMod_list, get_phos_residues))

# 每个样本：S/T/Y 频率
sample_loca_df <- df_phos_detect2 %>%
  tidyr::unnest(phos_residue) %>%
  dplyr::count(sample, phos_residue, name = "frequency") %>%
  arrange(sample, desc(frequency))

sample_loca_df


# 如果你想要和你原来 sample_count 那种“宽表矩阵”（列=1..max_count）
sample_loca_df_count <- sample_loca_df %>%
  tidyr::pivot_wider(
    names_from = phos_residue,
    values_from = frequency,
    values_fill = 0
  ) %>%
  as.data.frame()

# ============================================================
# H) “整合结果并输出”
# ============================================================
res1 <- left_join(sample_count, sample_loca_df_count, by = "sample")
res2 <- left_join(res1, total_peptides, by = "sample")
res3 <- left_join(res2, phos_peptides, by = "sample") %>% 
  mutate(
    `Phos%` = (Phos_Peptides / Total_Peptides)*100
  )
dput(names(res3))
res3 <- res3 %>% dplyr::select("sample", "Total_Peptides", "Phos_Peptides", 
                               "Phos%", "1", "2", "3", "S", "T", "Y")
#write.csv(res3, "../Results/磷酸化统计.csv", row.names = F)


# ============================================================
# I) “可视化：柱状图”
# ============================================================
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

# expr 是你读入后的矩阵数据框（包含 Protein.Group / Modified.Sequence / 各样本强度列）
# 自动识别样本列（把前面的注释列排除掉）
meta_cols <- c("Protein.Group","Genes","Modified.Sequence")
sample_cols <- df$sample_cols

# 1) 只取“真正的”磷酸化肽段行：精确匹配 (UniMod:21)
# 2) 可选：要求至少在任意一个样本里强度>0（更符合“检出”）
# 3) 同一 Protein.Group 内，对 Modified.Sequence 去重计数
df_protein_phos_n <- expr %>%
  filter(str_detect(Modified.Sequence, "\\(UniMod:21\\)")) %>%
  filter(if_any(all_of(sample_cols), ~ !is.na(.) & . > 0)) %>%   # 不想按强度筛就删掉这一行
  distinct(Protein.Group, Modified.Sequence) %>%                 # 同一肽段不同电荷/precursor 只算一次
  dplyr::count(Protein.Group, name = "n_phos_sites")                    # 这里“点位数”= unique(Modified.Sequence)数

# 分箱：1~15, 15+
df_plot <- df_protein_phos_n %>%
  mutate(bin = ifelse(n_phos_sites > 15, "15+", as.character(n_phos_sites))) %>%
  dplyr::count(bin, name = "n_protein") %>%
  mutate(bin = factor(bin, levels = c(as.character(1:15), "15+"))) %>%
  tidyr::complete(bin, fill = list(n_protein = 0))
if(F){
# 画柱状图
p <- ggplot(df_plot, aes(x = bin, y = n_protein)) +
  geom_col(fill = "#1496D4") +
  labs(
    title = "Phosphorylated Site Distribution",
    x = "Number of Phosphorylated Sites in a Protein",
    y = "The number of Phosphorylated Proteins"
  ) +
  theme_classic() +
  theme(
    # 标题居中 + 字体大小
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    # 坐标轴标题
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    # 坐标轴刻度
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14)
  )

p
ggsave("../Results/磷酸化修饰位点数量分布图.png", width = 8, height = 6)
}

################ 磷酸化蛋白 ##################
df_out <- df_phos_detect2
list_cols <- names(df_out)[sapply(df_out, is.list)]

df_out[list_cols] <- lapply(df_out[list_cols], function(x) {
  vapply(x, function(z) paste(z, collapse = ";"), character(1))
})

#write.csv(df_out, "./df_phos_detect2.csv", row.names = FALSE)
