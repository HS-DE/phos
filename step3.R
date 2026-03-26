library(dplyr)
library(stringr)

# 把 "LPS vs Control" -> "LPS_Control"，"LPS+IKK vs Control" -> "LPS_IKK_Control"
make_comp_name <- function(x) {
  x %>%
    str_replace_all("\\s*vs\\s*", "_") %>%          # vs -> _
    str_replace_all("\\+", "_") %>%                # + -> _
    str_replace_all("[^A-Za-z0-9_]", "_") %>%      # 其它非法字符 -> _
    str_replace_all("_+", "_") %>%                 # 压缩多个_
    str_replace_all("^_|_$", "")                   # 去首尾_
}

# 生成 KSEAapp 需要的 PX 格式（p 必须是 "NULL"，不能 NA）
make_px <- function(df_one_comp) {
  df_one_comp %>%
    transmute(
      Protein = Protein.Id,
      Gene    = Genes,
      Peptide = Stripped.Sequence,
      Residue.Both = str_replace_all(Residue.Both, "\\s*;\\s*", ";"),  # 去掉分号周围空格
      p  = "NULL",                       # ✅ 字符串，不是 NA
      FC = 2^log2FC                      # ✅ log2FC -> FC（B/A）
    ) %>%
    # KSEAapp 遇到 NA 会整行删除；这里我们直接剔除任何 NA/非有限值
    filter(
      !is.na(Protein), !is.na(Gene), !is.na(Peptide), !is.na(Residue.Both),
      !is.na(FC), is.finite(FC)
    ) %>%
    mutate(
      Protein = as.character(Protein),
      Gene = as.character(Gene),
      Peptide = as.character(Peptide),
      Residue.Both = as.character(Residue.Both),
      p = as.character(p)
    )
}

# 按 comparison 拆分
res_split <- split(res_fc_long, res_fc_long$comparison)

# 生成 PX 列表
px_list <- lapply(res_split, make_px)
names(px_list) <- make_comp_name(names(px_list))

# 放到全局变量里：LPS_Control, LPS_IKK_Control, ML162_LPS, ...
list2env(px_list, envir = .GlobalEnv)

# 例子：查看
head(LPS_Control)

# （可选）检查每个 comparison 有多少行被保留下来
sapply(px_list, nrow)

if(F) {
library(KSEAapp)

KSData_full <- read.csv("C:/Work/SH/Pub_database/phos/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv",
                        stringsAsFactors = FALSE)
# 一步到位：直接在工作目录输出 3 个文件（tiff 图 + 2 个 csv）
KSEA.Complete(KSData_full, PX = px_list$LPS_Control,
              NetworKIN = TRUE, NetworKIN.cutoff = 5,
              m.cutoff = 5, p.cutoff = 0.05)
res1 <- KSEA.Complete2(
  KSData_full,
  PX = px_list$LPS_Control,
  NetworKIN = TRUE,
  NetworKIN.cutoff = 5,
  m.cutoff = 5,
  p.cutoff = 0.05,
  plot_file = "../Results/KSEA/LPS_Control_barplot.tiff"
)
res2 <- KSEA.Complete2(
  KSData_full,
  PX = px_list$LPS_IKK_Control,
  NetworKIN = TRUE,
  NetworKIN.cutoff = 5,
  m.cutoff = 5,
  p.cutoff = 0.05,
  plot_file = "../Results/KSEA/LPS_IKK_Control_barplot.tiff"
)
res3 <- KSEA.Complete2(
  KSData_full,
  PX = px_list$LPS_IKK_LPS,
  NetworKIN = TRUE,
  NetworKIN.cutoff = 5,
  m.cutoff = 5,
  p.cutoff = 0.05,
  plot_file = "../Results/KSEA/LPS_IKK_LPS_barplot.tiff"
)
res4 <- KSEA.Complete2(
  KSData_full,
  PX = px_list$ML162_Control,
  NetworKIN = TRUE,
  NetworKIN.cutoff = 5,
  m.cutoff = 5,
  p.cutoff = 0.05,
  plot_file = "../Results/KSEA/ML162_Control_barplot.tiff"
)
res5 <- KSEA.Complete2(
  KSData_full,
  PX = px_list$ML162_LPS,
  NetworKIN = TRUE,
  NetworKIN.cutoff = 5,
  m.cutoff = 5,
  p.cutoff = 0.05,
  plot_file = "../Results/KSEA/ML162_LPS_barplot.tiff"
)
res6 <- KSEA.Complete2(
  KSData_full,
  PX = px_list$ML162_LPS_IKK,
  NetworKIN = TRUE,
  NetworKIN.cutoff = 5,
  m.cutoff = 5,
  p.cutoff = 0.05,
  plot_file = "../Results/KSEA/ML162_LPS_IKK_barplot.tiff"
)

scores1 = KSEA.Scores(KSData_full, px_list$LPS_Control, NetworKIN=TRUE, NetworKIN.cutoff=5)
scores2 = KSEA.Scores(KSData_full, px_list$LPS_IKK_Control, NetworKIN=TRUE, NetworKIN.cutoff=5)
scores3 = KSEA.Scores(KSData_full, px_list$LPS_IKK_LPS, NetworKIN=TRUE, NetworKIN.cutoff=5)
scores4 = KSEA.Scores(KSData_full, px_list$ML162_Control, NetworKIN=TRUE, NetworKIN.cutoff=5)
scores5 = KSEA.Scores(KSData_full, px_list$ML162_LPS, NetworKIN=TRUE, NetworKIN.cutoff=5)
scores6 = KSEA.Scores(KSData_full, px_list$ML162_LPS_IKK, NetworKIN=TRUE, NetworKIN.cutoff=5)

KSEA.Heatmap2(
  score.list = list(scores1, scores2, scores3, scores4, scores5, scores6),
  sample.labels = c("LPS vs Control", "LPS-IKK vs Control", "LPS-IKK vs LPS",
                    "ML162 vs Control", "ML162 vs LPS", "ML162 vs LPS-IKK"),
  stats = "FDR",
  m.cutoff = 5,
  p.cutoff = 0.05,
  sample.cluster = TRUE,
  png_file = "../Results/KSEA/KSEA_Merged_Heatmap.png",
  #plot_width_scale = 1,   # 画布更宽
  srtCol = 60,
  cexCol = 0.85
)
}
