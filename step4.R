source("C:/Work/SH/code/source/gene_enrichment.R")

gene_list <- lapply(names(res_list), function(nm) {
  df <- res_list[[nm]]
  df %>% 
    filter(is.finite(log2FC), abs(log2FC) > 1) %>%
    pull(Genes) %>%
    as.character()
})

names(gene_list) <- names(res_list)

# 把 gene_list 里每个向量都转成 “首字母大写，其余小写”
gene_list_mousecase <- gene_list
library(org.Mm.eg.db)
res_LPS_Control <- run_enrichment(gene_list_mousecase$LPS_Control, OrgDb = org.Mm.eg.db,
                      species_kegg = "mmu",
                      out_dir = "../Results/Enrichment/LPS_Control/",
                      file_tag = "")

res_LPS_IKK_Control <- run_enrichment(gene_list_mousecase$LPSIKK_Control, OrgDb = org.Mm.eg.db,
                                  species_kegg = "mmu",
                                  out_dir = "../Results/Enrichment/LPS_IKK_Control/",
                                  file_tag = "")

res_LPS_IKK_LPS <- run_enrichment(gene_list_mousecase$LPSIKK_LPS, OrgDb = org.Mm.eg.db,
                                      species_kegg = "mmu",
                                      out_dir = "../Results/Enrichment/LPS_IKK_LPS/",
                                      file_tag = "")

res_ML162_Control <- run_enrichment(gene_list_mousecase$ML162_Control, OrgDb = org.Mm.eg.db,
                                  species_kegg = "mmu",
                                  out_dir = "../Results/Enrichment/ML162_Control/",
                                  file_tag = "")

res_ML162_LPS <- run_enrichment(gene_list_mousecase$ML162_LPS, OrgDb = org.Mm.eg.db,
                                    species_kegg = "mmu",
                                    out_dir = "../Results/Enrichment/ML162_LPS/",
                                    file_tag = "")

res_ML162_LPS_IKK <- run_enrichment(gene_list_mousecase$ML162_LPSIKK, OrgDb = org.Mm.eg.db,
                                species_kegg = "mmu",
                                out_dir = "../Results/Enrichment/ML162_LPS_IKK/",
                                file_tag = "")


#--------------- 
# 加载所需的包
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(openxlsx)
library(pheatmap)
library(tibble)


df_pg <- data %>% dplyr::select(Protein.Group, sample_cols)
df_pg <- column_to_rownames(df_pg, var = "Protein.Group")

df_pg <- data.frame(lapply(df_pg, as.numeric))
colnames(df_pg) <- sample_metadata$样品信息
rownames(df_pg) <- data$Protein.Group

# 检查转换后的数据框
str(df_pg)


# 进行Log2转化并绘制箱线图
df_pg_log2 <- log2(df_pg + 1)  # 加1避免log(0)
df_pg_log2 <- rownames_to_column(df_pg_log2, var = "Protein.Group")
# 将数据转换为长格式，方便绘图
df_long <- df_pg_log2 %>%
  dplyr::select(Protein.Group, all_of(sample_metadata$样品信息)) %>%
  gather(key = "Sample", value = "Protein_Expression", -Protein.Group)

# 去除NA值
df_long <- na.omit(df_long)

# 绘制每个样本的蛋白质分布箱线图
ggplot(df_long, aes(x = Sample, y = Protein_Expression)) +
  geom_boxplot(aes(fill = Sample), alpha = 0.7) +
  theme_minimal() +
  labs(
    #title = "每个样本的蛋白质表达分布",
    x = "Sample",
    y = "Log2(intensity)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5)
  )
#ggsave("../Results/箱线图.pdf", plot = last_plot(), width = 6, height = 4)
#ggsave("../Results/箱线图.tiff", plot = last_plot(), width = 6, height = 4)
ggsave("../Results/箱线图.png", plot = last_plot(), width = 6, height = 4)



# 加载必要的包
library(pheatmap)

# 提取相关的样本列（去掉Protein.Group列）
df_pg_log2_samples <- df_pg_log2[, sample_metadata$样品信息]

# 计算Pearson相关性矩阵
cor_matrix <- cor(df_pg_log2_samples, method = "pearson", use = "pairwise.complete.obs")

# 打印相关性矩阵
print(cor_matrix)

# 绘制相关性热图
pheatmap(cor_matrix, 
         display_numbers = TRUE,  # 显示数字
         fontsize_number = 12,    # 设置数字字体大小
         cluster_rows = TRUE,     # 行聚类
         cluster_cols = TRUE,     # 列聚类
         color = colorRampPalette(c("#F7F7F7", "#FF0000"))(100),  # 热图颜色渐变
         main = "Pearson Correlation"#,  # 热图标题
         #filename = "../result/Correlation_Heatmap.png", # 保存图片路径
         #dpi = 300, width = 6, height = 6  # 图像保存设置
)
dev.off()


png("../Results/相关性分析.png", width = 6, height = 5, units = "in", res = 300)
pheatmap(cor_matrix, 
         display_numbers = TRUE,  # 显示数字
         fontsize_number = 12,    # 设置数字字体大小
         cluster_rows = TRUE,     # 行聚类
         cluster_cols = TRUE,     # 列聚类
         color = colorRampPalette(c("#F7F7F7", "#FF0000"))(100),  # 热图颜色渐变
         main = "Pearson Correlation"#,  # 热图标题
         #filename = "../result/Correlation_Heatmap.png", # 保存图片路径
         #dpi = 300, width = 6, height = 6  # 图像保存设置
)
dev.off()

