## =========================
## 1. 重新生成上游磷酸化矩阵
## =========================
source("./demo/R/01_build_phosphosite_matrix_from_diann.R")

prep <- prepare_phosphosite_matrix(
  pr_matrix_file = "./demo/data/4_Phos_Mouse_report.pr_matrix.tsv",
  fasta_file = "C:/Work/SH/Pub_database/添加GO注释/data/UP000000589_17155_Reviewed_20230517_uniprot-download_true_format_fasta_includeIsoform_true_query__28Mus_-2023.05.17-06.21.58.02.fasta",
  out_dir = "./demo/Results/prepare_phosphosite_matrix",
  sample_prefix_pattern = "YAS202601070072-3-",
  collapse_method = "max"
)

## =========================
## 2. 样本信息整理 + 磷酸化位点差异分析
## =========================

source("./demo/R/02_prepare_metadata_and_phosphosite_de.R")

de <- run_phosphosite_batch_de(
  final_mat = prep$final_mat,
  sample_cols = prep$sample_cols,
  sampleinfo_file = "./demo/data/sampleinfo.xlsx",
  samplegroup_file = "./demo/data/samplegroup.xlsx",
  compare_list_file = "./demo/data/compare_list.xlsx",
  out_dir = "./demo/Results/phosphosite_DE",
  samplegroup_skip = 10,
  sample_id_col = "样品编号",
  sample_name_col = "样品信息",
  samplegroup_cols = c("*样品名称", "*管上名称", "*样品类型", "*样品来源", "样本分组"),
  compare_exp_col = "实验",
  compare_ctrl_col = "对照",
  min_reps_for_limma = 3,
  low_rep_strategy = "group_fc_only"
)

## =========================
## 3. KSEA
## =========================

unlink("./demo/Results/phosphosite_KSEA", recursive = TRUE)

source("./demo/R/03_run_phosphosite_ksea.R")

ksea <- run_phosphosite_ksea(
  de_results = de$all_results,
  site_mat2 = prep$site_mat2,
  out_dir = "./demo/Results/phosphosite_KSEA",
  species = "mouse",
  keep_only_single_site = TRUE,
  min_sites = 10
)

#save.image("./demo/R/1.RData")

source("./demo/R/03b_plot_ksea_kinase_substrate_chord.R")

chord <- plot_all_ksea_kinase_substrate_chord(
  ksea_dir = "./demo/Results/phosphosite_KSEA",
  top_n_kinases = 8,
  max_sites_per_kinase = 10,
  min_abs_signed_score = 0
)

## =========================
## 4. motif
## =========================

source("./demo/R/04_phosphosite_motif_analysis.R")

motif <- run_phosphosite_motif_analysis(
  site_mat2 = prep$site_mat2,
  fasta_file = "C:/Work/SH/Pub_database/添加GO注释/data/UP000000589_17155_Reviewed_20230517_uniprot-download_true_format_fasta_includeIsoform_true_query__28Mus_-2023.05.17-06.21.58.02.fasta",
  out_dir = "./demo/Results/phosphosite_motif",
  flank = 6,
  keep_unique_site = TRUE,
  make_logo = TRUE
)

motif$summary_by_residue
motif$fail_summary

motif_full <- run_phosphosite_motif_analysis(
  site_mat2 = prep$site_mat2,
  fasta_file = "C:/Work/SH/Pub_database/添加GO注释/data/UP000000589_17155_Reviewed_20230517_uniprot-download_true_format_fasta_includeIsoform_true_query__28Mus_-2023.05.17-06.21.58.02.fasta",
  out_dir = "./demo/Results/phosphosite_motif_full_for_comparison",
  flank = 6,
  keep_unique_site = FALSE,
  make_logo = FALSE
)



motif_comp <- run_comparison_motif_analysis(
  motif_table = motif_full$motif_table,
  de_results = de$all_results,
  out_dir = "./demo/Results/phosphosite_motif/by_comparison",
  comparison_col = "comparison",
  logfc_col = "logFC",
  logfc_thresh = 1,
  residues = c("S", "T", "Y"),
  min_sequences_for_logo = 5,
  make_logo = TRUE
)

#motif_comp$match_summary
#head(motif_comp$trend_motif_table)
#
#
#motif_comp$trend_motif_table %>%
#  dplyr::count(comparison, trend_direction, abs_residue, name = "n_sites")


source("./demo/R/04b_motif_enrichment_Ronly.R")
motif_enrich <- run_r_motif_enrichment_pipeline(
  motif_table = motif$motif_table,
  trend_motif_table = motif_comp$trend_motif_table,
  out_dir = "./demo/Results/phosphosite_motif_enrichment_Ronly",
  flank = 6,
  run_all_sites = TRUE,
  run_all_sites_by_residue = TRUE,
  run_comparison_trends = TRUE,
  n_shuffle = 5,
  max_order = 2,
  top_seed = 30,
  min_count = 5,
  min_fold = 1.2,
  top_n = 20
)

motif_enrich$run_summary
motif_enrich$plot_files


## =========================
## 5. 亚细胞定位
## =========================
source("./demo/R/05a_get_uniprot_subcellular_annotation.R")

uniprot_loc <- make_uniprot_subcellular_annotation_for_phos(
 de_results = de$all_results,
 organism_id = 10090,
 reviewed_only = TRUE,
 use_only_proteins_in_de_results = TRUE,
 out_dir = "./demo/Results/subcellular_localization/uniprot_annotation",
 final_out_file = "./demo/data/subcellular_annotation_from_uniprot.csv"
)
head(uniprot_loc$annotation)

source("./demo/R/05_subcellular_localization_analysis.R")

subcell <- run_phosphosite_subcellular_localization(
  de_results = de$all_results,
  subcellular_annotation_file = "./demo/data/subcellular_annotation_from_uniprot.csv",
  out_dir = "./demo/Results/subcellular_localization_uniprot",
  diff_mode = "direction",
  keep_directions = c("Up", "Down")
)


## =========================
## 6. 结构域分析
## =========================
source("./demo/R/06a_get_uniprot_domain_annotation.R")

domain_anno <- make_uniprot_domain_annotation_for_phos(
  de_results = de$all_results,
  organism_id = 10090,
  reviewed_only = TRUE,
  use_only_proteins_in_de_results = TRUE,
  
  ## 建议：保留 UniProt FT 的 DOMAIN/REPEAT/MOTIF，再加 InterPro
  ## 不建议第一版保留 REGION，因为 REGION 里面大量是 Disordered，会污染 Top20
  include_feature_types = c("DOMAIN", "REPEAT", "MOTIF"),
  
  ## 关键：用 InterPro，不用 Pfam
  use_xref_interpro = TRUE,
  use_xref_pfam = FALSE,
  
  out_dir = "./demo/Results/domain_analysis/uniprot_annotation_interpro_only",
  final_out_file = "./demo/data/domain_annotation_from_uniprot_interpro_only.csv"
)


source("./demo/R/06b_map_interpro_entry_list_to_name.R")

domain_named <- map_domain_annotation_by_interpro_entry_list(
  domain_annotation_file = "./demo/data/domain_annotation_from_uniprot_interpro_only.csv",
  entry_list_file = "./demo/Results/domain_analysis/interpro_metadata/entry.list.txt",
  final_out_file = "./demo/data/domain_annotation_from_uniprot_interpro_only_named.csv",
  keep_entry_types = c("Domain", "Family", "Repeat", "Homologous_superfamily"),
  drop_unmapped_id_like_names = TRUE
)

source("./demo/R/06_domain_analysis.R")

domain_res <- run_phosphosite_domain_analysis(
  de_results = de$all_results,
  domain_annotation_file = "./demo/data/domain_annotation_from_uniprot_interpro_only_named.csv",
  out_dir = "./demo/Results/domain_analysis_uniprot_interpro_named",
  diff_mode = "direction",
  keep_directions = c("Up", "Down"),
  top_n = 20
)
