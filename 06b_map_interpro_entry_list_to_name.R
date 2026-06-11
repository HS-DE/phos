# ============================================================
# 06b. Map InterPro IDs to readable names using entry.list
#
# Input:
#   domain_annotation_from_uniprot_interpro_only.csv
#
# Mapping source:
#   https://ftp.ebi.ac.uk/pub/databases/interpro/releases/latest/entry.list
#
# Output:
#   domain_annotation_from_uniprot_interpro_only_named.csv
#
# Purpose:
#   Convert IPRxxxxx labels to readable InterPro entry names.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
})

.safe_dir_create <- function(path) {
  if (!is.null(path) && nzchar(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

is_ipr_id <- function(x) {
  stringr::str_detect(as.character(x), "^IPR[0-9]+$")
}

is_id_like_domain_name <- function(x) {
  stringr::str_detect(
    as.character(x),
    "^(IPR[0-9]+|PF[0-9]+|PS[0-9]+|SM[0-9]+|SSF[0-9]+|G3DSA:[0-9\\.]+)$"
  )
}

download_interpro_entry_list <- function(
    out_file = "./demo/Results/domain_analysis/interpro_metadata/entry.list"
) {
  .safe_dir_create(dirname(out_file))
  
  url <- "https://ftp.ebi.ac.uk/pub/databases/interpro/releases/latest/entry.list"
  
  message("Downloading InterPro entry.list...")
  message(url)
  
  entry_map <- readr::read_tsv(
    url,
    show_col_types = FALSE,
    progress = TRUE
  )
  
  readr::write_tsv(entry_map, out_file)
  
  message("Saved entry.list: ", out_file)
  message("Rows: ", nrow(entry_map))
  message("Columns: ", paste(colnames(entry_map), collapse = " | "))
  
  entry_map
}

read_or_download_interpro_entry_list <- function(
    entry_list_file = "./demo/Results/domain_analysis/interpro_metadata/entry.list"
) {
  if (file.exists(entry_list_file)) {
    message("Reading local entry.list: ", entry_list_file)
    
    readr::read_tsv(
      entry_list_file,
      show_col_types = FALSE,
      progress = FALSE
    )
  } else {
    download_interpro_entry_list(entry_list_file)
  }
}

map_domain_annotation_by_interpro_entry_list <- function(
    domain_annotation_file = "./demo/data/domain_annotation_from_uniprot_interpro_only.csv",
    entry_list_file = "./demo/Results/domain_analysis/interpro_metadata/entry.list",
    final_out_file = "./demo/data/domain_annotation_from_uniprot_interpro_only_named.csv",
    
    ## 推荐第一版保留这些类型
    keep_entry_types = c("Domain", "Family", "Repeat", "Homologous_superfamily"),
    
    ## 如果还有无法映射的 IPR/PF ID，是否直接丢掉
    drop_unmapped_id_like_names = TRUE,
    
    ## 是否把 IPR ID 替换成 entry.list 里的 ENTRY_NAME
    prefer_interpro_name = TRUE
) {
  .safe_dir_create(dirname(final_out_file))
  .safe_dir_create(dirname(entry_list_file))
  
  domain_anno <- read.csv(
    domain_annotation_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (!"Protein.Id" %in% colnames(domain_anno)) {
    stop("domain_annotation_file must contain Protein.Id")
  }
  
  if (!"Domain.Name" %in% colnames(domain_anno)) {
    stop("domain_annotation_file must contain Domain.Name")
  }
  
  if (!"Domain.ID" %in% colnames(domain_anno)) {
    domain_anno$Domain.ID <- NA_character_
  }
  
  entry_map <- read_or_download_interpro_entry_list(entry_list_file)
  
  required_cols <- c("ENTRY_AC", "ENTRY_TYPE", "ENTRY_NAME")
  missing_cols <- setdiff(required_cols, colnames(entry_map))
  
  if (length(missing_cols) > 0) {
    stop("entry.list missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  entry_map2 <- entry_map %>%
    transmute(
      InterPro.ID = as.character(ENTRY_AC),
      InterPro.Type = as.character(ENTRY_TYPE),
      InterPro.Name = as.character(ENTRY_NAME)
    ) %>%
    filter(!is.na(InterPro.ID), InterPro.ID != "") %>%
    filter(!is.na(InterPro.Name), InterPro.Name != "")
  
  if (!is.null(keep_entry_types)) {
    entry_map2 <- entry_map2 %>%
      filter(InterPro.Type %in% keep_entry_types)
  }
  
  out <- domain_anno %>%
    mutate(
      Domain.ID = as.character(Domain.ID),
      Domain.Name = as.character(Domain.Name),
      Domain.Name.Original = Domain.Name,
      
      ## 有些行 Domain.ID 是 IPR；
      ## 有些行 Domain.Name 本身是 IPR。
      InterPro.ID.For.Map = case_when(
        is_ipr_id(Domain.ID) ~ Domain.ID,
        is_ipr_id(Domain.Name) ~ Domain.Name,
        TRUE ~ NA_character_
      )
    ) %>%
    left_join(
      entry_map2,
      by = c("InterPro.ID.For.Map" = "InterPro.ID")
    ) %>%
    mutate(
      ## 只有原来的 Domain.Name 是 ID-like 时，才替换成 InterPro.Name；
      ## 如果原来已经是可读名称，比如 "IQ"、"Ig-like C2-type 1"，就保留。
      Domain.Name = case_when(
        isTRUE(prefer_interpro_name) &
          is_id_like_domain_name(Domain.Name.Original) &
          !is.na(InterPro.Name) &
          InterPro.Name != "" ~ InterPro.Name,
        
        TRUE ~ Domain.Name.Original
      ),
      
      Domain.ID = case_when(
        (is.na(Domain.ID) | Domain.ID == "") &
          is_ipr_id(Domain.Name.Original) ~ Domain.Name.Original,
        
        TRUE ~ Domain.ID
      ),
      
      Domain.Entry.Type = InterPro.Type
    )
  
  if (isTRUE(drop_unmapped_id_like_names)) {
    out <- out %>%
      filter(!is_id_like_domain_name(Domain.Name))
  }
  
  out <- out %>%
    select(
      any_of(c(
        "Protein.Id",
        "Genes",
        "Domain.ID",
        "Domain.Name",
        "Domain.Name.Original",
        "InterPro.ID.For.Map",
        "InterPro.Name",
        "Domain.Entry.Type",
        "Domain.Source",
        "Domain.Start",
        "Domain.End"
      ))
    ) %>%
    filter(!is.na(Domain.Name), Domain.Name != "") %>%
    distinct()
  
  write.csv(
    out,
    final_out_file,
    row.names = FALSE,
    fileEncoding = "GBK"
  )
  
  message("Saved named domain annotation: ", final_out_file)
  
  message("Before mapping, ID-like Domain.Name:")
  print(table(is_id_like_domain_name(domain_anno$Domain.Name), useNA = "ifany"))
  
  message("After mapping, ID-like Domain.Name:")
  print(table(is_id_like_domain_name(out$Domain.Name), useNA = "ifany"))
  
  message("Entry types used:")
  print(table(out$Domain.Entry.Type, useNA = "ifany"))
  
  message("Preview:")
  print(
    out %>%
      select(any_of(c(
        "Protein.Id",
        "Genes",
        "Domain.ID",
        "Domain.Name.Original",
        "Domain.Name",
        "Domain.Entry.Type",
        "Domain.Source"
      ))) %>%
      head(20)
  )
  
  invisible(out)
}

# ============================================================
# Example usage
# ============================================================
#
# source("./demo/R/06b_map_interpro_entry_list_to_name.R")
#
# domain_named <- map_domain_annotation_by_interpro_entry_list(
#   domain_annotation_file = "./demo/data/domain_annotation_from_uniprot_interpro_only.csv",
#   entry_list_file = "./demo/Results/domain_analysis/interpro_metadata/entry.list",
#   final_out_file = "./demo/data/domain_annotation_from_uniprot_interpro_only_named.csv",
#   keep_entry_types = c("Domain", "Family", "Repeat", "Homologous_superfamily"),
#   drop_unmapped_id_like_names = TRUE
# )