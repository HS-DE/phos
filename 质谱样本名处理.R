extract_sample_names <- function(expr, prefix_pattern = "_nanomics_") {
  
  col_names <- colnames(expr)
  
  # 自动检测样本列起始位置
  hit_cols <- grepl(prefix_pattern, col_names)
  
  if (!any(hit_cols)) {
    stop("未在列名中检测到 prefix_pattern：", prefix_pattern)
  }
  
  # 第一列命中即样本列起始
  start_col <- which(hit_cols)[1]
  
  sample_cols <- start_col:length(col_names)
  
  # 逐步处理
  x <- col_names[sample_cols]
  
  # 1. 去路径
  x <- basename(x)
  
  # 2. 去后缀
  x <- sub("\\.raw$", "", x)
  
  # 3. 去前缀
  x <- sub(paste0("^.*", prefix_pattern), "", x)
  
  # 回填
  col_names[sample_cols] <- x
  colnames(expr) <- col_names
  sample_name <- colnames(expr)[sample_cols]
  
  return(list(
    expr = expr,
    sample_cols = sample_name
  ))
  
}




if (F) {
  #使用方法：
  col_names_new <- extract_sample_names(
    expr = pg_data,
    #sample_cols = 6:length(col_names),
    prefix_pattern = "_nanomics_" # 默认为""，如果前缀更换，修改即可，如"_projectX_"
  )
  
}