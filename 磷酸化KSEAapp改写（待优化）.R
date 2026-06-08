KSEA.Complete2 <- function(KSData, PX,
                           NetworKIN = TRUE,
                           NetworKIN.cutoff = 3,
                           m.cutoff = 3,
                           p.cutoff = 0.05,
                           plot_file = "KSEA Bar Plot.tiff") {
  
  # ---------------------------
  # 1) PX 预处理：如果 Residue.Both 有分号，拆成一行一个位点
  # ---------------------------
  if (length(grep(";", PX$Residue.Both)) == 0) {
    new <- PX
    colnames(new)[c(2, 4)] <- c("SUB_GENE", "SUB_MOD_RSD")
    new$log2FC <- log2(abs(as.numeric(as.character(new$FC))))
    new <- new[complete.cases(new$log2FC), ]
  } else {
    double <- PX[grep(";", PX$Residue.Both), ]
    residues <- as.character(double$Residue.Both)
    residues <- as.matrix(residues, ncol = 1)
    split <- strsplit(residues, split = ";")
    x <- sapply(split, length)
    
    single <- data.frame(
      Protein = rep(double$Protein, x),
      Gene = rep(double$Gene, x),
      Peptide = rep(double$Peptide, x),
      Residue.Both = unlist(split),
      p = rep(double$p, x),
      FC = rep(double$FC, x),
      stringsAsFactors = FALSE
    )
    
    new <- PX[-grep(";", PX$Residue.Both), ]
    new <- rbind(new, single)
    colnames(new)[c(2, 4)] <- c("SUB_GENE", "SUB_MOD_RSD")
    new$log2FC <- log2(abs(as.numeric(as.character(new$FC))))
    new <- new[complete.cases(new$log2FC), ]
  }
  
  # ---------------------------
  # 2) 过滤 KSData：PSP 或 PSP+NetworKIN
  # ---------------------------
  if (NetworKIN == TRUE) {
    KSData.filtered <- KSData[grep("[a-z]", KSData$Source), ]
    KSData.filtered <- KSData.filtered[(KSData.filtered$networkin_score >= NetworKIN.cutoff), ]
  } else {
    KSData.filtered <- KSData[grep("PhosphoSitePlus", KSData$Source), ]
  }
  
  # ---------------------------
  # 3) 合并：得到 Kinase-Substrate links
  # ---------------------------
  KSData.dataset <- merge(KSData.filtered, new)
  KSData.dataset <- KSData.dataset[order(KSData.dataset$GENE), ]
  
  KSData.dataset$Uniprot.noIsoform <- sapply(
    KSData.dataset$KIN_ACC_ID,
    function(x) unlist(strsplit(as.character(x), split = "-"))[1]
  )
  
  KSData.dataset.abbrev <- KSData.dataset[, c(5, 1, 2, 16:19, 14)]
  colnames(KSData.dataset.abbrev) <- c(
    "Kinase.Gene", "Substrate.Gene", "Substrate.Mod",
    "Peptide", "p", "FC", "log2FC", "Source"
  )
  
  KSData.dataset.abbrev <- KSData.dataset.abbrev[order(
    KSData.dataset.abbrev$Kinase.Gene,
    KSData.dataset.abbrev$Substrate.Gene,
    KSData.dataset.abbrev$Substrate.Mod,
    KSData.dataset.abbrev$p
  ), ]
  
  KSData.dataset.abbrev <- aggregate(
    log2FC ~ Kinase.Gene + Substrate.Gene + Substrate.Mod + Source,
    data = KSData.dataset.abbrev,
    FUN = mean
  )
  KSData.dataset.abbrev <- KSData.dataset.abbrev[order(KSData.dataset.abbrev$Kinase.Gene), ]
  
  # ---------------------------
  # 4) 计算 kinase score（完全照原函数）
  # ---------------------------
  kinase.list <- as.vector(KSData.dataset.abbrev$Kinase.Gene)
  kinase.list <- as.matrix(table(kinase.list))
  
  Mean.FC <- aggregate(log2FC ~ Kinase.Gene, data = KSData.dataset.abbrev, FUN = mean)
  Mean.FC <- Mean.FC[order(Mean.FC[, 1]), ]
  Mean.FC$mS <- Mean.FC[, 2]
  Mean.FC$Enrichment <- Mean.FC$mS / abs(mean(new$log2FC, na.rm = TRUE))
  Mean.FC$m <- kinase.list
  Mean.FC$z.score <- ((Mean.FC$mS - mean(new$log2FC, na.rm = TRUE)) *
                        sqrt(Mean.FC$m)) / sd(new$log2FC, na.rm = TRUE)
  Mean.FC$p.value <- pnorm(-abs(Mean.FC$z.score))
  Mean.FC$FDR <- p.adjust(Mean.FC$p.value, method = "fdr")
  
  Mean.FC.filtered <- Mean.FC[(Mean.FC$m >= m.cutoff), -2]
  Mean.FC.filtered <- Mean.FC.filtered[order(Mean.FC.filtered$z.score), ]
  
  plot.height <- length(Mean.FC.filtered$z.score)^0.55
  Mean.FC.filtered$color <- "black"
  Mean.FC.filtered[(Mean.FC.filtered$p.value < p.cutoff) & (Mean.FC.filtered$z.score < 0), ncol(Mean.FC.filtered)] <- "blue"
  Mean.FC.filtered[(Mean.FC.filtered$p.value < p.cutoff) & (Mean.FC.filtered$z.score > 0), ncol(Mean.FC.filtered)] <- "red"
  
  # ---------------------------
  # 5) 画图：逻辑不改，只把文件名改为 plot_file
  # ---------------------------
  if (!is.null(plot_file) && !is.na(plot_file) && nzchar(plot_file)) {
    dir.create(dirname(plot_file), showWarnings = FALSE, recursive = TRUE)
    
    tiff(plot_file, width = 6 * 300, height = 300 * plot.height, res = 300, pointsize = 13)
    par(mai = c(1, 1, 0.4, 0.4))
    barplot(as.numeric(Mean.FC.filtered$z.score),
            col = Mean.FC.filtered$color,
            border = NA,
            xpd = FALSE,
            cex.names = 0.6,
            cex.axis = 0.8,
            xlab = "Kinase z-score",
            names.arg = Mean.FC.filtered$Kinase.Gene,
            horiz = TRUE,
            las = 1)
    dev.off()
  }
  
  # ---------------------------
  # 6) 返回结果：方便你在R里直接查看
  # ---------------------------
  return(list(
    scores_all = Mean.FC,
    scores_filtered = Mean.FC.filtered,
    links = KSData.dataset.abbrev,
    px_used = new,
    plot_file = plot_file
  ))
}
