library(gplots)
KSEA.Heatmap2 <- function(score.list, sample.labels, stats,
                          m.cutoff, p.cutoff, sample.cluster,
                          png_file = "KSEA.Merged.Heatmap.png",
                          # ✅ 不用 margins，靠画布宽度解决裁剪
                          srtCol = 60,
                          cexCol = 0.9,
                          plot_width_scale = 1) {
  
  filter.m <- function(dataset, m.cutoff) {
    filtered <- dataset[(dataset$m >= m.cutoff), ]
    return(filtered)
  }
  
  score.list.m <- lapply(score.list, function(...) filter.m(..., m.cutoff))
  
  for (i in 1:length(score.list.m)) {
    names <- colnames(score.list.m[[i]])[c(2:7)]
    colnames(score.list.m[[i]])[c(2:7)] <- paste(names, i, sep = ".")
  }
  
  master <- Reduce(function(...) merge(..., by = "Kinase.Gene", all = FALSE), score.list.m)
  row.names(master) <- master$Kinase.Gene
  
  columns <- as.character(colnames(master))
  merged.scores <- as.matrix(master[, grep("z.score", columns)])
  colnames(merged.scores) <- sample.labels
  merged.stats <- as.matrix(master[, grep(stats, columns)])
  
  asterisk <- function(matrix) {
    new <- data.frame()
    for (i in 1:nrow(matrix)) {
      for (j in 1:ncol(matrix)) {
        if (matrix[i, j] < p.cutoff) new[i, j] <- "*" else new[i, j] <- ""
      }
    }
    return(new)
  }
  merged.asterisk <- as.matrix(asterisk(merged.stats))
  
  create.breaks <- function(merged.scores) {
    if (min(merged.scores) < -1.6) {
      breaks.neg <- seq(-1.6, 0, length.out = 30)
      breaks.neg <- append(seq(min(merged.scores), -1.6, length.out = 10), breaks.neg)
      breaks.neg <- sort(unique(breaks.neg))
    } else {
      breaks.neg <- seq(-1.6, 0, length.out = 30)
    }
    if (max(merged.scores) > 1.6) {
      breaks.pos <- seq(0, 1.6, length.out = 30)
      breaks.pos <- append(breaks.pos, seq(1.6, max(merged.scores), length.out = 10))
      breaks.pos <- sort(unique(breaks.pos))
    } else {
      breaks.pos <- seq(0, 1.6, length.out = 30)
    }
    breaks.all <- unique(append(breaks.neg, breaks.pos))
    mycol.neg <- colorpanel(n = length(breaks.neg), low = "blue", high = "white")
    mycol.pos <- colorpanel(n = length(breaks.pos) - 1, low = "white", high = "red")
    mycol <- unique(append(mycol.neg, mycol.pos))
    list(breaks.all, mycol)
  }
  
  color.breaks <- create.breaks(merged.scores)
  
  plot.height <- nrow(merged.scores)^0.55
  plot.width  <- (ncol(merged.scores)^0.7) * plot_width_scale  # ✅ 画布变宽
  
  dir.create(dirname(png_file), showWarnings = FALSE, recursive = TRUE)
  png(png_file, width = plot.width * 300, height = plot.height * 300, res = 300, pointsize = 14)
  
  heatmap.2(
    merged.scores,
    Colv = sample.cluster,
    scale = "none",
    cellnote = merged.asterisk,
    notecol = "white",
    cexCol = cexCol,
    cexRow = 0.9,
    srtCol = srtCol,
    notecex = 1.4,
    col = color.breaks[[2]],
    density.info = "none",
    trace = "none",
    key = FALSE,
    breaks = color.breaks[[1]],
    lmat = rbind(c(0, 3), c(2, 1), c(0, 4)),
    lhei = c(0.4, 9.5, 0.6),
    lwid = c(0.5, 3)
    # ✅ 不再传 margins
  )
  
  dev.off()
  
  invisible(list(
    merged.scores = merged.scores,
    merged.stats = merged.stats,
    merged.asterisk = merged.asterisk,
    png_file = png_file
  ))
}
