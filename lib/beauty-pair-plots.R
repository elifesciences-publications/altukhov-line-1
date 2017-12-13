add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


# beatuy function for shadowed text
shadowtext <- function(x, y=NULL, labels, col='white', bg='black', 
                       theta= seq(0, 2*pi, length.out=50), r=0.1, ... ) {
  
  xy <- xy.coords(x,y)
  xo <- r*strwidth('A')
  yo <- r*strheight('A')
  
  # draw background text with small shift in x and y in background colour
  for (i in theta) {
    text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
  }
  # draw actual text in exact xy position in foreground colour
  text(xy$x, xy$y, labels, col=col, ... )
}

plotTwoExp <- function(m.affcol, ex1, ex2, plot.dens = F,
                       p1.sig.th = 0.7, p2.sig.th = 0.7, plot.empty = T) {
  # default colors
  col.pnt <- rgb(0, 0, 0, 1)
  col.png.sig <- "dodgerblue"
  col.hst <- "white"
  col.hst.sig <- "lightblue"
  col.hst.brd <- "black"
  col.lne <- "gray80"
  col.lne.sig <- "dodgerblue"
  col.lne.bsig <- "lightblue"
  col.lne.dsig <- "orange"
  col.dns <- "red"
  
  # chauvenet iterations
  n.it <- 1
  
  # subset matrix
  m.wk <- m.affcol[["amtx"]][, c(ex1, ex2)]
  m.wk <- m.wk[order(m.wk[, 1]), ]
  
  if (!plot.empty & sum(apply(m.wk, 1, function(x) sum(!is.na(x))) == 2) == 0) {
   return(FALSE) 
  }
  
  if (plot.dens) {
    x1 <- m.wk[, 1][!is.na(m.wk[, 1])]
    x2 <- m.wk[, 2][!is.na(m.wk[, 2])]
    
    for (i in 1:n.it) {
      x1 <- ChauvenetFilter(x1, th = 2)
      x2 <- ChauvenetFilter(x2, th = 2)
    }
    #x1 <- ChauvenetFilter(ChauvenetFilter(m.wk[, 1][!is.na(m.wk[, 1])], 2), 2)
    #x2 <- ChauvenetFilter(ChauvenetFilter(m.wk[, 2][!is.na(m.wk[, 2])], 2), 2)
    
    x1.d <- density(rnorm(3e3, mean(x1), sd(x1)))
    x2.d <- density(rnorm(3e3, mean(x2), sd(x2)))
    
    pv1 <- p.adjust(pnorm(q = m.wk[, 1], mean = mean(x1), sd = sd(x1), lower.tail = F), method = "BH")
    pv2 <- p.adjust(pnorm(q = m.wk[, 2], mean = mean(x2), sd = sd(x2), lower.tail = F), method = "BH")
    
    p1.sig.th <- m.wk[names(which.max(pv1[which(pv1 <= 0.05)])), 1]
    p2.sig.th <- m.wk[names(which.max(pv2[which(pv2 <= 0.05)])), 2]
    
    if(length(p1.sig.th) == 0) {
      p1.sig.th <- 1
    }
    if(length(p2.sig.th) == 0) {
      p2.sig.th <- 1
    }
  }
  
  # find significant points
  p1.sig <- m.wk[, 1] >= p1.sig.th
  p2.sig <- m.wk[, 2] >= p2.sig.th
  
  # create histogram objects
  hst.bns <- seq(0, 1, length.out = 51)
  h1 <- hist(m.wk[, 1], plot = F, 
             breaks = hst.bns)
  h2 <- hist(m.wk[, 2], plot = F,
             breaks = hst.bns)
  
  # prepare layout for figure
  layMat <- matrix(c(1, 3, 2), ncol = 1, byrow = TRUE)
  layout(layMat, heights=c(2/8, 4/8, 2/8))
  
  # plot histograms (as barplots)
  par(mar=c(0, 1, 3, 1))
  barplot(h1$counts, yaxt = "n",
          col = c(col.hst, col.hst.sig)[(hst.bns[-length(hst.bns)] >= p1.sig.th) + 1],
          border = col.hst.brd, space = 0)
  mtext(sprintf("%s", colnames(m.wk)[1]), side = 3, line = 0.5, cex = .6)
  mtext(sprintf("H/(H+L) signif: %s", round(p1.sig.th, 2)), side = 3,
        at = 45, line = 0.5, cex = .4)
  mtext(sprintf("Proteins: %s", sum(h1$counts)), side = 3,
        at = 5, line = 0.5, cex = .4)
  
  
  if (plot.dens) {
    lines(x = x1.d$x * length(h1$counts), y = x1.d$y * sum(h1$counts) * (1/50), col = col.dns,
          lwd = 1.5)
  }
  
  par(mar=c(3, 1, 0, 1))
  barplot(-h2$counts, yaxt = "n",
          col = c(col.hst, col.hst.sig)[(hst.bns[-length(hst.bns)] >= p2.sig.th) + 1],
          border = col.hst.brd, space = 0)
  mtext(sprintf("%s", colnames(m.wk)[2]), side = 1, line = 0.5, cex = .6)
  
  mtext(sprintf("H/(H+L) signif: %s", round(p2.sig.th, 2)), side = 1,
        at = 45, line = 0.5, cex = .4)
  mtext(sprintf("Proteins: %s", sum(h2$counts)), side = 1,
        at = 5, line = 0.5, cex = .4)
  
  
  
  if(plot.dens) {
    lines(x = x2.d$x * length(h2$counts), y = -x2.d$y * sum(h2$counts) * (1/50), col = col.dns,
          lwd = 1.5)
  }
  
  # plot middle figure with protein links
  par(mar=c(0, 1, 0, 1))
  plot(NA, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = F, ann = F)
  
  for (i in 1:nrow(m.wk)) {
    if (any(is.na(m.wk[i, ]))) next
    
    lines(x = m.wk[i, ], y = c(1, 0), col = col.lne)
  }
  
  for (i in 1:nrow(m.wk)) {
    if (any(is.na(m.wk[i, ]))) next
    
    is.txt <- F
    line.col <- col.lne
    if (p1.sig[i] & p2.sig[i]) {
      line.col <- col.lne.bsig
      is.txt <- T
      #dt.sig <<- rbind(dt.sig, data.frame(e1 = ex1, e2 = ex2, prot = names(p1.sig[i]), sig = 1))
    } else if (p1.sig[i] & !p2.sig[i]) {
      line.col <- col.lne.dsig
      is.txt <- T
      #dt.sig <<- rbind(dt.sig, data.frame(e1 = ex1, e2 = ex2, prot = names(p1.sig[i]), sig = 2))
    } else if (!p1.sig[i] & p2.sig[i]) {
      line.col <- col.lne.sig
      is.txt <- T
      #dt.sig <<- rbind(dt.sig, data.frame(e1 = ex1, e2 = ex2, prot = names(p1.sig[i]), sig = 3))
    }
    
    if (is.txt) {
      lines(x = m.wk[i, ], y = c(1, 0), col = line.col)
    }
  }
  
  text.line <- 0
  text.line.y <- seq(0.9, 0.1, by = -0.1)
  text.mline <- length(text.line.y)
  
  for (i in 1:nrow(m.wk)) {
    if (any(is.na(m.wk[i, ]))) next
    
    is.txt <- F
    line.col <- col.lne
    if (p1.sig[i] & p2.sig[i]) {
      line.col <- col.lne.bsig
      is.txt <- T
    } else if (p1.sig[i] & !p2.sig[i]) {
      line.col <- col.lne.dsig
      is.txt <- T
    } else if (!p1.sig[i] & p2.sig[i]) {
      line.col <- col.lne.sig
      is.txt <- T
    }
    
    if (is.txt) {
      x.dlt <- m.wk[i, 2] - m.wk[i, 1]
      shadowtext(rownames(m.wk)[i], x = m.wk[i, 1] + x.dlt * (1 - text.line.y[text.line + 1]),
                 y = text.line.y[text.line + 1],
                 col = "black", bg = line.col, 
                 cex = .6, srt = atan(m.wk[i, 2] - m.wk[i, 1]) * 180 / pi - 90)
      text.line <- (text.line + 1) %% text.mline
    }
  }
  
  
  p1.col <- c(col.pnt, col.png.sig)[p1.sig + 1]
  points(x = m.wk[, 1], y = rep(1, nrow(m.wk)), col = p1.col, pch = 22,
         cex = 1, bg = "gray")
  
  p2.col <- c(col.pnt, col.png.sig)[p2.sig + 1]
  points(x = m.wk[, 2], y = rep(0, nrow(m.wk)), col = p2.col, pch = 22,
         cex = 1, bg = "gray")
  
  legend(x = 0.05, y = 0.7, legend = c("Significant in both",
                                       "Significant in bottom",
                                       "Significant in top",
                                       "Non significant in both"),
         fill = c(col.lne.bsig, col.lne.sig, col.lne.dsig, col.lne),
         bty = "n", cex = .8)
  
  return(TRUE)
}

plotCompareTwoExp <- function(m.affcol, ex1, ex2){
  require(VennDiagram)
  # subset of matrix
  m.wk <- m.affcol[["amtx"]][, c(ex1, ex2)]
  m.wk <- m.wk[rowSums(is.na(m.wk)) != 2, ]
  m.wk[is.na(m.wk)] <- 0
  # correlation
  s.stat <- cor.test(m.wk[rowSums(m.wk == 0) == 0,ex1], 
                    m.wk[rowSums(m.wk == 0) == 0,ex2], 
                    method = "spearman")
  p.stat <- cor.test(m.wk[rowSums(m.wk == 0) == 0,ex1], 
                     m.wk[rowSums(m.wk == 0) == 0,ex2], 
                     method = "pearson")
  # plot
  layout(matrix(1:4, nrow = 2, byrow = T), widths = c(.5, .5))
  par(mar = c(4,4,4,0))
  plot(m.wk[rowSums(m.wk == 0) == 0,ex1], 
       m.wk[rowSums(m.wk == 0) == 0,ex2], 
       xlab = ex1, ylab = ex2, frame = F,
       xlim = c(0,1), ylim = c(0,1),
       col = "black", pch = 21, bg = "grey", cex = 1.5,
       main = sprintf("%s vs %s", ex1, ex2))
  points(m.wk[rowSums(m.wk == 0) != 0,ex1], 
         m.wk[rowSums(m.wk == 0) != 0,ex2],
         col = "grey", pch = 21, bg = "lightgrey", cex = 1)
  abline(a = 0, b = 1, col = "blue", lty = 2)
  text(x = .8, y = 1, "f(x) = x", col = 'blue')
  par(mar = c(2,2,2,2))
  plot.new()
  legend("topleft", 
         legend = c("Quantified in both\nI-DIRT experiments", 
                    "Quantified in only one\nI-DIRT experiment"),
         pch = 21, col = c("black", "grey"),
         pt.cex = c(1.5, 1), bty = "n", 
         pt.bg = c("grey", "lightgrey"),
         y.intersp = 2)
  text(0,.3, 
       sprintf(paste0("Spearman cor: %f\np-value: %f\n\n",
                      "Pearson cor: %f\np-value: %f"), 
               s.stat$estimate, s.stat$p.value,
               p.stat$estimate, p.stat$p.value),
       pos = 4, adj = 0)
  pr.list <- list(rownames(m.wk[m.wk[,ex1]!=0,]),
                  rownames(m.wk[m.wk[,ex2]!=0,]))
  names(pr.list) <- c(ex1,ex2)
  vd <- venn.diagram(pr.list, filename = NULL,
                     margin = .2, cex = 1, cat.dist = .07)
  gl <- grid.layout(nrow=2, ncol=2)
  
  # setup viewport
  vp <- viewport(layout.pos.col=1, layout.pos.row=2) 
  
  # init layout
  pushViewport(viewport(layout=gl))
  pushViewport(get("vp"))
  plot.new()
  grid.draw(vd)
  popViewport(1)
  return(TRUE)
}

plotScatterTwoExp <- function(m.affcol, ex1, ex2, cex = 1.5,
                              plot.blue.text = T){
  # subset of matrix
  m.wk <- m.affcol[["amtx"]][, c(ex1, ex2)]
  m.wk <- m.wk[rowSums(is.na(m.wk)) != 2, ]
  m.wk[is.na(m.wk)] <- 0
  
  plot(m.wk[rowSums(m.wk == 0) == 0,ex1], 
       m.wk[rowSums(m.wk == 0) == 0,ex2], 
       xlab = ex1, ylab = ex2, frame = F,
       xlim = c(0,1), ylim = c(0,1),
       col = "black", pch = 21, bg = "grey", cex = cex,
       main = sprintf("%s vs %s", ex1, ex2))
  points(m.wk[rowSums(m.wk == 0) != 0,ex1], 
         m.wk[rowSums(m.wk == 0) != 0,ex2],
         col = "grey", pch = 21, bg = "lightgrey", cex = 1)
  abline(a = 0, b = 1, col = "blue", lty = 2)
  if(plot.blue.text){
    text(x = .8, y = 1, "f(x) = x", col = 'blue')
  }
}

plotScatterIpTwoExp <- function(m.affcol, 
                                ex1, ex2, 
                                cex = 1.5,
                                main.title = NA,
                                plot.blue.text = T,
                                viz.proteins = c()){
  require(gtools)
  orf1.orf2 <- c("Orf2-untagged:", "sp|Q9UN81|LORF1_HUMAN")
  # subset of matrix
  m.wk <- m.affcol[["amtx"]][, c(ex1, ex2)]
  m.wk <- m.wk[rowSums(is.na(m.wk)) != 2, ]
  m.wk[is.na(m.wk)] <- 0
  
  # set layout
  layout(matrix(1:2, nrow = 1), widths = c(.7, .3))
  # plot main points
  par(mar = c(4,4,3,2))
  if(invalid(main.title)){
    main.title <- sprintf("%s vs %s", ex1, ex2)
  }
  plot(m.wk[rowSums(m.wk == 0) == 0,ex1], 
       m.wk[rowSums(m.wk == 0) == 0,ex2], 
       xlab = ex1, ylab = ex2, frame = F,
       xlim = c(0,1), ylim = c(0,1),
       col = add.alpha("grey50", .5), 
       pch = 19, 
       cex = cex,
       main = main.title)
  # plot y=x line
  abline(a = 0, b = 1, col = "blue", lty = 2)
  # plot zero points
  points(m.wk[rowSums(m.wk == 0) != 0, ex1], 
         m.wk[rowSums(m.wk == 0) != 0, ex2],
         col = "grey90", pch = 21, bg = "grey95", cex = cex)
  # highlight selected
  if(length(viz.proteins) > 0) {
    viz.rows <- unlist(sapply(viz.proteins[viz.proteins %nin% orf1.orf2], 
                              function(x){
      grep(x, rownames(m.wk), value = T)
    }))
    points(m.wk[viz.rows, ex1], 
           m.wk[viz.rows, ex2],
           col = "blue", pch = 21, bg = "lightblue", cex = cex)
  }
  # plot ORFs
  points(m.wk[orf1.orf2, ex1], 
         m.wk[orf1.orf2, ex2], 
         col = "black", 
         pch = 21, 
         bg = c("red", "lightblue"), cex = cex + .5)
  # add text to y=x line
  if(plot.blue.text){
    text(x = .8, y = 1, "f(x) = x", col = 'blue')
  }
  # plot legend
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = c("ORF1", "ORF2", 'From Cell article'),
         pch = 21, col = c("black", "black", "blue"),
         pt.bg = c("lightblue", "red", "lightblue"),
         pt.cex = c(cex + .5, cex + .5, cex),
         bty = "n")
}

plotHeatmap <- function(m.affcol, 
                        exps = NULL,
                        repl.list = NULL,
                        method = "spearman"){
  source("src/R/heatmap.3.R")
  
  # calculating distance as 1 - cor
  if(invalid(exps)) {
    dist.mrx <- 1 - cor(m.affcol[["amtx"]],
                        use = "pairwise.complete.obs", method = method) 
  } else {
    for(exp.excl in exps[exps %ni% colnames(m.affcol[["amtx"]])]){
      message(
        sprintf(
          "WARNING: experiment '%s' doesn't found in m.affcol. Excluded.", 
          exp.excl))
    }
    exps <- exps[exps %in% colnames(m.affcol[["amtx"]])]
    dist.mrx <- 1 - cor(m.affcol[["amtx"]][,exps],
                        use = "pairwise.complete.obs", method = method) 
  }
  
  if(!invalid(repl.list)){
    cols <- sapply(repl.list, function(x){
      c("white", "black")[as.numeric(rownames(dist.mrx) %in% x) + 1]
    })
    heatmap.3(dist.mrx, RowSideColors = t(cols), RowSideColorsSize = 1,
              KeyValueName = "Distance",
              side.height.fraction = 2,
              keysize = 1)
  } else {
    heatmap.3(dist.mrx, KeyValueName = "Distance", keysize = 1)
  }
  return(TRUE)
}

# plot 2 histograms
plotTwoExpDT <- function(dt, x1, x2, x1.name, x2.name, 
                         highlite = NULL, highlite.col = NULL,
                         highlite.names = NULL,
                         highlite.legend = NULL){
  # default colors
  col.pnt <- rgb(0, 0, 0, 1)
  col.hst <- "white"
  col.hst.sig <- "gray95"
  col.hst.brd <- "black"
  col.lne <- "gray50"
  orf1.col <- "#D55E00"
  orf2.col <- "#0072B2"
  if(is.null(highlite.col)){
    highlite.col <- rep('blueviolet', length(highlite))
  }
  
  # pg of ORF1 and ORF2 proteins
  orf1.orf2 <- c(orf1 = "57", 
                 orf2 = "59")
  
  # create histogram objects
  hst.bns <- seq(0, 1, length.out = 51)
  h1 <- hist(unlist(dt[, x1, with = F]), plot = F, breaks = hst.bns)
  h2 <- hist(unlist(dt[, x2, with = F]), plot = F, breaks = hst.bns)
  
  # prepare layout for figure
  layMat <- matrix(c(1, 3, 2), ncol = 1, byrow = TRUE)
  layout(layMat, heights=c(2/8, 4/8, 2/8))
  
  # plot histograms (as barplots)
  par(mar=c(.4, 4, 3, 1))
  hist(unlist(dt[, x1, with = F]), 
       breaks = hst.bns, main = "", las = 1,
       ylim = c(0, round(max(c(h1$counts, h2$counts) + 5), -1)),
       col = col.hst.sig, xaxt = "n")
  par(xpd = T)
  text(x = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), 
       y = 0, pos = 1, cex = .7, offset = .2,
       c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"),
       col = "gray20")
  par(xpd = F)
  mtext(sprintf("%s", x1.name), side = 3, line = 0.5, cex = .6)
  abline(v = unlist(dt[pg %in% orf1.orf2, x1, with = F]), 
         col = c(orf1.col, orf2.col), lwd = 2)
  
  par(mar=c(3, 4, .4, 1))
  hist(unlist(dt[, x2, with = F]), 
       breaks = hst.bns, main = "", las = 1, xaxt = "n",
       ylim = c(round(max(c(h1$counts, h2$counts) + 5), -1), 0),
       col = col.hst.sig)
  par(xpd = T)
  text(x = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), 
       y = 0, pos = 3, cex = .7, offset = .2,
       c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"),
       col = "gray20")
  par(xpd = F)
  mtext(x2.name, side = 1, line = 0.5, cex = .6)
  abline(v = unlist(dt[pg %in% orf1.orf2, x2, with = F]), 
         col = c(orf1.col, orf2.col), lwd = 2)
  
  # plot middle figure with protein links
  par(mar=c(0, 4, 0, 1))
  plot(NA, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = F, ann = F)
  
  for (i in 1:nrow(dt)) {
    lines(x = unlist(dt[i, c(x1, x2), with = F]), 
          y = c(1, 0), col = add.alpha(col.lne, .1), lwd = 1.5)
  }
  if(!is.null(highlite)){
    for (i in 1:length(highlite)) {
      lines(x = unlist(dt[pg == highlite[i], c(x1, x2), with = F]), 
            y = c(1, 0), col = highlite.col[i], lwd = 1.5)
    }
  }
  
  text.line <- 0
  text.line.y <- seq(0.9, 0.1, by = -0.2)
  text.mline <- length(text.line.y)
  if(length(highlite) > 0){
    for (i in 1:length(highlite)) {
        x.dlt <- dt[pg == highlite[i],][[x2]] - dt[pg == highlite[i],][[x1]]
        text(highlite.names[i], 
                   x = dt[pg == highlite[i],][[x1]] + x.dlt * (1 - text.line.y[text.line + 1]),
                   y = text.line.y[text.line + 1],
                   col = "black", bg = highlite.col[i], 
                   cex = .8, # theta = 0,
                   srt = atan(x.dlt) * 360 / pi - 90)
        text.line <- (text.line + 1) %% text.mline
    }
  }
  legend("right", legend = c("ORF1", "ORF2", "Link between proteins", highlite.legend$legend), 
         lwd = 2, col = c(orf1.col, orf2.col, "grey50", highlite.legend$col),
         bg = add.alpha("white", .5), box.lwd = 0)
}

