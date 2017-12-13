#' helpfull function for cumulative normal probability
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)

#' Check for significance for each protein in affinity object
#'
#' @param dt Affinity object
#' 
#' @export
signif_test <- function(dt, ...) UseMethod("signif_test")

#' 
#' @param apply.chauven Logical, apply or not Chauvenet filtration for
#'   affinity value vector
#' @param add.plot Plot histogram with significance
#' @param plot.sign Significance to mark values at histogram
#' 
#' @return vector with significance for each protein
#' @export
signif_test.affinity <- function(dt, apply.chauven,
                                 add.plot, plot.sign = 1e-2, lower.tail = F,
                                 p.col = "dodgerblue", scale = F,
                                 plot.main = "Affinity histogram", ...) {
  # get affinity vector
  tval <- dt[["aff"]]$aff
  # convert h/h+l to log h/l values
  tval <- log((tval^-1 - 1)^-1)
  
  # distribution quantilles
  r0 <- median(tval, na.rm = T)
  r1 <- sort(tval)[length(tval) * 0.8413]
  rm1 <- sort(tval)[length(tval) * 0.1587]
  
  # calc z
  if(lower.tail) {
    z <- (r0 - tval)/(r0 - rm1)
  } else {
    z <- (tval - r0)/(r1 - r0)
  }
  
  # calc significance
  sig.a <- 1/2 * erfc(z / sqrt(2))
  
  # adjusted p-value
  pv <- p.adjust(sig.a, method = "BH")
  names(pv) <- dt[["aff"]]$protein
  
  dt.mean <- mean(tval)
  dt.sd <- sd(tval)
  
  #pv <- p.adjust(pnorm(q = dt[["aff"]]$aff, mean = dt.mean, sd = dt.sd,
  #                     lower.tail = lower.tail), method = "BH")
  #names(pv) <- dt[["aff"]]$protein
  
  if (add.plot) {
    # add histogram
    pv.th <- plot.sign
    #hl.th <- qnorm(1 - pv.th, mean = dt.mean, sd = dt.sd)
    if (lower.tail) {
      hl.th <- max(dt[["aff"]]$aff[pv <= pv.th])
    } else {
      hl.th <- min(dt[["aff"]]$aff[pv <= pv.th])
    }
    
    if (scale) {
      if(abs(hl.th) == Inf) hl.th = 1
      dt[["aff"]]$aff <- (dt[["aff"]]$aff + (0.5 - dt.mean))*0.7/(hl.th + 0.5 - dt.mean)
      dt.mean <- mean(dt[["aff"]]$aff)
      dt.sd <- sd(dt[["aff"]]$aff)
      hl.th <- 0.7
    }
    
    v.brks <- seq(min(dt[["aff"]]$aff), max(dt[["aff"]]$aff), length.out = 50)
    v.clr <- rep("white", 50)
    if (lower.tail) {
      v.fdb <- max(which(v.brks <= hl.th))
      if (v.fdb == -Inf) { v.fdb <- 1}
      v.clr[1:(v.fdb - 1)] <- p.col
    } else {
      v.fdb <- min(which(v.brks >= hl.th))
      if (v.fdb == Inf) { v.fdb <- length(v.brks)}
      v.clr[(v.fdb):length(v.clr)] <- p.col
    }
    hist(dt[["aff"]]$aff, prob = T,
         main = plot.main, xlab = "Affinity",
         breaks = v.brks, col = v.clr)
    #lines(density(rnorm(n = 10*length(dt[["aff"]]$aff),
    #                    mean = dt.mean, sd = dt.sd)),
    #      col = "red", lwd = 2)
    abline(v = hl.th, lty = 2, col = "gray30")
  }
  
  return(pv)
}

#' @export
signif_test.st.single <- function(dt, apply.chauven,
                                  lower.tail = F, pv.th = 1, ...) {
  tval <- dt[["data"]]
  pv <- signif_test(tval, apply.chauven, lower.tail,
                    add.plot = F)
  return(tval[["aff"]][tval[["aff"]]$protein %in%
                         names(pv[pv <= pv.th])])
}


#' @export
signif_test.st.mix <- function(dt, apply.chauven,
                               lower.tail = F, pv.th = 1, ...) {
  pv1 <- signif_test(dt[["e1"]], apply.chauven, lower.tail,
                     add.plot = F)
  pv2 <- signif_test(dt[["e2"]], apply.chauven, lower.tail,
                     add.plot = F)
  pvm1 <- signif_test(dt[["m1"]], apply.chauven, lower.tail,
                      add.plot = F)
  pvm2 <- signif_test(dt[["m2"]], apply.chauven, lower.tail,
                      add.plot = F)
  dta <- merge(dt[["e1"]][["aff"]], dt[["e2"]][["aff"]], by = "protein")
  dta$sig1 <- dta$protein %in% names(pv1[pv1 <= pv.th])
  dta$sig2 <- dta$protein %in% names(pv2[pv2 <= pv.th])
  dta$exc1 <- dta$protein %in% names(pvm1[pvm1 <= pv.th])
  dta$exc2 <- dta$protein %in% names(pvm2[pvm2 <= pv.th])
  dta <- dta[(dta$sig1 & !dta$exc1) | (dta$sig2 & !dta$exc2), ]
  dta$aff <- (dta[["aff.x"]] + dta[["aff.y"]])/2
  
  return(dta[, c("protein", "aff"), with = F])
}


#' @export
signif_test.default <- function(dt, ...) {
  stop(sprintf("This class (%s) of experiment doesn't supported.",
               class(dt)))
}