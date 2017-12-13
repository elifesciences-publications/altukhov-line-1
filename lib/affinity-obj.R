# Affinity class introduced to work with 
# cleaned MaxQuant results, wihtout contaminants
# normalized, converted H/L to needed value e.t.c.
# Author: Dima Ischenko
# Date: 14.06.2016

#' Convert ms experiment object to affinity object
#' 
#' @param exp MS experiment object
#' @note Now supported only "mqexp" objects
#' 
#' @return Affinity object
#' @export
affinity <- function(exp, ...) UseMethod("affinity")

#'
#' @param minpep Minimal number of unique peptides to consider protein
#' @param avalue Which value use to calculate affinity:
#'   H/L, L/H, H/(H+L) or L/(H+L) (!not used now)
#' @param use.normalized Use or not Normalized ratios of H/L to
#'   calculate affinity
#' @param score.filt Filter proteins by Reverse search and score
#' @param protnames Specified names of selected proteins from
#'   mqexp object, if protnames is empty character ("") subset
#'   first protein for each protein group
#' @export
affinity.mqexp <- function(exp, minpep, avalue,
                           use.normalized,
                           score.filt,
                           protnames = c(""),
                           scale,
                           scale.sig = 5e-2,
                           scale.sig.aff = 0.7, 
                           select.first = F,...) {
  # check avalue
  posv <- c("H/L", "L/H", "H/(H+L)", "L/(H+L)")
  if (! avalue %in% posv) {
    stop("Wrong avalue type, possible: H/L, L/H, H/(H+L) or L/(H+L)")
  }
  # get data
  dt <- exp[["data"]]
  # if filter proteins by Reverse search and score
  if (score.filt) {
    # calc maximal decoy score
    mx.pep <- max(dt[["score"]][dt[["reverse"]] == "+"])
    # filter experiment data
    dt <- dt[dt$score > mx.pep, ]
  }
  dt <- dt[dt$rupeps >= minpep & dt$contaminant == "" & dt$reverse == "", ]
  # prepare affinity
  dt$val <- ifelse(rep(use.normalized, nrow(dt)),
                   dt$hln, dt$hl)
  if (avalue %in% c("L/H", "L/(H+L)")) dt$val <- 1 / dt$val
  if (avalue %in% c("H/(H+L)", "L/(H+L)")) dt$val <- (dt$val^-1 + 1)^-1
  dt$aff <- dt$val
  
  # subset proteins
  # wtf below?
  if(select.first){
    if (protnames == "") {
      dt <- dt[, .SD[1], by = list(pg)]
    } else {
      dt <- dt[dt$pid == 1, ]
    }
  }
  
  aff <- list(name = exp[["name"]],
              exp = exp, 
              aff = dt[, c("pg", "protein", "aff"), with = F],
              avalue = avalue,
              use.normalized = use.normalized,
              minpep = minpep,
              score.filt = score.filt)
  class(aff) <- c("affinity")
  
  aff[["scaled"]] = F
  aff[["scaled.par"]] <- c(0, 0)
  
  if (scale) {
    pv <- signif_test(aff, apply.chauven = T,
                lower.tail = F, add.plot = F)
    hl <- min(min(aff$aff$aff[pv <= scale.sig]), 1)
    
    tval <- aff$aff$aff
    # Chauvenet filtration
    n.it <- 0
    max.it = 10
    ch.th = 2
    while(n.it < max.it) {
      tval.n <- ChauvenetFilter(tval, th = ch.th)
      if (length(tval.n) == length(tval)) break
      tval <- tval.n
      n.it <- n.it + 1
    }
    afmn <- mean(tval)
    
    #afmn <- mean(aff$aff$aff)
    aff$aff$aff <- aff$aff$aff - afmn
    hl <- hl - afmn
    aff$aff$aff <- aff$aff$aff * (scale.sig.aff - 0.5)/hl
    aff$aff$aff <- aff$aff$aff + 0.5
    aff[["scaled"]] = T
    aff[["scaled.par"]] <- c(0.5, scale.sig.aff)
  }
  
  return(aff)
}

#' @export
affinity.default <- function(exp, ...) {
  stop(sprintf("This class (%s) of experiment doesn't supported.",
       class(exp)))
}