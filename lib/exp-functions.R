# New approach to experiment control
# experiments as separate objects
# Author: Dima Ischenko
# Date: 11.06.2016

# *******
# Usefull and helpfull functions
# *******

#' Simple boolean `not in` function.
#' @export
"%ni%" <- Negate("%in%")

#' Length of unique elements in vector.
#' 
#' @param x Vector with elements.
#' @return Length of unique elements.
#' @export
ulength <- function(x) length(unique(x))

#' Chauvenet filtration of distribution.
#'
#' @param Vector with data.
#' @param Constant for normalization.
#' @return Vector without filtered values.
#' @export
ChauvenetFilter <- function(x, th = 0.5) {
  x.mean <- mean(x)
  x.sd <- sd(x)
  
  xpr <- pnorm(x, x.mean, x.sd, lower.tail = F) * length(x)
  xpl <- pnorm(x, x.mean, x.sd, lower.tail = T) * length(x)
  return(x[xpr > th & xpl > th])
}

# *******
# Main artillery
# *******

#' Load experiment from file or data frame (data.table)
#' create mqexp object.
#' 
#' @param source Path to file or data.table
#' @param exp.name Experiment name.
#' @return mqexp object.
#' @export
mqexp <- function(source, exp.name, ...) UseMethod("mqexp")

#' Read file and then create mqexp object
#'
#' @export
mqexp.character <- function(source, exp.name) {
  require(data.table)
  
  # load file data
  dt.a <- fread(source, header = T)
  # colnames to upper case
  setnames(dt.a, toupper(colnames(dt.a)))
  
  # call mqexp for data table
  return(mqexp(dt.a, exp.name))
}

#'
#'
#' @export
mqexp.data.frame <- function(source, exp.name) {
  # experiment name to upper case
  exp.Name <- toupper(exp.name)
  
  old.names <- c("SCORE", "POTENTIAL CONTAMINANT")
  new.names <- c("PEP", "CONTAMINANT")
  for(i in 1:length(old.names)){
    if(old.names[i] %in% colnames(source)){
      setnames(source, old.names[i], new.names[i])
    }
  }
  
  # generate needed columns
  v.col <- c("MAJORITY PROTEIN IDS",
             sprintf(fmt = c("PEPTIDES %s", "RAZOR + UNIQUE PEPTIDES %s",
                             "RATIO H/L %s", "RATIO H/L NORMALIZED %s",
                             "RATIO H/L VARIABILITY [%%] %s"), exp.Name),
             "PEP", "CONTAMINANT", "REVERSE"
  )
  # report error if not all columns present in data file
  if (!all(v.col %in% colnames(source))) {
    stop(sprintf("Can not load %s experiment", exp.Name))
  }
  # subset needed columns
  dt.s <- source[, v.col, with = F]
  
  # split protein group to proteins
  # get list with splitted proteins for each protein group
  l.prot <- sapply(dt.s[["MAJORITY PROTEIN IDS"]], function(x) {
    x.u <- unlist(strsplit(x, split = ';'))
  })
  
  # generate final data
  dt.f <- dt.s[rep(1:nrow(dt.s), sapply(l.prot, length)),]
  dt.f[["PG"]] <- rep(1:nrow(dt.s), sapply(l.prot, length))
  # genrate id for proteins in PG
  dt.f[["PID"]] <- unlist(sapply(l.prot, function(x) seq(1, length(x), 1)))
  dt.f[["PROTEIN"]] <- unlist(l.prot)
  dt.f[["MAJORITY PROTEIN IDS"]] <- NULL
  
  # rename columns
  setnames(dt.f, c("peps", "rupeps", "hl", "hln", "hlv",
                   "score", "contaminant", "reverse", "pg", "pid",
                   "protein"))
  # manually convert needed columns to numeric
  for (col in c("peps", "rupeps", "hl", "hln", "hlv", "score")) {
    dt.f[[col]] <- as.numeric(dt.f[[col]])
  }
  
  # remove NA values
  dt.f <- dt.f[!is.na(dt.f$hl) & !is.na(dt.f$hln), ]
  
  # create and return mqexp object
  exp <- list(name = exp.name, data = dt.f)
  class(exp) <- "mqexp"
  return(exp)
}

#'
#'
#' @export
mqexp.default <- function(data, ...) {
  stop(sprintf("This class (%s) of source doesn't supported.",
               class(data)))
}

#' Function to load list of mqexps
#' 
#' @param source Path to file
#' @param exp.names Vector with names of experiments
#' @return List with "mqexp" objects
#' @export
get_list_mqexp <- function(source, exp.names) {
  require(data.table)
  
  # load file data
  dt.a <- fread(source, header = T)
  # colnames to upper case
  setnames(dt.a, toupper(colnames(dt.a)))
  
  return(lapply(exp.names, function(exp.name)
    mqexp(dt.a, exp.name)))
}
