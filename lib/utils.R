GetUniprotAnnotation <- function(x) {
  require(data.table)
  require(UniProt.ws)
  up <- UniProt.ws()
  up.kt <- "UNIPROTKB"
  up.columns <- c("UNIPROTKB", "GENES", "PROTEIN-NAMES")
  up.dt <- as.data.table(
    UniProt.ws::select(x = up, keys = x, columns = up.columns, keytype = up.kt))
  up.dt$`Protein` <- gsub("(.*?) \\(.*", "\\1", up.dt$`PROTEIN-NAMES`)
  up.dt$Gene <- gsub("(.*?) .*", "\\1", up.dt$GENES)
  setnames(up.dt, "UNIPROTKB", "Uniprot Symbol")
  return(up.dt[, .(`Uniprot Symbol`, Gene, Protein)])
}

'%nin%' <- Negate("%in%")

ulen <- function(x) length(unique(x))

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}
