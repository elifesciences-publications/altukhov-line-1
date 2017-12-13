# split-tandem analysis
library(data.table)
library(seqinr)

options(scipen=999)

# source "package"
source("lib/exp-functions.R")
source("lib/utils.R")
source("lib/affinity-obj.R")

# contaminants
cont <- names(read.fasta("data/contaminants.fasta"))

# colorblind palette
cbp <- list(grey = "#999999", 
            orange = "#E69F00", 
            lightblue = "#56B4E9", 
            green = "#009E73", 
            yellow = "#F0E442", 
            blue = "#0072B2", 
            red = "#D55E00", 
            pink = "#CC79A7")


# proteins from the cell paper
cell.pr <- fread("data/tbls/cell_article_37_proteins.txt")
cell.pr.l <- cell.pr$Gene
names(cell.pr.l) <- cell.pr$`Uniprot Symbol`

# names of the experiments
exps <- c("SILAC_F", "SILAC_MW")

# read exp template (from John)
dt.tmpl <- fread("data/tbls/experimentalDesignTemplate_ed2_210417.txt", header = T, sep = "\t")

# subset needed
dt.tmpl <- dt.tmpl[Aname %in% exps]

# get unique - H/L pair
dt.tmpl.e <- unique(dt.tmpl[, c("Experiment", "Aname", "heavy", "use.normalized",
                                "scale"), with = F])

# get affinity formula
dt.tmpl.e$avalue <- "H/(H+L)"
dt.tmpl.e$avalue[dt.tmpl.e$heavy == 0] <- "L/(H+L)"

# prepare datasets and experiments list
l.read <- list(data = "data/quantitation/062416_proteinGroups.txt", exps = dt.tmpl.e$Experiment)

# load all exp in list
l.exp <- get_list_mqexp(l.read$data, exps)
names(l.exp) <- exps

# get affinity object for each experiment
#  remove contaminants and decoy proteins
l.aff <- lapply(exps, function(i)
  affinity(l.exp[[i]], minpep = 2,
           avalue = dt.tmpl.e[Aname == i]$avalue,
           score.filt = T, scale = dt.tmpl.e[Aname == i]$scale,
           scale.sig =  5e-2, scale.sig.aff = 0.7,
           use.normalized = dt.tmpl.e[Aname == i]$use.normalized)
)
names(l.aff) <- exps

# selecting mov10
l.aff.mov10 <- lapply(exps, function(i)
  affinity(l.exp[[i]], minpep = 1,
           avalue = dt.tmpl.e[Aname == i]$avalue,
           score.filt = T, scale = dt.tmpl.e[Aname == i]$scale,
           scale.sig =  5e-2, scale.sig.aff = 0.7,
           use.normalized = dt.tmpl.e[Aname == i]$use.normalized)
)
names(l.aff.mov10) <- exps
# create empty matrix
aff.mov10.mtx <- matrix(NA, nrow = 1,
                        ncol = length(l.aff.mov10),
                        dimnames = list("mov10", names(l.aff)))
# fill with values
for (i in 1:length(l.aff.mov10)) {
  aff.mov10.mtx["mov10", names(l.aff.mov10)[i]] <-
    l.aff.mov10[[i]]$aff[grep("Q9HCE1", protein)]$aff
}

# get all proteins
prs <- unique(rbindlist(lapply(l.aff, function(x) x$aff[,.(pg, protein)])))
prs[["Uniprot Symbol"]] <- gsub(".*\\|(.*)\\|.*", "\\1", prs$protein)
v.pg <- unique(prs$pg)

# annotation
annotation <- GetUniprotAnnotation(prs$`Uniprot Symbol`)
annotation <- merge(annotation, prs, by = "Uniprot Symbol")
annotation[["protein"]] <- NULL
annotation[["is_in_cell_article"]] <- annotation$`Uniprot Symbol` %in% cell.pr$`Uniprot Symbol`
annotation[grep("Orf2", `Uniprot Symbol`)]$`Uniprot Symbol` <- "O00370"
annotation[grep("Orf1", `Uniprot Symbol`)]$Gene <- "ORF1"
annotation[grep("O00370", `Uniprot Symbol`)]$Gene <- "ORF2"
annotation <- annotation[order(pg)]
annotation <- cbind(annotation, annotation[, list(use = c(T, rep(F, length(`Uniprot Symbol`) - 1))) , by = list(pg)][, .(use)])
annotation[, is_multiple := ifelse(length(Gene) > 1, T, F) , by = list(pg)]
annotation$Gene_w <- annotation$Gene
annotation[is_multiple == T, Gene := paste0(Gene, " *"), ]
annotation <- rbind(annotation,
                    data.table("Uniprot Symbol" = "Q9HCE1",
                         Gene = 'MOV10†',
                         Gene_w = 'MOV10',
                         Protein = cell.pr[Gene == "MOV10"]$Protein,
                         pg = "mov10",
                         is_in_cell_article = T,
                         use = T,
                         is_multiple = F))
annotation <- annotation[`Uniprot Symbol` %nin% cont]
ann.agg <- annotation[, list(Genes = paste(Gene[!is.na(Gene)], collapse = ", "),
                             Uniprot = paste(`Uniprot Symbol`[!is.na(`Uniprot Symbol`)], collapse = ", "),
                             Protein = paste(Protein[!is.na(Protein)], collapse = ", "),
                             is_in_cell_article = c(F, T)[max(is_in_cell_article) + 1]), 
                      by = list(pg)]
ann.agg$pg <- as.character(ann.agg$pg)

# create empty matrix
aff.mtx <- matrix(NA, nrow = length(v.pg),
                  ncol = length(l.aff),
                  dimnames = list(v.pg, names(l.aff)))
# fill with values
for (i in 1:length(l.aff)) {
  aff.mtx[as.character(l.aff[[i]]$aff$pg), names(l.aff)[i]] <-
    l.aff[[i]]$aff$aff
}
aff.mtx <- rbind(aff.mtx, aff.mov10.mtx)
aff.mtx <- aff.mtx[rownames(aff.mtx) %in% annotation$pg, ]

# names of ORF1 and ORF2 proteins
orf1.orf2 <- c(orf1 = "sp|Q9UN81|LORF1_HUMAN", 
               orf2 = "Orf2-untagged(optimized):")

# matrix with affinity
mtx <- aff.mtx[, exps]

# protein groups
pg <- unique(rbindlist(lapply(l.aff, function(x) x$aff[,.(protein, pg)])))
pg <- pg[order(pg)]
pg <- cbind(pg, pg[, list(use = c(T, rep(F, length(protein) - 1))) , by = list(pg)][, .(use)])

# normalization
mtx.norm <- sapply(exps, function(exp) {
  exp.v <- mtx[, exp]
  exp.v.orf1 <- exp.v[as.character(pg[protein == orf1.orf2[1]]$pg)]
  exp.v.median <- median(exp.v, na.rm = T)
  exp.solv <- solve(matrix(c(exp.v.orf1, exp.v.median, c(1, 1)), nrow = 2, byrow = F), 
                    c(1, exp.v.median))
  exp.v.norm <- exp.solv[1] * exp.v + exp.solv[2]
  return(exp.v.norm)
})

# plot 2nd version
pdf("out/Split-tandem_normalization.pdf", width = 8, height = 4)
layout(matrix(1:2, nrow = 1))
par(mar = c(5,4,4,1))
h1.1 <- hist(mtx[,1], breaks = 30, plot = F)
h1.2 <- hist(mtx[,2], breaks = 20, plot = F)
h2.1 <- hist(mtx.norm[,1], breaks = 20, plot = F)
h2.2 <- hist(mtx.norm[,2], breaks = 20, plot = F)

plot(h1.1$mids, h1.1$counts, lwd = 2, type = "s",
     col = cbp$red, ylim = c(0, 30), xlim = c(0, 1),
     xlab = "ORF1p Co-Partitioning", ylab = "Number of proteins",
     frame = F)
par(new = T)
plot(h1.2$mids, h1.2$counts, lwd = 2, type = "s",
     col = cbp$blue, frame = F, xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 35))

plot(h2.1$mids, h2.1$counts, lwd = 2, type = "s",
     col = cbp$red, ylim = c(0, 30), xlim = c(0, 1),
     xlab = "Normalized ORF1p Co-Partitioning", ylab = "Number of proteins",
     frame = F)
par(new = T)
plot(h2.2$mids, h2.2$counts, lwd = 2, type = "s",
     col = cbp$blue, frame = F, xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 30))
dev.off()

# plot scatter plot
mtx.norm.2 <- mtx.norm
mtx.norm.2 <- mtx.norm.2[rowSums(is.na(mtx.norm.2)) == 0, ]
orf12.pt <- mtx.norm.2[rownames(mtx.norm.2) %in% 
                         annotation[`Uniprot Symbol` %in% cell.pr[co_IP_with == "ORF1/2"]$`Uniprot Symbol`]$pg, ]
orf2.pt <- mtx.norm.2[rownames(mtx.norm.2) %in% 
                        annotation[`Uniprot Symbol` %in% cell.pr[co_IP_with == "ORF2"]$`Uniprot Symbol`]$pg, ]
orf1.pt <- mtx.norm.2[rownames(mtx.norm.2) %in% 
                        annotation[`Uniprot Symbol` %in% cell.pr[co_IP_with == "ORF1"]$`Uniprot Symbol`]$pg, ]

pdf("out/Split-tandem_scatter_plot.pdf", width = 6, height = 6, useDingbats=FALSE)
plot(mtx.norm.2[rownames(mtx.norm.2) %nin% 
                  as.character(annotation[`Uniprot Symbol` %in% cell.pr$`Uniprot Symbol`]$pg), ], 
     frame = F, xlab = "ORF1p Co-Partitioning (Exp. 1)",
     ylab = "ORF1p Co-Partitioning (Exp. 2)", 
     pch = 20, col = cbp$grey, xlim = c(-.15, 1), ylim = c(-.15, 1))

points(orf12.pt, pch = 20, col = cbp$lightblue, cex = 1.5)
points(orf2.pt, pch = 20, col = cbp$pink, cex = 1.5)
points(orf1.pt, pch = 20, col = cbp$orange, cex = 1.5)
text(orf12.pt, col = cbp$lightblue, pos = 2,
     labels = annotation[pg %in% rownames(orf12.pt) & use]$Gene)
text(orf2.pt[1:5, ], col = cbp$pink, pos = 2, 
     labels = annotation[pg %in% rownames(orf2.pt[1:5, ]) &
                           is_in_cell_article]$Gene)
text(orf2.pt[6:nrow(orf2.pt), ], col = cbp$pink, pos = 4, 
     labels = annotation[pg %in% rownames(orf2.pt[6:nrow(orf2.pt), ]) &
                           is_in_cell_article]$Gene)

legend("top", legend = c(expression(paste("Significant in ", alpha, 
                                          "-ORF1p & ", alpha, "-ORF2p I-DIRT")),
                         expression(paste("Significant in ", alpha, 
                                          "-ORF2p I-DIRT only"))),
       pch = 20, col = c(cbp$lightblue, cbp$pink), pt.cex = 1.5, bty = "n")


dev.off()

# 5. Calculate the distances between node pairs and associated probabilities.
# PURA/PURB/PCNA
# PABPC1/PABPC4
# HSPA8/HSPA1A
# NAP1L1/IPO7

dt.norm <- as.data.table(mtx.norm.2)
dt.norm$pg <- rownames(mtx.norm.2)
dt.norm <- merge(dt.norm, ann.agg, by = "pg")
pr.pairs <- as.data.table(t(combn(dt.norm$pg, m = 2)))
pr.pairs$dist <- sapply(1:nrow(pr.pairs), function(x){
  a <- dt.norm[pg == pr.pairs[x]$V1]
  b <- dt.norm[pg == pr.pairs[x]$V2]
  ab.dist <- sqrt((a$SILAC_F - b$SILAC_F) ^ 2 + (a$SILAC_MW - b$SILAC_MW) ^ 2)
  ab.dist
})
pr.pairs.list <- list(c("PABPC1", "PABPC4"),
                      c("HSPA8", "HSPA1A"),
                      c("NAP1L1", "IPO7"))
pr.pairs.pgs <- lapply(pr.pairs.list, function(x){
  annotation[Gene_w %in% x]$pg
})
pr.pairs.dist <- sapply(pr.pairs.pgs, function(x){
  pr.pairs[V1 %in% x & V2 %in% x]$dist
})
pr.pairs.prob <- sapply(pr.pairs.dist, function(x){
  sum(pr.pairs$dist[pr.pairs$dist <= x]) / sum(pr.pairs$dist)
})
pr.triples <- as.data.table(t(combn(dt.norm$pg, m = 3)))
pr.triples <- merge(pr.triples, pr.pairs, by = c("V1", "V2"))
pr.triples <- merge(pr.triples, pr.pairs, 
                    by.x = c("V1", "V3"),
                    by.y = c("V1", "V2"))
pr.triples <- merge(pr.triples, pr.pairs, 
                    by.x = c("V2", "V3"),
                    by.y = c("V1", "V2"))
pr.triples$triple.dist <- rowMeans(pr.triples[,.(dist, dist.x, dist.y)])

# plot hist
tr.hist <- hist(pr.triples$triple.dist, breaks = 200)
plot(tr.hist$mids, tr.hist$counts, type = "s", frame = F,
     lwd = 2, main = "Histogram of the average distance\nbetween nodes of all three-node groups",
     xlab = "Average distance", ylab = "Number of node groups")
pr.triples.list <- c("PURA", "PURB", "PCNA")
pr.triples.pgs <- annotation[Gene %in% pr.triples.list]$pg
abline(v = pr.triples[V1 %in% pr.triples.pgs &
                        V2 %in% pr.triples.pgs &
                        V3 %in% pr.triples.pgs]$triple.dist,
       col = "red", lwd = 2, lty = 2)
pr.triples.prob <- pr.triples[V1 %in% pr.triples.pgs &
                                V2 %in% pr.triples.pgs &
                                V3 %in% pr.triples.pgs]$triple.dist / 
  sum(pr.triples$triple.dist)


dup.hist <- hist(pr.pairs$dist, breaks = 100, plot = F)
plot(dup.hist$mids, dup.hist$counts, type = "s", frame = F,
     lwd = 2, main = "Histogram of the distance\nbetween nodes of all two-node groups",
     xlab = "Distance", ylab = "Number of node groups")

# write output table
dt.out <- dt.norm[,.(Genes, Uniprot, Protein, SILAC_F, SILAC_MW, pg)]
dt.out$Avg_affinity <- rowMeans(dt.out[,.(SILAC_F, SILAC_MW)])
dt.out$color <- sapply(dt.out$pg, function(x){
  if(x %in% rownames(orf12.pt)) {
    return("blue")
  } else if(x %in% rownames(orf2.pt)){
    return("pink")
  } else {
    return("gray")
  }
})
dt.out <- dt.out[order(Avg_affinity, decreasing = T)]
dt.out <- rbind(dt.out[color == "blue"],
                dt.out[color == "pink"],
                dt.out[color == "gray"])
dt.out$pg <- NULL
setnames(dt.out, c("SILAC_F", "SILAC_MW", "color"),
         c("Exp1_affinity", "Exp2_affinity", "Color"))
dt.out <- dt.out[,.(Genes, Uniprot, Protein, Color, 
                    Exp1_affinity, Exp2_affinity, Avg_affinity)]
dt.out$Genes <- gsub(" \\*", "", dt.out$Genes)
dt.out[, c("Exp1_affinity", "Exp2_affinity", "Avg_affinity")] <-
  round(dt.out[, .(Exp1_affinity, Exp2_affinity, Avg_affinity)], 4)
write.table(dt.out, "out/Split-tandem_supplementary_table.txt", sep = "\t", row.names = F,
            quote = F)
