# Time series
library(data.table)
library(seqinr)
library(plotrix)
library(proxy)
library(IRanges)

options(scipen=999)
# source "package"
source("lib/exp-functions.R")
source("lib/affinity-obj.R")
source("lib/utils.R")

# colorblind palette
cbp <- list(grey = "#999999", 
            orange = "#E69F00", 
            lightblue = "#56B4E9", 
            green = "#009E73", 
            yellow = "#F0E442", 
            blue = "#0072B2", 
            red = "#D55E00", 
            pink = "#CC79A7")

# contaminants
cont <- names(read.fasta("data/contaminants.fasta"))

# proteins from the cell article
cell.pr <- fread("data/tbls/cell_article_37_proteins.txt")
cell.pr.l <- cell.pr$Gene
names(cell.pr.l) <- cell.pr$`Uniprot Symbol`

# names of the experiments
exps.qe <- c("mt302_t0", "mt302_t30sec", "mt302_t5min", "mt302_30min" )

# read exp template (from John)
dt.tmpl <- fread("data/tbls/experimentalDesignTemplate_ed2_210417.txt", header = T, sep = "\t")

# subset needed
dt.tmpl <- dt.tmpl[Aname %in% exps.qe]

# get unique - H/L pair
dt.tmpl.e <- unique(dt.tmpl[, c("Experiment", "Aname", "heavy", "use.normalized",
                                "scale"), with = F])
dt.tmpl.e$avalue <- "H/(H+L)"
dt.tmpl.e$avalue[dt.tmpl.e$heavy == 0] <- "L/(H+L)"

# prepare datasets and experiments list
l.read <- list(data = "data/quantitation/062416_proteinGroups.txt", exps = dt.tmpl.e$Aname)

# load all exp in list
l.exp <- get_list_mqexp(l.read$data, exps.qe)
names(l.exp) <- exps.qe

# get affinity object for each experiment
#  remove contaminants and decoy proteins
l.aff <- lapply(exps.qe, function(i)
  affinity(l.exp[[i]], minpep = 2,
           avalue = dt.tmpl.e[Experiment == i]$avalue,
           score.filt = T, scale = dt.tmpl.e[Experiment == i]$scale,
           scale.sig =  5e-2, scale.sig.aff = 0.7,
           use.normalized = dt.tmpl.e[Experiment == i]$use.normalized)
)
names(l.aff) <- exps.qe

# get all proteins
prs <- unique(rbindlist(lapply(l.aff, function(x) x$aff[,.(pg, protein)])))
prs[["Uniprot Symbol"]] <- gsub(".*\\|(.*)\\|.*", "\\1", prs$protein)
v.pg <- unique(prs$pg)

# annotation
annotation <- GetUniprotAnnotation(prs$`Uniprot Symbol`)
annotation <- merge(annotation, prs, by = "Uniprot Symbol")
annotation[["protein"]] <- NULL
annotation[grep("Orf2", `Uniprot Symbol`)]$`Uniprot Symbol` <- "O00370"
annotation[grep("Orf1", `Uniprot Symbol`)]$Gene <- "ORF1"
annotation[grep("O00370", `Uniprot Symbol`)]$Gene <- "ORF2"
annotation[["is_in_cell_article"]] <- annotation$`Uniprot Symbol` %in% cell.pr$`Uniprot Symbol`
annotation <- annotation[order(pg)]
annotation <- cbind(annotation, annotation[, list(use = c(T, rep(F, length(`Uniprot Symbol`) - 1))) , by = list(pg)][, .(use)])
annotation[, is_multiple := ifelse(length(Gene) > 1, T, F) , by = list(pg)]
annotation[is_multiple == T, Gene := paste0(Gene, "*"), ]
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

# names of ORF1 and ORF2 proteins
orf1.orf2 <- c(orf1 = "sp|Q9UN81|LORF1_HUMAN", 
               orf2 = "Orf2-untagged(optimized):")

# matrix with affinity
mtx.qe <- aff.mtx[, exps.qe]
mtx.qe.comp <- mtx.qe[rowSums(is.na(mtx.qe)) < 2 & rowSums(is.na(mtx.qe[,c(1,4)])) == 0, ]
mtx.cell <- mtx.qe[rownames(mtx.qe) %in% ann.agg[is_in_cell_article == T]$pg,]
mtx.cell.dt <- as.data.table(mtx.cell)
mtx.cell.dt$pg <- rownames(mtx.cell)
mtx.cell.dt <- merge(mtx.cell.dt, ann.agg)
# protein groups
pg <- unique(rbindlist(lapply(l.aff, function(x) x$aff[,.(protein, pg)])))
pg <- pg[order(pg)]
pg <- cbind(pg, pg[, list(use = c(T, rep(F, length(protein) - 1))) , 
                   by = list(pg)][, .(use)])

# plot densities
pdf("out/Exchange_affinity_distribution.pdf", width = 7, height = 6)
plot(density(mtx.qe[,1], na.rm = T), lwd = 2, col = cbp$red, frame = F,
     main = "Distribution of H/(H+L) value",
     xlab = "H/(H+L)")
lines(density(mtx.qe[,2], na.rm = T), lwd = 2, col = cbp$yellow)
lines(density(mtx.qe[,3], na.rm = T), lwd = 2, col = cbp$green)
lines(density(mtx.qe[,4], na.rm = T), lwd = 2, col = cbp$blue)
legend("topleft", col = c(cbp$red, cbp$yellow, cbp$green, cbp$blue),
       legend = c("0 point", "30 seconds", "5 minutes", "30 minutes"),
       bty = "n", lwd = 2)
dev.off()

# plot
mtx.1 <- mtx.qe.comp[rownames(mtx.qe.comp) %nin% ann.agg[is_in_cell_article == T]$pg &
                    complete.cases(mtx.qe.comp), ]
mtx.2 <- mtx.qe.comp[rownames(mtx.qe.comp) %in% ann.agg[is_in_cell_article == T]$pg, ]

# find distance
pdf("out/Exchange_dendogram.pdf", width = 5, height = 7)
par(mar = c(4,4,2,6))
mtx.2.dist <- dist(mtx.2, method = 'cosine')
hc <- hclust(mtx.2.dist)
clusters <- cutree(hc, 3)
hc$labels <- annotation[match(hc$labels, pg)]$Gene
plot(as.dendrogram(hc), horiz = T)
dev.off()

# plot clusters
pdf("out/Exchange_lineplot.pdf", width = 7, height = 6, useDingbats = F)
par(mar = c(6,5,4,4))
plot(NA, xlim = c(0, 30 * 60), 
     ylim = c(0,1), frame = F,
     xlab = "", ylab = "H/(H+L)", xaxt = "n")
time.intervals <- c(0, 30, 5 * 60, 30 * 60)
for(i in 1:nrow(mtx.1)){
  line.col = "lightgray"
  line.type = "l"
  lines(x = time.intervals, y = mtx.1[i, ],
        col = line.col, type = line.type, lwd = .5)
}
cl.cols <- c(cbp$lightblue, cbp$pink, cbp$orange)
b.pch <- c(21, 24)
b.col = c(cbp$green, cbp$yellow)
for(i in 1:nrow(mtx.2)){
  line.col = cl.cols[clusters[rownames(mtx.2)[i]]]
  line.type = "l"
  lines(x = time.intervals, y = mtx.2[i, ],
        col = line.col, type = line.type, lwd = 2, cex = .3)
}
par(xpd = T)
text(x = 30 * 60, y = mtx.2[, 4], pos = 4,
     labels = annotation[match(rownames(mtx.2), pg)]$Gene, 
     col = cl.cols[clusters[rownames(mtx.2)]], cex = .7)
axis(1, at = c(0, 30, 5 * 60, 30 * 60), las = 2, cex.axis = 0.7,
     labels = c("0", "30 seconds", "5 minutes", "30 minutes"))
legend("topright", lwd = 2, col = cl.cols, 
       legend = c("Cluster 1", "Cluster 2", "Cluster 3"),
       bty = "n")
dev.off()

# write output
mtx.dt <- as.data.table(rbind(mtx.1, mtx.2))
mtx.dt$pg <- rownames(rbind(mtx.1, mtx.2))
mtx.dt <- merge(ann.agg[,.(Genes, Uniprot, Protein, pg, is_in_cell_article)], 
                mtx.dt, by = "pg")
mtx.dt$Genes <- gsub("\\*", "", mtx.dt$Genes )
mtx.dt$Color <- "grey"
mtx.dt[match(names(clusters), pg)]$Color <- c("blue", "pink", "orange")[clusters]
setnames(mtx.dt, c("mt302_t0", "mt302_t30sec", "mt302_t5min", "mt302_30min"),
         c("0_point", "30_seconds", "5_minutes", "30_minutes"))
mtx.dt <- mtx.dt[order(is_in_cell_article, Color, decreasing = T)]
mtx.dt <- mtx.dt[,.(Genes, Uniprot, Protein, Color,
                    `0_point`, `30_seconds`, `5_minutes`, `30_minutes`)]
write.table(mtx.dt, "out/Exchange_supplementary_table.txt",
            sep = "\t", row.names = F, quote = F)

# write output distance matrix
out.dist.mtx <- as.matrix(mtx.2.dist)
out.dist.dt <- as.data.table(out.dist.mtx)
colnames(out.dist.dt) <- gsub("\\*", "", ann.agg[match(colnames(out.dist.mtx), pg)]$Genes)
out.dist.dt$Genes <- gsub("\\*", "", ann.agg[match(rownames(out.dist.mtx), pg)]$Genes)
write.table(out.dist.dt, "out/Exchange_distance_matrix.txt",
            sep = "\t", row.names = F, quote = F)
