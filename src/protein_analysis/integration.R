# Integration
library(data.table)
library(seqinr)
library(plotrix)
library(VennDiagram)
library(proxy)
library(IRanges)
library(ade4)

options(scipen=999)
# source "package"
source("lib/exp-functions.R")
source("lib/affinity-obj.R")
source("lib/utils.R")
source("lib/constants.R")

# names of the experiments
exps <- list(tandem = c("SILAC_F", "SILAC_MW"), 
             rnase = c("RNAseSILAC_3", "SILACRNAse_4"), 
             rt_mut = c("401_624_Fusion", "MAP_3_624", "MAP_4_624"),
             en_mut = c("401_567_Fusion", "MAP_1_567", "MAP_2_567"),
             exchange = c("mt302_t0", "mt302_t30sec", "mt302_t5min", "mt302_30min" ))

# subset needed
dt.tmpl <- dt.tmpl[Aname %in% unlist(exps)]

# get unique - H/L pair
dt.tmpl.e <- unique(dt.tmpl[, c("Experiment", "Aname", "heavy", "use.normalized",
                                "scale"), with = F])
dt.tmpl.e$avalue <- "H/(H+L)"
dt.tmpl.e$avalue[dt.tmpl.e$heavy == 0] <- "L/(H+L)"

# prepare datasets and experiments list
l.read <- list(data = "data/quantitation/062416_proteinGroups.txt", exps = dt.tmpl.e$Experiment)

# load all exp in list
l.exp <- get_list_mqexp(l.read$data, l.read$exps)
names(l.exp) <- l.read$exps

# get affinity object for each experiment
#  remove contaminants and decoy proteins
l.aff <- lapply(l.read$exps, function(i) {
  affinity(l.exp[[i]], minpep = 2,
           avalue = dt.tmpl.e[Experiment == i]$avalue,
           score.filt = T, scale = dt.tmpl.e[Experiment == i]$scale,
           scale.sig =  5e-2, scale.sig.aff = 0.7,
           use.normalized = dt.tmpl.e[Experiment == i]$use.normalized)
})
names(l.aff) <- unique(l.read$exps)

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

# matrix with affinity
mtx <- aff.mtx[, l.read$exps]
mtx <- mtx[, dt.tmpl.e[match(unlist(exps), Aname)]$Experiment]
mtx.cell <- mtx[rownames(mtx) %in% ann.agg[is_in_cell_article == T]$pg, ]

# calculate distance
int.dist <- dist(mtx.cell, method = 'euclidean')
int.dist <- scales::rescale(int.dist, to = c(0, 0.9))
int.dist[is.na(int.dist)] <- 1
int.dist.all <- dist(mtx, method = 'euclidean')
int.dist.all[is.na(int.dist.all)] <- 1
int.hc <- hclust(int.dist)
int.hc$labels <- annotation[match(int.hc$labels, pg)]$Gene
pdf("out/integration_dendo_230517.pdf", width = 5, height = 7)
par(mar = c(4, 4, 4, 8))
plot(as.dendrogram(int.hc), horiz = T)
dev.off()
int.clusters <- cutree(int.hc, k = 6)
par(mar = c(6, 5, 5, 3))
plot(NA, xlim = c(1, 14), ylim = c(0, 1), frame = F, xaxt = "n",
     ylab = "Affinity", xlab = "", las = 2)
axis(side = 1, at = 1:length(l.read$exps), labels = l.read$exps, las = 2, 
     cex.axis = .8)
for(i in 1:nrow(mtx)){
  lines(x = 1:14, mtx[i, ], col = "gray90")
}
for(i in 1:nrow(mtx.cell)){
  i.cl <- int.clusters[annotation[match(rownames(mtx.cell[i, , drop = F]), pg)]$Gene]
  lines(x = 1:14, mtx.cell[i, ], col = cbp[[i.cl + 1]],
        lwd = 2)
}

# clusters
pdf("out/integration_clusters_230517.pdf", width = 7, height = 5)
layout(matrix(1:6, nrow = 2))
for(i in 1:max(int.clusters)){
  plot(NA, xlim = c(1, 14), ylim = c(0, 1), frame = F, xaxt = "n",
       ylab = "Affinity", xlab = "Experiments", las = 2)
  axis(side = 1, at = 1:length(l.read$exps), # labels = l.read$exps, las = 2, 
       cex.axis = .8, labels = rep("", 14))
  for(j in 1:nrow(mtx)){
    lines(x = 1:14, mtx[j, ], col = "gray90")
  }
  mtx.cell.i <- mtx.cell[rownames(mtx.cell) %in% 
                           annotation[Gene %in% names(int.clusters[int.clusters == i])]$pg, , drop = F]
  for(k in 1:nrow(mtx.cell.i)){
    lines(x = 1:14, mtx.cell.i[k, ], col = cbp[[i + 1]],
          lwd = 2)
  }
}
dev.off()

# mds
int.fit <- cmdscale(int.dist.all, eig=TRUE, k=2) # k is the number of dim
int.col <- unlist(cbp[int.clusters[annotation[match(attr(int.dist, 'Labels'), pg)]$Gene] + 1])
pdf("out/integration_mds_230517_v2.pdf", width = 6, height = 6, useDingbats = F)
par(xpd = F)
plot(int.fit$points[, 1], int.fit$points[, 2], 
     xlim = c(-1, 1.5), ylim = c(-1, 1),
     xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric MDS",	type="n", frame = F)
points(int.fit$points[, 1], int.fit$points[, 2], 
       col = "lightgray", pch = 21)
points(int.fit$points[rownames(mtx.cell), 1],
       int.fit$points[rownames(mtx.cell), 2],
       col = int.col, pch = 19)
for(i in 1:max(int.clusters)){
  pts.i <- int.fit$points[rownames(int.fit$points) %in% 
                            annotation[match(names(int.clusters[int.clusters == i]), Gene)]$pg, ]
  pts.i <- lapply(1:2, function(x) pts.i[, x])
  names(pts.i) <- c("x", "y")
  chuld <- lapply(pts.i, "[", chull(pts.i))
  polygon(chuld, border = cbp[[i + 1]], col = add.alpha(cbp[[i + 1]], .5))
}
text(int.fit$points[rownames(mtx.cell), 1], 
     int.fit$points[rownames(mtx.cell), 2],
     labels = annotation[match(rownames(mtx.cell), pg)]$Gene,
     col = int.col, 
     pos = 3, cex=.7)
dev.off()

# Why some proteins have small distance?
pr.why <- c("MOV10", "ORF1*", "PABPC4", "PABPC1", "ZCCHC3", "UPF1")
pr.col <- c(rep("red", 4), rep("blue", 2))
names(pr.col) <- annotation[match(pr.why, Gene)]$pg
plot(NA, xlim = c(1, 14), ylim = c(0, 1), frame = F, xaxt = "n",
     ylab = "Affinity", xlab = "", las = 2)
axis(side = 1, at = 1:length(l.read$exps), 
     labels = l.read$exps, las = 2, 
     cex.axis = .8)
mtx.cell.why <- mtx.cell[rownames(mtx.cell) %in% 
                         annotation[Gene %in% pr.why]$pg, , drop = F]
for(k in 1:nrow(mtx.cell.why)){
  lines(x = 1:14, mtx.cell.why[k, ], 
        col = pr.col[rownames(mtx.cell.why)[k]],
        lwd = 2)
}

# trying to work with multiple distances
d1 <- dist(mtx.cell, method = 'cosine')
d2 <- dist(mtx.cell, method = 'euclidean')
d1[is.na(d1)] <- 1
d2 <- scales::rescale(d2, to = c(0, 0.9))
d2[is.na(d2)] <- 1
d12 <- log(d1 * d2)
d12 <- scales::rescale(d12, to = c(0, 1))
d12[is.infinite(d12)] <- 0 

d12.hc <- hclust(as.dist(d12))
d12.hc$labels <- annotation[match(d12.hc$labels, pg)]$Gene
pdf("out/Integration_dendogram.pdf")
par(mar = c(4, 4, 4, 8))
plot(as.dendrogram(d12.hc), horiz = T)
dev.off()

# clusters
d12.clusters <- cutree(d12.hc, k = 5)
d12.clusters[c("HSP90AB1", "HSP90AA1", "UBB*", "TUBB4B", "TUBB", "ORF2")] <- 6
pdf("out/Integration_clusters.pdf", width = 7, height = 5, useDingbats = F)
layout(matrix(1:6, nrow = 2))
for(i in 1:max(d12.clusters)){
  plot(NA, xlim = c(1, 14), ylim = c(0, 1), frame = F, xaxt = "n",
       ylab = "Affinity", xlab = "Experiments", las = 2)
  for(j in 1:nrow(mtx)){
    lines(x = 1:14, mtx[j, ], col = "gray90")
  }
  mtx.cell.i <- mtx.cell[rownames(mtx.cell) %in% 
                           annotation[Gene %in% names(d12.clusters[d12.clusters == i])]$pg, , drop = F]
  for(k in 1:nrow(mtx.cell.i)){
    lines(x = 1:14, mtx.cell.i[k, ], col = cbp[[i + 1]],
          lwd = 2)
    for(clm in 1:ncol(mtx.cell.i)){
      if(clm == 1) {
        if(!is.na(mtx.cell.i[k, clm]) & is.na(mtx.cell.i[k, clm + 1])){
          points(x = clm, mtx.cell.i[k, clm], col = cbp[[i + 1]],
                 pch = 20)
        } 
      } else if(clm != 1 & clm != ncol(mtx.cell.i)){
        if(!is.na(mtx.cell.i[k, clm]) & is.na(mtx.cell.i[k, clm + 1]) & is.na(mtx.cell.i[k, clm - 1])) {
          points(x = clm, mtx.cell.i[k, clm], col = cbp[[i + 1]],
                 pch = 20)
        }
      } else if(clm == ncol(mtx.cell.i)) {
        if(!is.na(mtx.cell.i[k, clm]) & is.na(mtx.cell.i[k, clm - 1])){
          points(x = clm, mtx.cell.i[k, clm], col = cbp[[i + 1]],
                 pch = 20)
        }
      } 
    }
  }
  abline(v = c(2, 4, 7, 10) + .5, lty = 2)
}
dev.off()

# write output table
mtx.dt <- as.data.table(mtx)
mtx.dt$pg <- rownames(mtx)
mtx.dt <- merge(ann.agg, mtx.dt, by = 'pg')
mtx.dt <- merge(mtx.dt,
                data.table(pg = as.character(annotation[match(names(d12.clusters), Gene)]$pg),
                           Cluster = d12.clusters,
                           Color = names(d12.col)), by = "pg", all.x = T)
mtx.dt <- mtx.dt[order(Cluster, decreasing = F), 
                 .(Genes, Uniprot, Protein, SILAC_F, SILAC_MW, 
                   RNAseSILAC_3, SILACRNAse_4, `401_624_Fusion`,
                   MAP_3_567, MAP_4_567, `401_567_Fusion`, 
                   MAP_1_567, MAP_2_567, mt302_t0, mt302_t30sec, 
                   mt302_t5min, mt302_30min, Cluster, Color)]
setnames(mtx.dt, c("SILAC_F", "SILAC_MW", 
                   "RNAseSILAC_3", "SILACRNAse_4", "401_624_Fusion",
                   "MAP_3_567", "MAP_4_567", "401_567_Fusion", 
                   "MAP_1_567", "MAP_2_567", "mt302_t0", "mt302_t30sec", 
                   "mt302_t5min", "mt302_30min"),
         c("Tandem_Exp1_affinity", "Tandem_Exp2_affinity",
           "RNAse_Exp1_affinity", "RNAse_Exp2_affinity",
           "RT_Exp1_affinity", "RT_Exp2_affinity", "RT_Exp3_affinity",
           "EN_Exp1_affinity", "EN_Exp2_affinity", "EN_Exp3_affinity",
           "Exchange_0_point", "Exchange_30_seconds", "Exchange_5_minutes",
           "Exchange_30_minutes"))
mtx.dt$Genes <- gsub("\\*", "", mtx.dt$Genes)
mtx.dt$Genes <- gsub("NA, ", "", mtx.dt$Genes)
write.table(mtx.dt, 
            "out/Integration_supplementary_table.txt", 
            sep = "\t", row.names = F, quote = F)

# write output distance matrix
d12.dt <- as.data.table(round(as.matrix(d12), 5))
colnames(d12.dt) <- gsub("\\*", "", ann.agg[match(colnames(as.matrix(d12)), pg)]$Genes)
colnames(d12.dt) <- gsub("NA, ", "", colnames(d12.dt))
d12.dt$Genes <- gsub("\\*", "", ann.agg[match(rownames(as.matrix(d12)), pg)]$Genes)
d12.dt$Genes <- gsub("NA, ", "", d12.dt$Genes)
write.table(d12.dt, "out/Integration_distance_matrix.txt",
            sep = "\t", row.names = F, quote = F)
