# RNAse
library(data.table)
library(seqinr)
library(plotrix)

options(scipen=999)
# source "package"
source("lib/exp-functions.R")
source("lib/affinity-obj.R")
source("lib/utils.R")
source("lib/constants.R")

# names of the experiments
exps <- c("RNAseSILAC_3", "SILACRNAse_4")

# subset needed
dt.tmpl <- dt.tmpl[Aname %in% exps]
# get unique - H/L pair
dt.tmpl.e <- unique(dt.tmpl[, c("Experiment", "Aname", "heavy", "use.normalized",
                                "scale"), with = F])
dt.tmpl.e$avalue <- "H/(H+L)"
dt.tmpl.e$avalue[dt.tmpl.e$heavy == 0] <- "L/(H+L)"

# prepare datasets and experiments list
l.read <- list(data = "data/quantitation/proteinGroups.txt", exps = dt.tmpl.e$Experiment)

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
mtx <- aff.mtx[, exps]

# protein groups
pg <- unique(rbindlist(lapply(l.aff, function(x) x$aff[,.(protein, pg)])))
pg <- pg[order(pg)]
pg <- cbind(pg, pg[, list(use = c(T, rep(F, length(protein) - 1))) , by = list(pg)][, .(use)])

# normalization ----
mtx.norm <- 2 * t(t(mtx) - colMeans(mtx, na.rm = T))
mtx.norm <- mtx.norm[gsub("sp\\|(.*)\\|.*", "\\1", rownames(mtx.norm)) %nin% 
                       as.character(annotation[`Uniprot Symbol` %in% cont]$pg), ]

# plot 1st version ----
pdf("out/RNAse_normalization.pdf", width = 8, height = 5)
layout(matrix(1:2, nrow = 1))
par(mar = c(5,5,8,1))
h1.1 <- hist(mtx[,1], breaks = 30, ylim = c(0,35),
     col = add.alpha(cbp$red,.3), xlim = c(0, 1),
     xlab = "RNase sensitivity", ylab = "Number of proteins",
     main = "")
h1.2 <- hist(mtx[,2], breaks = 40, add = T,
     col = add.alpha(cbp$blue, .3))
abline(v = colMeans(mtx, na.rm = T), lwd = 2, col = c(cbp$red, cbp$blue))
legend("topright", bty = "n", inset = c(.55, -.4), xpd = T,
       fill = c(add.alpha(cbp$red, .3), add.alpha(cbp$blue, .3)),
       legend = c("1st replicate", "2nd replicate"))
legend("topright", bty = "n", lwd = 2, col = c(cbp$red, cbp$blue),
       legend = c("Mean sensitivity for 1st replicate",
                  "Mean sensitivity for 2nd replicate"), 
       inset = c(.2, -.25), xpd = T)

h2.1 <- hist(mtx.norm[,1], breaks = 20, ylim = c(0,40),
     col = add.alpha(cbp$red,.3), xlim = c(-1, 1),
     xlab = "Normalized RNase sensitivity", ylab = "Number of proteins",
     main = "")
h2.2 <- hist(mtx.norm[,2], breaks = 20, add = T,
     col = add.alpha(cbp$blue, .3))
dev.off()


# plot distribution of affinity-value ----
mtx.norm.2 <- mtx.norm
mtx.norm.2 <- mtx.norm.2[rowSums(is.na(mtx.norm.2)) == 0, ]

# trying to select threshold
sens <- sqrt(mtx.norm.2[, 1] ^ 2 + mtx.norm.2[, 2] ^ 2)
sens <- sens[sens < median(sens) * 2]
h0.1 <- hist(sens, breaks = 20, plot = F)
h0.2 <- hist(rnorm(n = 10000, sd = sd(sens), 
                   mean = mean(sens)), breaks = 20, plot = F)
pdf("out/RNAse_normality_test.pdf", width = 7, height = 4)
layout(matrix(1:2, nrow = 1))
# distribution of RNAse sensitivivty
plot(h0.1$mids, h0.1$counts, type = "s", lwd = 2, col = cbp$lightblue, 
     frame = F, xlim = c(0, .3), xlab = "RNase sensitivity",
     ylab = "Frequency", main = "Distribution of normalized\nRNAse sensitivity")
lines(h0.2$mids, h0.2$counts / 65, type = "s", lwd = 2, col = "black")
# q-q plot
qqnorm(sort(ChauvenetFilter(mtx.norm.2[, 1], 2)), frame = F,
       pch = 21, bg = "gray80", col = "black")
qqline(ChauvenetFilter(mtx.norm.2[, 1], 2), lty = 1, col = cbp$green, 
       lwd = 1)
dev.off()
# shapiro.test
shapiro.test(sens)
# value thershold according to p-value = 0.001
v.th <- qnorm(1 - 0.001, mean = mean(sens), sd = sd(sens))


tmp.mtx <- mtx.norm[complete.cases(mtx.norm), ]
tmp.v <- (tmp.mtx[, 1] - tmp.mtx[, 2]) / sqrt(2)
mtx.norm.2 <- mtx.norm.2[rownames(mtx.norm.2) %in% names(ChauvenetFilter(tmp.v, 2)), ]

# plot scatter plot ----
orf12.pt <- mtx.norm.2[rownames(mtx.norm.2) %in% 
                        annotation[`Uniprot Symbol` %in% cell.pr[co_IP_with == "ORF1/2"]$`Uniprot Symbol`]$pg, ]
orf2.pt <- mtx.norm.2[rownames(mtx.norm.2) %in% 
                        annotation[`Uniprot Symbol` %in% cell.pr[co_IP_with == "ORF2"]$`Uniprot Symbol`]$pg, ]
orf1.pt <- mtx.norm.2[rownames(mtx.norm.2) %in% 
                        annotation[`Uniprot Symbol` %in% cell.pr[co_IP_with == "ORF1"]$`Uniprot Symbol`]$pg, ]
v.th.pt <- mtx.norm.2[
  apply(mtx.norm.2, 1, function(x) { sqrt(x[1] ^ 2 + x[2] ^ 2) > v.th }) & 
    rownames(mtx.norm.2) %nin% annotation[`Uniprot Symbol` %in% cell.pr$`Uniprot Symbol`]$pg, ]

pdf("out/RNAse_scatter_plot.pdf", width = 6, height = 6, useDingbats=FALSE)
par(mar = c(5,5,5,5))
plot(mtx.norm.2[rownames(mtx.norm.2) %nin% 
                  as.character(annotation[`Uniprot Symbol` %in% cell.pr$`Uniprot Symbol`]$pg), ], 
     frame = F, xlab = "RNAse Sensitivity (Exp. 1)",
     ylab = "RNAse Sensitivity (Exp. 2)", 
     pch = 20, col = cbp$grey, xlim = c(-.5, 1), ylim = c(-.5, 1))

points(orf12.pt, pch = 20, col = cbp$lightblue, cex = 1.5)
points(orf2.pt, pch = 20, col = cbp$pink, cex = 1.5)
points(orf1.pt, pch = 20, col = cbp$orange, cex = 1.5)
points(v.th.pt, pch = 20, col = "black", cex = 1.5)
text(orf12.pt, col = cbp$lightblue, pos = 4,
     labels = annotation[pg %in% rownames(orf12.pt) & use][order(pg)]$Gene)
text(orf2.pt[1:5, ], col = cbp$pink, pos = 2, 
     labels = annotation[pg %in% rownames(orf2.pt[1:5, ]) &
                           is_in_cell_article][order(pg)]$Gene)
text(orf2.pt[6:nrow(orf2.pt), ], col = cbp$pink, pos = 4, 
     labels = annotation[pg %in% rownames(orf2.pt[6:nrow(orf2.pt), ]) &
                           is_in_cell_article][order(pg)]$Gene)
text(v.th.pt, col = "black", pos = 2, 
   labels = annotation[pg %in% rownames(v.th.pt) & use][order(pg)]$Gene)
draw.circle(0, 0, radius = v.th, border = cbp$grey)

legend("top", legend = c(expression(paste("Significant in ", alpha, 
                                          "-ORF1p & ", alpha, "-ORF2p I-DIRT")),
                         expression(paste("Significant in ", alpha, 
                                          "-ORF2p I-DIRT only"))),
       pch = 20, col = c(cbp$lightblue, cbp$pink), pt.cex = 1.5, bty = "n")


dev.off()

# write output table
mtx.norm.2.dt <- as.data.table(mtx.norm.2)
mtx.norm.2.dt$pg <- rownames(mtx.norm.2)
mtx.norm.2.dt$Avg_affinity <- rowMeans(mtx.norm.2.dt[,.(RNAseSILAC_3, SILACRNAse_4)])
dt.norm <- ann.agg[pg %in% rownames(mtx.norm.2), .(pg, Genes, Uniprot, Protein)]
dt.norm$Genes <- gsub("\\*", "", dt.norm$Genes)
dt.norm$color <- sapply(dt.norm$pg, function(x){
  if(x %in% rownames(orf12.pt)) {
    return("blue")
  } else if(x %in% rownames(orf2.pt)){
    return("pink")
  } else if(x %in% rownames(v.th.pt)){
    return("black")
  } else {
    return("gray")
  }
})
dt.norm[Genes == "RPL5"]$color <- "gray"
dt.norm$pvalue <- sapply(dt.norm$pg, function(x){
   pnorm(sqrt(mtx.norm.2[x, 2] ^ 2 + mtx.norm.2[x, 1] ^ 2), mean = mean(sens), sd = sd(sens), lower.tail = F)
})
dt.norm[pvalue > 0.001 & color != "black"]$pvalue <- NA
dt.norm <- merge(dt.norm, mtx.norm.2.dt, by = "pg")
dt.norm$pg <- NULL
dt.norm <- rbind(dt.norm[color == "blue"][order(pvalue, decreasing = F)],
                 dt.norm[color == "black"][order(pvalue, decreasing = F)],
                 dt.norm[color == "pink"],
                 dt.norm[color == "gray"])
setnames(dt.norm, c("RNAseSILAC_3", "SILACRNAse_4", "pvalue", "color"), 
         c("Exp1_affinity", "Exp2_affinity", "P-value", "Color"))
dt.norm[, c("Exp1_affinity", "Exp2_affinity", "Avg_affinity")] <- round(dt.norm[,.(Exp1_affinity, Exp2_affinity, Avg_affinity)], 4)
write.table(dt.norm, "out/RNAse_supplementary_table.txt", sep = "\t", row.names = F,
            quote = F)
