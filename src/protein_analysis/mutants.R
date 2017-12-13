# ORF2p catalytic point mutations alter the proteomes purified L1 RNPs
library(data.table)
library(seqinr)

options(scipen=999)

# source "package"
source("lib/exp-functions.R")
source("lib/affinity-obj.R")
source("lib/beauty-pair-plots.R")
source("lib/utils.R")
source("lib/constants.R")


# names of the experiments
exps <- list(
  RT_mut = c("401_624_Fusion", "MAP_3_624", "MAP_4_624"),
  ENDO_mut = c("401_567_Fusion", "MAP_1_567", "MAP_2_567")
)

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
l.exp <- get_list_mqexp(l.read$data, dt.tmpl.e$Experiment)
names(l.exp) <- dt.tmpl.e$Experiment

# get affinity object for each experiment
#  remove contaminants and decoy proteins
l.aff <- lapply(dt.tmpl.e$Experiment, function(i)
  affinity(l.exp[[i]], minpep = 2,
           avalue = dt.tmpl.e[Experiment == i]$avalue,
           score.filt = T, scale = dt.tmpl.e[Experiment == i]$scale,
           scale.sig =  5e-2, scale.sig.aff = 0.7,
           use.normalized = dt.tmpl.e[Experiment == i]$use.normalized)
)
names(l.aff) <- dt.tmpl.e$Aname

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
annotation <- annotation[`Uniprot Symbol` %nin% cont]
annotation$pg <- as.character(annotation$pg)
ann.agg <- annotation[, list(Genes = paste(Gene[!is.na(Gene)], collapse = ", "),
                             Uniprot = paste(`Uniprot Symbol`[!is.na(`Uniprot Symbol`)], collapse = ", "),
                             Protein = paste(Protein[!is.na(Protein)], collapse = ", "),
                             is_in_cell_article = c(F, T)[max(is_in_cell_article) + 1]), 
                      by = list(pg)]

# create empty matrix
aff.mtx <- matrix(NA, nrow = length(v.pg),
                  ncol = length(l.aff),
                  dimnames = list(v.pg, names(l.aff)))
# fill with values
for (i in 1:length(l.aff)) {
  aff.mtx[as.character(l.aff[[i]]$aff$pg), names(l.aff)[i]] <-
    l.aff[[i]]$aff$aff
}
aff.mtx <- aff.mtx[rownames(aff.mtx) %in% annotation$pg, ]

# matrix with affinity
mtx <- aff.mtx[, unlist(exps)]

# normalization
mtx.norm <- sapply(unlist(exps), function(exp) {
  exp.v <- mtx[, exp]
  exp.v.orfs <- exp.v[as.character(prs[protein %in% orf1.orf2]$pg)]
  orf2.delta <- 0.5 - exp.v.orfs[as.character(prs[protein == orf1.orf2[2]]$pg)]
  exp.v.norm <- exp.v + orf2.delta
  return(exp.v.norm)
})
mtx.norm <- mtx.norm[apply(mtx.norm, 1, function(x){
  if(sum(is.na(x[grep("RT", names(x))])) <= 1 & 
     sum(is.na(x[grep("ENDO", names(x))])) <= 1) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}), ]


# creating data.table
dt.norm <- as.data.table(mtx.norm)
dt.norm$pg <- rownames(mtx.norm)
dt.norm.cell <- dt.norm[pg %in% as.character(prs[`Uniprot Symbol` %in% 
                                      c(cell.pr$`Uniprot Symbol`, orf1.orf2)]$pg)]
dt.norm.cell <- merge(dt.norm.cell, ann.agg, by = "pg")

# calculating mean value for replicates
dt.norm.cell$RT_mut_mean <- 1 - rowMeans(
  dt.norm.cell[, grep("RT", colnames(dt.norm.cell), value = T), with = F],
  na.rm = T
  )
dt.norm.cell$ENDO_mut_mean <- 1 - rowMeans(
  dt.norm.cell[, grep("ENDO", colnames(dt.norm.cell), value = T), with = F],
  na.rm = T
)
dt.norm.cell$delta <- dt.norm.cell$RT_mut_mean - dt.norm.cell$ENDO_mut_mean

# calculate se
dt.norm.cell$RT_sd <- apply(dt.norm.cell, 1, function(x){
  rt.mut <- as.numeric(x[grep("RT_mut[0-9]", names(x))])
  sd(rt.mut, na.rm = T)
})
dt.norm.cell$ENDO_sd <- apply(dt.norm.cell, 1, function(x){
  endo.mut <- as.numeric(x[grep("ENDO_mut[0-9]", names(x))])
  sd(endo.mut, na.rm = T)
})
dt.norm.cell <- dt.norm.cell[order(delta, decreasing = T)]

# calculating significance
dt.norm.cell$p.value <- apply(dt.norm.cell, 1, function(x){
  rt.mut <- as.numeric(x[grep("RT_mut[0-9]", names(x))])
  endo.mut <- as.numeric(x[grep("ENDO_mut[0-9]", names(x))])
  if(ulength(rt.mut) == 1 | ulength(endo.mut) == 1){
    return(1)
  } else {
    t.test(rt.mut, endo.mut)$p.value
  }
})
dt.norm.cell$p.value.adj <- p.adjust(dt.norm.cell$p.value, method = "BH")

# calculate significance 1-sample t-test
dt.norm <- merge(ann.agg[,.(Genes, Uniprot, Protein, pg, is_in_cell_article)], 
                 dt.norm, by = "pg")
dt.norm.cell$RT_Avg_affinity <- apply(1 - mtx.norm[dt.norm.cell$pg, ], 1, function(x){
  mean(x[grep("RT", names(x))], na.rm = T)
})
dt.norm.cell$ENDO_Avg_affinity <- apply(1 - mtx.norm[dt.norm.cell$pg, ], 1, function(x){
  mean(x[grep("ENDO", names(x))], na.rm = T)
})
dt.norm$RT_Avg_affinity <- apply(1 - mtx.norm[dt.norm$pg, ], 1, function(x){
  mean(x[grep("RT", names(x))], na.rm = T)
})
dt.norm$ENDO_Avg_affinity <- apply(1 - mtx.norm[dt.norm$pg, ], 1, function(x){
  mean(x[grep("ENDO", names(x))], na.rm = T)
})

dt.norm.cell[, RT_pvalue :=
  round(apply(1 - mtx.norm[pg,, drop = F], 1, function(x){
    z <- x[grep("RT", names(x))]
    if(length(unique(z)) == 1) {
      return(NaN)
    } else {
      return(t.test(z, mu = 0.5)$p.value)
    }
}), 5)]
dt.norm.cell[, ENDO_pvalue :=
          round(apply(1 - mtx.norm[pg, , drop = F], 1, function(x){
            z <- x[grep("ENDO", names(x))]
            if(length(unique(z)) == 1) {
              return(NaN)
            } else {
              return(t.test(z, mu = 0.5)$p.value)
            }
          }), 5)]

dt.norm.cell[, RT_pvalue.adj := round(p.adjust(RT_pvalue, method = "BH"), 6)]
dt.norm.cell[, ENDO_pvalue.adj := round(p.adjust(ENDO_pvalue, method = "BH"), 6)]
dt.norm.cell[is.na(RT_pvalue.adj)]$RT_pvalue.adj <- 1
dt.norm.cell[is.na(ENDO_pvalue.adj)]$ENDO_pvalue.adj <- 1

# plot barplots
pdf("out/Mutants_barplot.pdf", width = 8, height = 5, useDingbats = F)
par(mar = c(6, 5, 4, 1))
bp <- barplot(matrix(c(dt.norm.cell$ENDO_mut_mean, dt.norm.cell$RT_mut_mean), 
                 nrow = 2, byrow = T), beside = T, cex.names = 0.7,
        names.arg = annotation[pg %in% dt.norm.cell$pg & 
                                 (is_in_cell_article == T | Gene == "ORF2")][
                                   match(dt.norm.cell$pg, pg)]$Gene, 
        las = 2, ylim = c(0,1), axes = F,
        ylab = "", col = c(cbp$pink, cbp$lightblue))
abline(h = c(2/3, 1/3, 1/3 + 1/6), lwd = c(2,2,2), col = "gray", lty = c(1,1,2))
yaxis.ticks <- c(1/9, 0.2, 1/3, 1/2, 2/3, 4/5, 8/9)
axis(2, yaxis.ticks, las = 2,
     labels = as.character(round(yaxis.ticks / (1 - yaxis.ticks), 2)))
par(xpd = T)
text(x = apply(bp, 2, mean)[dt.norm.cell$p.value.adj < 0.05 & 
                              dt.norm.cell$p.value.adj >= 0.01], 
     y = -0.05, labels = "*")
text(x = apply(bp, 2, mean)[dt.norm.cell$p.value.adj < 0.01], 
     y = -0.05, labels = "**")
legend("topleft", legend = c("EN vs WT", "RT vs WT"), 
       fill = c(cbp$pink, cbp$lightblue), bty = "n")

segments(c(bp[1, ], bp[2, ]), 
         c(dt.norm.cell$ENDO_mut_mean, dt.norm.cell$RT_mut_mean) - 
           c(dt.norm.cell$ENDO_sd, dt.norm.cell$RT_sd), 
         c(bp[1, ], bp[2, ]), 
         c(dt.norm.cell$ENDO_mut_mean, dt.norm.cell$RT_mut_mean) + 
           c(dt.norm.cell$ENDO_sd, dt.norm.cell$RT_sd),
         lwd = 1)

arrows(c(bp[1, ], bp[2, ]), 
       c(dt.norm.cell$ENDO_mut_mean, dt.norm.cell$RT_mut_mean) - 
         c(dt.norm.cell$ENDO_sd, dt.norm.cell$RT_sd), 
       c(bp[1, ], bp[2, ]), 
       c(dt.norm.cell$ENDO_mut_mean, dt.norm.cell$RT_mut_mean) + 
         c(dt.norm.cell$ENDO_sd, dt.norm.cell$RT_sd),
       lwd = 1, angle = 90,
       code = 3, length = 0.01)

points(x = bp[1, dt.norm.cell$ENDO_pvalue.adj > .05],
       y = (dt.norm.cell$ENDO_mut_mean + dt.norm.cell$ENDO_sd + .03)[
         dt.norm.cell$ENDO_pvalue.adj > .05],
       pch = 24, cex = .5)
points(x = bp[2, dt.norm.cell$RT_pvalue.adj > .05],
       y = (dt.norm.cell$RT_mut_mean + dt.norm.cell$RT_sd + .03)[
         dt.norm.cell$RT_pvalue.adj > .05],
       pch = 24, cex = .5)

dev.off()

# plot distributions
pdf("out/Mutants_affinities_histogram.pdf", width = 8, height = 5)
layout(matrix(1:ncol(mtx.norm), nrow = 2, byrow = T))
for(i in 1:ncol(mtx.norm)){
  hist(mtx.norm[, i], breaks = 20, col = "grey95", xlim = c(0, 1),
       main = colnames(mtx.norm)[i], xlab = "Normalized affinity")
  abline(v = mtx.norm[as.character(prs[protein %in% orf1.orf2]$pg), i], col = c('red', "blue"), 
         lwd = 2)
  legend("topright", bty = "n", legend = c("ORF1", "ORF2"), lwd = 2, 
         col = c("red", "blue"))
}
dev.off()


# highlight proteins from the cell article
cell.pr.pg <- dt.norm[pg %in% annotation[is_in_cell_article == T]$pg 
                       & RT_Avg_affinity <= ENDO_Avg_affinity
                      ]$pg
highlite.col = c(cbp$green, cbp$pink)[
  cell.pr.pg %in% dt.norm.cell[p.value.adj < .05]$pg + 1]
pdf("out/Mutants_paired_histogram_r2l.pdf", width = 6, height = 5)
plotTwoExpDT(dt.norm, x1 = "ENDO_Avg_affinity", x2 = "RT_Avg_affinity",
             x1.name = "ENDO", x2.name = "RT",
             highlite = cell.pr.pg,
             highlite.col = highlite.col,
             highlite.names = annotation[is_in_cell_article == T][match(cell.pr.pg, pg)]$Gene, 
             highlite.legend = data.table(legend = c("Not significantly different",
                                                     "Significantly different"),
                                          col = c(cbp$green, cbp$pink)))
dev.off()
cell.pr.pg <- dt.norm[pg %in% annotation[is_in_cell_article == T]$pg 
                      & RT_Avg_affinity >= ENDO_Avg_affinity
                      ]$pg
highlite.col = c(cbp$green, cbp$pink)[
  cell.pr.pg %in% dt.norm.cell[p.value.adj < .05]$pg + 1]
pdf("out/Mutants_paired_histogram_l2r.pdf", width = 6, height = 5)
plotTwoExpDT(dt.norm, x1 = "ENDO_Avg_affinity", x2 = "RT_Avg_affinity",
             x1.name = "ENDO", x2.name = "RT",
             highlite = cell.pr.pg,
             highlite.col = highlite.col,
             highlite.names = annotation[is_in_cell_article == T][match(cell.pr.pg, pg)]$Gene, 
             highlite.legend = data.table(legend = c("Not significantly different",
                                                     "Significantly different"),
                                          col = c(cbp$green, cbp$pink)))
dev.off()

cell.pr.pg <- dt.norm[pg %in% annotation[is_in_cell_article == T]$pg]$pg
highlite.col = c(cbp$green, cbp$pink)[
  cell.pr.pg %in% dt.norm.cell[p.value.adj < .05]$pg + 1]
pdf("out/Mutants_paired_histogram.pdf", width = 6, height = 5)
plotTwoExpDT(dt.norm, x1 = "ENDO_Avg_affinity", x2 = "RT_Avg_affinity",
             x1.name = "ENDO", x2.name = "RT",
             highlite = cell.pr.pg,
             highlite.col = highlite.col,
             highlite.legend = data.table(legend = c("Not significantly different",
                                                     "Significantly different"),
                                          col = c(cbp$green, cbp$pink)))
dev.off()

# write output
dt.norm$RT_Avg_affinity <- NULL
dt.norm$ENDO_Avg_affinity <- NULL
dt.norm <- merge(dt.norm, 
                 dt.norm.cell[,.(pg, p.value.adj, RT_Avg_affinity, ENDO_Avg_affinity, RT_pvalue.adj, ENDO_pvalue.adj)], 
                 by = "pg", all.x = T)
dt.norm$Color <- "grey"
dt.norm[p.value.adj < .05]$Color <- "pink"
dt.norm[p.value.adj >= .05 & is_in_cell_article == T]$Color <- "green"
dt.norm[is_in_cell_article == F]$p.value.adj <- NA
dt.norm <- dt.norm[order(p.value.adj, decreasing = T)]
dt.norm <- dt.norm[,.(Genes, Uniprot, Protein, Color, RT_mut1, RT_mut2, RT_mut3,
                      ENDO_mut1, ENDO_mut2, ENDO_mut3,
                      RT_Avg_affinity, ENDO_Avg_affinity, p.value.adj, 
                      RT_pvalue.adj, ENDO_pvalue.adj)]
setnames(dt.norm, c("p.value.adj", "RT_mut1", "RT_mut2", "RT_mut3",
                    "ENDO_mut1", "ENDO_mut2", "ENDO_mut3", 
                    "RT_pvalue.adj", "ENDO_pvalue.adj"), 
         c("Adjusted_P-value", "RT_Exp1_affinity", "RT_Exp2_affinity", "RT_Exp3_affinity",
           "EN_Exp1_affinity", "EN_Exp2_affinity", "EN_Exp3_affinity",
           "RT_vs_WT_Adjusted_P-value", "EN_vs_WT_Adjusted_P-value"))
dt.norm$Genes <- gsub("\\*", "", dt.norm$Genes)
dt.norm[,c("RT_Exp1_affinity", "RT_Exp2_affinity", "RT_Exp3_affinity",
           "EN_Exp1_affinity", "EN_Exp2_affinity", "EN_Exp3_affinity",
           "RT_Avg_affinity", "ENDO_Avg_affinity")] <- round(dt.norm[,.(RT_Exp1_affinity, RT_Exp2_affinity, RT_Exp3_affinity,
                                                                        EN_Exp1_affinity, EN_Exp2_affinity, EN_Exp3_affinity, 
                                                                        RT_Avg_affinity, ENDO_Avg_affinity)], 6)
dt.norm$Genes <- gsub("NA, ", "", dt.norm$Genes)

write.table(dt.norm, "out/Mutants_supplementary_table.txt", 
            row.names = F, sep = "\t", quote = F)
