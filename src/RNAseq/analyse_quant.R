library(magrittr)
library(data.table)
library(DESeq2)
library(vsn)
library(ggplot2)
library(hexbin)
library(biomaRt)
library(gplots)
library(GenomicAlignments)
options(stringsAsFactors = F)

source("lib/utils.R")
source("lib/constants.R")

# LD401 gtf
ld401 <- fread("data/rnaseq/LD401.gtf")

# bam files
list.bam.files <- list.files("data/rnaseq/mapping/", pattern = "*.bam$", 
                             recursive = T, full.names = T)

# read counts
results.path <- "data/rnaseq/mapping/"
read.counts <- data.table(V1 = character())
for(file.i in list.files(results.path, pattern = 'ReadsPerGene', full.names = T, 
                       recursive = T)){
  read.counts.i <- fread(file.i, 
                       header = F, skip = 4)
  read.counts.i <- read.counts.i[, .(V1, V2)]
  setnames(read.counts.i, "V2", gsub(paste0(results.path, "/(.*)/.*"), "\\1", 
                                     file.i))
  read.counts <- merge(read.counts, read.counts.i, by = "V1", all = T, fill = 0)
}
read.counts <- read.counts[rowSums(read.counts[,2:ncol(read.counts)] > 10) >= 3]

# getting gene names
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
genes <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name', 
                            'start_position', 'end_position', 'strand', 'gene_biotype'), 
               filters = 'ensembl_gene_id', values = read.counts$V1, mart = ensembl)
read.counts.g <- as.data.table(merge(genes, read.counts, by.x = "ensembl_gene_id", by.y = "V1", all.y = T))
read.counts.g <- read.counts.g[read.counts.g$ensembl_gene_id %in% readLines(out.genes), ]
write.table(read.counts.g, "out/RNAseq_read_counts.txt", sep = "\t", row.names = F)

read.counts.dt <- as.data.frame(read.counts.g[,.(ensembl_gene_id, `2_401_Inp_10_11_S1`, 
                                               `2_567_Elu_10_16_S8`, `3_401_Sup_10_11_S2`, 
                                               `4_401_Elu_10_11_S3`, `4_624_Elu_10_16_S9`, 
                                               `6_401_Inp_10_12_S4`, `6_567_Elu_10_10_S7`, 
                                               `7_401_Sup_10_12_S5`, `8_401_Elu_10_12_S6`)])

read.counts.dt <- unique(read.counts.dt)
rownames(read.counts.dt) <- read.counts.dt$ensembl_gene_id
read.counts.dt$ensembl_gene_id <- NULL
# corr matrix
pdf("out/RNAseq_correlation_matrix.pdf", width = 8, height = 8)
heatmap.2(cor(read.counts.dt[,c("2_401_Inp_10_11_S1","2_567_Elu_10_16_S8",
                             "3_401_Sup_10_11_S2","4_401_Elu_10_11_S3",
                             "4_624_Elu_10_16_S9","6_401_Inp_10_12_S4",
                             "6_567_Elu_10_10_S7","7_401_Sup_10_12_S5",
                             "8_401_Elu_10_12_S6")], method = 'spearman'),
          cellnote = round(cor(read.counts.dt, method = 'spearman'), 2),
          main = "Spearman correlation",
          notecol="black",
          density.info="none",
          trace="none", cexRow = .9, cexCol = .9,
          margins =c(12,12), key = F,
          dendrogram = "col",
          col= colorRampPalette(c("grey", "#ff9999",  "#ff3333", "#ff0000", "#cc0000"))(n = 299))  
dev.off()

# sample information
sample.info <- read.table("data/rnaseq/sample_info", sep = "\t", header = T)
rownames(sample.info) <- sample.info$sample_name
sample.info$sample_name <- NULL

# generate the DESeqDataSet
DESeq.ds <- DESeqDataSetFromMatrix(countData = read.counts.dt[,c("2_401_Inp_10_11_S1","2_567_Elu_10_16_S8",
                                                                 "3_401_Sup_10_11_S2","4_401_Elu_10_11_S3",
                                                                 "4_624_Elu_10_16_S9","6_401_Inp_10_12_S4",
                                                                 "6_567_Elu_10_10_S7","7_401_Sup_10_12_S5",
                                                                 "8_401_Elu_10_12_S6")],
                                   colData = sample.info[c("2_401_Inp_10_11_S1","2_567_Elu_10_16_S8",
                                                           "3_401_Sup_10_11_S2","4_401_Elu_10_11_S3",
                                                           "4_624_Elu_10_16_S9","6_401_Inp_10_12_S4",
                                                           "6_567_Elu_10_10_S7","7_401_Sup_10_12_S5",
                                                           "8_401_Elu_10_12_S6"), ,drop = F],
                                   design = ~condition)

# remove genes without any counts
DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds)) > 0, ]

# calculate the size factor and add it to the data set
DESeq.ds <- estimateSizeFactors(DESeq.ds)
size.factors <- stack(sizeFactors(DESeq.ds))
write.table(size.factors, "out/RNAseq_size_factors.txt",
            sep = "\t", row.names = F)

#  normalized read counts
counts.sf_normalized <- t(t(counts(DESeq.ds)) / size.factors$values)

# transform size-factor normalized read counts to log2 scale using a pseudocount of 1
log.norm.counts <- log2(counts.sf_normalized + 1)
norm.counts.dt <- as.data.table(log.norm.counts)
norm.counts.dt$ensembl_gene_id <- rownames(log.norm.counts)
norm.counts.dt <- merge(genes, norm.counts.dt, by = "ensembl_gene_id", all.y = T)
write.table(norm.counts.dt, "out/log2_normalized_read_counts.txt",
            sep = "\t", row.names = F)

pdf("out/Boxplots_rnaseq_normalization.pdf", width = 6, height = 5)
layout(matrix(1:2, nrow = 1))
par(mar = c(14, 4, 4, 2), cex = .8)
boxplot(log2(assay(DESeq.ds) + 1), lwd = 2, frame = F, las = 2,
        main = "Log2-transformed\nread counts", ylab = "log2(read counts)")

# box plots of log2-transformed read counts
boxplot(log.norm.counts, lwd = 2, frame = F, las = 2,
        main = "Log2-transformed\nnormalized read counts", 
        ylab = "log2(normalized read counts)")
dev.off()

# presence of gene types
read.counts.m <- melt(read.counts.g[, c(7:ncol(read.counts.g)), with = F], id.vars = 'gene_biotype')
read.counts.m[is.na(gene_biotype)]$gene_biotype <- 
  paste0('Synthetic L1 pLD', gsub("[0-9]{1-2}_(.*?)_.*", "\\1", read.counts.m[is.na(gene_biotype)]$variable))
read.counts.m[gene_biotype != 'protein_coding' & 
                !grepl('Synthetic', gene_biotype)]$gene_biotype <- "Other"
read.counts.m.agg <- read.counts.m[, list(sum_reads = sum(value)), by = .(gene_biotype, variable)]
read.counts.m.agg <- rbindlist(lapply(unique(read.counts.m.agg$variable), function(v) {
  tmp <- read.counts.m.agg[variable == v]
  tmp$perc <- round(tmp$sum_r / sum(tmp$sum_r) * 100, 4)
  tmp
}))
types <- unique(read.counts.m.agg$gene_biotype)
names(types) <- types
types[grep("Synthetic", types)] <- "Synthetic L1"
colors.types <- unlist(cbp[as.numeric(as.factor(types)) + 3])
names(colors.types) <- names(types)
types.u <- unique(names(types))
types.u <- types.u[c(4,5,7,1,3,6,2)]
types.u <- factor(types.u, levels = types.u)
names(types.u) <- types.u
pdf("out/RNAseq_pie_charts.pdf", width = 11, height = 11)
par(mar = c(4,10,4,10))
layout(matrix(1:4, nrow = 2))
for(col.i in unique(read.counts.m.agg$variable)) {
  pie.i <- read.counts.m.agg[variable == col.i]
  pie.i <- pie.i[match(as.character(sort(types.u[pie.i$gene_biotype])), pie.i$gene_biotype), ]
  pie(pie.i$perc, labels = paste0(pie.i$gene_biotype, " ", round(pie.i$perc, 1), "%"),
      main = col.i, col = add.alpha(colors.types[pie.i$gene_biotype], .7), clockwise = T, init.angle =180 * 3)
}
dev.off()


# Replicate average 
abundance.avg.norm <- read.counts.g
for(sam in rownames(sample.info)) {
  abundance.avg.norm[, sam] <-
    abundance.avg.norm[, sam, with = F] / 
           size.factors[ind == sam]$values
}
for(gr in unique(sample.info$condition)) {
  abundance.avg.norm[, gr] <- rowMeans(
    abundance.avg.norm[, rownames(sample.info[sample.info$condition == gr,,drop = F]), 
                       with = F])
}
abundance.avg.norm <- abundance.avg.norm[, c(colnames(abundance.avg.norm)[1:7],
                                             unique(sam.to.avg$condition)),
                                         with = F]
write.table(abundance.avg.norm, "out/RNASeq_avg_normalized_read_counts.txt",
            sep = "\t", row.names = F)

# LD401 Coverage
gene_coords <- unlist(ld401[grep("gene_id \"LD401\"", V9), .(V4, V5)])
total.range <- GRanges("LD401", IRanges(min(gene_coords), max(gene_coords)))
range.cov.to.avg <- 
  lapply(list.bam.files[sapply(rownames(sample.info), grep, list.bam.files)], 
         function(bam.file.name){
           ga <- readGAlignments(bam.file.name, 
                                 index = paste0(bam.file.name, ".bai"),
                                 param = ScanBamParam(which = total.range))
           ga.cov <- coverage(ga)
           as.numeric(ga.cov[total.range][[1]])
         })
names(range.cov.to.avg) <- rownames(sample.info)
range.cov.to.avg.norm <- lapply(names(range.cov.to.avg), function(x){
  range.cov.to.avg[[x]] / size.factors[ind == x]$values
})
names(range.cov.to.avg.norm) <- names(range.cov.to.avg)
range.cov.avg <- lapply(unique(sample.info$condition), function(cond) {
  cond.sams <- rownames(sample.info[sample.info$condition == cond,, drop = F])
  rowMeans(matrix(
    unlist(range.cov.to.avg.norm[cond.sams]),
    ncol = length(range.cov.to.avg.norm[cond.sams]), 
    byrow = F))
})

names(range.cov.avg) <- unique(sample.info$condition)

# Plot coverage
gr1 <- c("401_Elu", "401_Inp", "401_Sup")
gr2 <- c("567_Elu", "624_Elu", "401_Inp")
list.coords <- list(ORF1 = c(35, 1080),
                    ORF2 = c(1131, 4968),
                    Tags = c(4969, 5197),
                    "3'UTR" = c(5198, 5482))
pdf("out/RNAseq_coverage_plot.pdf", width = 9, height = 6)
for(gr in list(gr1, gr2)) {
  y.min.pos <- -max(unlist(range.cov.avg))/10
  y.max.pos <- max(unlist(range.cov.avg))
  
  plot(NA, 
       xlim = c(0, 6000),
       ylim = c(y.min.pos, y.max.pos),
       xlab = 'pLD401 Coordinate', ylab = 'Normalized Coverage', frame = F,
       main = "")
  colors <- c(cbp$green, cbp$blue, cbp$pink)
  names(colors) <- gr
  for(gr.i in gr) {
    lines(x = min(total.range@ranges[1]):max(total.range@ranges[1]),
          y = range.cov.avg[[gr.i]], col = colors[gr.i], lwd = 2)
  }
  legend("topleft", legend = gr, bg = add.alpha("white", .5), box.lty = 0,
         fill = colors, cex = .8)
  abline(h = 0, lwd = 2, lty = 2, col = "lightgrey")
  for(i in 1:length(list.coords)) {
    polygon(x = c(list.coords[[i]][1], list.coords[[i]][1], 
                  list.coords[[i]][2], list.coords[[i]][2]),
            y = c(y.min.pos/2, y.min.pos/6, y.min.pos/6, y.min.pos/2), 
            col = "lightgrey")
    text(x = mean(list.coords[[i]]), y = y.min.pos/2, 
         labels = names(list.coords)[i], pos = 1, cex = .8)
  }
}
dev.off()

# plot avg pie charts
abundance.gr.m <- melt(abundance.avg.norm[, c(7:ncol(abundance.avg.norm)), with = F], 
                       id.vars = 'gene_biotype')
abundance.gr.m[is.na(gene_biotype)]$gene_biotype <- 
  paste0('Synthetic L1 pLD', gsub("([0-9]{3})_.*", "\\1", 
                                  abundance.gr.m[is.na(gene_biotype)]$variable))
abundance.gr.m[gene_biotype != "protein_coding" & 
                 !grepl('Synthetic', gene_biotype)]$gene_biotype <- 'Other'
abundance.gr.m.agg <- abundance.gr.m[, list(sum_reads = sum(value)), by = .(gene_biotype, variable)]
abundance.gr.m.agg <- rbindlist(lapply(unique(abundance.gr.m.agg$variable), function(v) {
  tmp <- abundance.gr.m.agg[variable == v]
  tmp$perc <- round(tmp$sum_r / sum(tmp$sum_r) * 100, 4)
  tmp
}))
types <- unique(abundance.gr.m.agg$gene_biotype)
names(types) <- types
types[grep("Synthetic", types)] <- "Synthetic L1"
colors.types <- unlist(cbp[as.numeric(as.factor(types)) + 2])
names(colors.types) <- names(types)
types.u <- unique(names(types))
types.u <- types.u[c(4,5,7,1,3,6,2)]
types.u <- factor(types.u, levels = types.u)
names(types.u) <- types.u
pdf("out/RNAseq_avg_pie_charts.pdf", width = 11, height = 11)
par(mar = c(4,10,4,10))
layout(matrix(1:4, nrow = 2))
for(col.i in unique(abundance.gr.m.agg$variable)) {
  pie.i <- abundance.gr.m.agg[variable == col.i]
  pie.i <- pie.i[match(as.character(sort(types.u[pie.i$gene_biotype])), pie.i$gene_biotype), ]
  pie(pie.i$perc, labels = paste0(pie.i$gene_biotype, " ", round(pie.i$perc, 1), "%"),
      main = col.i, col = add.alpha(colors.types[pie.i$gene_biotype], .7), clockwise = T, 
      init.angle = 180 * 3)
}
dev.off()

