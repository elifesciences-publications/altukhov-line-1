# run this scrip after running src/RNAseq/analyse_quant.R

options(stringsAsFactors = F)
library(data.table)

# mapping statistics
read.stats <- data.table()
results.path <- "data/rnaseq/mapping/"
for(file.i in list.files(results.path, pattern = "Log.final.out", 
                         full.names = T, recursive = T)){
  read.stats.i <- readLines(file.i)
  read.stats.i.count <- strsplit(read.stats.i[6], "\t")[[1]][2]
  read.stats.i.sample <- gsub(paste0(results.path, "/(.*)/.*"), "\\1", 
                              file.i)
  read.stats.i.unmapped <- as.numeric(gsub("(.*)%", "\\1", strsplit(read.stats.i[30], "\t")[[1]][2]))
  read.stats.i.unmapped2 <- as.numeric(gsub("(.*)%", "\\1", strsplit(read.stats.i[31], "\t")[[1]][2]))
  read.stats.i.unmapped3 <- as.numeric(gsub("(.*)%", "\\1", strsplit(read.stats.i[29], "\t")[[1]][2]))
  read.stats.i.mapped <- 100 - read.stats.i.unmapped - read.stats.i.unmapped2 - read.stats.i.unmapped3
  read.stats <- rbind(read.stats, data.table(sample = read.stats.i.sample,
                                             count = read.stats.i.count,
                                             count_mapped = as.numeric(read.stats.i.count) * read.stats.i.mapped / 100,
                                             mapped = paste0(read.stats.i.mapped, "%")))
}

# read counts
read.counts <- fread("out/RNAseq_read_counts.txt")
read.stats <- merge(read.stats, stack(colSums(read.counts[,8:ncol(read.counts)])),
                    by.x = "sample", by.y = "ind")
setnames(read.stats, "values", "mapped_on_genes")
read.stats$percent_mapped_on_genes <- paste0(round(
  read.stats$mapped_on_genes / as.numeric(read.stats$count) * 100, 2), "%")

# counts by categories
pc.c <- stack(colSums(read.counts[gene_biotype == 'protein_coding', read.stats$sample, with = F]))
ld.c <- stack(colSums(read.counts[ensembl_gene_id == 'LD401', read.stats$sample, with = F]))
o.c <- stack(colSums(read.counts[ensembl_gene_id != 'LD401' & gene_biotype != 'protein_coding', 
                                 read.stats$sample, with = F]))
read.stats$count <- as.numeric(read.stats$count)
read.stats <- merge(read.stats, pc.c, by.x = "sample", by.y = "ind")
setnames(read.stats, "values", "protein_coding")
read.stats$percent_protein_coding <- read.stats$protein_coding / read.stats$mapped_on_genes * 100
read.stats <- merge(read.stats, ld.c, by.x = "sample", by.y = "ind")
setnames(read.stats, "values", "LD401")
read.stats$percent_LD401 <- read.stats$LD401 / read.stats$mapped_on_genes * 100
read.stats <- merge(read.stats, o.c, by.x = "sample", by.y = "ind")
setnames(read.stats, "values", "other")
read.stats$percent_other <- read.stats$other / read.stats$mapped_on_genes * 100

# getting size factors
size.factors <- fread("out/RNAseq_size_factors.txt")
read.stats <- merge(read.stats, size.factors, by.x = "sample", by.y = "ind")
read.stats$mapped_on_genes_norm <- round(read.stats$mapped_on_genes / read.stats$values, 0)
read.stats$LD401_norm <- round(read.stats$LD401 / read.stats$values, 0)
read.stats$protein_coding_norm <- round(read.stats$protein_coding / read.stats$values, 0)
read.stats$Other_norm <- round(read.stats$other / read.stats$values, 0)

write.table(read.stats, "out/RNAseq_mapping_statistics.txt", sep ="\t", 
            row.names = F)
