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
# read exp template
dt.tmpl <- fread("data/tbls/experimentalDesignTemplate.txt", header = T, sep = "\t")
# names of ORF1 and ORF2 proteins
orf1.orf2 <- c(orf1 = "sp|Q9UN81|LORF1_HUMAN", 
               orf2 = "Orf2-untagged(optimized):")
# out genes rnaseq
out.genes <- "out/rnaseq_genes.txt"
