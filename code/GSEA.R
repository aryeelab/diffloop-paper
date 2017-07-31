#source("https://bioconductor.org/biocLite.R")
#biocLite("piano")

library(diffloop)
library(piano)

# Table 0
# Make table of differential loops with annotation for text

km_full <- readRDS("../output/qcPOL2comp.rds")
h3k27ac.k <- rmchr(padGRanges(bedToGRanges("../data/annotation/K562_H3K27ac.bed"), pad = 1000))
h3k27ac.m <- rmchr(padGRanges(bedToGRanges("../data/annotation/MCF7_H3K27ac.bed"), pad = 1000))
enhancer <- union(h3k27ac.m, h3k27ac.k)
promoter <- padGRanges(getHumanTSS(), pad = 1000)
km_full@colData$groups <- c("k562", "k562", "mcf7", "mcf7")
km_res <- quickAssoc(km_full)
km_assoc <- annotateLoops(km_res, enhancer = enhancer, promoter = promoter)
km_assoc <- loopGenes(km_assoc)
df <- summary(km_assoc)
k <- df[df$logFC < 0 & df$FDR < 0.01, ]
m <- df[df$logFC > 0 & df$FDR < 0.01, ]

kr <- cbind(sum(k$loop.type == "e-p"), sum(k$loop.type == "p-p"), sum(k$loop.type == "e-e"), sum(k$loop.type == "none"), length(k$loop.type))
mr <- cbind(sum(m$loop.type == "e-p"), sum(m$loop.type == "p-p"), sum(m$loop.type == "e-e"), sum(m$loop.type == "none"), length(m$loop.type))
tr <- kr + mr

t1out <- data.frame(rbind(kr, mr, tr))
rownames(t1out) <- c("K562", "MCF7", "Total")
colnames(t1out) <- c("E-P", "P-P", "E-E", "None", "Total")

write.table(t1out, file = "../tables/Table0.txt", sep = "\t", quote = FALSE)

# Keep EP loops
km_res.ep <- keepEPloops(km_res, enhancer, promoter)
df2 <- summary(km_res.ep)

# Table 1
# Strongest Associations
sdf <- rbind(head(df2[order(df2$logFC),] , 5), head(df2[order(df2$logFC, decreasing = TRUE),] , 5))
write.table(sdf[,c(1:11,14,17,18,20)], file = "../tables/Table1.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Table 2 - Gene set analysis
gsc <- loadGSC("../data/msigdb/h.all.v5.2.symbols.gmt")
s <- strsplit(df2[,20], split = ",")
km <- data.frame(gene = unlist(s), pval = rep(df2[,17], sapply(s, length)), logfc = rep(df2[,14], sapply(s, length)))
# Rank by fold change-signed p-value
km$stat <- (1-km$pval) * sign(km$logfc)

# K562
k <- km[order(km$stat, km$pval), ]
k <- k[!duplicated(k$gene), ]
stat <- -k$stat
names(stat) <- k$gene
head(k)
head(stat)
system.time(gsa <- runGSA(stat, geneSetStat="wilcoxon", signifMethod="nullDist", gsc=gsc, gsSizeLim=c(5,300)))
gsa <- GSAsummaryTable(gsa)
gsa <- gsa[order(gsa[,"p (dist.dir.up)"]),]
k_gsa <- gsa[gsa[,"p adj (dist.dir.up)"] < 0.1, ]
k_gsa[, c(1,5)]
write.table(k_gsa[, c(1,5)], file = "../tables/Table_k562.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# MCF-7
m <- km[rev(order(km$stat, km$pval)), ]
m <- m[!duplicated(m$gene), ]
stat <- m$stat
names(stat) <- m$gene
head(m)
head(stat)
system.time(gsa <- runGSA(stat, geneSetStat="wilcoxon", signifMethod="nullDist", gsc=gsc, gsSizeLim=c(5,300)))
gsa <- GSAsummaryTable(gsa)
o <- order(gsa[,"p (dist.dir.up)"])
gsa <- gsa[o,]
m_gsa <- gsa[gsa[,"p adj (dist.dir.up)"] < 0.1, ]
m_gsa[,c(1,5)]
write.table(m_gsa[, c(1,5)], file = "../tables/Table_mcf7.txt", sep = "\t", quote = FALSE, row.names = FALSE)
intersect(m$gene[1:20], gsc$gsc[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]])
