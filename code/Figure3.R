library(DESeq2)
library(diffloop)
library(ggplot2)
library(data.table)
library(GenomicRanges)
library(reshape2)
library(cowplot)

if (basename(getwd()) != "code") setwd("code")
source("theme_Publication.R")

#' # Epigenetics Plot
######## A
#+ cache = TRUE, message = FALSE
km_full <- readRDS("../output/qcPOL2comp.rds")
km_full@colData$groups <- c("k562", "k562", "mcf7", "mcf7")
km_res <- quickAssoc(km_full)
km_res <- annotateAnchors.bigwig(km_res, "../data/bigwigs/K562-DNase.bw")
km_res <- annotateAnchors.bigwig(km_res, "../data/bigwigs/MCF7-DNase.bw")

# Annotate anchors; link anchors together and take the max
log2FC.dnase <- log2(mcols(km_res@anchors)$MCF7.DNase/mcols(km_res@anchors)$K562.DNase) 
log2FC.dnase <- log2FC.dnase - mean(log2FC.dnase, na.rm = TRUE)

mcols(km_res@anchors) <- as.data.frame(cbind(mcols(km_res@anchors), data.frame(log2FC.dnase)))
big.km.summary <- summary(km_res)

# Keep with largest absolute value
idx <- as.logical(abs(big.km.summary$log2FC.dnase_1) < abs(big.km.summary$log2FC.dnase_2))
idx[is.na(idx)] <- FALSE
diff.dnase.max <- as.matrix((sapply(1:length(idx), function(i) {
    if (idx[i]) {
        big.km.summary[i, ]$log2FC.dnase_2
    } else {
        big.km.summary[i, ]$log2FC.dnase_1
    }
})))
big.km.summary$diff.dnase.max <- as.numeric(diff.dnase.max)
summary(lm(data = big.km.summary, logFC ~ diff.dnase.max))

# Bin data for boxplots
big.km.summary$diff.max.bin <- cut(big.km.summary$diff.dnase.max, seq(-12, 12, 3))

f31 <- ggplot(big.km.summary[complete.cases(big.km.summary), ]) + geom_hline(yintercept = 0) + theme_bw() + 
    ggtitle("Differential loops stratified by differential chromatin accessibility") + 
    geom_violin(aes(diff.max.bin, logFC), fill = "dodgerblue", trim = TRUE, scale = "width") +
    xlab("Max change in anchor log FC DNase binding") + ylab("log FC Loops")  + theme_Publication()

######## B

#' # DNA Methylation SubFigure
km_res <- annotateAnchors.bed(km_res, "../data/450k-methyl/K562-450k.bedgraph")
km_res <- annotateAnchors.bed(km_res, "../data/450k-methyl/MCF7-450k.bedgraph")

# Annotate anchors; link anchors together and take the max
mcf7methyl <- mcols(km_res@anchors)$MCF7.450k
k562methyl <- mcols(km_res@anchors)$K562.450k
change.methyl <- mcf7methyl - k562methyl
log2FC.methyl <- log2(mcf7methyl / k562methyl)

t.test(mcf7methyl, k562methyl, paired = TRUE)
mcols(km_res@anchors) <- as.data.frame(cbind(mcols(km_res@anchors), data.frame(cbind(log2FC.methyl, change.methyl))))
big.meth.summary <- summary(km_res)

# Keep with largest absolute value
idx <- as.logical(abs(big.meth.summary$change.methyl_1) < abs(big.meth.summary$change.methyl_2))
idx[is.na(idx) | is.nan(idx)] <- FALSE
diff.methyl.max <- as.matrix((sapply(1:length(idx), function(i) {
    if (idx[i]) {
        big.meth.summary[i, ]$change.methyl_2
    } else {
        big.meth.summary[i, ]$change.methyl_1
    }
})))
big.meth.summary$diff.methyl.max <- as.numeric(diff.methyl.max)
methyl.df <- big.meth.summary[, c("logFC", "diff.methyl.max", "FDR")]
methyl.df <- methyl.df[complete.cases(methyl.df) & !is.infinite(methyl.df$diff.methyl.max), ]

# Bin data for boxplots
methyl.df$diff.methyl.max.bin <- cut(methyl.df$diff.methyl.max, seq(-1, 1, 0.25))


f32 <- ggplot(methyl.df[complete.cases(methyl.df), ], aes(diff.methyl.max.bin, logFC)) +
    geom_violin(aes(diff.methyl.max.bin, logFC), fill = "firebrick", trim = TRUE, scale = "width") + 
    geom_hline(yintercept = 0) + theme_bw() + ggtitle("Differential loops stratified by differential methylation") + 
    xlab("Max change in anchor methylation") + ylab("log FC Loops") + theme_Publication()


####### C
km_full <- readRDS("../output/qcPOL2comp.rds")
km_full@colData$groups <- c("k562", "k562", "mcf7", "mcf7")
km_res <- quickAssoc(km_full)
km_res <- annotateAnchors.bigwig(km_res, "../data/bigwigs/K562-RAD21.bw")
km_res <- annotateAnchors.bigwig(km_res, "../data/bigwigs/MCF7-RAD21.bw")

# Annotate anchors; link anchors together and take the max
log2FCrad21 <- log2(mcols(km_res@anchors)$MCF7.RAD21/mcols(km_res@anchors)$K562.RAD21)
log2FCrad21 <- log2FCrad21 - mean(log2FCrad21, na.rm = TRUE)
mcols(km_res@anchors) <- as.data.frame(cbind(mcols(km_res@anchors), data.frame(log2FCrad21)))
big.km.summary <- summary(km_res)

# Keep with largest absolute value
idx <- as.logical(abs(big.km.summary$log2FCrad21_1) < abs(big.km.summary$log2FCrad21_2))
idx[is.na(idx)] <- FALSE
diff.rad21.max <- as.matrix((sapply(1:length(idx), function(i) {
    if (idx[i]) {
        big.km.summary[i, ]$log2FCrad21_2
    } else {
        big.km.summary[i, ]$log2FCrad21_1
    }
})))
big.km.summary$diff.rad21.max <- as.numeric(diff.rad21.max)
summary(lm(data = big.km.summary, logFC ~ diff.rad21.max))

# Bin data for boxplots
big.km.summary$diff.max.bin <- cut(big.km.summary$diff.rad21.max, seq(-6, 6, 2))

f33 <- ggplot(big.km.summary[complete.cases(big.km.summary), ], aes(diff.max.bin, logFC)) +
    geom_violin(aes(diff.max.bin, logFC), fill = "dodgerblue", trim = TRUE, scale = "width") + 
    geom_hline(yintercept = 0) + theme_bw() + ggtitle("Differential loops stratified by differential cohesin localization") + 
    xlab("Max change in anchor log FC RAD21 localization") + ylab("log FC Loops") + theme_Publication()


##### D
km_full <- readRDS("../output/qcPOL2comp.rds")
km_full@colData$groups <- c("k562", "k562", "mcf7", "mcf7")
km_res <- quickAssoc(km_full)
km_res <- annotateAnchors.bigwig(km_res, "../data/bigwigs/K562-H3K27ac.bw")
km_res <- annotateAnchors.bigwig(km_res, "../data/bigwigs/MCF7-H3K27ac.bw")

# Annotate anchors; link anchors together and take the max
log2FC.H3K27ac <- log2(mcols(km_res@anchors)$MCF7.H3K27ac/mcols(km_res@anchors)$K562.H3K27ac)
log2FC.H3K27ac <- log2FC.H3K27ac - mean(log2FC.H3K27ac)
mcols(km_res@anchors) <- as.data.frame(cbind(mcols(km_res@anchors), data.frame(log2FC.H3K27ac)))
big.km.summary <- summary(km_res)

# Keep with largest absolute value
idx <- as.logical(abs(big.km.summary$log2FC.H3K27ac_1) < abs(big.km.summary$log2FC.H3K27ac_2))
idx[is.na(idx)] <- FALSE
diff.H3K27ac.max <- as.matrix((sapply(1:length(idx), function(i) {
    if (idx[i]) {
        big.km.summary[i, ]$log2FC.H3K27ac_2
    } else {
        big.km.summary[i, ]$log2FC.H3K27ac_1
    }
})))
big.km.summary$diff.H3K27ac.max <- as.numeric(diff.H3K27ac.max)
summary(lm(data = big.km.summary, logFC ~ diff.H3K27ac.max))

# Bin data for boxplots
big.km.summary$diff.max.bin <- cut(big.km.summary$diff.H3K27ac.max, seq(-6, 8, 2))

f34 <- ggplot(big.km.summary[complete.cases(big.km.summary), ], aes(diff.max.bin, logFC)) +
    geom_violin(aes(diff.max.bin, logFC), fill = "dodgerblue", trim = TRUE, scale = "width") + 
    geom_hline(yintercept = 0) + theme_bw() + ggtitle("Differential loops stratified by differential H3K27ac") + 
    xlab("Max change in anchor log FC H3K27ac") + ylab("log FC Loops") + theme_Publication()

png("../figures/Figure3.png", width = 4800, height = 4800, res = 300)
plot_grid(f31, f32, f33, f34, ncol = 2, labels = c('A', 'B', 'C', 'D'), label_size = 20)
dev.off()
