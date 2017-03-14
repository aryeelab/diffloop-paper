library(DESeq2)
library(diffloop)
library(ggplot2)
library(data.table)
library(GenomicRanges)
library(reshape2)
library(cowplot)
source("theme_Publication.R")

if (basename(getwd()) != "code") setwd("code")

rdsRes <- "../output/DESeq-Res.rds"
if (file.exists(rdsRes)) {
    res <- readRDS(rdsRes)
} else {
    files <- grep("counts", list.files("../data/rna-seq/"), value = TRUE)
    files <- paste("../data/rna-seq/", files, sep = "")
    condition <- c("k562", "k562", "k562", "mcf7", "mcf7", "mcf7")
    names <- c("k1", "k2", "k3", "m1", "m2", "m3")
    sampleTable <- data.frame(sampleName = names, fileName = files, condition = condition)
    dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = getwd(), design = ~condition)
    
    dds <- dds[rowSums(counts(dds)) > 1, ]  #remove zero counts
    dds$condition <- factor(dds$condition, levels = c("k562", "mcf7"))
    dds <- DESeq(dds)
    res <- results(dds)
    res <- res[complete.cases(res), ]
    res <- res[order(res$pvalue), ]
    saveRDS(res, file = rdsRes)
}

plotMA(res)

#' ## Load ChIA-PET Data/Perform diffloop
km_full <- readRDS("../output/qcPOL2comp.rds")
h3k27ac.k <- rmchr(padGRanges(bedToGRanges("../data/annotation/K562_H3K27ac.bed"), pad = 1000))
h3k27ac.m <- rmchr(padGRanges(bedToGRanges("../data/annotation/MCF7_H3K27ac.bed"), pad = 1000))
enhancer <- union(h3k27ac.m, h3k27ac.k)
promoter <- padGRanges(getHumanTSS(), pad = 1000)
km_full@colData$groups <- c("k562", "k562", "mcf7", "mcf7")

km_full_anno <- annotateLoops(km_full, enhancer = enhancer, promoter = promoter)
sum(km_full_anno@rowData$loop.type == "e-p")


km_res <- quickAssoc(km_full)
km_res.ep <- keepEPloops(km_res, enhancer, promoter)
km_res.ep <- annotateAnchors.bigwig(km_res.ep, "../data/bigwigs/K562-H3K27ac.bw")
km_res.ep <- annotateAnchors.bigwig(km_res.ep, "../data/bigwigs/MCF7-H3K27ac.bw")
km_res.ep <- annotateAnchors.bed(km_res.ep, "../data/450k-methyl/K562-450k.bedgraph")
km_res.ep <- annotateAnchors.bed(km_res.ep, "../data/450k-methyl/MCF7-450k.bedgraph")

mcols(km_res.ep@anchors) <- cbind(log2(mcols(km_res.ep@anchors)$MCF7.H3K27ac/mcols(km_res.ep@anchors)$K562.H3K27ac), 
    mcols(km_res.ep@anchors)$MCF7.450k - mcols(km_res.ep@anchors)$K562.450k)


#' ## Link the datasets
#+ cache = TRUE
km.linked <- annotateLoops.dge(km_res.ep, res, multiple = FALSE)
dim(km.linked)
length(unique(km.linked@rowData$gene.tss))

df <- summary(km.linked)
df$logFC_bin <- cut(df$logFC, seq(-9, 9, 3))
df$logFC_bin2 <- cut(df$logFC, seq(-9, 9, 3))

df$log2FoldChange_bin <- cut(df$log2FoldChange, seq(-20, 20, 5))

distal_enhancer <- rep(0, dim(df)[1])
distal_methyl <- rep(0, dim(df)[1])
distal_enhancer[df$anchor.tss == 1] <- df$V1_2[df$anchor.tss == 1]
distal_enhancer[df$anchor.tss == 2] <- df$V1_1[df$anchor.tss == 2]
df$distal_enhancer <- distal_enhancer

distal_methyl[df$anchor.tss == 1] <- df$V2_2[df$anchor.tss == 1]
distal_methyl[df$anchor.tss == 2] <- df$V2_1[df$anchor.tss == 2]
df$distal_methyl <- distal_methyl

proximal_enhancer <- rep(0, dim(df)[1])
proximal_enhancer[df$anchor.tss == 1] <- df$V1_1[df$anchor.tss == 1]
proximal_enhancer[df$anchor.tss == 2] <- df$V1_2[df$anchor.tss == 2]
df$proximal_enhancer <- proximal_enhancer

proximal_methyl <- rep(0, dim(df)[1])
proximal_methyl[df$anchor.tss == 1] <- df$V2_1[df$anchor.tss == 1]
proximal_methyl[df$anchor.tss == 2] <- df$V2_2[df$anchor.tss == 2]
df$proximal_methyl <- proximal_methyl

df$Differential_Loop <- df$FDR < 0.01

df.ssig <- df[df$padj < 0.01 & df$FDR < 0.01, ]



f4A <- ggplot(df[complete.cases(df) & df$padj < 0.01, ], aes(logFC_bin, log2FoldChange)) +
    geom_violin(aes(logFC_bin, log2FoldChange), fill = "dodgerblue", trim = TRUE, scale = "width") + 
    geom_hline(yintercept = 0) +  theme_bw() +  theme_Publication() +
    labs(title = "Differential expression stratified by E-P Loops", x = "log FC Loops", 
    y = "log FC Gene Expression")

f4B <- qplot(df$distal_enhancer, df$log2FoldChange) + theme_bw() + labs(title = "Differential Expression by enhancer H3K27ac", 
    x = "log FC H3K27ac at enhancer anchor", y = "log FC of Transcript Expression") + geom_point(aes(colour = df$Differential_Loop)) + theme_Publication() +
     scale_colour_manual(name="",  values =c("dimgrey", "red"))+ theme(legend.position="none") + geom_smooth(method = "lm", se = FALSE, color = "black")

f4C <- qplot(df$proximal_enhancer, df$log2FoldChange) + theme_bw() + labs(title = "Differential Expression by promoter H3K27ac", 
    x = "log FC H3K27ac at promoter anchor", y = "log FC of Transcript Expression") + geom_point(aes(colour = df$Differential_Loop)) + theme_Publication() +
    scale_colour_manual(name="",  values =c("dimgray", "red"))+ theme(legend.position="none") + geom_smooth(method = "lm", se = FALSE, color = "black")

f4D <- qplot(df$distal_methyl, df$log2FoldChange) + theme_bw() + labs(title = "Differential Expression by enhancer DNA Methylation", 
    x = "Change in DNA methylation at enhancer anchor", y = "log FC of Transcript Expression") + geom_point(aes(colour = df$Differential_Loop)) +
    theme_Publication() + scale_colour_manual(name="",  values =c("dimgray", "red"))+ theme(legend.position="none") + geom_smooth(method = "lm", se = FALSE, color = "black")

f4E <- qplot(df$proximal_methyl, df$log2FoldChange) + theme_bw() + labs(title = "Differential Expression by promoter DNA Methylation", 
    x = "Change in DNA methylation at promoter anchor", y = "log FC of Transcript Expression") + geom_point(aes(colour = df$Differential_Loop)) +
    theme_Publication() + scale_colour_manual(name="",  values =c("dimgray", "red")) + theme(legend.position="none")+ geom_smooth(method = "lm", se = FALSE, color = "black")


png("../figures/Figure4.png", width = 4800, height = 7200, res = 300)
plot_grid(f4A, plot_grid(f4B, f4C, f4D, f4E, ncol = 2, labels = c('B','C','D','E'), label_size = 20),
          ncol = 1, rel_heights = c(0.5, 1), label_size = 20, labels = c('A', NULL))
dev.off()


############################# 


#' ## Top biology
# Highly expressed w/ loop in breast cancer
df.ssig[df.ssig$logFC < -5 & df.ssig$log2FoldChange < -5, ]

# Highly expresssed w/ loop in leukemia
df.ssig[df.ssig$logFC > 9 & df.ssig$log2FoldChange > 9, ]

#' ## Weird outliers
# Highly expressed in Breast cancer regardless of loops
a <- dim(df.ssig[df.ssig$logFC > 0 & df.ssig$log2FoldChange < 0, ])[1]
# Highly expressed in leukemia cancer regardless of loops
b <- dim(df.ssig[df.ssig$logFC < 0 & df.ssig$log2FoldChange > 0, ])[1]

1 - (a + b) / dim(df.ssig)[1]

dim(df.ssig)