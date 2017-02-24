##################
# This script makes sets the stage for analysis including loading data into memory.
# It also filters loops and saves the .rds file for differential
# loop calling, which is used in downstream analyses/figures/tables.
##################

library(diffloop)
library(ggplot2)
library(GenomicRanges)
library(reshape2)


options(scipen=999)
if (basename(getwd()) != "code") setwd("code")
source("theme_Publication.R")

# Process Full data
full <- loopsMake.mango("../data/raw_mango")
dim(full)
full <- updateLDGroups(full, c("K562", "K562", "MCF7", "MCF7"))
samples <- c("HCT116", "K562_r1", "K562_r2", "MCF7_r1", "MCF_r2")
full <- subsetLoops(full, full@rowData$loopWidth >= 5000) # remove and loops that merged together from import

# Remove regions of CNV
k562.cnv <- makeGRangesFromDataFrame(setNames(read.table("../data/cnv/K562-CNV.bedLogR")[,1:4],
                                              c("chr", "start", "end", "type")), keep.extra.columns = TRUE)
mcf7.cnv <- makeGRangesFromDataFrame(setNames(read.table("../data/cnv/MCF7-CNV.bedLogR")[,1:4],
                                              c("chr", "start", "end", "type")), keep.extra.columns = TRUE)

cnv.regions <- union(mcf7.cnv[mcols(mcf7.cnv)$type != "normal"], k562.cnv[mcols(k562.cnv)$type != "normal"])
noCNV <- removeRegion(full, rmchr(cnv.regions))
dim(noCNV)

# Filter Mango interactions
mangoSig <- mangoCorrection(full, FDR = 0.01)
dim(mangoSig)

# save RDS
realdat22 <- filterLoops(mangoSig, width = 5000, nreplicates = 2, nsamples = 2)
saveRDS(realdat22, "../output/qcPOL2comp.rds")
dim(realdat22)

# Identify differential loops 
km_res <- quickAssoc(realdat22)
sum(km_res@rowData$FDR < 0.01)

# Annotate Loops
h3k27ac.k <- rmchr(padGRanges(bedToGRanges("../data/annotation/K562_H3K27ac.bed"), pad = 1000))
h3k27ac.m <- rmchr(padGRanges(bedToGRanges("../data/annotation/MCF7_H3K27ac.bed"), pad = 1000))
enhancer <- union(h3k27ac.m, h3k27ac.k)
promoter <- padGRanges(getHumanTSS(), pad = 1000)

km_full_anno <- annotateLoops(km_res, enhancer = enhancer, promoter = promoter)
sum(km_full_anno@rowData$loop.type == "e-p" & km_full_anno@rowData$FDR < 0.01)

# Look at top hits to see what's interesting
# Then we actually made the figure at our Shiny implementation
# At DNAlandscapeR.aryeelab.org