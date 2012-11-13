library(GenomicRanges)
library(rtracklayer)

wdd <- '/Users/cwhelan/Documents/gene_rearrange/svpipeline/venter_chr2_100bp_dip'
testname <- 'newdiscfilt_as175_bt2_nomapq_nofilt_readGroupsBt2scoreMinL01_None_None_2500_25_3_sam_2000_fbbecea_175'

l1 <- read.table(paste(wd, '/', testname, '_l1.wig.gz', sep=""), skip=2, col.names=c("loc", "l1"))
w <- read.table(paste(wd, '/', testname, '_w0.wig.gz', sep=""), skip=2, col.names=c("loc", "w"))
mu <- read.table(paste(wd, '/', testname, '_mu1.wig.gz', sep=""), skip=2, col.names=c("loc", "mu"))
lrHet <- read.table(paste(wd, '/', testname, '_lrHet.wig.gz', sep=""), skip=2, col.names=c("loc", "lrHet"))
lrHom <- read.table(paste(wd, '/', testname, '_lrHom.wig.gz', sep=""), skip=2, col.names=c("loc", "lrHom"))

features <- merge(l1, w, mu, lrHet, lrHom, by="loc")

windows <- GRanges(seqnames="2", ranges=IRanges(start=seq(0,242751150,by=25), end=seq(24, 242751174, by=25)), strand="*")

hap1 <- import('~/Documents/gene_rearrange/svpipeline/venter_chr2_100bp_dip/HuRef.homozygous_indels_inversion.061109.chr2.Deletions.sim_hap1.bed.gz', asRangedData=FALSE)
hap2 <- import('~/Documents/gene_rearrange/svpipeline/venter_chr2_100bp_dip/HuRef.homozygous_indels_inversion.061109.chr2.Deletions.sim_hap2.bed.gz', asRangedData=FALSE)

elementMetadata(windows)[,"genotype"] <- 0

hap1gt50 <- hap1[end(hap1) - start(hap1) >= 50]
hap2gt50 <- hap2[end(hap2) - start(hap2) >= 50]

elementMetadata(windows)[as.matrix(findOverlaps(hap1gt50, windows))[,2],"genotype"] <- 1
elementMetadata(windows)[as.matrix(findOverlaps(hap2gt50, windows))[,2],"genotype"] <- 2

genotype <- data.frame(loc=start(windows), genotype=elementMetadata(windows)[,"genotype"])

features <- merge(features, genotype, by="loc")

rm(w)
rm(mu)
rm(l1)
rm(lrHet)
rm(lrHom)
rm(windows)

gc()