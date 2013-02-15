library(GenomicRanges)
library(rtracklayer)
library(plyr)
library(RColorBrewer)
library(VariantAnnotation)
library(xtable)

drawPerfLine <- function(perf, col, lineType, maxFP) {
  calls <- c(perf$Calls, min(perf$Calls) - min(perf$TP), 0)
  tp <- c(perf$TP, 0, 0)
  fp <- calls - tp
  lines(fp[fp <= maxFP], tp[fp <= maxFP], col=col, lty=lineType, lwd=3)
}

plotROC <- function(perfs, perfNames, totalDels, main, sim=TRUE, maxTP=NA, legendLoc="bottomright") {
  maxFP <- totalDels
  if (is.na(maxTP)) {
    maxTP <- totalDels
  }
  plot(0, type="n", ylim=c(0, maxTP), xlim=c(0,maxFP), xlab=ifelse(sim, "False Positives", "Novel Predictions"), ylab="True Positives", main=main)
  if (length(perfs) <= 9) {
    perfCols <- brewer.pal(length(perfs), "Set1")
  } else {
    perfCols <- rainbow(length(perfs))
  }
  lineTypes <- seq(1,length(perfs))
  mapply(drawPerfLine, perfs, perfCols, lineTypes, MoreArgs=list(maxFP=maxFP))  
  legend(legendLoc, legend=perfNames, col=perfCols, lty=lineTypes, lwd=3, cex=.9)
}

overlap <- function(a,b,x,y) {
  pmin(b,y) - pmax(a,x)
}

sizeclasses <- function(starts, ends) {
  sizelens <- ends - starts + 1
  sizeclasses <- rep(NA, length(sizelens))
  sizeclasses[sizelens <= 100] <- "40-100bp"
  sizeclasses[sizelens > 100 & sizelens <= 250] <- "101-250bp"
  sizeclasses[sizelens > 250 & sizelens <= 500] <- "251-500bp"
  sizeclasses[sizelens > 500 & sizelens <= 1000] <- "501-1000bp"
  sizeclasses[sizelens > 1000] <- "> 1000bp"  
  factor(sizeclasses, levels=c("40-100bp","101-250bp","251-500bp","501-1000bp","> 1000bp"),ordered=TRUE)
}

processPredictions <- function(name, predfile, trueDels) {
  predswithhits <- read.table(predfile)
  names(predswithhits) <- c('predchrom', 'predstart', 'predend', 'truechrom', 'truestart', 'trueend')
  
  overlaps <- with(predswithhits, overlap(predstart,predend,truestart,trueend))
  
  predswithhits$overlaps <- overlaps
  
  maxhits <- ddply(predswithhits,.(predchrom,predstart,predend),function(x){x[which.max(x$overlaps),]})
      
  maxhits$tok <- paste(maxhits$truechrom, ':', maxhits$truestart, '-', maxhits$trueend, sep="")
  
  maxhits <- merge(maxhits,mcols(trueDels), by="tok")
  maxhits$name <- name
  maxhits
}

#chr2 100bp diploid deletions
totalDels <- 400
cloudbreak <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/cloudbreak_venterchr2allindels100bpdip_readgroup1_mt0.8_del_perf.txt', header=TRUE, sep="\t")
breakdancer <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort_35_2_3_del.perf.txt', header=TRUE)
gasv <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/gasvpro.perf.txt', header=TRUE)
delly <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort.delly_q0_c5_del.perf.txt', header=TRUE)
pindel <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort_pindel_D.perf.txt', header=TRUE)

perfsList <- list(cloudbreak=cloudbreak, breakdancer=breakdancer, pindel=pindel, gasv=gasv, delly=delly, )
pdf('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/CHR2SIM_DELS_ROC.pdf', width=10)
par(xpd=T, mar=par()$mar+c(0,0,0,7))
plotROC(perfsList, 
        c("Cloudbreak", "Breakdancer","Pindel", "GASVPro", "DELLY"), 
        totalDels, "Deletions in Venter diploid chr2 simulation",legendLoc=xy.coords(425,200), maxTP=350)
dev.off()

#chr2 100bp diploid insertions
totalInsertions <- 80
cloudbreak <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/cloudbreak_venterchr2allindels100bpdip_readgroup1_ins_perf.txt', header=TRUE, sep="\t")
breakdancer <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/breakdancer_insertions.perf.txt', header=TRUE)
pindel <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/pindel_insertions.perf.txt', header=TRUE)

perfsList <- list(cloudbreak=cloudbreak, breakdancer=breakdancer, pindel=pindel)
pdf('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/CHR2SIM_INS_ROC.pdf', width=10)
par(xpd=T, mar=par()$mar+c(0,0,0,7))
plotROC(perfsList, 
        c("Cloudbreak", "Breakdancer","Pindel"), 
        totalInsertions, "Insertions in Venter diploid chr2 simulation",legendLoc=xy.coords(85,50), maxTP=80)
dev.off()

trueDels <- import.bed('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/HuRef.homozygous_indels.061109.chr2.deletions.gff.bed', asRangedData=FALSE)
trueDelshap2 <- import.bed('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/HuRef.homozygous_indels.061109.chr2.sim_hap2.deletions.gff.bed', asRangedData=FALSE)

mcols(trueDels)$tok <- paste(seqnames(trueDels), ':', start(trueDels)-1, '-', end(trueDels), sep="")
mcols(trueDels)$haps <- 1

hap2toks <- paste(seqnames(trueDelshap2), ':', start(trueDelshap2)-1, '-', end(trueDelshap2), sep="")
mcols(trueDels[mcols(trueDels)$tok %in% hap2toks])$haps <- 2

mcols(trueDels)$sizeclasses <- sizeclasses(start(trueDels), end(trueDels))

segdups <- import.bed('~/genomes/1kg/hg18/segmental_duplications.bed', asRangedData=FALSE)

mcols(trueDels)$segdup <- 0
mcols(trueDels[as.matrix(findOverlaps(segdups,trueDels))[,2]])$segdup <- 1

repMask <- import.bed('~/genomes/1kg/hg18/repeatmasker_b36_merged.bed.gz', asRangedData=FALSE)

mcols(trueDels)$repMask <- 0
mcols(trueDels[as.matrix(findOverlaps(repMask,trueDels))[,2]])$repMask <- 1

cb <- processPredictions("Cloudbreak", '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/cloudbreak_venterchr2allindels100bpdip_readgroup1_deletions_210.hits.txt', trueDels)

bd <- processPredictions('Breakdancer', '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort_35_2_3.bd_out.dels.pairgte4.hits.txt', trueDels)

gasv <- processPredictions('GASVPro', '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/BAMToGASV.gasvpro.in.ALL.MCMCThreshold.clusters.pruned.clusters.gte12.hits.txt', trueDels)

delly <- processPredictions('DELLY', '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort.delly_q0_c5_del.gte7.hits.txt', trueDels)

pindel <- processPredictions('Pindel', '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort_pindel_D_gte33.hits.txt', trueDels)

allpreds <- rbind(as.data.frame(cb),as.data.frame(bd),as.data.frame(gasv),as.data.frame(delly),as.data.frame(pindel))

allpreds <- ddply(allpreds, .(tok), function(x) { cbind(x, list(discoveries=rep(dim(x)[1],dim(x)[1]))) })

trueDelsGte40 <- trueDels[end(trueDels) - start(trueDels) >= 40]
trueDelSizes <- table(trueDelsGte40$sizeclasses)
trueDelHaps <- table(trueDelsGte40$haps)
trueDelRepmask <- table(trueDelsGte40$repMask)

predsBySize <- table(allpreds[,c('name','sizeclasses')])
exBySize <- table(allpreds[which(allpreds$discoveries == 1),c('name','sizeclasses')])

predsByRepMask <- table(allpreds[,c('name','repMask')])
exPredsByRepMask <- table(allpreds[which(allpreds$discoveries == 1),c('name','repMask')])

predsByHap <- table(allpreds[,c('name','haps')])
exPredsByHap <- table(allpreds[which(allpreds$discoveries == 1),c('name','haps')])

predsBySegDup <- table(allpreds[,c('name','segdup')])
exPredsBySegDup <- table(allpreds[which(allpreds$discoveries == 1),c('name','segdup')])

print(xtable(trueDelSizes))
print(xtable(trueDelHaps))
print(xtable(trueDelRepmask))

print(xtable(predsBySize))
print(xtable(exBySize))

print(xtable(predsByRepMask))
print(xtable(exPredsByRepMask))

# genotyping with extra data
extraData <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/cloudbreak_venterchr2allindels100bpdip_readgroup1_deletions_210.bed', skip=1)
names(extraData) <- c('predchrom', 'predstart', 'predend', 'peaknum', 'maxscore', 'svtype', 'avgMu', 'minMu', 'maxMu', 'avgW', 'minW', 'maxW')

cbpredsWithExtraData <- merge(as.data.frame(cb), extraData)

hapCM <- table(cbpredsWithExtraData$haps, cbpredsWithExtraData$avgW < .2)
print(hapCM)
print((hapCM[1,1] + hapCM[2,2]) / sum(hapCM))

# NA18507

trueDelsNA18507 <- import.bed('~/Documents/gene_rearrange/svpipeline/NA18507/MillsDIP1KG_merged.bed', asRangedData=FALSE)

mcols(trueDelsNA18507)$tok <- paste(seqnames(trueDelsNA18507), ':', start(trueDelsNA18507)-1, '-', end(trueDelsNA18507), sep="")

mcols(trueDelsNA18507)$sizeclasses <- sizeclasses(start(trueDelsNA18507), end(trueDelsNA18507))

repMaskHg19 <- import.bed('~/genomes/1kg/hg19/repeatmasker.bed', asRangedData=FALSE)

mcols(trueDelsNA18507)$repMask <- 0
mcols(trueDelsNA18507[as.matrix(findOverlaps(repMaskHg19,trueDelsNA18507))[,2]])$repMask <- 1

trueDelsNA18507Gte40 <- trueDelsNA18507[end(trueDelsNA18507) - start(trueDelsNA18507) >= 40]
trueDelsNA18507Sizes <- table(trueDelsNA18507Gte40$sizeclasses)

trueDelRepmask <- table(trueDelsNA18507Gte40$repMask)

na185071KGVariants <- readVcf('~/Documents/gene_rearrange/svpipeline/NA18507/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.NA18507.vcf.gz', "hg19")

mcols(rowData(na185071KGVariants))$gt <- as.factor(geno(na185071KGVariants)$GT)
mcols(rowData(na185071KGVariants))$hap <- ifelse(mcols(rowData(na185071KGVariants))$gt == "1|1", 2, 1)

na185071KGDels <- rowData(na185071KGVariants)[nchar(fixed(na185071KGVariants)$REF) > nchar(unlist(fixed(na185071KGVariants)$ALT))]

na185071KGDels <- na185071KGDels[width(na185071KGDels) >= 40]

mcols(trueDelsNA18507)$hap <- NA

ol <- findOverlaps(trueDelsNA18507, na185071KGDels) 

mcols(trueDelsNA18507[as.matrix(ol)[,1]])$hap <- mcols(na185071KGDels[as.matrix(ol)[,2]])$hap

# ROC curve
totalDels <- 10000
cloudbreak <- read.table('~/Documents/gene_rearrange/svpipeline/NA18507/new_gold_na18507final_readGroupsRazerS3_i94_s99_m1000_None_None_25000_25_sam_6_67eedd7_-1_-1_3_lrHeterozygous.allsvs.short40.perf.txt', header=TRUE, sep="\t")
breakdancer <- read.table('~/Documents/gene_rearrange/svpipeline/NA18507/NA18507_35_2_3_singlecpu_perf.txt', header=TRUE)
gasv <- read.table('~/Documents/gene_rearrange/svpipeline/NA18507/gasvpro.perf.txt', header=TRUE)
delly <- read.table('~/Documents/gene_rearrange/svpipeline/NA18507/NA18507.delly_q0_c5_del.perf.txt', header=TRUE)
pindel <- read.table('~/Documents/gene_rearrange/svpipeline/NA18507/bwa_pindel_D.perf.txt', header=TRUE)

perfsList <- list(cloudbreak=cloudbreak, breakdancer=breakdancer, pindel=pindel, gasv=gasv, delly=delly)
pdf('~/Documents/svpipeline/manuscript/NA18507_ROC.pdf', width=10)
par(xpd=T, mar=par()$mar+c(0,0,0,7))
plotROC(perfsList, c("Cloudbreak", "Breakdancer","Pindel","GASVPro", "DELLY", "Cloudbreak"), totalDels, "NA18507", legendLoc=xy.coords(10750,700), maxTP=1300, sim=FALSE)
dev.off()

# break down predictons

cbHitsNA18507 <- processPredictions("Cloudbreak", '~/Documents/gene_rearrange/svpipeline/NA18507/na18507final_readGroupsRazerS3_i94_s99_m1000_None_None_25000_25_sam_6_67eedd7_-1_-1_3_deletions_260.hits.txt', trueDelsNA18507)

bdHitsNA18507 <- processPredictions('Breakdancer', '~/Documents/gene_rearrange/svpipeline/NA18507/NA18507_35_2_3.bd_out.dels.pairgte4pt9.hits.txt', trueDelsNA18507)

gasvHitsNA18507 <- processPredictions('GASVPro', '~/Documents/gene_rearrange/svpipeline/NA18507/BAMToGASV.gasvpro.in.ALL.MCMCThreshold.clusters.pruned.clusters.gte14pt8.hits.txt', trueDelsNA18507)

dellyHitsNA18507 <- processPredictions('DELLY', '~/Documents/gene_rearrange/svpipeline/NA18507/NA18507.delly_q0_c5_del.gte8pt6.hits.txt', trueDelsNA18507)

pindelHitsNA18507 <- processPredictions('Pindel', '~/Documents/gene_rearrange/svpipeline/NA18507/bwa_pindel.noshort.gte40pt7.hits.txt', trueDelsNA18507)

allpredsNA18507 <- rbind(as.data.frame(cbHitsNA18507),as.data.frame(bdHitsNA18507),as.data.frame(gasvHitsNA18507),as.data.frame(dellyHitsNA18507),as.data.frame(pindelHitsNA18507))
allpredsNA18507 <- ddply(allpredsNA18507, .(tok), function(x) { cbind(x, list(discoveries=rep(dim(x)[1],dim(x)[1]))) })

predsBySizeNA18507 <- table(allpredsNA18507[,c('name','sizeclasses')])
exBySizeNA18507 <- table(allpredsNA18507[which(allpredsNA18507$discoveries == 1),c('name','sizeclasses')])

predsByRepMaskNA18507 <- table(allpredsNA18507[,c('name','repMask')])
exByRepMaskNA18507 <- table(allpredsNA18507[which(allpredsNA18507$discoveries == 1),c('name','repMask')])

# genotyping 

extraDataNA18507 <- read.table('~/Documents/gene_rearrange/svpipeline/NA18507/na18507final_readGroupsRazerS3_i94_s99_m1000_None_None_25000_25_sam_6_67eedd7_-1_-1_3_deletions_260.bed', skip=1)
names(extraDataNA18507) <- c('predchrom', 'predstart', 'predend', 'peaknum', 'maxscore', 'svType', 'avgMu', 'minMu', 'maxMu', 'avgW', 'minW', 'maxW')

cbpredsWithExtraDataNA18507 <- merge(as.data.frame(cbHitsNA18507), extraDataNA18507)

cbpredsWithExtraDataNA18507Genotyped <- cbpredsWithExtraDataNA18507[!is.na(cbpredsWithExtraDataNA18507$hap),]

hapCMNA18507 <- table(cbpredsWithExtraDataNA18507$avgW < .2, cbpredsWithExtraDataNA18507$hap)
print(hapCMNA18507)
print((hapCMNA18507[1,1] + hapCMNA18507[2,2]) / sum(hapCMNA18507))

millsDelsWithGenotypes <- read.table('~/Documents/Papers2/Articles/2011/Mills/Supplemental/dels_gt40_with_genotypes.csv', header=TRUE, sep="\t")
millsGenotypeData <- read.table('~/Documents/Papers2/Articles/2011/Mills/Supplemental/Mills_genotypes_NA18507_with_ids_gt40.txt', header=TRUE, sep="\t")

millsFullGenotypedDels <- merge(millsDelsWithGenotypes,millsGenotypeData)

millsGenotypedRanges <- GRanges(seqnames=substr(millsFullGenotypedDels$CHR,4,length(millsFullGenotypedDels$CHR)), ranges=IRanges(start=millsFullGenotypedDels$START, end=millsFullGenotypedDels$STOP), mcols=DataFrame(hap=millsFullGenotypedDels$SED00002_NA18507.CEL))

# gem param grid search
#chr2 100bp diploid deletions
totalDels <- 400
l <- list()
for (e in c(4,5,6,7,8,9,10)) {
  for (q in c(1,3,5,7,9,1000)) {
    name <- paste0("gem_e", e, "_q", q)
    fname <- paste0('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/', name, "_del_perf.txt")
    tryCatch(l <- c(l, list(read.table(fname, header=TRUE, sep="\t"))),
             error = function(e) print(paste(name, e)))
    names(l)[length(l)] <- name
  }
}

l <- c(l, list(read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/cloudbreak_venterchr2allindels100bpdip_readgroup1_del_perf.txt', header=TRUE, sep="\t")))

plotROC(head(l, 30), head(names(l), 30), 
        totalDels, "Deletions in Venter diploid chr2 simulation")

lapply(l, function(x) {max(x[x$TPR > .9, 'TP'])})