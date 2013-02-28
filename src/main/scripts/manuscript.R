library(GenomicRanges)
library(rtracklayer)
library(plyr)
library(RColorBrewer)
library(VariantAnnotation)
library(brew)

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

sizeclasses <- function(sizelens) {
  sc <- rep("< 40bp", length(sizelens))
  sc[sizelens >= 40 & sizelens <= 100] <- "40-100bp"
  sc[sizelens > 100 & sizelens <= 250] <- "101-250bp"
  sc[sizelens > 250 & sizelens <= 500] <- "251-500bp"
  sc[sizelens > 500 & sizelens <= 1000] <- "501-1000bp"
  sc[sizelens > 1000] <- "> 1000bp"  
  factor(sc, levels=c("< 40bp", "40-100bp","101-250bp","251-500bp","501-1000bp","> 1000bp"),ordered=TRUE)
}

processDeletionPredictions <- function(name, predfile, trueDels) {
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

processInsertionPredictions <- function(name, predfile, trueInsertions) {
  predswithhits <- read.table(predfile)
  names(predswithhits) <- c('predchrom', 'predstart', 'predend', 'truechrom', 'truestart', 'trueend','trueLength')
  
  overlaps <- with(predswithhits, overlap(predstart,predend,truestart,trueend+trueLength))
  
  predswithhits$overlaps <- overlaps
  
  maxhits <- ddply(predswithhits,.(predchrom,predstart,predend),function(x){x[which.max(x$overlaps),]})
  
  maxhits$tok <- paste(maxhits$truechrom, ':', maxhits$truestart, '-', maxhits$trueend, sep="")
  maxhits <- merge(maxhits,mcols(trueInsertions), by="tok")
  maxhits$name <- name
  maxhits
}

printPredTableValue <- function(t, x, y) {
  d <- t[x,y]
  if (t[x,y] == max(t[,y])) {
    return(paste0("\\textbf{",formatC(d,digits=3,format="fg"),"}"))
  } else {
    return(formatC(d,digits=3,format="fg"))
  }
}

printPredTableValueFromList <- function(x, l) {
  d <- l[x]
  if (d == max(l)) {
    return(paste0("\\textbf{",formatC(d,digits=3,format="fg"),"}"))
  } else {
    return(formatC(d,digits=3,format="fg"))
  }
}

printPredTableCell <- function(preds, ex, x, y) {
  return(paste0(printPredTableValue(preds, x, y), " (", printPredTableValue(ex, x, y), ")"))
}

# input files

# chr2

# deletions
cloudbreakChr2DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/cloudbreak_venterchr2allindels100bpdip_readgroup1_mt0.8_del_perf.txt'
breakdancerChr2DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort_35_2_3_del.perf.txt'
gasvProChr2DeletionPerfFile <-  '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/gasvpro.perf.txt'
dellyChr2DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort.delly_q0_c5_del.perf.txt'
pindelChr2DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort_pindel_D.perf.txt'

cloudbreakChr2DeletionPredictionsFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/cloudbreak_venterchr2allindels100bpdip_readgroup1_deletions_210.bed'

cloudbreakChr2DeletionHitsFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/cloudbreak_venterchr2allindels100bpdip_readgroup1_deletions_210.hits.txt'
breakdancerChr2DeletionHitsFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort_35_2_3.bd_out.dels.pairgte4.hits.txt'
gasvProChr2DeletionHitsFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/BAMToGASV.gasvpro.in.ALL.MCMCThreshold.clusters.pruned.clusters.gte12.hits.txt'
dellyChr2DeletionHitsFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort.delly_q0_c5_del.gte7.hits.txt'
pindelChr2DeletionHitsFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort_pindel_D_gte33.hits.txt'

# insertions
cloudbreakChr2InsertionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/cloudbreak_venterchr2allindels100bpdip_readgroup1_ins_mfw5_gt40_nosplit.perf.txt'
breakdancerChr2InsertionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/breakdancer_insertions.perf.txt'
pindelChr2InsertionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/pindel_insertions.perf.txt'

cloudbreakChr2InsertionPredictionsFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/cloudbreak_venterchr2allindels100bpdip_readgroup1_insertions_mfw5_gt40_nosplit_gt0p5.bed'

cloudbreakChr2InsertionHitsFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/cloudbreak_venterchr2allindels100bpdip_readgroup1_ins_mfw5_gt40_nosplit.gte0p5.hits.txt'
breakdancerChr2InsertionHitsFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort_35_2_3.bd_out.ins.gte2.hits.txt'
pindelChr2InsertionHitsFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/pindel_insertions.hits.txt'

# NA18507

# Deletions
cloudbreakNA18507DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/new_gold_na18507final_readGroupsRazerS3_i94_s99_m1000_None_None_25000_25_sam_6_67eedd7_-1_-1_3_lrHeterozygous.allsvs.short40.perf.txt'
breakdancerNA18507DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/NA18507_35_2_3_singlecpu_perf.txt'
gasvProNA18507DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/gasvpro.perf.txt'
dellyNA18507DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/NA18507.delly_q0_c5_del.perf.txt'
pindelNA18507DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/bwa_pindel_D.perf.txt'

cloudbreakNA18507DelPredictionsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/na18507final_readGroupsRazerS3_i94_s99_m1000_None_None_25000_25_sam_6_67eedd7_-1_-1_3_deletions_260.bed'
breakdancerNA18507DelPredictionsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/final_predictions/breakdancer/NA18507_35_2_3.bd_out.dels.pairgte4pt9.bed'
gasvProNA18507DelPredictionsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/final_predictions/GASVPro/BAMToGASV.gasvpro.in.ALL.MCMCThreshold.clusters.pruned.clusters.gte14pt8.bed'
dellyNA18507DelPredictionsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/final_predictions/delly/NA18507.delly_q0_c5_del.gte8pt6.bed'
pindelNA18507DelPredictionsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/final_predictions/pindel/bwa_pindel.noshort.gte40pt7.bed'

cloudbreakNA18507DeletionHitsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/na18507final_readGroupsRazerS3_i94_s99_m1000_None_None_25000_25_sam_6_67eedd7_-1_-1_3_deletions_260.hits.txt'
breakdancerNA18507DeletionHitsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/NA18507_35_2_3.bd_out.dels.pairgte4pt9.hits.txt'
gasvProNA18507DeletionHitsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/BAMToGASV.gasvpro.in.ALL.MCMCThreshold.clusters.pruned.clusters.gte14pt8.hits.txt'
dellyNA18507DeletionHitsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/NA18507.delly_q0_c5_del.gte8pt6.hits.txt'
pindelNA18507DeletionHitsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/bwa_pindel.noshort.gte40pt7.hits.txt'

# insertions
cloudbreakNA18507InsertionsPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/na18507final_readGroupsRazerS3_i94_s99_m1000_None_None_25000_25_sam_6_67eedd7_-1_-1_3_ins_nosplit_mfw5_perf.txt'
breakdancerNA18507InsertionsPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/NA18507_35_2_3_insertions.perf.txt'
pindelNA18507InsertionsPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/pindel_insertions.perf.txt'

cloudbreakNA18507InsertionPredictionsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/na18507final_readGroupsRazerS3_i94_s99_m1000_None_None_25000_25_sam_6_67eedd7_-1_-1_3_insertions_mfw5_nosplit_gt0p5.bed'
breakdancerNA18507InsertionPredictionsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/NA18507_35_2_3.bd_out.insertions.gte2.bed'
pindelNA18507InsertionPredictionsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/pindel_insertions.bed'

cloudbreakNA18507InsertionHitsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/na18507final_readGroupsRazerS3_i94_s99_m1000_None_None_25000_25_sam_6_67eedd7_-1_-1_3_insertions_nosplit_mfw5_gt0p5.hits.txt'
breakdancerNA18507InsertionHitsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/NA18507_35_2_3.bd_out.ins.gte2.hits.txt'
pindelNA18507InsertionHitsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/pindel_insertions.gte4.hits.txt'

# reference files
repMask <- import.bed('~/genomes/1kg/hg18/repeatmasker_b36_merged.bed.gz', asRangedData=FALSE)

#chr2 100bp diploid deletions
totalDels <- 400
cloudbreakChr2DeletionPerf <- read.table(cloudbreakChr2DeletionPerfFile, header=TRUE, sep="\t")
breakdancerChr2DeletionPerf <- read.table(breakdancerChr2DeletionPerfFile, header=TRUE)
gasvChr2DeletionPerf <- read.table(gasvProChr2DeletionPerfFile, header=TRUE)
dellyChr2DeletionPerf <- read.table(dellyChr2DeletionPerfFile, header=TRUE)
pindelChr2DeletionPerf <- read.table(pindelChr2DeletionPerfFile, header=TRUE)

perfsList <- list(cloudbreak=cloudbreakChr2DeletionPerf, breakdancer=breakdancerChr2DeletionPerf, pindel=pindelChr2DeletionPerf, gasv=gasvChr2DeletionPerf, delly=dellyChr2DeletionPerf)
pdf('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/CHR2SIM_DELS_ROC.pdf', width=10)
par(xpd=T, mar=par()$mar+c(0,0,0,7))
plotROC(perfsList, 
        c("Cloudbreak", "Breakdancer","Pindel", "GASVPro", "DELLY"), 
        totalDels, "Deletions in Venter diploid chr2 simulation",legendLoc=xy.coords(425,200), maxTP=350)
dev.off()

# analyze predictions at a given threshold

# first load up information about the deletions
trueDels <- import.bed('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/HuRef.homozygous_indels.061109.chr2.deletions.gff.bed', asRangedData=FALSE)
trueDelshap2 <- import.bed('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/HuRef.homozygous_indels.061109.chr2.sim_hap2.deletions.gff.bed', asRangedData=FALSE)

mcols(trueDels)$tok <- paste(seqnames(trueDels), ':', start(trueDels)-1, '-', end(trueDels), sep="")
mcols(trueDels)$haps <- 1

hap2toks <- paste(seqnames(trueDelshap2), ':', start(trueDelshap2)-1, '-', end(trueDelshap2), sep="")
mcols(trueDels[mcols(trueDels)$tok %in% hap2toks])$haps <- 2

mcols(trueDels)$sizeclasses <- sizeclasses(end(trueDels) - start(trueDels) + 1)

segdups <- import.bed('~/genomes/1kg/hg18/segmental_duplications.bed', asRangedData=FALSE)

mcols(trueDels)$segdup <- 0
mcols(trueDels[as.matrix(findOverlaps(segdups,trueDels))[,2]])$segdup <- 1

mcols(trueDels)$repMask <- 0
mcols(trueDels[as.matrix(findOverlaps(repMask,trueDels))[,2]])$repMask <- 1

# load true positive predictions
cb <- processDeletionPredictions("Cloudbreak", cloudbreakChr2DeletionHitsFile, trueDels)
bd <- processDeletionPredictions('Breakdancer', breakdancerChr2DeletionHitsFile, trueDels)
gasv <- processDeletionPredictions('GASVPro', gasvProChr2DeletionHitsFile, trueDels)
delly <- processDeletionPredictions('DELLY', dellyChr2DeletionHitsFile, trueDels)
pindel <- processDeletionPredictions('Pindel', pindelChr2DeletionHitsFile, trueDels)

allpreds <- rbind(as.data.frame(cb),as.data.frame(bd),as.data.frame(gasv),as.data.frame(delly),as.data.frame(pindel))

allpreds <- ddply(allpreds, .(tok), function(x) { cbind(x, list(discoveries=rep(dim(x)[1],dim(x)[1]))) })

# tabulate against true deletion data
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

# genotyping with extra data
extraData <- read.table(cloudbreakChr2DeletionPredictionsFile, skip=1)
names(extraData) <- c('predchrom', 'predstart', 'predend', 'peaknum', 'maxscore', 'svtype', 'avgMu', 'minMu', 'maxMu', 'avgW', 'minW', 'maxW')

cbpredsWithExtraData <- merge(as.data.frame(cb), extraData)

chr2HapCM <- table(cbpredsWithExtraData$haps, cbpredsWithExtraData$avgW < .2)

#chr2 100bp diploid insertions
totalInsertions <- 80
cloudbreakChr2InsertionPerf <- read.table(cloudbreakChr2InsertionPerfFile, header=TRUE, sep="\t")
cloudbreakChr2InsertionPerf <- cloudbreakChr2InsertionPerf[seq(1,dim(cloudbreakChr2InsertionPerf)[1],by=4),]
breakdancerChr2InsertionPerf <- read.table(breakdancerChr2InsertionPerfFile, header=TRUE)
pindelChr2InsertionPerf <- read.table(pindelChr2InsertionPerfFile, header=TRUE)

perfsListInsertions <- list(cloudbreak=cloudbreakChr2InsertionPerf, breakdancer=breakdancerChr2InsertionPerf, pindel=pindelChr2InsertionPerf)
pdf('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/CHR2SIM_INS_ROC.pdf', width=10)
par(xpd=T, mar=par()$mar+c(0,0,0,7))
plotROC(perfsListInsertions, 
        c("Cloudbreak", "Breakdancer","Pindel"), 
        totalInsertions, "Insertions in Venter diploid chr2 simulation",legendLoc=xy.coords(85,50), maxTP=80)
dev.off()

# analyze predictions at a given threshold

# first load up information about the insertions
trueInsertions <- import.bed('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/HuRef.homozygous_indels.061109.chr2.insertions.gff.bed', asRangedData=FALSE)
trueInsertionsHap2 <- import.bed('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/HuRef.homozygous_indels.061109.chr2.sim_hap2.insertions.gff.bed', asRangedData=FALSE)

mcols(trueInsertions)$length <- as.numeric(mcols(trueInsertions)$name)

mcols(trueInsertions)$tok <- paste(seqnames(trueInsertions), ':', start(trueInsertions)-1, '-', end(trueInsertions), sep="")
mcols(trueInsertions)$haps <- 1

hap2toksInsertions <- paste(seqnames(trueInsertionsHap2), ':', start(trueInsertionsHap2)-1, '-', end(trueInsertionsHap2), sep="")
mcols(trueInsertions[mcols(trueInsertions)$tok %in% hap2toksInsertions])$haps <- 2

mcols(trueInsertions)$sizeclasses <- factor(as.character(sizeclasses(mcols(trueInsertions)$length)))

segdups <- import.bed('~/genomes/1kg/hg18/segmental_duplications.bed', asRangedData=FALSE)

mcols(trueInsertions)$segdup <- 0
mcols(trueInsertions[as.matrix(findOverlaps(segdups,trueInsertions))[,2]])$segdup <- 1

mcols(trueInsertions)$repMask <- 0
mcols(trueInsertions[as.matrix(findOverlaps(repMask,trueInsertions))[,2]])$repMask <- 1

# load true positive predictions
cbInsertions <- processInsertionPredictions("Cloudbreak", cloudbreakChr2InsertionHitsFile, trueInsertions)

bdInsertions <- processInsertionPredictions('Breakdancer', breakdancerChr2InsertionHitsFile, trueInsertions)

pindelInsertions <- processInsertionPredictions('Pindel', pindelChr2InsertionHitsFile, trueInsertions)

allpredsInsertions <- rbind(as.data.frame(cbInsertions),as.data.frame(bdInsertions),as.data.frame(pindelInsertions))

allpredsInsertions <- ddply(allpredsInsertions, .(tok), function(x) { cbind(x, list(discoveries=rep(dim(x)[1],dim(x)[1]))) })

# tabulate against true insertion data
trueInsertionsGte40 <- trueInsertions[mcols(trueInsertions)$length >= 40]
trueInsertionsSizes <- table(trueInsertionsGte40$sizeclasses)
trueInsertionsHaps <- table(trueInsertionsGte40$haps)
trueInsertionsRepmask <- table(trueInsertionsGte40$repMask)

insertionPredsBySize <- table(allpredsInsertions[,c('name','sizeclasses')])
insertionExBySize <- table(allpredsInsertions[which(allpredsInsertions$discoveries == 1),c('name','sizeclasses')])

insertionPredsByRepMask <- table(allpredsInsertions[,c('name','repMask')])
insertionExPredsByRepMask <- table(allpredsInsertions[which(allpredsInsertions$discoveries == 1),c('name','repMask')])

insertionPredsByHap <- table(allpredsInsertions[,c('name','haps')])
insertionExPredsByHap <- table(allpredsInsertions[which(allpredsInsertions$discoveries == 1),c('name','haps')])

insertionPredsBySegDup <- table(allpredsInsertions[,c('name','segdup')])
insertionExPredsBySegDup <- table(allpredsInsertions[which(allpredsInsertions$discoveries == 1),c('name','segdup')])

# genotyping insertions
extraInsertionData <- read.table(cloudbreakChr2InsertionPredictionsFile, skip=1)
names(extraInsertionData) <- c('predchrom', 'predstart', 'predend', 'peaknum', 'maxscore', 'svtype', 'avgMu', 'minMu', 'maxMu', 'avgW', 'minW', 'maxW')

cbInsertionPredsWithExtraData <- merge(as.data.frame(cbInsertions), extraInsertionData, by=c("predchrom", "predstart"))

insertionHapCM <- table(cbInsertionPredsWithExtraData$haps, cbInsertionPredsWithExtraData$minW < .3)
print(insertionHapCM)
print((insertionHapCM[1,1] + insertionHapCM[2,2]) / sum(insertionHapCM))


# NA18507

trueDelsNA18507 <- import.bed('~/Documents/gene_rearrange/svpipeline/NA18507/MillsDIP1KG_merged.bed', asRangedData=FALSE)

mcols(trueDelsNA18507)$tok <- paste(seqnames(trueDelsNA18507), ':', start(trueDelsNA18507)-1, '-', end(trueDelsNA18507), sep="")

mcols(trueDelsNA18507)$sizeclasses <- sizeclasses(end(trueDelsNA18507) - start(trueDelsNA18507) + 1)

repMaskHg19 <- import.bed('~/genomes/1kg/hg19/repeatmasker.bed', asRangedData=FALSE)

mcols(trueDelsNA18507)$repMask <- 0
mcols(trueDelsNA18507[as.matrix(findOverlaps(repMaskHg19,trueDelsNA18507))[,2]])$repMask <- 1

trueDelsNA18507Gte40 <- trueDelsNA18507[end(trueDelsNA18507) - start(trueDelsNA18507) >= 40]
trueDelsNA18507Sizes <- table(trueDelsNA18507Gte40$sizeclasses)

trueDelRepmaskNA18507 <- table(trueDelsNA18507Gte40$repMask)

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
NA18507cloudbreakDelsPerf <- read.table(cloudbreakNA18507DeletionPerfFile, header=TRUE, sep="\t")
NA18507breakdancerDelsPerf <- read.table(breakdancerNA18507DeletionPerfFile, header=TRUE)
NA18507gasvDelsPerf <- read.table(gasvProNA18507DeletionPerfFile, header=TRUE)
NA18507dellyDelsPerf <- read.table(dellyNA18507DeletionPerfFile, header=TRUE)
NA18507pindelDelsPerf <- read.table(pindelNA18507DeletionPerfFile, header=TRUE)

perfsList <- list(cloudbreak=NA18507cloudbreakDelsPerf, breakdancer=NA18507breakdancerDelsPerf, pindel=NA18507pindelDelsPerf, gasv=NA18507gasvDelsPerf, delly=NA18507dellyDelsPerf)
pdf('~/Documents/svpipeline/manuscript/NA18507_DELS_ROC.pdf', width=10)
par(xpd=T, mar=par()$mar+c(0,0,0,7))
plotROC(perfsList, c("Cloudbreak", "Breakdancer","Pindel","GASVPro", "DELLY", "Cloudbreak"), totalDels, "NA18507", legendLoc=xy.coords(10750,700), maxTP=1300, sim=FALSE)
dev.off()

# break down predictons

NA18507delCutoffs=list(cloudbreak=2.60, breakdancer=4.9, GASVPro=14.8, DELLY=8.6, Pindel=40.7)

totalCloudbreakNA18507DelPredictions <- dim(read.table(cloudbreakNA18507DelPredictionsFile, skip=1))[1]
totalBreakdancerNA18507DelPredictions <- dim(read.table(breakdancerNA18507DelPredictionsFile))[1]
totalGASVProNA18507DelPredictions <- dim(read.table(gasvProNA18507DelPredictionsFile))[1]
totalDELLYNA18507DelPredictions <- dim(read.table(dellyNA18507DelPredictionsFile))[1]
totalPindelNA18507DelPredictions <- dim(read.table(pindelNA18507DelPredictionsFile))[1]

totalNA18507DelPredictions <- list(Cloudbreak=totalCloudbreakNA18507DelPredictions, Breakdancer=totalBreakdancerNA18507DelPredictions, GASVPro=totalGASVProNA18507DelPredictions, DELLY=totalDELLYNA18507DelPredictions, Pindel=totalPindelNA18507DelPredictions)

cbHitsNA18507 <- processDeletionPredictions("Cloudbreak", cloudbreakNA18507DeletionHitsFile, trueDelsNA18507)
bdHitsNA18507 <- processDeletionPredictions('Breakdancer', breakdancerNA18507DeletionHitsFile, trueDelsNA18507)
gasvHitsNA18507 <- processDeletionPredictions('GASVPro', gasvProNA18507DeletionHitsFile, trueDelsNA18507)
dellyHitsNA18507 <- processDeletionPredictions('DELLY', dellyNA18507DeletionHitsFile, trueDelsNA18507)
pindelHitsNA18507 <- processDeletionPredictions('Pindel', pindelNA18507DeletionHitsFile, trueDelsNA18507)

totalNumNA18507DeletionHits <- list(list(Cloudbreak=dim(cbHitsNA18507)[1], Breakdancer=dim(bdHitsNA18507)[1], GASVPro=dim(gasvHitsNA18507)[1], DELLY=dim(dellyHitsNA18507)[1], Pindel=dim(pindelHitsNA18507)[1]))

NA18507DelPrecision <- unlist(totalNumNA18507DeletionHits) / unlist(totalNA18507DelPredictions)
NA18507DelRecall <- unlist(totalNumNA18507DeletionHits) / length(trueDelsNA18507Gte40)

allpredsNA18507 <- rbind(as.data.frame(cbHitsNA18507),as.data.frame(bdHitsNA18507),as.data.frame(gasvHitsNA18507),as.data.frame(dellyHitsNA18507),as.data.frame(pindelHitsNA18507))
allpredsNA18507 <- ddply(allpredsNA18507, .(tok), function(x) { cbind(x, list(discoveries=rep(dim(x)[1],dim(x)[1]))) })

predsBySizeNA18507 <- table(allpredsNA18507[,c('name','sizeclasses')])
exBySizeNA18507 <- table(allpredsNA18507[which(allpredsNA18507$discoveries == 1),c('name','sizeclasses')])

predsByRepMaskNA18507 <- table(allpredsNA18507[,c('name','repMask')])
exByRepMaskNA18507 <- table(allpredsNA18507[which(allpredsNA18507$discoveries == 1),c('name','repMask')])

# genotyping 

extraDataNA18507 <- read.table(cloudbreakNA18507DelPredictionsFile, skip=1)
names(extraDataNA18507) <- c('predchrom', 'predstart', 'predend', 'peaknum', 'maxscore', 'svType', 'avgMu', 'minMu', 'maxMu', 'avgW', 'minW', 'maxW')

cbpredsWithExtraDataNA18507 <- merge(as.data.frame(cbHitsNA18507), extraDataNA18507)

cbpredsWithExtraDataNA18507Genotyped <- cbpredsWithExtraDataNA18507[!is.na(cbpredsWithExtraDataNA18507$hap),]

NA18507HapCM <- table(cbpredsWithExtraDataNA18507$hap, cbpredsWithExtraDataNA18507$avgW < .2)

millsDelsWithGenotypes <- read.table('~/Documents/Papers2/Articles/2011/Mills/Supplemental/dels_gt40_with_genotypes.csv', header=TRUE, sep=",")
millsGenotypeData <- read.table('~/Documents/Papers2/Articles/2011/Mills/Supplemental/Mills_genotypes_NA18507_with_ids_gt40.txt', header=TRUE, sep="\t")

millsFullGenotypedDels <- merge(millsDelsWithGenotypes,millsGenotypeData)

millsGenotypedRanges <- GRanges(seqnames=substr(millsFullGenotypedDels$CHR,4,length(millsFullGenotypedDels$CHR)), ranges=IRanges(start=millsFullGenotypedDels$START, end=millsFullGenotypedDels$STOP), mcols=DataFrame(hap=millsFullGenotypedDels$SED00002_NA18507.CEL))

ol <- findOverlaps(trueDelsNA18507, millsGenotypedRanges)  

mcols(trueDelsNA18507[as.matrix(ol)[,1]])$hap <- mcols(na185071KGDels[as.matrix(ol)[,2]])$hap

#chr2 100bp diploid insertions
totalInsertions <- 20000
cloudbreakNA18507Insertions <- read.table(cloudbreakNA18507InsertionsPerfFile, header=TRUE, sep="\t")
cloudbreakNA18507Insertions <- cloudbreakNA18507Insertions[seq(1,dim(cloudbreakNA18507Insertions)[1],by=4),]
breakdancerNA18507Insertions <- read.table(breakdancerNA18507InsertionsPerfFile, header=TRUE)
pindelNA18507Insertions <- read.table(pindelNA18507InsertionsPerfFile, header=TRUE)

perfsListInsertionsNA18507 <- list(cloudbreak=cloudbreakNA18507Insertions, breakdancer=breakdancerNA18507Insertions, pindel=pindelNA18507Insertions)
pdf('~/Documents/gene_rearrange/svpipeline/NA18507/NA18507_INS_ROC.pdf', width=10)
par(xpd=T, mar=par()$mar+c(0,0,0,7))
plotROC(perfsListInsertionsNA18507, 
        c("Cloudbreak", "Breakdancer","Pindel"), 
        totalInsertions, "Insertions in NA18507",legendLoc=xy.coords(22000,100), maxTP=300, sim=FALSE)
dev.off()

# analyze predictions at a given threshold

# first load up information about the insertions
trueInsertionsNA18507 <- import.bed('~/Documents/gene_rearrange/svpipeline/NA18507/MillsDIP1KG_insertions_merged.bed', asRangedData=FALSE)
mcols(trueInsertionsNA18507)$length <- as.numeric(mcols(trueInsertionsNA18507)$name)

mcols(trueInsertionsNA18507)$tok <- paste(seqnames(trueInsertionsNA18507), ':', start(trueInsertionsNA18507)-1, '-', end(trueInsertionsNA18507), sep="")

mcols(trueInsertionsNA18507)$shorttok <- paste(seqnames(trueInsertionsNA18507), ':', start(trueInsertionsNA18507)-1)

# why do I need to do this?? ARGH
# otherwise:  length of 'dimnames' [2] not equal to array extent
mcols(trueInsertionsNA18507)$sizeclasses <- factor(as.character(sizeclasses(mcols(trueInsertionsNA18507)$length)))

mcols(trueInsertionsNA18507)$repMask <- 0
mcols(trueInsertionsNA18507[as.matrix(findOverlaps(repMaskHg19,trueInsertionsNA18507))[,2]])$repMask <- 1

trueInsertionsNA18507Sizes <- table(trueInsertionsNA18507Gte40$sizeclasses)

trueInsertionsNA18507Repmask <- table(trueInsertionsNA18507Gte40$repMask)

na185071KGInsLengths <- nchar(unlist(fixed(na185071KGVariants)$ALT)) - nchar(fixed(na185071KGVariants)$REF)
na185071KGIns <- rowData(na185071KGVariants)[na185071KGInsLengths >= 40]


mcols(trueInsertionsNA18507)$hap <- NA

# right now not getting any genotypes from the 1KG insertions
ol <- findOverlaps(trueInsertionsNA18507, na185071KGIns) 
mcols(trueInsertionsNA18507[as.matrix(ol)[,1]])$hap <- mcols(trueInsertionsNA18507[as.matrix(ol)[,2]])$hap

millsInsertionsWithIds <- read.table('/Users/cwhelan/Documents/Papers2/Articles/2011/Mills/Supplemental/Mills_2011_NA18507_INS.b37.bed')
names(millsInsertionsWithIds) <- c("chrom", "start", "end", "length", "probeset_id")

millsGenotypes <- read.table('/Users/cwhelan/Documents/Papers2/Articles/2011/Mills/Supplemental/Mills_genotypes_NA18507_with_ids.txt', header=TRUE, sep="\t")

genotypedMillsInsertions <- merge(millsInsertionsWithIds, millsGenotypes, by="probeset_id")
 genotypedMillsInsertions$shorttok <- paste(genotypedMillsInsertions$chrom, ':', genotypedMillsInsertions$start)

totalCloudbreakNA18507InsertionPredictions <- dim(read.table(cloudbreakNA18507InsertionPredictionsFile, skip=1))[1]
totalBreakdancerNA18507InsertionPredictions <- dim(read.table(breakdancerNA18507InsertionPredictionsFile))[1]
totalPindelNA18507InsertionPredictions <- dim(read.table(pindelNA18507InsertionPredictionsFile))[1]

# load true positive predictions
cbInsertionsNA18507 <- processInsertionPredictions("Cloudbreak", cloudbreakNA18507InsertionHitsFile, trueInsertionsNA18507)
bdInsertionsNA18507 <- processInsertionPredictions('Breakdancer', breakdancerNA18507InsertionHitsFile, trueInsertionsNA18507)
pindelInsertionsNA18507 <- processInsertionPredictions('Pindel', pindelNA18507InsertionHitsFile, trueInsertionsNA18507)

allpredsInsertionsNA18507 <- rbind(as.data.frame(cbInsertionsNA18507),as.data.frame(bdInsertionsNA18507),as.data.frame(pindelInsertionsNA18507))
allpredsInsertionsNA18507 <- ddply(allpredsInsertionsNA18507, .(tok), function(x) { cbind(x, list(discoveries=rep(dim(x)[1],dim(x)[1]))) })

totalNA18507InsertionPredictions <- list(Cloudbreak=totalCloudbreakNA18507InsertionPredictions, Breakdancer=totalBreakdancerNA18507InsertionPredictions, Pindel=totalPindelNA18507InsertionPredictions)
totalNumNA18507InsertionHits <- list(list(Cloudbreak=dim(cbInsertionsNA18507)[1], Breakdancer=dim(bdInsertionsNA18507)[1], Pindel=dim(pindelInsertionsNA18507)[1]))

# tabulate against true insertion data
trueInsertionsNA18507Gte40 <- trueInsertionsNA18507[mcols(trueInsertionsNA18507)$length >= 40]
trueInsertionsNA18507Sizes <- table(trueInsertionsNA18507Gte40$sizeclasses)
trueInsertionsNA18507Haps <- table(trueInsertionsNA18507Gte40$haps)
trueInsertionsNA18507Repmask <- table(trueInsertionsNA18507Gte40$repMask)

NA18507InsertionPrecision <- unlist(totalNumNA18507InsertionHits) / unlist(totalNA18507InsertionPredictions)
NA18507InsertionRecall <- unlist(totalNumNA18507InsertionHits) / length(trueInsertionsNA18507Gte40)

NA18507insertionPredsBySize <- table(allpredsInsertionsNA18507[,c('name','sizeclasses')])
NA18507insertionExPredsBySize <- table(allpredsInsertionsNA18507[which(allpredsInsertionsNA18507$discoveries == 1),c('name','sizeclasses')])

NA18507insertionPredsByRepMask <- table(allpredsInsertionsNA18507[,c('name','repMask')])
NA18507insertionExPredsByRepMask <- table(allpredsInsertionsNA18507[which(allpredsInsertionsNA18507$discoveries == 1),c('name','repMask')])

NA18507insertionPredsByHap <- table(allpredsInsertionsNA18507[,c('name','hap')])
NA18507insertionExPredsByHap <- table(allpredsInsertions[which(allpredsInsertionsNA18507$discoveries == 1),c('name','haps')])

#NA18507insertionPredsBySegDup <- table(allpredsInsertionsNA18507[,c('name','segdup')])
#NA18507insertionExPredsBySegDup <- table(allpredsInsertions[which(allpredsInsertionsNA18507$discoveries == 1),c('name','segdup')])

# genotyping insertions
#cloudbreakNA18507InsertionsWithExtraData <- read.table(cloudbreakNA18507InsertionPredictionsFile, skip=1)
#names(cloudbreakNA18507InsertionsWithExtraData) <- c('predchrom', 'predstart', 'predend', 'peaknum', 'maxscore', 'svtype', 'avgMu', 'minMu', 'maxMu', 'avgW', 'minW', 'maxW')

#NA18507cbInsertionPredsWithExtraData <- merge(as.data.frame(cbInsertionsNA18507), cloudbreakNA18507InsertionsWithExtraData)

#NA18507insertionHapCM <- table(NA18507cbInsertionPredsWithExtraData$hap, NA18507cbInsertionPredsWithExtraData$avgW < .2)
#print(NA18507insertionHapCM)
#print((NA18507insertionHapCM[1,1] + NA18507insertionHapCM[2,2]) / sum(NA18507insertionHapCM))

# LaTeX tables

# chr 2 Deletions by size
tableEnv = new.env()
assign("totals", trueDelSizes, envir=tableEnv)
assign("preds", predsBySize, envir=tableEnv)
assign("exPreds", exBySize, envir=tableEnv)
brew('~/Documents/svpipeline/manuscript/chr2DeletionPreds.table.brew.tex', env=tableEnv)

# chr 2 Insertions by size
tableEnv = new.env()
assign("totals", trueInsertionsSizes, envir=tableEnv)
assign("preds", insertionPredsBySize, envir=tableEnv)
assign("exPreds", insertionExBySize, envir=tableEnv)
brew('~/Documents/svpipeline/manuscript/chr2InsertionPredsBySize.table.brew.tex', env=tableEnv)

# NA18507 Deletions by size
tableEnv = new.env()
assign("totals", trueDelsNA18507Sizes, envir=tableEnv)
assign("precision", NA18507DelPrecision, envir=tableEnv)
assign("recall", NA18507DelRecall, envir=tableEnv)
assign("preds", predsBySizeNA18507, envir=tableEnv)
assign("exPreds", exBySizeNA18507, envir=tableEnv)
brew('~/Documents/svpipeline/manuscript/NA18507DeletionPredsBySize.table.brew.tex', env=tableEnv)

# NA18507 Insertions by size
tableEnv = new.env()
assign("totals", trueInsertionsNA18507Sizes, envir=tableEnv)
assign("precision", NA18507InsertionPrecision, envir=tableEnv)
assign("recall", NA18507InsertionRecall, envir=tableEnv)
assign("preds", NA18507insertionPredsBySize, envir=tableEnv)
assign("exPreds", NA18507insertionExPredsBySize, envir=tableEnv)
brew('~/Documents/svpipeline/manuscript/NA18507InsertionPredsBySize.table.brew.tex', env=tableEnv)

# deletions in repetitive regions
tableEnv = new.env()
assign("chr2Totals", trueDelRepmask, envir=tableEnv)
assign("chr2Preds", predsByRepMask, envir=tableEnv)
assign("chr2ExPreds", exPredsByRepMask, envir=tableEnv)
assign("NA18507Totals", trueDelRepmaskNA18507, envir=tableEnv)
assign("NA18507Preds", predsByRepMaskNA18507, envir=tableEnv)
assign("NA18507ExPreds", exByRepMaskNA18507, envir=tableEnv)
brew('~/Documents/svpipeline/manuscript/deletionRepmask.table.brew.tex', env=tableEnv)

# insertions in repetitive regions
tableEnv = new.env()
assign("chr2Totals", trueInsertionsRepmask, envir=tableEnv)
assign("chr2Preds", insertionPredsByRepMask, envir=tableEnv)
assign("chr2ExPreds", insertionExPredsByRepMask, envir=tableEnv)
assign("NA18507Totals", trueInsertionsNA18507Repmask, envir=tableEnv)
assign("NA18507Preds", NA18507insertionPredsByRepMask, envir=tableEnv)
assign("NA18507ExPreds", NA18507insertionExPredsByRepMask, envir=tableEnv)
brew('~/Documents/svpipeline/manuscript/insertionRepmask.table.brew.tex', env=tableEnv)

# genotyping
tableEnv = new.env()
assign("chr2HapCM", chr2HapCM, envir=tableEnv)
assign("NA18507HapCM", NA18507HapCM, envir=tableEnv)
brew('~/Documents/svpipeline/manuscript/deletionGenotypeAccuracy.table.brew.tex', env=tableEnv)

