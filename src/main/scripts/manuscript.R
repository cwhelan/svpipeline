library(GenomicRanges)
library(rtracklayer)
library(plyr)
library(RColorBrewer)
library(VariantAnnotation)
library(xtable)
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
  print(head(maxhits))
  maxhits <- merge(maxhits,mcols(trueInsertions), by="tok")
  print(head(maxhits))  
  maxhits$name <- name
  maxhits
}

printPredTableValue <- function(t, x, y) {
  d <- t[x,y]
  if (t[x,y] == max(t[,y])) {
    return(paste0("\\textbf{",formatC(d,digits=3,format="f"),"}"))
  } else {
    return(format(d,digits=3,format="f"))
  }
}

printPredTableValueFromList <- function(x, l) {
  d <- l[x]
  if (d == max(l)) {
    return(paste0("\\textbf{",formatC(d,digits=3,format="f"),"}"))
  } else {
    return(formatC(d,digits=3,format="f"))
  }
}

printPredTableCell <- function(preds, ex, x, y) {
  return(paste0(printPredTableValue(preds, x, y), " (", printPredTableValue(ex, x, y), ")"))
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

repMask <- import.bed('~/genomes/1kg/hg18/repeatmasker_b36_merged.bed.gz', asRangedData=FALSE)

mcols(trueDels)$repMask <- 0
mcols(trueDels[as.matrix(findOverlaps(repMask,trueDels))[,2]])$repMask <- 1

# load true positive predictions
cb <- processDeletionPredictions("Cloudbreak", '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/cloudbreak_venterchr2allindels100bpdip_readgroup1_deletions_210.hits.txt', trueDels)

bd <- processDeletionPredictions('Breakdancer', '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort_35_2_3.bd_out.dels.pairgte4.hits.txt', trueDels)

gasv <- processDeletionPredictions('GASVPro', '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/BAMToGASV.gasvpro.in.ALL.MCMCThreshold.clusters.pruned.clusters.gte12.hits.txt', trueDels)

delly <- processDeletionPredictions('DELLY', '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort.delly_q0_c5_del.gte7.hits.txt', trueDels)

pindel <- processDeletionPredictions('Pindel', '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort_pindel_D_gte33.hits.txt', trueDels)

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
extraData <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/cloudbreak_venterchr2allindels100bpdip_readgroup1_deletions_210.bed', skip=1)
names(extraData) <- c('predchrom', 'predstart', 'predend', 'peaknum', 'maxscore', 'svtype', 'avgMu', 'minMu', 'maxMu', 'avgW', 'minW', 'maxW')

cbpredsWithExtraData <- merge(as.data.frame(cb), extraData)

chr2HapCM <- table(cbpredsWithExtraData$haps, cbpredsWithExtraData$avgW < .2)

#chr2 100bp diploid insertions
totalInsertions <- 80
cloudbreakInsertions <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/cloudbreak_venterchr2allindels100bpdip_readgroup1_ins_mfw3_gt40_nosplit.perf.txt', header=TRUE, sep="\t")
cloudbreakInsertions <- cloudbreakInsertions[seq(1,dim(cloudbreakInsertions)[1],by=4),]
breakdancerInsertions <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/breakdancer_insertions.perf.txt', header=TRUE)
pindelInsertions <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/pindel_insertions.perf.txt', header=TRUE)

perfsListInsertions <- list(cloudbreak=cloudbreakInsertions, breakdancer=breakdancerInsertions, pindel=pindelInsertions)
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

repMask <- import.bed('~/genomes/1kg/hg18/repeatmasker_b36_merged.bed.gz', asRangedData=FALSE)

mcols(trueInsertions)$repMask <- 0
mcols(trueInsertions[as.matrix(findOverlaps(repMask,trueInsertions))[,2]])$repMask <- 1

# load true positive predictions
cbInsertions <- processInsertionPredictions("Cloudbreak", '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/cloudbreak_venterchr2allindels100bpdip_readgroup1_ins_mfw3_gt40_nosplit.gt0p5.hits.txt', trueInsertions)

bdInsertions <- processInsertionPredictions('Breakdancer', '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort_35_2_3.bd_out.ins.gte3.hits.txt', trueInsertions)

pindelInsertions <- processInsertionPredictions('Pindel', '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/pindel_insertions.hits.txt', trueInsertions)

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
extraInsertionData <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/cloudbreak_venterchr2allindels100bpdip_readgroup1_insertions_mfw3_gt40_nosplit_gt0p5.bed', skip=1)
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
NA18507cloudbreakDelsPerfDelsPerf <- read.table('~/Documents/gene_rearrange/svpipeline/NA18507/new_gold_na18507final_readGroupsRazerS3_i94_s99_m1000_None_None_25000_25_sam_6_67eedd7_-1_-1_3_lrHeterozygous.allsvs.short40.perf.txt', header=TRUE, sep="\t")
NA18507breakdancerDelsPerf <- read.table('~/Documents/gene_rearrange/svpipeline/NA18507/NA18507_35_2_3_singlecpu_perf.txt', header=TRUE)
NA18507gasvDelsPerf <- read.table('~/Documents/gene_rearrange/svpipeline/NA18507/gasvpro.perf.txt', header=TRUE)
NA18507dellyDelsPerf <- read.table('~/Documents/gene_rearrange/svpipeline/NA18507/NA18507.delly_q0_c5_del.perf.txt', header=TRUE)
NA18507pindelDelsPerf <- read.table('~/Documents/gene_rearrange/svpipeline/NA18507/bwa_pindel_D.perf.txt', header=TRUE)

perfsList <- list(cloudbreak=NA18507cloudbreakDelsPerf, breakdancer=NA18507breakdancerDelsPerf, pindel=NA18507pindelDelsPerf, gasv=NA18507gasvDelsPerf, delly=NA18507dellyDelsPerf)
pdf('~/Documents/svpipeline/manuscript/NA18507_ROC.pdf', width=10)
par(xpd=T, mar=par()$mar+c(0,0,0,7))
plotROC(perfsList, c("Cloudbreak", "Breakdancer","Pindel","GASVPro", "DELLY", "Cloudbreak"), totalDels, "NA18507", legendLoc=xy.coords(10750,700), maxTP=1300, sim=FALSE)
dev.off()

# break down predictons

NA18507delCutoffs=list(cloudbreak=2.60, breakdancer=4.9, GASVPro=14.8, DELLY=8.6, Pindel=40.7)

totalCloudbreakNA18507DelPredictions <- dim(read.table('~/Documents/gene_rearrange/svpipeline/NA18507/na18507final_readGroupsRazerS3_i94_s99_m1000_None_None_25000_25_sam_6_67eedd7_-1_-1_3_deletions_260.bed', skip=1))[1]
totalBreakdancerNA18507DelPredictions <- dim(read.table('~/Documents/gene_rearrange/svpipeline/NA18507/final_predictions/breakdancer/NA18507_35_2_3.bd_out.dels.pairgte4pt9.bed'))[1]
totalGASVProNA18507DelPredictions <- dim(read.table('~/Documents/gene_rearrange/svpipeline/NA18507/final_predictions/GASVPro/BAMToGASV.gasvpro.in.ALL.MCMCThreshold.clusters.pruned.clusters.gte14pt8.bed'))[1]
totalDELLYNA18507DelPredictions <- dim(read.table('~/Documents/gene_rearrange/svpipeline/NA18507/final_predictions/delly/NA18507.delly_q0_c5_del.gte8pt6.bed'))[1]
totalPindelNA18507DelPredictions <- dim(read.table('~/Documents/gene_rearrange/svpipeline/NA18507/final_predictions/pindel/bwa_pindel.noshort.gte40pt7.bed'))[1]

totalNA18507DelPredictions <- list(Cloudbreak=totalCloudbreakNA18507DelPredictions, Breakdancer=totalBreakdancerNA18507DelPredictions, GASVPro=totalGASVProNA18507DelPredictions, DELLY=totalDELLYNA18507DelPredictions, Pindel=totalPindelNA18507DelPredictions)

cbHitsNA18507 <- processDeletionPredictions("Cloudbreak", '~/Documents/gene_rearrange/svpipeline/NA18507/na18507final_readGroupsRazerS3_i94_s99_m1000_None_None_25000_25_sam_6_67eedd7_-1_-1_3_deletions_260.hits.txt', trueDelsNA18507)
bdHitsNA18507 <- processDeletionPredictions('Breakdancer', '~/Documents/gene_rearrange/svpipeline/NA18507/NA18507_35_2_3.bd_out.dels.pairgte4pt9.hits.txt', trueDelsNA18507)
gasvHitsNA18507 <- processDeletionPredictions('GASVPro', '~/Documents/gene_rearrange/svpipeline/NA18507/BAMToGASV.gasvpro.in.ALL.MCMCThreshold.clusters.pruned.clusters.gte14pt8.hits.txt', trueDelsNA18507)
dellyHitsNA18507 <- processDeletionPredictions('DELLY', '~/Documents/gene_rearrange/svpipeline/NA18507/NA18507.delly_q0_c5_del.gte8pt6.hits.txt', trueDelsNA18507)
pindelHitsNA18507 <- processDeletionPredictions('Pindel', '~/Documents/gene_rearrange/svpipeline/NA18507/bwa_pindel.noshort.gte40pt7.hits.txt', trueDelsNA18507)

totalNumNA18507Hits <- list(list(Cloudbreak=dim(cbHitsNA18507)[1], Breakdancer=dim(bdHitsNA18507)[1], GASVPro=dim(gasvHitsNA18507)[1], DELLY=dim(dellyHitsNA18507)[1], Pindel=dim(pindelHitsNA18507)[1]))

NA18507DelPrecision <- unlist(totalNumNA18507Hits) / unlist(totalNA18507DelPredictions)
NA18507DelRecall <- unlist(totalNumNA18507Hits) / length(trueDelsNA18507Gte40)

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

NA18507HapCM <- table(cbpredsWithExtraDataNA18507$hap, cbpredsWithExtraDataNA18507$avgW < .2)

millsDelsWithGenotypes <- read.table('~/Documents/Papers2/Articles/2011/Mills/Supplemental/dels_gt40_with_genotypes.csv', header=TRUE, sep=",")
millsGenotypeData <- read.table('~/Documents/Papers2/Articles/2011/Mills/Supplemental/Mills_genotypes_NA18507_with_ids_gt40.txt', header=TRUE, sep="\t")

millsFullGenotypedDels <- merge(millsDelsWithGenotypes,millsGenotypeData)

millsGenotypedRanges <- GRanges(seqnames=substr(millsFullGenotypedDels$CHR,4,length(millsFullGenotypedDels$CHR)), ranges=IRanges(start=millsFullGenotypedDels$START, end=millsFullGenotypedDels$STOP), mcols=DataFrame(hap=millsFullGenotypedDels$SED00002_NA18507.CEL))

ol <- findOverlaps(trueDelsNA18507, millsGenotypedRanges)  

mcols(trueDelsNA18507[as.matrix(ol)[,1]])$hap <- mcols(na185071KGDels[as.matrix(ol)[,2]])$hap

#chr2 100bp diploid insertions
totalInsertions <- 20000
cloudbreakNA18507Insertions <- read.table('~/Documents/gene_rearrange/svpipeline/NA18507/na18507final_readGroupsRazerS3_i94_s99_m1000_None_None_25000_25_sam_6_67eedd7_-1_-1_3_ins_nosplit_mfw3_perf.txt', header=TRUE, sep="\t")
cloudbreakNA18507Insertions <- cloudbreakNA18507Insertions[seq(1,dim(cloudbreakNA18507Insertions)[1],by=4),]
breakdancerNA18507Insertions <- read.table('~/Documents/gene_rearrange/svpipeline/NA18507/NA18507_35_2_3_insertions.perf.txt', header=TRUE)
pindelNA18507Insertions <- read.table('~/Documents/gene_rearrange/svpipeline/NA18507/pindel_insertions.perf.txt', header=TRUE)

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

trueInsertionsNA18507Gte40 <- trueInsertionsNA18507[mcols(trueInsertionsNA18507)$length >= 40]
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

# load true positive predictions
cbInsertionsNA18507 <- processInsertionPredictions("Cloudbreak", '~/Documents/gene_rearrange/svpipeline/NA18507/', trueInsertionsNA18507)

bdInsertionsNA18507 <- processInsertionPredictions('Breakdancer', '~/Documents/gene_rearrange/svpipeline/NA18507/NA18507_35_2_3.bd_out.ins.gte2.hits.txt', trueInsertionsNA18507)

pindelInsertionsNA18507 <- processInsertionPredictions('Pindel', '~/Documents/gene_rearrange/svpipeline/NA18507//pindel_insertions.hits.txt', trueInsertions)

allpredsInsertionsNA18507 <- rbind(as.data.frame(cbInsertionsNA18507),as.data.frame(bdInsertionsNA18507),as.data.frame(pindelInsertionsNA18507))

allpredsInsertionsNA18507 <- ddply(allpredsInsertionsNA18507, .(tok), function(x) { cbind(x, list(discoveries=rep(dim(x)[1],dim(x)[1]))) })

# tabulate against true insertion data
trueInsertionsNA18507Gte40 <- trueInsertionsNA18507[end(trueInsertions) - start(trueInsertions) >= 40]
trueInsertionsNA18507Sizes <- table(trueInsertionsNA18507Gte40$sizeclasses)
trueInsertionsNA18507Haps <- table(trueInsertionsNA18507Gte40$haps)
trueInsertionsNA18507Repmask <- table(trueInsertionsNA18507Gte40$repMask)

insertionPredsBySize <- table(allpredsInsertions[,c('name','sizeclasses')])
insertionExBySize <- table(allpredsInsertions[which(allpredsInsertions$discoveries == 1),c('name','sizeclasses')])

insertionPredsByRepMask <- table(allpredsInsertions[,c('name','repMask')])
insertionExPredsByRepMask <- table(allpredsInsertions[which(allpredsInsertions$discoveries == 1),c('name','repMask')])

insertionPredsByHap <- table(allpredsInsertions[,c('name','haps')])
insertionExPredsByHap <- table(allpredsInsertions[which(allpredsInsertions$discoveries == 1),c('name','haps')])

insertionPredsBySegDup <- table(allpredsInsertions[,c('name','segdup')])
insertionExPredsBySegDup <- table(allpredsInsertions[which(allpredsInsertions$discoveries == 1),c('name','segdup')])

# genotyping insertions
extraInsertionData <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/cloudbreak_venterchr2allindels100bpdip_readgroup1_insertions_gt0p5.bed', skip=1)
names(extraInsertionData) <- c('predchrom', 'predstart', 'predend', 'peaknum', 'maxscore', 'svtype', 'avgMu', 'minMu', 'maxMu', 'avgW', 'minW', 'maxW')

cbInsertionPredsWithExtraData <- merge(as.data.frame(cbInsertions), extraInsertionData)

insertionHapCM <- table(cbInsertionPredsWithExtraData$haps, cbInsertionPredsWithExtraData$avgW < .2)
print(insertionHapCM)
print((insertionHapCM[1,1] + insertionHapCM[2,2]) / sum(insertionHapCM))

# LaTeX tables

# chr 2 by size
tableEnv = new.env()
assign("totals", trueDelSizes, envir=tableEnv)
assign("preds", predsBySize, envir=tableEnv)
assign("exPreds", exBySize, envir=tableEnv)
brew('~/Documents/svpipeline/manuscript/chr2preds.table.brew.tex', env=tableEnv)

# NA18507 by size
tableEnv = new.env()
assign("totals", trueDelsNA18507Sizes, envir=tableEnv)
assign("precision", NA18507DelPrecision, envir=tableEnv)
assign("recall", NA18507DelRecall, envir=tableEnv)
assign("totals", trueDelsNA18507Sizes, envir=tableEnv)
assign("preds", predsBySizeNA18507, envir=tableEnv)
assign("exPreds", exBySizeNA18507, envir=tableEnv)
brew('~/Documents/svpipeline/manuscript/NA18507preds.table.brew.tex', env=tableEnv)

# repetitive regions
tableEnv = new.env()
assign("chr2Totals", trueDelRepmask, envir=tableEnv)
assign("chr2Preds", predsByRepMask, envir=tableEnv)
assign("chr2ExPreds", exPredsByRepMask, envir=tableEnv)
assign("NA18507Totals", trueDelRepmaskNA18507, envir=tableEnv)
assign("NA18507Preds", predsByRepMaskNA18507, envir=tableEnv)
assign("NA18507ExPreds", exByRepMaskNA18507, envir=tableEnv)
brew('~/Documents/svpipeline/manuscript/repmask.table.brew.tex', env=tableEnv)

# genotyping
tableEnv = new.env()
assign("chr2HapCM", chr2HapCM, envir=tableEnv)
assign("NA18507HapCM", NA18507HapCM, envir=tableEnv)
brew('~/Documents/svpipeline/manuscript/genotypeaccuracy.table.brew.tex', env=tableEnv)


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