library(GenomicRanges)
library(rtracklayer)
library(plyr)

drawPerfLine <- function(perf, col, maxFP) {
  calls <- c(perf$Calls, min(perf$Calls) - min(perf$TP), 0)
  tp <- c(perf$TP, 0, 0)
  fp <- calls - tp
  lines(fp[fp <= maxFP], tp[fp <= maxFP], col=col, lwd=3)
}

plotROC <- function(perfs, perfNames, totalDels, main, sim=TRUE) {
  maxFP <- totalDels
  plot(0, type="n", ylim=c(0, totalDels), xlim=c(0,maxFP), xlab=ifelse(sim, "False Positives", "Novel Predictions"), ylab="True Positives", main=main)
  perfCols <- rainbow(length(perfs))
  mapply(drawPerfLine, perfs, perfCols, MoreArgs=list(maxFP=maxFP))  
  legend("bottomright", legend=perfNames, col=perfCols, lwd=3, cex=.9)
}

overlap <- function(a,b,x,y) {
  pmin(b,y) - pmax(a,x)
}

sizeclasses <- function(starts, ends) {
  sizelens <- ends - starts + 1
  sizeclasses <- rep(NA, length(sizelens))
  sizeclasses[sizelens <= 50] <- "25-50bp"
  sizeclasses[sizelens > 50 & sizelens <= 100] <- "51-100bp"
  sizeclasses[sizelens > 100 & sizelens <= 250] <- "101-250bp"
  sizeclasses[sizelens > 250 & sizelens <= 500] <- "251-500bp"
  sizeclasses[sizelens > 500 & sizelens <= 1000] <- "501-1000bp"
  sizeclasses[sizelens > 1000] <- "> 1000bp"  
  factor(sizeclasses, levels=c("25-50bp","51-100bp","101-250bp","251-500bp","501-1000bp","> 1000bp"),ordered=TRUE)
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

#chr2 gt 50 100bp diploid
totalDels <- 500
cloudbreak <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_100bp_dip/chr2mwf5_readGroupsRazerS3_new3_i94_s99_m1000_None_None_25000_25_sam_6_32fc04d_-1_-1_3_lrHeterozygous.e1df167.perf.txt', header=TRUE, sep="\t")
breakdancer <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_100bp_dip/final_predictions/breakdancer/bwa_new3_sort_35_2_3.bd_out.perf.txt', header=TRUE)
gasv <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_100bp_dip/final_predictions/GASVPro/BAMToGASV.gasvpro.in.ALL.MCMCThreshold.clusters.pruned.clusters.perf.txt', header=TRUE)
delly <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_100bp_dip/final_predictions/delly/delly_q0_c5_del_bwa_new3.perf9c902d5.txt', header=TRUE)
pindel <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_100bp_dip/final_predictions/pindel/bwa_new3.pindel.perf.txt', header=TRUE)

cloudbreak <- cloudbreak[seq(1,dim(cloudbreak)[1],by=4),]

perfsList <- list(breakdancer=breakdancer, gasv=gasv, delly=delly, pindel=pindel, cloudbreak=cloudbreak)
#pdf('~/Documents/gene_rearrange/svpipeline/CHR2SIM_ROC_NEW.pdf')
plotROC(perfsList, c("Breakdancer","GASVPro", "DELLY", "Pindel", "Cloudbreak"), totalDels, "diploid chr2 Simulated (30X)")
#dev.off()

trueDels <- import.bed('~/Documents/gene_rearrange/svpipeline/venter_chr2_100bp_dip/HuRef.homozygous_indels_inversion.061109.chr2.Deletions.sim_hap1.bed.gz', asRangedData=FALSE)
trueDelshap2 <- import.bed('~/Documents/gene_rearrange/svpipeline/venter_chr2_100bp_dip/HuRef.homozygous_indels_inversion.061109.chr2.Deletions.sim_hap2.bed.gz', asRangedData=FALSE)

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

cb <- processPredictions("Cloudbreak", '~/Documents/gene_rearrange/svpipeline/venter_chr2_100bp_dip/test_mufilt_hits2.txt', trueDels)

bd <- processPredictions('Breakdancer', '~/Documents/gene_rearrange/svpipeline/venter_chr2_100bp_dip/final_predictions/breakdancer/bwa_new3_sort_35_2_3.bd_out.gte44.hits.txt', trueDels)

gasv <- processPredictions('GASVPro', '~/Documents/gene_rearrange/svpipeline/venter_chr2_100bp_dip/final_predictions/GASVPro/BAMToGASV.gasvpro.in.ALL.MCMCThreshold.clusters.pruned.clusters.gte8.hits.txt', trueDels)

delly <- processPredictions('DELLY', '~/Documents/gene_rearrange/svpipeline/venter_chr2_100bp_dip/final_predictions/delly/bwa_new3.delly_q0_c5_del.gte4.hits.txt', trueDels)

pindel <- processPredictions('Pindel', '~/Documents/gene_rearrange/svpipeline/venter_chr2_100bp_dip/final_predictions/pindel/bwa_new3.pindel.gte21.hits.txt', trueDels)

allpreds <- rbind(as.data.frame(cb),as.data.frame(bd),as.data.frame(gasv),as.data.frame(delly),as.data.frame(pindel))

allpreds <- ddply(allpreds, .(tok), function(x) { cbind(x, list(discoveries=rep(dim(x)[1],dim(x)[1]))) })

table(trueDels$sizeclasses)
table(allpreds[,c('name','sizeclasses')])
table(allpreds[which(allpreds$discoveries == 1),c('name','sizeclasses')])
