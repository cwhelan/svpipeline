library(GenomicRanges)
library(rtracklayer)
library(plyr)
library(RColorBrewer)
library(VariantAnnotation)
library(brew)
library(ggplot2)

drawPerfLine <- function(perf, col, lineType, maxFP) {
  calls <- c(perf$Calls, min(perf$Calls) - min(perf$TP), 0)
  tp <- c(perf$TP, 0, 0)
  fp <- calls - tp
  lines(fp[fp <= maxFP], tp[fp <= maxFP], col=col, lty=lineType, lwd=3)
}

rocColors <- function(numLines) {
  if (numLines <= 8) { 
    perfCols <- brewer.pal(numLines, "Dark2")
  } else {
    perfCols <- rainbow(numLines)
  }
  return(perfCols)
}

rocLineTypes <- function(numLines) {
  seq(1,numLines)
}

getColLtyMapping <- function(lineNames) {
  myColors <- rocColors(length(lineNames))
  myLtys <- rocLineTypes(length(lineNames))
  colLtyMapping <- list()
  for (n in 1:length(lineNames)) {
    colLtyMapping[[lineNames[n]]] <- list(myColors[n], myLtys[n])
  }
  return(colLtyMapping)
}

getCols <- function(keys, colLtyMapping) {
  unlist(lapply(colLtyMapping[keys], function(x) {x[1]}))
}

getLtys <- function(keys, colLtyMapping) {
  unlist(lapply(colLtyMapping[keys], function(x) {x[2]}))
}

plotROC <- function(perfs, perfNames, totalDels, main, sim=TRUE, maxTP=NA, legendLoc="bottomright",legendCex=.9, 
                    colLtyMapping=getColLtyMapping(perfNames), ...) {
  maxFP <- totalDels
  if (is.na(maxTP)) {
    maxTP <- totalDels
  }
  plot(0, type="n", ylim=c(0, maxTP), xlim=c(0,maxFP), xlab=ifelse(sim, "False Positives", "Novel Predictions"), ylab="True Positives", main=main, ...)
  mapply(drawPerfLine, perfs, getCols(perfNames, colLtyMapping), getLtys(perfNames, colLtyMapping), MoreArgs=list(maxFP=maxFP))  
  if (! is.null(legendLoc)) {
    legend(legendLoc, legend=perfNames, col=getCols(perfNames, colLtyMapping), 
           lty=getLtys(perfNames, colLtyMapping), lwd=3, cex=legendCex)
  }
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

numPredictionsAtLevel <- function(perfTable, threshold) {
  return(perfTable[perfTable$Thresh == min(perfTable[perfTable$Thresh >= threshold,'Thresh']),'Calls'])
}

formatNumber <- function(d) {
  formatC(d,digits=3,format="fg",width=-1,drop0trailing=TRUE)
}

printPredTableValue <- function(t, x, y) {
  if (! y %in% colnames(t)) {
    return(0)
  }
  d <- t[x,y]
  if (d == max(t[,y])) {
    return(paste0("\\textbf{",formatNumber(d),"}"))
  } else {
    return(formatNumber(d))
  }
}

printPredTableValueFromList <- function(x, l) {
  d <- ifelse(is.na(l[x]), 0, l[x])
  if (d == max(l)) {
    return(paste0("\\textbf{",formatNumber(d),"}"))
  } else {
    return(formatNumber(d))
  }
}

printPredTableCell <- function(preds, ex, x, y) {
  return(paste0(printPredTableValue(preds, x, y), " (", printPredTableValue(ex, x, y), ")"))
}

# input files

# chr2

# deletions
cloudbreakGEMChr2DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/gem_e6_s2_del.perf.txt'
cloudbreakBWAChr2DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/bwa_del.perf.txt'
cloudbreakBWAXAChr2DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/bwa_xa2multi_del.perf.txt'
breakdancerChr2DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort_35_2_3_del.perf.txt'
gasvProChr2DeletionPerfFile <-  '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/gasvpro_gte40.perf.txt'
dellyChr2DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort.delly_q0_c5_gte40.perf.txt'
dellyBRChr2DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort.delly_q0_c5_br.perf.txt'
pindelChr2DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort_pindel_D_gte40.perf.txt'

# bwa mem meappings
breakdancerBWAMEMChr2DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/bwa_mem_human_b36_male_chr2_venterindels_c15_i100_s30_rl100_breakdancer.dels.perf.txt'

cloudbreakGEMChr2DeletionPredictionsFileFDR10 <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/gem_del_calls_gt2p29.bed'
cloudbreakBWAChr2DeletionPredictionsFileFDR10 <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/bwa_del_1p98.bed'

cloudbreakGEMChr2DeletionHitsFileFDR10 <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/gem_del_calls_gt2p29.hits.txt'
cloudbreakBWAChr2DeletionHitsFileFDR10 <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/bwa_dels_1p98.hits.txt'
breakdancerChr2DeletionHitsFileFDR10 <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort_35_2_3_del.gt4_gte40.hits.txt'
gasvProChr2DeletionHitsFileFDR10 <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/gasvpro_gte40_gt11.hits.txt'
dellyChr2DeletionHitsFileFDR10 <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort.delly_q0_c5_gte40_gt7.hits.txt'
dellyBRChr2DeletionHitsFileFDR10 <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort.delly_q0_c5_br_gt11.hits.txt'
pindelChr2DeletionHitsFileFDR10 <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort_pindel_D_gte40_gt13.hits.txt'

chr2MaxSensitivityDelCutoffs=list(cloudbreakGEM=0.85, breakdancer=2, GASVPro=4, DELLY=2, DELLYSR=2, Pindel=4, cloudbreakBWA=0.8)

cloudbreakGEMChr2DeletionPredictionsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/gem_del_calls_gt0p85.bed'
cloudbreakBWAChr2DeletionPredictionsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/bwa_dels_0p8.bed'
modilChr2DeletionPredictionsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/modil_del_gt1.bed'

cloudbreakGEMChr2DeletionHitsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/gem_del_calls_gt0p85.hits.txt'
cloudbreakBWAChr2DeletionHitsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/bwa_dels_0p8.hits.txt'
breakdancerChr2DeletionHitsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort_35_2_3_del.gt2_gte40.hits.txt'
gasvProChr2DeletionHitsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/gasvpro_gte40_gt4.hits.txt'
dellyChr2DeletionHitsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort.delly_q0_c5_gte40_gt2.hits.txt'
dellyBRChr2DeletionHitsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort.delly_q0_c5_br_gt2.hits.txt'
pindelChr2DeletionHitsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort_pindel_D_gte40_gt4.hits.txt'
modilChr2DeletionHitsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/modil_del_gt1.hits.txt'

cloudbreakGEMChr2Perf1ReportPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/gem_e6_s2_del_d1.perf.txt'
cloudbreakGEMChr2Perf10ReportPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/gem_e6_s2_del_d10.perf.txt'
cloudbreakGEMChr2Perf100ReportPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/gem_e6_s2_del_d100.perf.txt'
cloudbreakGEMChr2Perf1000ReportPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/gem_e6_s2_del_d1000.perf.txt'
cloudbreakGEMChr2Perf5000ReportPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/gem_e6_s2_del_d5000.perf.txt'
cloudbreakGEMChr2Perf50000ReportPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/gem_e6_s2_del_d50000.perf.txt'
cloudbreakGEMChr2Perf100000ReportPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/gem_e6_s2_del_d100000.perf.txt'

# insertions
cloudbreakGEMChr2InsertionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/gem_e6_s2_ins.perf.txt'
cloudbreakBWAChr2InsertionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/bwa_ins.perf.txt'
cloudbreakBWATweakChr2InsertionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/bwa_instweak.perf.txt'
cloudbreakBWAXATweakChr2InsertionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/bwa_xa2multi_instweak.perf.txt'
breakdancerChr2InsertionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort_35_2_3.bd_out.ins_gte40.perf.txt'
pindelChr2InsertionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/pindel_insertions.gte40.perf.txt'

# length adjusted insertions
cloudbreakBWAChr2LAInsertionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/bwa_ins.length_adjusted.perf.txt'
cloudbreakBWATweakChr2LAInsertionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/bwa_instweak.length_adjusted.perf.txt'
breakdancerChr2LAInsertionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/breakdancer_ins.length_adjusted.perf.txt'
pindelChr2LAInsertionPerfFile <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/pindel_ins.length_adjusted.perf.txt'

cloudbreakGEMChr2InsertionPredictionsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/gem_ins_calls_gt0p26.bed'
cloudbreakBWAChr2InsertionPredictionsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/bwa_ins_0p26.bed'
#cloudbreakBWATweakChr2InsertionPredictionsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/bwa_ins_tweak_0p26.bed'
cloudbreakBWATweakChr2InsertionPredictionsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/ins_test.bed'
modilChr2InsertionPredictionsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/modil_ins_gte1.bed'

chr2MaxSensitivityInsCutoffs=list(cloudbreakGEM=0.26, breakdancer=2, Pindel=4, MoDIL=1, cloudbreakBWA=0.26, cloudbreakBWATweak=0.26)

cloudbreakGEMChr2InsertionHitsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/gem_ins_calls_gt0p26.hits.txt'
cloudbreakBWAChr2InsertionHitsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/bwa_ins_0p26.hits.txt'
#cloudbreakBWATweakChr2InsertionHitsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/bwa_ins_tweak_0p26.hits.txt'
cloudbreakBWATweakChr2InsertionHitsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/ins_test.hits.txt'
breakdancerChr2InsertionHitsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort_35_2_3.bd_out.ins_gte40.gt2.hits.txt'
pindelChr2InsertionHitsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/pindel_insertions.gte40.gt4.hits.txt'
modilChr2InsertionHitsFileMaxSensitivity <- '~/Documents/gene_rearrange/svpipeline/venter_chr2_allindels_100bp_dip/modil_ins_gte1.hits.txt'

# NA18507

# Deletions
cloudbreakGEMNA18507DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/1KG_imprecise_sv_dels.perf.txt'
cloudbreakBWANA18507DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/bwa_del.perf.txt'
breakdancerNA18507DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/NA18507_35_2_3_singlecpu_gte40.new1KGCalls.perf.txt'
gasvProNA18507DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/gasvpro.gte40.new1KGCalls.perf.txt'
dellyNA18507DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/NA18507.delly_q0_c5_del.new1KGCalls.perf.txt'
dellyBRNA18507DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/NA18507.delly_q0_c5_del_br.new1KGCalls.perf.txt'
pindelNA18507DeletionPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/bwa_pindel_D_st40_gte40.new1KGCalls.perf.txt'

cloudbreakGEMNA18507DelPredictionsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/gem_del_calls_gt2p82.bed'
#cloudbreakBWANA18507DelPredictionsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/bwa_del_calls_gt2p44.bed'
cloudbreakBWANA18507DelPredictionsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/bwa_del_calls_sparseves_gt2p44.bed'

NA18507delCutoffs=list(cloudbreakGEM=2.82, breakdancer=4.9, GASVPro=13.6, DELLY=8.6, DELLYSR=13.6, Pindel=16.03, cloudbreakBWA=2.44)

cloudbreakGEMNA18507DeletionHitsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/gem_del_calls_gt2p82.new1KGCalls.hits.txt'
#cloudbreakBWANA18507DeletionHitsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/bwa_del_calls_gt2p44.hits.txt'
cloudbreakBWANA18507DeletionHitsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/bwa_del_calls_sparseves_gt2p44.hits.txt'
breakdancerNA18507DeletionHitsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/NA18507_35_2_3_singlecpu_gte40.gt4p9.new1KGCalls.hits.txt'
gasvProNA18507DeletionHitsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/gasvpro.gte40.gt13p6.new1KGCalls.hits.txt'
dellyNA18507DeletionHitsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/NA18507.delly_q0_c5_del_gte40.gt8p6.new1KGCalls.hits.txt'
dellyBRNA18507DeletionHitsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/NA18507.delly_q0_c5_br_del_gte40.gt13p6.new1KGCalls.hits.txt'
pindelNA18507DeletionHitsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/bwa_pindel_D_st40_gte40.gt16.new1KGCalls.hits.txt'

# insertions
cloudbreakGEMNA18507InsertionsPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/gem_gmm_300reducers_ins_perf.txt'
cloudbreakBWANA18507InsertionsPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/bwa_ins.perf.txt'
cloudbreakBWATweakNA18507InsertionsPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/bwa_instweak.perf.txt'
breakdancerNA18507InsertionsPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/NA18507_35_2_3_insertions.gte40.perf.txt'
pindelNA18507InsertionsPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/bwa_pindel_ins_gte40.perf.txt'

# length adjusted insertions
cloudbreakBWANA18507LAInsertionsPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/bwa_ins.length_adjusted.perf.txt'
cloudbreakBWATweakNA18507LAInsertionsPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/bwa_instweak.length_adjusted.perf.txt'
breakdancerNA18507LAInsertionsPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/breakdancer_ins.length_adjusted.perf.txt'
pindelNA18507LAInsertionsPerfFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/pindel_ins.length_adjusted.perf.txt'

NA18507insertionCutoffs=list(cloudbreakGEM=1.08, breakdancer=2, Pindel=4, cloudbreakBWATweak=1.16)

cloudbreakGEMNA18507InsertionPredictionsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/gem_ins_calls_gt1p08.hits.txt'
# cloudbreakBWATweakNA18507InsertionPredictionsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/bwa_instweak_calls_gt1p16.bed'
cloudbreakBWATweakNA18507InsertionPredictionsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/bwa_instweak_calls_gt1p16.bed'

cloudbreakGEMNA18507InsertionHitsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/gem_ins_calls_gt1p08.hits.txt'
breakdancerNA18507InsertionHitsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/NA18507_35_2_3.bd_out.ins.gte2.hits.txt'
pindelNA18507InsertionHitsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/pindel_insertions.gte4.hits.txt'
#cloudbreakBWATweakNA18507InsertionHitsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/bwa_instweak_calls_gt1p16.hits.txt'
cloudbreakBWATweakNA18507InsertionHitsFile <- '~/Documents/gene_rearrange/svpipeline/NA18507/bwa_instweak_calls_sparseve_gt1p16.hits.txt'

genotypingAlphaCutoff <- .35

# reference files
repMask <- import.bed('~/genomes/1kg/hg18/repeatmasker_b36_merged.bed.gz', asRangedData=FALSE)

# output files

chr2DeletionsROCOutputFile <- '~/Documents/svpipeline/figures/CHR2SIM_DELS_ROC.pdf'
chr2DeletionsBWAMEMROCOutputFile <- '~/Documents/svpipeline/figures/CHR2SIM_BWAMEM_DELS_ROC.pdf'
chr2InsertionsROCOutputFile <- '~/Documents/svpipeline/figures/CHR2SIM_INS_ROC.pdf'
chr2LAInsertionsROCOutputFile <- '~/Documents/svpipeline/figures/CHR2SIM_INS_ROC_LA.pdf'
chr2ROCsOutputFile <- '~/Documents/svpipeline/figures/CHR2SIM_ROC_COMBINED_ROCS.pdf'
chr2ROCsPosterOutputFile <- '~/Documents/svpipeline/figures/CHR2SIM_ROC_COMBINED_ROCS_POSTER.pdf'
chr2ROCsMultipleMappingsFile <- '~/Documents/svpipeline/figures/CHR2SIM_ROCS_MULTIPLE_MAPPINGS.pdf'
NA18507DeletionsROCOutputFile <- '~/Documents/svpipeline/figures/NA18507_DELS_ROC.pdf'
NA18507InsertionsROCOutputFile <- '~/Documents/svpipeline/figures/NA18507_INS_ROC.pdf'
NA18507LAInsertionsROCOutputFile <- '~/Documents/svpipeline/figures/NA18507_INS_ROC_LA.pdf'
NA18507ROCsOutputFile <- '~/Documents/svpipeline/figures/NA18507_COMBINED_ROCS.pdf'
NA18507ROCsPosterOutputFile <- '~/Documents/svpipeline/figures/NA18507_COMBINED_ROCS_POSTER.pdf'
chr2DelReportsComparisonROCOutputFile <- '~/Documents/svpipeline/figures/CHR2_SIM_DEL_REPORTSCOMPARISON_ROC.pdf'
breakpointResolutionOutputFile <- '~/Documents/svpipeline/figures/breakpoint_resolution.pdf'
runtimeByNumberOfNodesOutputFile <- '~/Documents/svpipeline/figures/runtimeByNodes.pdf'

#chr2 100bp diploid deletions
totalDelsChr2 <- 400
cloudbreakGEMChr2DeletionPerf <- read.table(cloudbreakGEMChr2DeletionPerfFile, header=TRUE, sep="\t")
cloudbreakBWAChr2DeletionPerf <- read.table(cloudbreakBWAChr2DeletionPerfFile, header=TRUE, sep="\t")
cloudbreakBWAXAChr2DeletionPerf <- read.table(cloudbreakBWAXAChr2DeletionPerfFile, header=TRUE, sep="\t")
breakdancerChr2DeletionPerf <- read.table(breakdancerChr2DeletionPerfFile, header=TRUE)
breakdancerBWAMEMChr2DeletionPerf <- read.table(breakdancerBWAMEMChr2DeletionPerfFile, header=TRUE)
gasvChr2DeletionPerf <- read.table(gasvProChr2DeletionPerfFile, header=TRUE)
dellyChr2DeletionPerf <- read.table(dellyChr2DeletionPerfFile, header=TRUE)
dellyBRChr2DeletionPerf <- read.table(dellyBRChr2DeletionPerfFile, header=TRUE)
pindelChr2DeletionPerf <- read.table(pindelChr2DeletionPerfFile, header=TRUE)

perfsListDelsChr2 <- list(cloudbreak=cloudbreakBWAChr2DeletionPerf, breakdancer=breakdancerChr2DeletionPerf, pindel=pindelChr2DeletionPerf, gasv=gasvChr2DeletionPerf, delly=dellyChr2DeletionPerf, dellyBR=dellyBRChr2DeletionPerf)

perfsListDelsMultimapChr2 <- list(breakdancer=breakdancerChr2DeletionPerf, pindel=pindelChr2DeletionPerf, gasv=gasvChr2DeletionPerf, delly=dellyChr2DeletionPerf, dellyBR=dellyBRChr2DeletionPerf, cloudbreakBWA=cloudbreakBWAChr2DeletionPerf, cloudbreakBWAM=cloudbreakBWAXAChr2DeletionPerf, cloudbreakGEM=cloudbreakGEMChr2DeletionPerf)

myColLtyMapping <- getColLtyMapping(c("Cloudbreak", "Breakdancer","Pindel", "GASVPro", "DELLY-RP", "DELLY-SR", "Cloudbreak-BWA-PE-MM", "Cloudbreak-GEM-SE-MM", "Breakdancer-BWAMEM"))

pdf(chr2DeletionsROCOutputFile, width=10)
par(xpd=T, mar=par()$mar+c(0,0,0,7))
plotROC(perfsListDelsChr2, 
        c("Cloudbreak", "Breakdancer","Pindel", "GASVPro", "DELLY-RP", "DELLY-SR"), 
        totalDelsChr2, "Deletions in Venter diploid chr2 simulation",legendLoc=xy.coords(425,200), maxTP=350, 
        colLtyMapping=myColLtyMapping)
dev.off()

perfsListDelsChr2 <- list(cloudbreak=cloudbreakBWAChr2DeletionPerf, breakdancer=breakdancerChr2DeletionPerf, pindel=pindelChr2DeletionPerf, gasv=gasvChr2DeletionPerf, delly=dellyChr2DeletionPerf, dellyBR=dellyBRChr2DeletionPerf)
cairo_ps("~/Documents/svpipeline/figures/chr2_del_roc_cairo.ps", width=10, height=7)
par(xpd=T, mar=par()$mar+c(.5,.5,.5,13))
plotROC(perfsListDelsChr2, 
        c("Cloudbreak", "Breakdancer","Pindel", "GASVPro", "DELLY-RP", "DELLY-SR"), 
        totalDelsChr2, "Deletions",legendLoc=xy.coords(425,250), maxTP=350,
        cex.main=1.5, cex.axis=1.5, cex.lab=1.5, legendCex=1.5, colLtyMapping=myColLtyMapping)
dev.off()

# reports comparison
cloudbreakGEMChr2DelPerf1Report <- read.table(cloudbreakGEMChr2Perf1ReportPerfFile, header=TRUE, sep="\t")
cloudbreakGEMChr2DelPerf10Report <- read.table(cloudbreakGEMChr2Perf10ReportPerfFile, header=TRUE, sep="\t")
cloudbreakGEMChr2DelPerf100Report <- read.table(cloudbreakGEMChr2Perf100ReportPerfFile, header=TRUE, sep="\t")
cloudbreakGEMChr2DelPerf1000Report <- read.table(cloudbreakGEMChr2Perf1000ReportPerfFile, header=TRUE, sep="\t")
cloudbreakGEMChr2DelPerf5000Report <- read.table(cloudbreakGEMChr2Perf5000ReportPerfFile, header=TRUE, sep="\t")
cloudbreakGEMChr2DelPerf50000Report <- read.table(cloudbreakGEMChr2Perf50000ReportPerfFile, header=TRUE, sep="\t")
cloudbreakGEMChr2DelPerf100000Report <- read.table(cloudbreakGEMChr2Perf100000ReportPerfFile, header=TRUE, sep="\t")

perfsListDelsChr2CompareReports <- list(cloudbreak1=cloudbreakGEMChr2DelPerf1Report, 
                                        cloudbreak10=cloudbreakGEMChr2DelPerf10Report, 
                                        cloudbreak100=cloudbreakGEMChr2DelPerf100Report, 
                                        cloudbreak1000=cloudbreakGEMChr2DelPerf1000Report, 
                                        cloudbreak5000=cloudbreakGEMChr2DelPerf5000Report, 
                                        cloudbreak50000=cloudbreakGEMChr2DelPerf50000Report, 
                                        cloudbreak100000=cloudbreakGEMChr2DelPerf100000Report)
                                        
pdf(chr2DelReportsComparisonROCOutputFile, width=10)
par(xpd=T, mar=par()$mar+c(0,0,0,7))
plotROC(perfsListDelsChr2CompareReports, 
        c("1 Mapping", "10 Mappings","100 Mappings", "1000 Mappings", "5000 Mappings", "50000 Mappings", "100000 Mappings"), 
        totalDelsChr2, "Performance on Venter diploid chr2 simulation by number of mappings considered",legendLoc=xy.coords(425,200), maxTP=350)
dev.off()

# breakdancer BWA vs BWA-MEM
perfsListDelsChr2BWAMEM <- list(breakdancer=breakdancerChr2DeletionPerf, breakdancerBWAMEM=breakdancerBWAMEMChr2DeletionPerf)
pdf(chr2DeletionsBWAMEMROCOutputFile, width=10)
par(xpd=T, mar=par()$mar+c(0,0,0,10))
plotROC(perfsListDelsChr2BWAMEM, 
        c("Breakdancer","Breakdancer-BWAMEM"), 
        totalDelsChr2, "Deletions in Venter diploid chr2 simulation",legendLoc=xy.coords(425,200), maxTP=350, 
        colLtyMapping=myColLtyMapping)
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

# tabulate info about true dels
trueDelsGte40 <- trueDels[end(trueDels) - start(trueDels) >= 40]
trueDelSizes <- table(trueDelsGte40$sizeclasses)
trueDelHaps <- table(trueDelsGte40$haps)
trueDelRepmask <- table(trueDelsGte40$repMask)

# test with FDR 10

# load true positive predictions
cbFDR10 <- processDeletionPredictions("Cloudbreak", cloudbreakBWAChr2DeletionHitsFileFDR10, trueDels)
bdFDR10 <- processDeletionPredictions('Breakdancer', breakdancerChr2DeletionHitsFileFDR10, trueDels)
gasvFDR10 <- processDeletionPredictions('GASVPro', gasvProChr2DeletionHitsFileFDR10, trueDels)
dellyFDR10 <- processDeletionPredictions('DELLY-RP', dellyChr2DeletionHitsFileFDR10, trueDels)
dellyBRFDR10 <- processDeletionPredictions('DELLY-SR', dellyBRChr2DeletionHitsFileFDR10, trueDels)
pindelFDR10 <- processDeletionPredictions('Pindel', pindelChr2DeletionHitsFileFDR10, trueDels)
allpredsFDR10 <- rbind(as.data.frame(cbFDR10),as.data.frame(bdFDR10),as.data.frame(gasvFDR10),as.data.frame(dellyFDR10),as.data.frame(dellyBRFDR10),as.data.frame(pindelFDR10))
allpredsFDR10 <- ddply(allpredsFDR10, .(tok), function(x) { cbind(x, list(discoveries=rep(dim(x)[1],dim(x)[1]))) })

# tabulate against true deletion data
predsBySizeFDR10 <- table(allpredsFDR10[,c('name','sizeclasses')])
exBySizeFDR10 <- table(allpredsFDR10[which(allpredsFDR10$discoveries == 1),c('name','sizeclasses')])

predsByRepMaskFDR10 <- table(allpredsFDR10[,c('name','repMask')])
exPredsByRepMaskFDR10 <- table(allpredsFDR10[which(allpredsFDR10$discoveries == 1),c('name','repMask')])

predsByHapFDR10 <- table(allpredsFDR10[,c('name','haps')])
exPredsByHapFDR10 <- table(allpredsFDR10[which(allpredsFDR10$discoveries == 1),c('name','haps')])

predsBySegDupFDR10 <- table(allpredsFDR10[,c('name','segdup')])
exPredsBySegDupFDR10 <- table(allpredsFDR10[which(allpredsFDR10$discoveries == 1),c('name','segdup')])

# genotyping with extra data
extraDataFDR10 <- read.table(cloudbreakGEMChr2DeletionPredictionsFileFDR10, skip=1)
names(extraDataFDR10) <- c('predchrom', 'predstart', 'predend', 'peaknum', 'maxscore', 'svtype', 'avgMu', 'minMu', 'maxMu', 'avgW', 'minW', 'maxW')

cbpredsWithExtraDataFDR10 <- merge(as.data.frame(cbFDR10), extraDataFDR10)
chr2HapCMFDR10 <- table(cbpredsWithExtraDataFDR10$haps, cbpredsWithExtraDataFDR10$avgW < genotypingAlphaCutoff)

# test with Max sensitivity

# load true positive predictions
cbMaxSensitivity <- processDeletionPredictions("Cloudbreak", cloudbreakBWAChr2DeletionHitsFileMaxSensitivity, trueDels)
bdMaxSensitivity <- processDeletionPredictions('Breakdancer', breakdancerChr2DeletionHitsFileMaxSensitivity, trueDels)
gasvMaxSensitivity <- processDeletionPredictions('GASVPro', gasvProChr2DeletionHitsFileMaxSensitivity, trueDels)
dellyMaxSensitivity <- processDeletionPredictions('DELLY-RP', dellyChr2DeletionHitsFileMaxSensitivity, trueDels)
dellyBRMaxSensitivity <- processDeletionPredictions('DELLY-SR', dellyBRChr2DeletionHitsFileMaxSensitivity, trueDels)
pindelMaxSensitivity <- processDeletionPredictions('Pindel', pindelChr2DeletionHitsFileMaxSensitivity, trueDels)
modilMaxSensitivity <- processDeletionPredictions('MoDIL', modilChr2DeletionHitsFileMaxSensitivity, trueDels)
allpredsMaxSensitivity <- rbind(as.data.frame(cbMaxSensitivity),as.data.frame(bdMaxSensitivity),as.data.frame(gasvMaxSensitivity),as.data.frame(dellyMaxSensitivity),as.data.frame(dellyBRMaxSensitivity),as.data.frame(pindelMaxSensitivity),as.data.frame(modilMaxSensitivity))
allpredsMaxSensitivity <- ddply(allpredsMaxSensitivity, .(tok), function(x) { cbind(x, list(discoveries=rep(dim(x)[1],dim(x)[1]))) })

modilChr2DeletionPredictionsMaxSensitivity <- read.table(modilChr2DeletionPredictionsFileMaxSensitivity)

totalCloudbreakChr2MaxSensitivityDelPredictions <- numPredictionsAtLevel(cloudbreakBWAChr2DeletionPerf, chr2MaxSensitivityDelCutoffs['cloudbreakBWA'])
totalBreakdancerchr2MaxSensitivityDelPredictions <- numPredictionsAtLevel(breakdancerChr2DeletionPerf, chr2MaxSensitivityDelCutoffs['breakdancer'])
totalGASVProchr2MaxSensitivityDelPredictions <- numPredictionsAtLevel(gasvChr2DeletionPerf, chr2MaxSensitivityDelCutoffs['GASVPro'])
totalDELLYchr2MaxSensitivityDelPredictions <- numPredictionsAtLevel(dellyChr2DeletionPerf, chr2MaxSensitivityDelCutoffs['DELLY'])
totalDELLYBRchr2MaxSensitivityDelPredictions <- numPredictionsAtLevel(dellyBRChr2DeletionPerf, chr2MaxSensitivityDelCutoffs['DELLYSR'])
totalPindelchr2MaxSensitivityDelPredictions <- numPredictionsAtLevel(pindelChr2DeletionPerf, chr2MaxSensitivityDelCutoffs['Pindel'])
totalMoDILchr2MaxSensitivityDelPredictions <- dim(modilChr2DeletionPredictionsMaxSensitivity)[1]

totalchr2MaxSensitivityDelPredictions <- list(Cloudbreak=totalCloudbreakChr2MaxSensitivityDelPredictions, 
                                              Breakdancer=totalBreakdancerchr2MaxSensitivityDelPredictions, 
                                              GASVPro=totalGASVProchr2MaxSensitivityDelPredictions, 
                                              DELLY=totalDELLYchr2MaxSensitivityDelPredictions, 
                                              DELLYSR=totalDELLYBRchr2MaxSensitivityDelPredictions, 
                                              Pindel=totalPindelchr2MaxSensitivityDelPredictions, 
                                              MoDIL=totalMoDILchr2MaxSensitivityDelPredictions)
totalchr2MaxSensitivityDelHits <- list(Cloudbreak=dim(cbMaxSensitivity)[1], Breakdancer=dim(bdMaxSensitivity)[1], GASVPro=dim(gasvMaxSensitivity)[1], DELLY=dim(dellyMaxSensitivity)[1], DELLYSR=dim(dellyBRMaxSensitivity)[1], Pindel=dim(pindelMaxSensitivity)[1], MoDIL=dim(modilMaxSensitivity)[1])

chr2MaxSensitivityDelPrecision <- unlist(totalchr2MaxSensitivityDelHits) / unlist(totalchr2MaxSensitivityDelPredictions)
chr2MaxSensitivityDelRecall <- unlist(totalchr2MaxSensitivityDelHits) / length(trueDelsGte40)

# tabulate against true deletion data
predsBySizeMaxSensitivity <- table(allpredsMaxSensitivity[,c('name','sizeclasses')])
exBySizeMaxSensitivity <- table(allpredsMaxSensitivity[which(allpredsMaxSensitivity$discoveries == 1),c('name','sizeclasses')])

predsByRepMaskMaxSensitivity <- table(allpredsMaxSensitivity[,c('name','repMask')])
exPredsByRepMaskMaxSensitivity <- table(allpredsMaxSensitivity[which(allpredsMaxSensitivity$discoveries == 1),c('name','repMask')])

predsByHapMaxSensitivity <- table(allpredsMaxSensitivity[,c('name','haps')])
exPredsByHapMaxSensitivity <- table(allpredsMaxSensitivity[which(allpredsMaxSensitivity$discoveries == 1),c('name','haps')])

predsBySegDupMaxSensitivity <- table(allpredsMaxSensitivity[,c('name','segdup')])
exPredsBySegDupMaxSensitivity <- table(allpredsMaxSensitivity[which(allpredsMaxSensitivity$discoveries == 1),c('name','segdup')])

# genotyping with extra data
extraDataMaxSensitivity <- read.table(cloudbreakBWAChr2DeletionPredictionsFileMaxSensitivity, skip=1)
names(extraDataMaxSensitivity) <- c('predchrom', 'predstart', 'predend', 'peaknum', 'maxscore', 'svtype', 'avgMu', 'minMu', 'maxMu', 'avgW', 'minW', 'maxW')

cbpredsWithExtraDataMaxSensitivity <- merge(as.data.frame(cbMaxSensitivity), extraDataMaxSensitivity)

chr2HapCMMaxSensitivity <- table(cbpredsWithExtraDataMaxSensitivity$haps, cbpredsWithExtraDataMaxSensitivity$avgW < genotypingAlphaCutoff)
chr2DeletionHapAccuracy <- (chr2HapCMMaxSensitivity['1','FALSE'] + chr2HapCMMaxSensitivity['2','TRUE']) / sum(chr2HapCMMaxSensitivity)

#chr2 100bp diploid insertions
totalInsertionsChr2 <- 140
cloudbreakGEMChr2InsertionPerf <- read.table(cloudbreakGEMChr2InsertionPerfFile, header=TRUE, sep="\t")
cloudbreakGEMChr2InsertionPerfSmoothed <- cloudbreakGEMChr2InsertionPerf[seq(1,dim(cloudbreakGEMChr2InsertionPerf)[1],by=10),]

cloudbreakBWAChr2InsertionPerf <- read.table(cloudbreakBWAChr2InsertionPerfFile, header=TRUE, sep="\t")
cloudbreakBWAChr2InsertionPerfSmoothed <- cloudbreakBWAChr2InsertionPerf[seq(1,dim(cloudbreakBWAChr2InsertionPerf)[1],by=10),]

cloudbreakBWATweakChr2InsertionPerf <- read.table(cloudbreakBWATweakChr2InsertionPerfFile, header=TRUE, sep="\t")
cloudbreakBWATweakChr2InsertionPerfSmoothed <- cloudbreakBWATweakChr2InsertionPerf[seq(1,dim(cloudbreakBWATweakChr2InsertionPerf)[1],by=11),]

cloudbreakBWAXATweakChr2InsertionPerf <- read.table(cloudbreakBWAXATweakChr2InsertionPerfFile, header=TRUE, sep="\t")
cloudbreakBWAXATweakChr2InsertionPerfSmoothed <- cloudbreakBWAXATweakChr2InsertionPerf[seq(1,dim(cloudbreakBWAXATweakChr2InsertionPerf)[1],by=13),]

breakdancerChr2InsertionPerf <- read.table(breakdancerChr2InsertionPerfFile, header=TRUE)
pindelChr2InsertionPerf <- read.table(pindelChr2InsertionPerfFile, header=TRUE)

perfsListInsertionsChr2 <- list(cloudbreak=cloudbreakBWATweakChr2InsertionPerfSmoothed, breakdancer=breakdancerChr2InsertionPerf, pindel=pindelChr2InsertionPerf)

perfsListInsertionsMultimapChr2 <- list(breakdancer=breakdancerChr2InsertionPerf, pindel=pindelChr2InsertionPerf, 
                                        cloudbreakBWA=cloudbreakBWATweakChr2InsertionPerfSmoothed, cloudbreakBWAM=cloudbreakBWAXATweakChr2InsertionPerfSmoothed, cloudbreakGEM=cloudbreakGEMChr2InsertionPerfSmoothed)

pdf(chr2InsertionsROCOutputFile, width=10)
par(xpd=T, mar=par()$mar+c(0,0,0,8))
plotROC(perfsListInsertionsChr2, 
        c("Cloudbreak", "Breakdancer","Pindel"), 
        totalInsertionsChr2, "Insertions in Venter diploid chr2 simulation",
        legendLoc=xy.coords(totalInsertionsChr2 * 1.05,50), maxTP=totalInsertionsChr2, colLtyMapping=myColLtyMapping)
dev.off()

cloudbreakBWAChr2LAInsertionPerf <- read.table(cloudbreakBWAChr2LAInsertionPerfFile, header=TRUE, sep="\t")
cloudbreakBWAChr2LAInsertionPerfSmoothed <- cloudbreakBWAChr2LAInsertionPerf[seq(1,dim(cloudbreakBWAChr2LAInsertionPerf)[1],by=10),]

cloudbreakBWATweakChr2LAInsertionPerf <- read.table(cloudbreakBWATweakChr2LAInsertionPerfFile, header=TRUE, sep="\t")
cloudbreakBWATweakChr2LAInsertionPerfSmoothed <- cloudbreakBWATweakChr2LAInsertionPerf[seq(1,dim(cloudbreakBWATweakChr2LAInsertionPerf)[1],by=10),]

breakdancerChr2LAInsertionPerf <- read.table(breakdancerChr2LAInsertionPerfFile, header=TRUE)
pindelChr2LAInsertionPerf <- read.table(pindelChr2LAInsertionPerfFile, header=TRUE)

perfsListInsertionsChr2LA <- list(breakdancer=breakdancerChr2LAInsertionPerf, pindel=pindelChr2LAInsertionPerf, cloudbreakBWA=cloudbreakBWAChr2LAInsertionPerfSmoothed, cloudbreakBWATweak=cloudbreakBWATweakChr2LAInsertionPerfSmoothed)
pdf(chr2LAInsertionsROCOutputFile, width=10)
par(xpd=T, mar=par()$mar+c(0,0,0,8))
plotROC(perfsListInsertionsChr2LA, 
        c("Breakdancer","Pindel", "Cloudbreak-BWA", "Cloudbreak-BWA-2"), 
        totalInsertionsChr2, "Length-Adjusted Insertions in Venter diploid chr2 simulation",
        legendLoc=xy.coords(totalInsertionsChr2 * 1.05,50), maxTP=totalInsertionsChr2)
dev.off()

perfsListInsertionsChr2 <- list(cloudbreak=cloudbreakBWATweakChr2InsertionPerfSmoothed, breakdancer=breakdancerChr2InsertionPerf, pindel=pindelChr2InsertionPerf)
cairo_ps("~/Documents/svpipeline/figures/chr2_ins_roc_cairo.ps", width=10, height=7)
par(xpd=T, mar=par()$mar+c(.5,.5,.5,13))
plotROC(perfsListInsertionsChr2, 
        c("Cloudbreak", "Breakdancer","Pindel"), 
        totalInsertionsChr2, "Insertions",legendLoc=xy.coords(totalInsertionsChr2 * 1.05,50), maxTP=totalInsertionsChr2,
        cex.main=1.5,cex.axis=1.5,cex.lab=1.5,legendCex=1.5, colLtyMapping=myColLtyMapping)
dev.off()

pdf(chr2ROCsOutputFile, width=15)
par(xpd=T, mfrow=(c(1,2)), oma=par()$oma+c(0,0,0,5))
plotROC(perfsListDelsChr2, 
        c("Cloudbreak", "Breakdancer","Pindel", "GASVPro", "DELLY-RP", "DELLY-SR"), 
        totalDelsChr2, "Deletions in Venter diploid chr2 simulation",legendLoc=NULL, maxTP=350, colLtyMapping=myColLtyMapping)
par(mar=par()$mar+c(0,0,0,6))
plotROC(perfsListInsertionsChr2, 
        c("Cloudbreak", "Breakdancer","Pindel"), 
        totalInsertionsChr2, "Insertions in Venter diploid chr2 simulation",legendLoc=NA, maxTP=totalInsertionsChr2)
perfNames <- c("Cloudbreak", "Breakdancer","Pindel", "GASVPro", "DELLY-RP", "DELLY-SR")
legend(xy.coords(85,50), legend=perfNames, col=getCols(perfNames, myColLtyMapping), lty=getLtys(perfNames, myColLtyMapping), lwd=3, cex=.9)
dev.off()

pdf(chr2ROCsPosterOutputFile, width=18)
par(xpd=T, mfrow=(c(1,2)), oma=par()$oma+c(0,0,0,3), cex=1.8)
plotROC(perfsListDelsChr2, 
        c("Cloudbreak", "Breakdancer","Pindel", "GASVPro", "DELLY-RP", "DELLY-SR"), 
        totalDelsChr2, "Deletions - Chr2 Simulation",legendLoc=NULL, maxTP=350, colLtyMapping=myColLtyMapping)
par(mar=par()$mar+c(0,0,0,8))
plotROC(perfsListInsertionsChr2, 
        c("Cloudbreak", "Breakdancer","Pindel"), 
        totalInsertionsChr2, "Insertions - Chr2 Simulation",legendLoc=NA, maxTP=totalInsertionsChr2, colLtyMapping=myColLtyMapping)
perfNames <- c("Cloudbreak", "Breakdancer","Pindel", "GASVPro", "DELLY-RP", "DELLY-SR")       
legend(xy.coords(95,70), legend=perfNames, col=getCols(perfNames, myColLtyMapping), lty=getLtys(perfNames, myColLtyMapping), lwd=3)
dev.off()

pdf(chr2ROCsMultipleMappingsFile, width=20)
par(xpd=T, mfrow=(c(1,2)), oma=par()$oma+c(0,0,0,6))
plotROC(perfsListDelsMultimapChr2, 
        c("Breakdancer","Pindel", "GASVPro", "DELLY-RP", "DELLY-SR", 
          "Cloudbreak", "Cloudbreak-BWA-PE-MM", "Cloudbreak-GEM-SE-MM"), 
        totalDelsChr2, "Deletions",legendLoc=NULL, maxTP=350, colLtyMapping=myColLtyMapping)
par(mar=par()$mar+c(0,0,0,12))
plotROC(perfsListInsertionsMultimapChr2, 
        c("Breakdancer","Pindel", "Cloudbreak", "Cloudbreak-BWA-PE-MM", "Cloudbreak-GEM-SE-MM"), 
        totalInsertionsChr2, "Insertions",legendLoc=NA,  maxTP=totalInsertionsChr2, colLtyMapping=myColLtyMapping)
perfNames <- c("Breakdancer","Pindel", "GASVPro", "DELLY-RP", "DELLY-SR", "Cloudbreak", "Cloudbreak-BWA-PE-MM", "Cloudbreak-GEM-SE-MM")
legend(xy.coords(totalInsertionsChr2 * 1.1,100), legend=perfNames, col=getCols(perfNames, myColLtyMapping), lty=getLtys(perfNames, myColLtyMapping), lwd=3, cex=.9)
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

trueInsertionsGte40 <- trueInsertions[mcols(trueInsertions)$length >= 40]
trueInsertionsSizes <- table(trueInsertionsGte40$sizeclasses)
trueInsertionsHaps <- table(trueInsertionsGte40$haps)
trueInsertionsRepmask <- table(trueInsertionsGte40$repMask)

# load true positive predictions
cbInsertionsMaxSensitivity <- processInsertionPredictions("Cloudbreak", cloudbreakBWAChr2InsertionHitsFileMaxSensitivity, trueInsertions)
bdInsertionsMaxSensitivity <- processInsertionPredictions('Breakdancer', breakdancerChr2InsertionHitsFileMaxSensitivity, trueInsertions)
pindelInsertionsMaxSensitivity <- processInsertionPredictions('Pindel', pindelChr2InsertionHitsFileMaxSensitivity, trueInsertions)
modilInsertionsMaxSensitivity <- processInsertionPredictions('MoDIL', modilChr2InsertionHitsFileMaxSensitivity, trueInsertions)
allpredsInsertionsMaxSensitivity <- rbind(as.data.frame(cbInsertionsMaxSensitivity),as.data.frame(bdInsertionsMaxSensitivity),as.data.frame(pindelInsertionsMaxSensitivity),as.data.frame(modilInsertionsMaxSensitivity))
allpredsInsertionsMaxSensitivity <- ddply(allpredsInsertionsMaxSensitivity, .(tok), function(x) { cbind(x, list(discoveries=rep(dim(x)[1],dim(x)[1]))) })

modilChr2InsertionPredictionsMaxSensitivity <- read.table(modilChr2InsertionPredictionsFileMaxSensitivity)

totalCloudbreakChr2MaxSensitivityInsPredictions <- numPredictionsAtLevel(cloudbreakBWAChr2InsertionPerf, chr2MaxSensitivityInsCutoffs['cloudbreakBWA'])
totalBreakdancerchr2MaxSensitivityInsPredictions <- numPredictionsAtLevel(breakdancerChr2InsertionPerf, chr2MaxSensitivityInsCutoffs['breakdancer'])
totalPindelchr2MaxSensitivityInsPredictions <- numPredictionsAtLevel(pindelChr2InsertionPerf, chr2MaxSensitivityInsCutoffs['Pindel'])
totalMoDILchr2MaxSensitivityInsPredictions <- dim(modilChr2InsertionPredictionsMaxSensitivity)[1]

totalchr2MaxSensitivityInsPredictions <- list(Cloudbreak=totalCloudbreakChr2MaxSensitivityInsPredictions, Breakdancer=totalBreakdancerchr2MaxSensitivityInsPredictions, Pindel=totalPindelchr2MaxSensitivityInsPredictions, MoDIL=totalMoDILchr2MaxSensitivityInsPredictions)
totalchr2MaxSensitivityInsHits <- list(Cloudbreak=dim(cbInsertionsMaxSensitivity)[1], Breakdancer=dim(bdInsertionsMaxSensitivity)[1], Pindel=dim(pindelInsertionsMaxSensitivity)[1], MoDIL=dim(modilInsertionsMaxSensitivity)[1])

chr2MaxSensitivityInsPrecision <- unlist(totalchr2MaxSensitivityInsHits) / unlist(totalchr2MaxSensitivityInsPredictions)
chr2MaxSensitivityInsRecall <- unlist(totalchr2MaxSensitivityInsHits) / length(trueInsertionsGte40)

# tabulate against true insertion data
insertionPredsBySizeMaxSensitivity <- table(allpredsInsertionsMaxSensitivity[,c('name','sizeclasses')])
insertionExBySizeMaxSensitivity <- table(allpredsInsertionsMaxSensitivity[which(allpredsInsertionsMaxSensitivity$discoveries == 1),c('name','sizeclasses')])

insertionPredsByRepMaskMaxSensitivity <- table(allpredsInsertionsMaxSensitivity[,c('name','repMask')])
insertionExPredsByRepMaskMaxSensitivity <- table(allpredsInsertionsMaxSensitivity[which(allpredsInsertionsMaxSensitivity$discoveries == 1),c('name','repMask')])

insertionPredsByHapMaxSensitivity <- table(allpredsInsertionsMaxSensitivity[,c('name','haps')])
insertionExPredsByHapMaxSensitivity <- table(allpredsInsertionsMaxSensitivity[which(allpredsInsertionsMaxSensitivity$discoveries == 1),c('name','haps')])

insertionPredsBySegDupMaxSensitivity <- table(allpredsInsertionsMaxSensitivity[,c('name','segdup')])
insertionExPredsBySegDupMaxSensitivity <- table(allpredsInsertionsMaxSensitivity[which(allpredsInsertionsMaxSensitivity$discoveries == 1),c('name','segdup')])

# genotyping insertions
extraInsertionDataMaxSensitivity <- read.table(cloudbreakBWAChr2InsertionPredictionsFileMaxSensitivity, skip=1)
names(extraInsertionDataMaxSensitivity) <- c('predchrom', 'predstart', 'predend', 'peaknum', 'maxscore', 'svtype', 'avgMu', 'minMu', 'maxMu', 'avgW', 'minW', 'maxW')

cbInsertionPredsWithExtraDataMaxSensitivity <- merge(as.data.frame(cbInsertionsMaxSensitivity), extraInsertionDataMaxSensitivity, by=c("predchrom", "predstart"))

insertionHapCMMaxSensitivity <- table(cbInsertionPredsWithExtraDataMaxSensitivity$haps, cbInsertionPredsWithExtraDataMaxSensitivity$minW < genotypingAlphaCutoff)
print((insertionHapCMMaxSensitivity[1,1] + insertionHapCMMaxSensitivity[2,2]) / sum(insertionHapCMMaxSensitivity))

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

# na185071KGVariants <- readVcf('~/Documents/gene_rearrange/svpipeline/NA18507/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.NA18507.vcf.gz', "hg19")
# 
# mcols(rowData(na185071KGVariants))$gt <- as.factor(geno(na185071KGVariants)$GT)
# mcols(rowData(na185071KGVariants))$hap <- ifelse(mcols(rowData(na185071KGVariants))$gt == "1|1", 2, 1)
# 
# na185071KGDels <- rowData(na185071KGVariants)[nchar(fixed(na185071KGVariants)$REF) > nchar(unlist(fixed(na185071KGVariants)$ALT))]
# 
# na185071KGDels <- na185071KGDels[width(na185071KGDels) >= 40]

mcols(trueDelsNA18507)$hap <- NA

newNA185071KGData <- read.table('/Users/cwhelan/Documents/gene_rearrange/svpipeline/NA18507/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.NA18507.indels_and_svs.txt.gz',
                                    col.names=c('chrom', 'pos', 'ref', 'alt', 'gtr', 'gt', 'type', 'svtype', 'imprecise', 'end', 'cipos', 'ciend'))

newNA185071KGDelsData <- newNA185071KGData[newNA185071KGData$type == 'SV' 
                                           | (newNA185071KGData$type == 'INDEL' & 
                                                 (unlist(lapply(as.character(newNA185071KGData$ref), nchar)) > unlist(lapply(as.character(newNA185071KGData$alt), nchar)))),]

newNA185071KGDels <- GRanges(seqnames=newNA185071KGDelsData$chrom,
                             ranges=IRanges(start=newNA185071KGDelsData$pos,
                                            end=ifelse(newNA185071KGDelsData$type == "INDEL",
                                                       newNA185071KGDelsData$pos + unlist(lapply(as.character(newNA185071KGDelsData$ref), nchar)),
                                                       as.numeric(ifelse(newNA185071KGDelsData$end == ".", "0", as.character(newNA185071KGDelsData$end))))),
                            hap=ifelse(newNA185071KGDelsData$gt == "1|1", 2, 1))
newNA185071KGDels <- newNA185071KGDels[width(newNA185071KGDels) >= 40]

ol <- findOverlaps(trueDelsNA18507, newNA185071KGDels)  

mcols(trueDelsNA18507[as.matrix(ol)[,1]])$hap <- mcols(newNA185071KGDels[as.matrix(ol)[,2]])$hap

# ROC curve
totalDelsNA18507 <- 15000
NA18507cloudbreakGEMDelsPerf <- read.table(cloudbreakGEMNA18507DeletionPerfFile, header=TRUE, sep="\t")
NA18507cloudbreakBWADelsPerf <- read.table(cloudbreakBWANA18507DeletionPerfFile, header=TRUE, sep="\t")
NA18507breakdancerDelsPerf <- read.table(breakdancerNA18507DeletionPerfFile, header=TRUE)
NA18507gasvDelsPerf <- read.table(gasvProNA18507DeletionPerfFile, header=TRUE)
NA18507dellyDelsPerf <- read.table(dellyNA18507DeletionPerfFile, header=TRUE)
NA18507dellyBRDelsPerf <- read.table(dellyBRNA18507DeletionPerfFile, header=TRUE)
NA18507pindelDelsPerf <- read.table(pindelNA18507DeletionPerfFile, header=TRUE)

perfsListDelsNA18507 <- list(cloudbreak=NA18507cloudbreakBWADelsPerf, breakdancer=NA18507breakdancerDelsPerf, pindel=NA18507pindelDelsPerf, gasv=NA18507gasvDelsPerf, delly=NA18507dellyDelsPerf, dellyBR=NA18507dellyBRDelsPerf)
pdf(NA18507DeletionsROCOutputFile, width=10)
par(xpd=T, mar=par()$mar+c(0,0,0,7))
plotROC(perfsListDelsNA18507, c("Cloudbreak", "Breakdancer","Pindel","GASVPro", "DELLY-RP", "DELLY-SR"), totalDelsNA18507, 
        "Deletions in NA18507", legendLoc=xy.coords(15750,1250), maxTP=2000, sim=FALSE, colLtyMapping=myColLtyMapping)
dev.off()

cairo_ps("~/Documents/svpipeline/figures/NA18507_del_roc_cairo.ps", width=10, height=7)
par(xpd=T, mar=par()$mar+c(.5,.5,.5,13))
plotROC(perfsListDelsNA18507, 
        c("Cloudbreak", "Breakdancer","Pindel", "GASVPro", "DELLY-RP", "DELLY-SR"), 
        totalDelsNA18507, "Deletions",legendLoc=xy.coords(15750,1250), maxTP=2000, sim=FALSE,
        cex.main=1.5, cex.axis=1.5, cex.lab=1.5, legendCex=1.5, colLtyMapping=myColLtyMapping)
dev.off()


# break down predictons

totalCloudbreakNA18507DelPredictions <- dim(read.table(cloudbreakBWANA18507DelPredictionsFile, skip=1))[1]

totalBreakdancerNA18507DelPredictions <- numPredictionsAtLevel(NA18507breakdancerDelsPerf, NA18507delCutoffs['breakdancer'])
totalGASVProNA18507DelPredictions <- numPredictionsAtLevel(NA18507gasvDelsPerf, NA18507delCutoffs['GASVPro'])
totalDELLYNA18507DelPredictions <- numPredictionsAtLevel(NA18507dellyDelsPerf, NA18507delCutoffs['DELLY'])
totalDELLYBRNA18507DelPredictions <- numPredictionsAtLevel(NA18507dellyBRDelsPerf, NA18507delCutoffs['DELLYSR'])
totalPindelNA18507DelPredictions <- numPredictionsAtLevel(NA18507pindelDelsPerf, NA18507delCutoffs['Pindel'])

totalNA18507DelPredictions <- list(Cloudbreak=totalCloudbreakNA18507DelPredictions, 
                                   Breakdancer=totalBreakdancerNA18507DelPredictions, 
                                   GASVPro=totalGASVProNA18507DelPredictions, 
                                   DELLY=totalDELLYNA18507DelPredictions, 
                                   DELLYBR=totalDELLYBRNA18507DelPredictions, 
                                   Pindel=totalPindelNA18507DelPredictions)

cbHitsNA18507 <- processDeletionPredictions("Cloudbreak", cloudbreakBWANA18507DeletionHitsFile, trueDelsNA18507)
bdHitsNA18507 <- processDeletionPredictions('Breakdancer', breakdancerNA18507DeletionHitsFile, trueDelsNA18507)
gasvHitsNA18507 <- processDeletionPredictions('GASVPro', gasvProNA18507DeletionHitsFile, trueDelsNA18507)
dellyHitsNA18507 <- processDeletionPredictions('DELLY-RP', dellyNA18507DeletionHitsFile, trueDelsNA18507)
dellyBRHitsNA18507 <- processDeletionPredictions('DELLY-SR', dellyBRNA18507DeletionHitsFile, trueDelsNA18507)
pindelHitsNA18507 <- processDeletionPredictions('Pindel', pindelNA18507DeletionHitsFile, trueDelsNA18507)

totalNumNA18507DeletionHits <- list(list(Cloudbreak=dim(cbHitsNA18507)[1], 
                                         Breakdancer=dim(bdHitsNA18507)[1], 
                                         GASVPro=dim(gasvHitsNA18507)[1], 
                                         DELLY=dim(dellyHitsNA18507)[1], 
                                         DELLYSR=dim(dellyBRHitsNA18507)[1], 
                                         Pindel=dim(pindelHitsNA18507)[1]))

NA18507DelPrecision <- unlist(totalNumNA18507DeletionHits) / unlist(totalNA18507DelPredictions)
NA18507DelRecall <- unlist(totalNumNA18507DeletionHits) / length(trueDelsNA18507Gte40)

allpredsNA18507 <- rbind(as.data.frame(cbHitsNA18507),
                         as.data.frame(bdHitsNA18507),
                         as.data.frame(gasvHitsNA18507),
                         as.data.frame(dellyHitsNA18507),
                         as.data.frame(dellyBRHitsNA18507),
                         as.data.frame(pindelHitsNA18507))
allpredsNA18507 <- ddply(allpredsNA18507, .(tok), function(x) { cbind(x, list(discoveries=rep(dim(x)[1],dim(x)[1]))) })

predsBySizeNA18507 <- table(allpredsNA18507[,c('name','sizeclasses')])
exBySizeNA18507 <- table(allpredsNA18507[which(allpredsNA18507$discoveries == 1),c('name','sizeclasses')])

predsByRepMaskNA18507 <- table(allpredsNA18507[,c('name','repMask')])
exByRepMaskNA18507 <- table(allpredsNA18507[which(allpredsNA18507$discoveries == 1),c('name','repMask')])

predsByHapNA18507 <- table(allpredsNA18507[,c('name','hap')])
exByHapNA18507 <- table(allpredsNA18507[which(allpredsNA18507$discoveries == 1),c('name','hap')])

# genotyping 

extraDataNA18507 <- read.table(cloudbreakBWANA18507DelPredictionsFile, skip=1)
names(extraDataNA18507) <- c('predchrom', 'predstart', 'predend', 'peaknum', 'maxscore', 'svType', 'avgMu', 'minMu', 'maxMu', 'avgW', 'minW', 'maxW')

cbpredsWithExtraDataNA18507 <- merge(as.data.frame(cbHitsNA18507), extraDataNA18507)

cbpredsWithExtraDataNA18507Genotyped <- cbpredsWithExtraDataNA18507[!is.na(cbpredsWithExtraDataNA18507$hap),]

NA18507HapCM <- table(cbpredsWithExtraDataNA18507$hap, cbpredsWithExtraDataNA18507$avgW < genotypingAlphaCutoff)
NA18507DeletionHapAccuracy <- (NA18507HapCM['1','FALSE'] + NA18507HapCM['2','TRUE']) / sum(NA18507HapCM)

millsDelsWithGenotypes <- read.table('~/Documents/Papers2/Articles/2011/Mills/Supplemental/dels_gt40_with_genotypes.csv', header=TRUE, sep=",")
millsGenotypeData <- read.table('~/Documents/Papers2/Articles/2011/Mills/Supplemental/Mills_genotypes_NA18507_with_ids_gt40.txt', header=TRUE, sep="\t")

millsFullGenotypedDels <- merge(millsDelsWithGenotypes,millsGenotypeData)

millsGenotypedRanges <- GRanges(seqnames=substr(millsFullGenotypedDels$CHR,4,length(millsFullGenotypedDels$CHR)), ranges=IRanges(start=millsFullGenotypedDels$START, end=millsFullGenotypedDels$STOP), mcols=DataFrame(hap=millsFullGenotypedDels$SED00002_NA18507.CEL))

#ol <- findOverlaps(trueDelsNA18507, millsGenotypedRanges)  

# mcols(trueDelsNA18507[as.matrix(ol)[,1]])$hap <- mcols(na185071KGDels[as.matrix(ol)[,2]])$hap

#NA18507 100bp diploid insertions
totalInsertionsNA18507 <- 8000
cloudbreakGEMNA18507InsertionsPerf <- read.table(cloudbreakGEMNA18507InsertionsPerfFile, header=TRUE, sep="\t")
cloudbreakGEMNA18507InsertionsPerfSmoothed <- cloudbreakGEMNA18507InsertionsPerf[seq(1,dim(cloudbreakGEMNA18507InsertionsPerf)[1],by=4),]

cloudbreakBWANA18507InsertionsPerf <- read.table(cloudbreakBWANA18507InsertionsPerfFile, header=TRUE, sep="\t")
cloudbreakBWANA18507InsertionsPerfSmoothed <- cloudbreakBWANA18507InsertionsPerf[seq(1,dim(cloudbreakBWANA18507InsertionsPerf)[1],by=4),]

cloudbreakBWATweakNA18507InsertionsPerf <- read.table(cloudbreakBWATweakNA18507InsertionsPerfFile, header=TRUE, sep="\t")
cloudbreakBWATweakNA18507InsertionsPerfSmoothed <- cloudbreakBWATweakNA18507InsertionsPerf[seq(1,dim(cloudbreakBWATweakNA18507InsertionsPerf)[1],by=4),]

breakdancerNA18507InsertionsPerf <- read.table(breakdancerNA18507InsertionsPerfFile, header=TRUE)
pindelNA18507InsertionsPerf <- read.table(pindelNA18507InsertionsPerfFile, header=TRUE)

perfsListInsertionsNA18507 <- list(cloudbreak=cloudbreakBWATweakNA18507InsertionsPerfSmoothed, breakdancer=breakdancerNA18507InsertionsPerf, pindel=pindelNA18507InsertionsPerf)
pdf(NA18507InsertionsROCOutputFile, width=10)
par(xpd=T, mar=par()$mar+c(0,0,0,9))
plotROC(perfsListInsertionsNA18507, 
        c("Cloudbreak", "Breakdancer","Pindel"), 
        totalInsertionsNA18507, "Insertions in NA18507",legendLoc=xy.coords(8500,100), maxTP=300, sim=FALSE, colLtyMapping=myColLtyMapping)
dev.off()

#Length adjusted NA18507 100bp diploid insertions
totalInsertionsNA18507 <- 8000
cloudbreakBWANA18507LAInsertionsPerf <- read.table(cloudbreakBWANA18507LAInsertionsPerfFile, header=TRUE, sep="\t")
cloudbreakBWANA18507LAInsertionsPerfSmoothed <- cloudbreakBWANA18507LAInsertionsPerf[seq(1,dim(cloudbreakBWANA18507LAInsertionsPerf)[1],by=4),]

cloudbreakBWATweakNA18507LAInsertionsPerf <- read.table(cloudbreakBWATweakNA18507LAInsertionsPerfFile, header=TRUE, sep="\t")
cloudbreakBWATweakNA18507LAInsertionsPerfSmoothed <- cloudbreakBWATweakNA18507LAInsertionsPerf[seq(1,dim(cloudbreakBWATweakNA18507LAInsertionsPerf)[1],by=4),]

breakdancerNA18507LAInsertionsPerf <- read.table(breakdancerNA18507LAInsertionsPerfFile, header=TRUE)
pindelNA18507LAInsertionsPerf <- read.table(pindelNA18507LAInsertionsPerfFile, header=TRUE)

perfsListLAInsertionsNA18507 <- list(breakdancer=breakdancerNA18507InsertionsPerf, pindel=pindelNA18507LAInsertionsPerf, cloudbreakBWA=cloudbreakBWANA18507LAInsertionsPerfSmoothed, cloudbreakBWA2=cloudbreakBWATweakNA18507LAInsertionsPerfSmoothed)
pdf(NA18507LAInsertionsROCOutputFile, width=10)
par(xpd=T, mar=par()$mar+c(0,0,0,17))
plotROC(perfsListLAInsertionsNA18507, 
        c("Breakdancer","Pindel", "Cloudbreak-BWA", "Cloudbreak-BWA-2"), 
        totalInsertionsNA18507, "Length-adjusted Insertions in NA18507",legendLoc=xy.coords(8500,200), maxTP=400, sim=FALSE, colLtyMapping=myColLtyMapping)
dev.off()


cairo_ps("~/Documents/svpipeline/figures/NA18507_ins_roc_cairo.ps", width=10, height=7)
par(xpd=T, mar=par()$mar+c(.5,.5,.5,13))
plotROC(perfsListInsertionsNA18507, 
        c("Cloudbreak", "Breakdancer","Pindel"), 
        totalInsertionsNA18507, "Insertions",legendLoc=xy.coords(8500,100), maxTP=300, sim=FALSE,
        cex.main=1.5,cex.axis=1.5,cex.lab=1.5,legendCex=1.5, colLtyMapping=myColLtyMapping)
dev.off()


pdf(NA18507ROCsOutputFile, width=15)
par(xpd=T, mfrow=(c(1,2)), oma=par()$oma+c(0,0,0,5))
plotROC(perfsListDelsNA18507, 
        c("Cloudbreak", "Breakdancer","Pindel", "GASVPro", "DELLY-RP", "DELLY-SR"), 
        totalDelsNA18507, "Deletions in NA18507",legendLoc=NULL, maxTP=2000, sim=FALSE, colLtyMapping=myColLtyMapping)
par(mar=par()$mar+c(0,0,0,6))
plotROC(perfsListInsertionsNA18507, 
        c("Cloudbreak", "Breakdancer","Pindel"), 
        totalInsertionsNA18507, "Insertions in NA18507",legendLoc=NA, maxTP=300, sim=FALSE, colLtyMapping=myColLtyMapping)
perfNames <- c("Cloudbreak", "Breakdancer","Pindel", "GASVPro", "DELLY-RP", "DELLY-SR")
legend(xy.coords(8500,175), legend=perfNames, col=getCols(perfNames, myColLtyMapping), lty=getLtys(perfNames, myColLtyMapping), lwd=3, cex=.9)
dev.off()

pdf(NA18507ROCsPosterOutputFile, width=18)
par(xpd=T, mfrow=(c(1,2)), oma=par()$oma+c(0,0,0,3), cex=1.8)
plotROC(perfsListDelsNA18507, 
        c("Cloudbreak", "Breakdancer","Pindel", "GASVPro", "DELLY-RP", "DELLY-SR"), 
        totalDelsNA18507, "Deletions - NA18507",legendLoc=NULL, maxTP=2000, sim=FALSE)
par(mar=par()$mar+c(0,0,0,8))
plotROC(perfsListInsertionsNA18507, 
        c("Cloudbreak", "Breakdancer","Pindel"), 
        totalInsertionsNA18507, "Insertions - NA18507",legendLoc=NA, maxTP=300, sim=FALSE)
legend(xy.coords(9750,250), legend=perfNames, col=getCols(perfNames, myColLtyMapping), lty=getLtys(perfNames, myColLtyMapping), lwd=3)
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

mcols(trueInsertionsNA18507)$hap <- NA

# # right now not getting any genotypes from the 1KG insertions
# 
# na185071KGInsLengths <- nchar(unlist(fixed(na185071KGVariants)$ALT)) - nchar(fixed(na185071KGVariants)$REF)
# na185071KGIns <- rowData(na185071KGVariants)[na185071KGInsLengths >= 40]
# 
# ol <- findOverlaps(trueInsertionsNA18507, na185071KGIns) 
# mcols(trueInsertionsNA18507[as.matrix(ol)[,1]])$hap <- mcols(trueInsertionsNA18507[as.matrix(ol)[,2]])$hap

millsInsertionsWithIds <- read.table('/Users/cwhelan/Documents/Papers2/Articles/2011/Mills/Supplemental/Mills_2011_NA18507_INS.b37.bed')
names(millsInsertionsWithIds) <- c("chrom", "start", "end", "length", "probeset_id")

millsGenotypes <- read.table('/Users/cwhelan/Documents/Papers2/Articles/2011/Mills/Supplemental/Mills_genotypes_NA18507_with_ids.txt', header=TRUE, sep="\t")

genotypedMillsInsertions <- merge(millsInsertionsWithIds, millsGenotypes, by="probeset_id")
 genotypedMillsInsertions$shorttok <- paste(genotypedMillsInsertions$chrom, ':', genotypedMillsInsertions$start)

totalCloudbreakNA18507InsertionPredictions <- numPredictionsAtLevel(cloudbreakBWATweakNA18507InsertionsPerf, NA18507insertionCutoffs['cloudbreakBWATweak']) 
totalBreakdancerNA18507InsertionPredictions <- numPredictionsAtLevel(breakdancerNA18507InsertionsPerf, NA18507insertionCutoffs['breakdancer']) 
totalPindelNA18507InsertionPredictions <- numPredictionsAtLevel(pindelNA18507InsertionsPerf, NA18507insertionCutoffs['Pindel']) 

# load true positive predictions
cbInsertionsNA18507 <- processInsertionPredictions("Cloudbreak", cloudbreakBWATweakNA18507InsertionHitsFile, trueInsertionsNA18507)
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
NA18507insertionExPredsByHap <- table(allpredsInsertionsNA18507[which(allpredsInsertionsNA18507$discoveries == 1),c('name','hap')])

#NA18507insertionPredsBySegDup <- table(allpredsInsertionsNA18507[,c('name','segdup')])
#NA18507insertionExPredsBySegDup <- table(allpredsInsertions[which(allpredsInsertionsNA18507$discoveries == 1),c('name','segdup')])

# genotyping insertions
#cloudbreakNA18507InsertionsWithExtraData <- read.table(cloudbreakBWANA18507InsertionPredictionsFile, skip=1)
#names(cloudbreakNA18507InsertionsWithExtraData) <- c('predchrom', 'predstart', 'predend', 'peaknum', 'maxscore', 'svtype', 'avgMu', 'minMu', 'maxMu', 'avgW', 'minW', 'maxW')

#NA18507cbInsertionPredsWithExtraData <- merge(as.data.frame(cbInsertionsNA18507), cloudbreakNA18507InsertionsWithExtraData)

#NA18507insertionHapCM <- table(NA18507cbInsertionPredsWithExtraData$hap, NA18507cbInsertionPredsWithExtraData$avgW < .2)
#print(NA18507insertionHapCM)
#print((NA18507insertionHapCM[1,1] + NA18507insertionHapCM[2,2]) / sum(NA18507insertionHapCM))

# LaTeX tables

# chr 2 Deletions by size at FDR10
tableEnv = new.env()
assign("totals", trueDelSizes, envir=tableEnv)
assign("preds", predsBySizeFDR10, envir=tableEnv)
assign("exPreds", exBySizeFDR10, envir=tableEnv)
brew('~/Documents/svpipeline/manuscript/chr2DeletionPredsFDR10.table.brew.tex', env=tableEnv)

# chr 2 Deletions by size at max sensitivity
tableEnv = new.env()
assign("totals", trueDelSizes, envir=tableEnv)
assign("preds", predsBySizeMaxSensitivity, envir=tableEnv)
assign("exPreds", exBySizeMaxSensitivity, envir=tableEnv)
assign("precision", chr2MaxSensitivityDelPrecision, envir=tableEnv)
assign("recall", chr2MaxSensitivityDelRecall, envir=tableEnv)
brew('~/Documents/svpipeline/manuscript/chr2DeletionPredsMaxSensitivity.table.brew.tex', env=tableEnv)

# chr 2 Insertions by size
tableEnv = new.env()
assign("totals", trueInsertionsSizes, envir=tableEnv)
assign("preds", insertionPredsBySizeMaxSensitivity, envir=tableEnv)
assign("exPreds", insertionExBySizeMaxSensitivity, envir=tableEnv)
assign("precision", chr2MaxSensitivityInsPrecision, envir=tableEnv)
assign("recall", chr2MaxSensitivityInsRecall, envir=tableEnv)
brew('~/Documents/svpipeline/manuscript/chr2InsertionPredsBySize.table.brew.tex', env=tableEnv)

# chr 2 Deletions and insertions by size at max sensitivity (combined table)
tableEnv = new.env()
assign("totalDels", trueDelSizes, envir=tableEnv)
assign("delPreds", predsBySizeMaxSensitivity, envir=tableEnv)
assign("exDelPreds", exBySizeMaxSensitivity, envir=tableEnv)
assign("delPrecision", chr2MaxSensitivityDelPrecision, envir=tableEnv)
assign("delRecall", chr2MaxSensitivityDelRecall, envir=tableEnv)
assign("totalIns", trueInsertionsSizes, envir=tableEnv)
assign("insPreds", insertionPredsBySizeMaxSensitivity, envir=tableEnv)
assign("exInsPreds", insertionExBySizeMaxSensitivity, envir=tableEnv)
assign("insPrecision", chr2MaxSensitivityInsPrecision, envir=tableEnv)
assign("insRecall", chr2MaxSensitivityInsRecall, envir=tableEnv)
brew('~/Documents/svpipeline/manuscript/chr2DeletionAndInsertionPredsMaxSensitivity.table.brew.tex', env=tableEnv)

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

# NA18507 Deletions and Insertions combined table
tableEnv = new.env()
assign("totalDels", trueDelsNA18507Sizes, envir=tableEnv)
assign("delPrecision", NA18507DelPrecision, envir=tableEnv)
assign("delRecall", NA18507DelRecall, envir=tableEnv)
assign("delPreds", predsBySizeNA18507, envir=tableEnv)
assign("exDelPreds", exBySizeNA18507, envir=tableEnv)
assign("totalIns", trueInsertionsNA18507Sizes, envir=tableEnv)
assign("insPrecision", NA18507InsertionPrecision, envir=tableEnv)
assign("insRecall", NA18507InsertionRecall, envir=tableEnv)
assign("insPreds", NA18507insertionPredsBySize, envir=tableEnv)
assign("exInsPreds", NA18507insertionExPredsBySize, envir=tableEnv)
brew('~/Documents/svpipeline/manuscript/NA18507DeletionAndInsertionPredsBySize.table.brew.tex', env=tableEnv)

# deletions in repetitive regions
tableEnv = new.env()
assign("chr2Totals", trueDelRepmask, envir=tableEnv)
assign("chr2Preds", predsByRepMaskFDR10, envir=tableEnv)
assign("chr2ExPreds", exPredsByRepMaskFDR10, envir=tableEnv)
assign("NA18507Totals", trueDelRepmaskNA18507, envir=tableEnv)
assign("NA18507Preds", predsByRepMaskNA18507, envir=tableEnv)
assign("NA18507ExPreds", exByRepMaskNA18507, envir=tableEnv)
brew('~/Documents/svpipeline/manuscript/deletionRepmask.table.brew.tex', env=tableEnv)

# insertions in repetitive regions
tableEnv = new.env()
assign("chr2Totals", trueInsertionsRepmask, envir=tableEnv)
assign("chr2Preds", insertionPredsByRepMaskMaxSensitivity, envir=tableEnv)
assign("chr2ExPreds", insertionExPredsByRepMaskMaxSensitivity, envir=tableEnv)
assign("NA18507Totals", trueInsertionsNA18507Repmask, envir=tableEnv)
assign("NA18507Preds", NA18507insertionPredsByRepMask, envir=tableEnv)
assign("NA18507ExPreds", NA18507insertionExPredsByRepMask, envir=tableEnv)
brew('~/Documents/svpipeline/manuscript/insertionRepmask.table.brew.tex', env=tableEnv)

# genotyping
tableEnv = new.env()
assign("chr2HapCM", chr2HapCMFDR10, envir=tableEnv)
assign("NA18507HapCM", NA18507HapCM, envir=tableEnv)
brew('~/Documents/svpipeline/manuscript/deletionGenotypeAccuracy.table.brew.tex', env=tableEnv)

# breakpoint resolution
pdf(breakpointResolutionOutputFile)
allpredsFDR10$lendiff <- with(allpredsFDR10, abs((predend - predstart) - (trueend - truestart)))
ggplot(aes(y=lendiff,x=name),data=allpredsFDR10) + geom_boxplot() + xlab("") + ylab("Difference in length")
dev.off()

# plot EC2 runtimes on chr2 features with different numbers of reducers
nodes <- c(5,10,20,30,40)
featureCompTimes <- c(1491,1049,552,412,458)
runtimeDF <- data.frame(nodes=nodes, featureCompTimes=featureCompTimes)
pdf(runtimeByNumberOfNodesOutputFile)
ggplot(aes(x=nodes,y=featureCompTimes), data=runtimeDF) + geom_line() + xlab("Nodes") + ylab("Feature Computation Time (seconds)") + labs(title="Runtime by Number of Nodes in Cluster")
dev.off()

bestChr2Runtimes <- data.frame(name=c('Cloudbreak', 'Breakdancer', 'Pindel','GASVPro','DELLY' ), runtimes=c(290,653,4885,3339,1964))
palette <- rocColors(5)

renameColorNames <- function(names, nameMap=list()) {
  lapply(as.character(bestChr2Runtimes$name), function(x) { if (x %in% names(nm)) { nm[[x]] } else { x } })
}

colorsByRuntime <- function(runtimes, colLtyMapping, nameMap=list()) {
  runtimes$name <- renameColorNames(runtimes$name, nameMap)
  unname(unlist(lapply(colLtyMapping[as.character(runtimes[order(runtimes$runtimes, runtimes$name), 'name'])], function(x) {x[1]})))
}

bestChr2Runtimes$colorNames <- renameColorNames(bestChr2Runtimes$name, nameMap=list(DELLY="DELLY-RP"))

cairo_pdf('~/Documents/svpipeline/figures/chr2BestRuntimes.pdf', height=6.5, width=3)
ggplot(bestChr2Runtimes, aes(x=reorder(name,runtimes), y=runtimes)) + geom_point(aes(colour=reorder(name,runtimes)), size=5) +
  scale_colour_manual(values=colorsByRuntime(bestChr2Runtimes, myColLtyMapping, nameMap=list(DELLY="DELLY-RP"))) + 
  xlab("") + 
  ylab("Runtime (s)") + theme_bw()  + theme(axis.text.x=element_text(angle=-90)) + theme(text=element_text(size=18)) +
  theme(panel.grid.major.x=element_line(linetype = c("28"), color="black")) + theme(legend.position="none") +
  theme(panel.border=element_rect(colour="black"))
dev.off()

bestNA18507Runtimes <- data.frame(name=c('Cloudbreak', 'Breakdancer', 'Pindel','GASVPro','DELLY' ), runtimes=c(824,5586,28587,52385,20224))
bestNA18507Runtimes$colorNames <- renameColorNames(bestNA18507Runtimes$name, nameMap=list(DELLY="DELLY-RP"))

cairo_pdf('~/Documents/svpipeline/figures/NA18507BestRuntimes.pdf', height=6.5, width=3)
ggplot(bestNA18507Runtimes, aes(x=reorder(name,runtimes), y=runtimes)) + geom_point(aes(colour=reorder(name,runtimes)), size=5) +
  scale_colour_manual(values=colorsByRuntime(bestChr2Runtimes, myColLtyMapping, nameMap=list(DELLY="DELLY-RP"))) + xlab("") + 
  ylab("Runtime (s)") + theme_bw()  + theme(axis.text.x=element_text(angle=-90)) + theme(text=element_text(size=18)) +
  theme(panel.grid.major.x=element_line(linetype = c("28"), color="black")) + theme(legend.position="none") +
  theme(panel.border=element_rect(colour="black"))
dev.off()

library(grid)
setEPS()
postscript('~/Documents/svpipeline/figures/test_chr2.eps', width=20)
layout(c(1,2,3), c(4,4,2))
par(xpd=T, mfrow=(c(1,2)), oma=par()$oma+c(0,0,0,6))
plotROC(perfsListDelsMultimapChr2, 
        c("Breakdancer","Pindel", "GASVPro", "DELLY-RP", "DELLY-SR", 
          "Cloudbreak", "Cloudbreak-BWA-PE-MM", "Cloudbreak-GEM-SE-MM"), 
        totalDelsChr2, "Deletions",legendLoc=NULL, maxTP=350, colLtyMapping=myColLtyMapping)
par(mar=par()$mar+c(0,0,0,12))
plotROC(perfsListInsertionsMultimapChr2, 
        c("Breakdancer","Pindel", "Cloudbreak", "Cloudbreak-BWA-PE-MM", "Cloudbreak-GEM-SE-MM"), 
        totalInsertionsChr2, "Insertions",legendLoc=NA,  maxTP=totalInsertionsChr2, colLtyMapping=myColLtyMapping)
perfNames <- c("Breakdancer","Pindel", "GASVPro", "DELLY-RP", "DELLY-SR", "Cloudbreak", "Cloudbreak-BWA-PE-MM", "Cloudbreak-GEM-SE-MM")
legend(xy.coords(totalInsertionsChr2 * 1.1,100), legend=perfNames, col=getCols(perfNames, myColLtyMapping), lty=getLtys(perfNames, myColLtyMapping), lwd=3, cex=.9)
ggplot(bestChr2Runtimes, aes(x=reorder(name,runtimes), y=runtimes)) + geom_point(aes(colour=reorder(name,runtimes)), size=5) +
  scale_colour_manual(values=palette[order(bestChr2Runtimes$runtimes)]) + xlab("") + 
  ylab("Runtime (s)") + theme_bw()  + theme(axis.text.x=element_text(angle=-90)) + theme(text=element_text(size=18)) +
  theme(panel.grid.major.x=element_line(linetype = c("28"), color="black")) + theme(legend.position="none") +
  theme(panel.border=element_rect(colour="black"))
dev.off()
