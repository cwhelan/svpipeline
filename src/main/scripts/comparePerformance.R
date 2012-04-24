drawPerfLine <- function(perf, col, maxFP) {
  calls <- c(max(max(perf$Calls), maxFP + max(perf$TP)), perf$Calls, min(perf$Calls) - min(perf$TP), 0)
  tp <- c(max(perf$TP), perf$TP, 0, 0)
  lines(calls - tp, tp, col=col, lwd=3)
}

plotROC <- function(perfs, perfNames, totalDels) {
  maxFP <- totalDels
  plot(0, type="n", ylim=c(0, totalDels), xlim=c(0,maxFP), xlab="False Positives", ylab="True Positives")
  perfCols <- rainbow(length(perfs))
  mapply(drawPerfLine, perfs, perfCols, MoreArgs=list(maxFP=maxFP))  
  legend("bottomright", legend=perfNames, col=perfCols, lwd=3)
}

totalDels <- 197
hydraPerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim/hydra_perf_chr2_sim_gt125.txt', header=TRUE)
svpipelinePerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim/jcvi_chr2_sim_t180_nosd_or_piledup_deletion_scores_perf_mw3_gt125.txt', header=TRUE)
breakdancerPerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim/breakdancer_bedpe_perf_gt125.txt', header=TRUE)

perfsList <- list(svpipeline=svpipelinePerf, hydra=hydraPerf, breakdancer=breakdancerPerf)
plotROC(perfsList, c("SVPipeline", "Hydra", "Breakdancer"), totalDels)

totalDels <- 372
hydraPerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim/hydra_perf_chr2_sim_gt50.txt', header=TRUE)
svpipelinePerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim/jcvi_chr2_sim_t180_nosd_or_piledup_deletion_scores_perf_mw3_gt50.txt', header=TRUE)
breakdancerPerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim/breakdancer_bedpe_perf_gt50.txt', header=TRUE)

perfsList <- list(svpipeline=svpipelinePerf, hydra=hydraPerf, breakdancer=breakdancerPerf)
plotROC(perfsList, c("SVPipeline", "Hydra", "Breakdancer"), totalDels)
