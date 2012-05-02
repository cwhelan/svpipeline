drawPerfLine <- function(perf, col, maxFP) {
  calls <- c(max(max(perf$Calls), maxFP + max(perf$TP)), perf$Calls, min(perf$Calls) - min(perf$TP), 0)
  tp <- c(max(perf$TP), perf$TP, 0, 0)
  fp <- calls - tp
  lines(fp[fp <= maxFP], tp[fp <= maxFP], col=col, lwd=3)
}

plotROC <- function(perfs, perfNames, totalDels, main, sim=TRUE) {
  maxFP <- totalDels
  plot(0, type="n", ylim=c(0, totalDels), xlim=c(0,maxFP), xlab=ifelse(sim, "False Positives", "Novel Predictions"), ylab="True Positives", main=main)
  perfCols <- rainbow(length(perfs))
  mapply(drawPerfLine, perfs, perfCols, MoreArgs=list(maxFP=maxFP))  
  legend("top", legend=perfNames, col=perfCols, lwd=3, cex=.5)
}

#chr2 gt 125
totalDels <- 197
hydraPerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim/hydra_perf_chr2_sim_gt125.txt', header=TRUE)
svpipelinePerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim/jcvi_chr2_sim_t180_nosd_or_piledup_deletion_scores_perf_mw3_gt125.txt', header=TRUE)
breakdancerPerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim/breakdancer_bedpe_perf_gt125.txt', header=TRUE)

perfsList <- list(svpipeline=svpipelinePerf, hydra=hydraPerf, breakdancer=breakdancerPerf)
plotROC(perfsList, c("SVPipeline", "Hydra", "Breakdancer"), totalDels)

#chr2 gt 50
totalDels <- 372
hydraPerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim/hydra2_perf_gt50.txt', header=TRUE)
mo3Perf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim/mo3_r25_gt50_quant_perf.txt', header=TRUE)
mo5Perf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim/mo5_r25_gt50_quant_perf.txt', header=TRUE)
breakdancerPerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim/breakdancer2_perf_gt50.txt', header=TRUE)
gasvPerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim/gasv_perf_gt50.txt', header=TRUE)

perfsList <- list(hydra=hydraPerf, breakdancer=breakdancerPerf, gasv=gasvPerf, mo3=mo3Perf)
plotROC(perfsList, c("Hydra", "Breakdancer", "GASV", "SVP"), totalDels, "chr2 Simulated (30X)")

# chr17
totalDels <- 199
hydra <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr17_sim/hydra2_perf_gt50.txt', header=TRUE)
breakdancer <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr17_sim/breakdancer2_perf_gt50.txt', header=TRUE)
gasv <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr17_sim/gasv_perf_gt50.txt', header=TRUE)
mo3 <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr17_sim/mo3_gt50_quant_perf.txt', header=TRUE)
mo3blgn <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr17_sim/mo3_biglognorm_mw3_perf.txt', header=TRUE)

aa3 <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr17_sim/test_mapability_weighting_avg_cutoff_nosd_and_mw3_perf.txt', header=TRUE)
ma3 <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr17_sim/test_mapability_weighting_min_cutoff_nosd_and_mw3_perf.txt', header=TRUE)
aa5 <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr17_sim/test_mapability_weighting_avg_cutoff_nosd_and_mw5_perf.txt', header=TRUE)
ma5 <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr17_sim/test_mapability_weighting_min_cutoff_nosd_and_mw5_perf.txt', header=TRUE)

ao3 <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr17_sim/test_mapability_weighting_avg_cutoff_nosd_or_mw3_perf.txt', header=TRUE)
#mo3 <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr17_sim/test_mapability_weighting_min_cutoff_nosd_or_mw3_perf.txt', header=TRUE)
ao5 <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr17_sim/test_mapability_weighting_avg_cutoff_nosd_or_mw5_perf.txt', header=TRUE)
mo5 <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr17_sim/test_mapability_weighting_min_cutoff_nosd_or_mw5_perf.txt', header=TRUE)

perfsList <- list(hydra, breakdancer, aa3, ma3, aa5, ma5, ao3, mo3, ao5, mo5)
#mo3short <- mo3[mo3$FP < totalDels,]
perfsList <- list(hydra, breakdancer, gasv, mo3)
plotROC(perfsList, c("Hydra", "Breakdancer", "GASV", "SVP"), totalDels, "chr17 Simulated (100X)")

# NA07051
totalDels <- 2548
hydra <- read.table('~/Documents/gene_rearrange/svpipeline/NA07051/hydra_perf.txt', header=TRUE)
breakdancer <- read.table('~/Documents/gene_rearrange/svpipeline/NA07051/breakdancer_perf.txt', header=TRUE)
#gasv <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr17_sim/gasv_perf_gt50.txt', header=TRUE)
mo3 <- read.table('~/Documents/gene_rearrange/svpipeline/NA07051/mo3_mw3_perf.txt', header=TRUE)

perfsList <- list(hydra, breakdancer, mo3)
plotROC(perfsList, c("Hydra", "Breakdancer", "SVP"), totalDels, "NA07051 (3X)")



