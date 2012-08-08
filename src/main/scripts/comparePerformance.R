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
  legend("bottomright", legend=perfNames, col=perfCols, lwd=3, cex=.75)
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
cloudbreak <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim/b292f1f35e6fc997581e19b773ae4dd33beaf58e_r25_perf.txt', header=TRUE)
breakdancerPerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim/breakdancer2_perf_gt50.txt', header=TRUE)
gasvPerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim/gasv_perf_gt50.txt', header=TRUE)
cloudbreak.new <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim/test_max_insert_25000_sdseqfilter_t180_perf.txt', header=TRUE)

perfsList <- list(hydra=hydraPerf, breakdancer=breakdancerPerf, gasv=gasvPerf, mo3=cloudbreak, cbnew=cloudbreak.new)
pdf('~/Documents/svpipeline/CHR2SIM_ROC_NEW.pdf')
plotROC(perfsList, c("Hydra", "Breakdancer","GASV", "Cloudbreak", "Cloudbreak (NEW)"), totalDels, "chr2 Simulated (30X)")
dev.off()

#chr2 gt 50 LOW COVERAGE
totalDels <- 372
hydraPerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim_lc/hydra_perf.txt', header=TRUE)
cloudbreak <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim_lc/b292f1f35e6fc997581e19b773ae4dd33beaf58e_r25_perf.txt', header=TRUE)
breakdancerPerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim_lc/breakdancer_rmdup_perf.txt', header=TRUE)
breakdancerSensitive <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim_lc/breakdancer_rmdup_sensitive_perf.txt', header=TRUE)
cloudbreak.new <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim_lc/test_max_insert_25000_sdseqfilter_t180.perf.txt', header=TRUE)
gasvPerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim_lc/gasv_perf.txt', header=TRUE)

perfsList <- list(hydra=hydraPerf, breakdancer=breakdancerPerf, gasv=gasvPerf, cloudbreak=cloudbreak, cloudbreak.new=cloudbreak.new)
pdf('~/Documents/svpipeline/CHR2SIMLC_ROC_NEW.pdf')
plotROC(perfsList, c("Hydra", "Breakdancer", "GASV", "Cloudbreak", "Cloudbreak (NEW)"), totalDels, "chr2 Simulated (5X)")
dev.off()

#chr2 gt 50 V. LOW COVERAGE
totalDels <- 372
hydraPerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim_vlc/hydra_perf.txt', header=TRUE)
cloudbreak <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim_vlc/b292f1f35e6fc997581e19b773ae4dd33beaf58e_r25_perf.txt', header=TRUE)
breakdancerPerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim_vlc/breakdancer_rmdup_perf.txt', header=TRUE)
breakdancerSensitive <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim_vlc/breakdancer_rmdup_sensitive_perf.txt', header=TRUE)
gasvPerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim_vlc/gasv_perf.txt', header=TRUE)

perfsList <- list(hydra=hydraPerf, breakdancer=breakdancerPerf, gasv=gasvPerf, cloudbreak=cloudbreak)
pdf('~/Documents/svpipeline/CHR2SIM_VLC_ROC.pdf')
plotROC(perfsList, c("Hydra", "Breakdancer", "GASV", "Cloudbreak"), totalDels, "chr2 Simulated (3X)")
dev.off()


#chr2 gt 50 V. LOW COVERAGE w / mutations
totalDels <- 372
hydraPerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim_vlc_mut/hydra_perf.txt', header=TRUE)
cloudbreak <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim_vlc_mut/b292f1f35e6fc997581e19b773ae4dd33beaf58e_r25_perf.txt', header=TRUE)
breakdancerPerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim_vlc_mut/breakdancer_rmdup_perf.txt', header=TRUE)
breakdancerSensitive <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim_vlc_mut/breakdancer_rmdup_sensitive_perf.txt', header=TRUE)
gasvPerf <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr2_sim_vlc_mut/gasv_perf.txt', header=TRUE)

perfsList <- list(hydra=hydraPerf, breakdancer=breakdancerPerf, gasv=gasvPerf, cloudbreak=cloudbreak)
pdf('~/Documents/svpipeline/CHR2SIM_VLC_MUT_ROC.pdf')
plotROC(perfsList, c("Hydra", "Breakdancer", "GASV", "Cloudbreak"), totalDels, "chr2 Simulated (3X)")
dev.off()

# chr17
totalDels <- 199
hydra <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr17_sim/hydra2_perf_gt50.txt', header=TRUE)
breakdancer <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr17_sim/breakdancer2_perf_gt50.txt', header=TRUE)
gasv <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr17_sim/gasv_perf_gt50.txt', header=TRUE)

cloudBreak <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr17_sim/b292f1f35e6fc997581e19b773ae4dd33beaf58e_r25_perf.txt', header=TRUE)
cloudBreak.new <- read.table('~/Documents/gene_rearrange/svpipeline/venter_chr17_sim/test_max_insert_25000_sdseqfilter_t180.perf.txt', header=TRUE)

perfsList <- list(hydra, breakdancer, gasv, cloudbreak, cloudbreak.new)
pdf('~/Documents/svpipeline/CHR17SIM_ROC.pdf')
plotROC(perfsList, c("Hydra", "Breakdancer", "GASV", "Cloudbreak", "Cloudbreak (NEW)"), totalDels, "chr17 Simulated (100X)")
dev.off()

# NA07051
totalDels <- 2548
hydra <- read.table('~/Documents/gene_rearrange/svpipeline/NA07051/hydra_perf.txt', header=TRUE)
breakdancer <- read.table('~/Documents/gene_rearrange/svpipeline/NA07051/breakdancer_perf.txt', header=TRUE)
gasv <- read.table('~/Documents/gene_rearrange/svpipeline/NA07051/gasv_perf.txt', header=TRUE)
cloudbreak <- read.table('~/Documents/gene_rearrange/svpipeline/NA07051/b292f1f35e6fc997581e19b773ae4dd33beaf58e_r25_perf.txt', header=TRUE)
cloudbreak.new <- read.table('~/Documents/gene_rearrange/svpipeline/NA07051/test_max_insert_25000_sdhsd_filter_t120_perf.txt', header=TRUE)

perfsList <- list(hydra, breakdancer, gasv, cloudbreak, cloudbreak.new
pdf('~/Documents/svpipeline/NA07051_ROC_NEW.pdf')
plotROC(perfsList, c("Hydra", "Breakdancer", "GASV", "Cloudbreak", "Cloudbreak (NEW)"), totalDels, "NA07051 (3X)", sim=FALSE)
dev.off()


