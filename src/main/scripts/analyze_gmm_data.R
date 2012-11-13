library(GenomicRanges)
library(rtracklayer)

wdd <- '/Users/cwhelan/Documents/gene_rearrange/svpipeline/venter_chr2_100bp_dip'
testname <- 'newdiscfilt_as175_bt2_nomapq_nofilt_readGroupsBt2scoreMinL01_None_None_2500_25_3_sam_2000_fbbecea_175'

l1 <- read.table(paste(wd, '/', testname, '_l1.wig.gz', sep=""), skip=2, col.names=c("loc", "l1"))
w <- read.table(paste(wd, '/', testname, '_w0.wig.gz', sep=""), skip=2, col.names=c("loc", "w"))
mu <- read.table(paste(wd, '/', testname, '_mu1.wig.gz', sep=""), skip=2, col.names=c("loc", "mu"))
lrHet <- read.table(paste(wd, '/', testname, '_lrHet.wig.gz', sep=""), skip=2, col.names=c("loc", "lrHet"))
lrHom <- read.table(paste(wd, '/', testname, '_lrHom.wig.gz', sep=""), skip=2, col.names=c("loc", "lrHom"))

features <- merge(l1, w, by="loc")
features <- merge(features, mu, by="loc")
features <- merge(features, lrHet, by="loc")
features <- merge(features, lrHom, by="loc")

windows <- GRanges(seqnames="2", ranges=IRanges(start=seq(0,242751150,by=25), end=seq(24, 242751174, by=25)), strand="*")

hap1file <- '/l2/users/whelanch/gene_rearrange/data/jcvi/HuRef.homozygous_indels_inversion.061109.chr2.Deletions.sim_hap1.bed.gz'
hap2file <- '/l2/users/whelanch/gene_rearrange/data/jcvi/HuRef.homozygous_indels_inversion.061109.chr2.Deletions.sim_hap2.bed.gz'

hap1 <- import(hap1file, asRangedData=FALSE)
hap2 <- import(hap2file, asRangedData=FALSE)

elementMetadata(windows)[,"genotype"] <- 0

hap1gt50 <- hap1[end(hap1) - start(hap1) >= 50]
hap2gt50 <- hap2[end(hap2) - start(hap2) >= 50]

elementMetadata(windows)[as.matrix(findOverlaps(hap1gt50, windows))[,2],"genotype"] <- 1
elementMetadata(windows)[as.matrix(findOverlaps(hap2gt50, windows))[,2],"genotype"] <- 2

genotype <- data.frame(loc=start(windows), genotype=elementMetadata(windows)[,"genotype"])


features <- merge(features, genotype, by="loc")

features$sv <- features$genotype != 0

tryit <- function(x,y) {

    #testindex <- 1:(2 * nrow(x)/3)
    testindex=sample(1:dim(x)[1],dim(x)[1]/3)

    xtrain <- x[-testindex,]
    ytrain <- y[-testindex]

    xtrain <- scale(xtrain,center=TRUE,scale=TRUE)

    l <- LiblineaR(as.matrix(xtrain), as.matrix(ytrain), type=1, cost=1.6, verbose=TRUE)

    xtest <- x[testindex,]
    ytest <- y[testindex]
    xtest=scale(xtest,attr(xtrain,"scaled:center"),attr(xtrain,"scaled:scale"))

        pr=FALSE
    p=predict(l,xtest,proba=pr,decisionValues=TRUE)

    res=table(p$predictions,ytest)
    print(res)

    l
}

l <- tryit(features[,2:6], features[,7])

l <- tryit(features[,2:6], features[,8])


library(LiblineaR)


rm(w)
rm(mu)
rm(l1)
rm(lrHet)
rm(lrHom)
rm(windows)

attach(features)

index <- 1:nrow(features)
testindex <- sample(index, trunc(length(index)/3))
testset <- features[testindex,]
trainset <- features[-testindex,]

svm.model <- svm(genotype ~ l1 + w + mu + lrHet + lrHom, data = trainset, cost = 100, gamma = 1)

gc()