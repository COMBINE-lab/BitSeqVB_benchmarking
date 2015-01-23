#source("http://bioconductor.org/biocLite.R")
#biocLite("casper")
library(GenomicFeatures)
genome='hg19'
library(Rsamtools)
library(casper)
library(parallel)
mcr <- 4




#To generate a transcriptome from a gtf le use the following code, changing the file name for your gtf file (including the full path).
library('rtracklayer')

genDB<-makeTranscriptDbFromUCSC(genome=genome, tablename="refGene")
hg19DB <- procGenome(genDB=genDB, genome=genome, mc.cores=mcr)

setwd("../tophat/A1.tophat.out")
bamFile <- 'accepted_hits.sorted.bam'
rep1 <- wrapKnown(bamFile, verbose=TRUE, genomeDB=hg19DB, readLength=76, rpkm=FALSE, priorq=2, mc.cores.int=mcr, mc.cores=mcr)
save(rep1,file = "expCasper.RData")

setwd("../A2.tophat.out")
bamFile <- 'accepted_hits.sorted.bam'
rep1 <- wrapKnown(bamFile, verbose=TRUE, genomeDB=hg19DB, readLength=76, rpkm=FALSE, priorq=2, mc.cores.int=mcr, mc.cores=mcr)
save(rep1,file = "expCasper.RData")

setwd("../A3.tophat.out")
bamFile <- 'accepted_hits.sorted.bam'
rep1 <- wrapKnown(bamFile, verbose=TRUE, genomeDB=hg19DB, readLength=76, rpkm=FALSE, priorq=2, mc.cores.int=mcr, mc.cores=mcr)
save(rep1,file = "expCasper.RData")

setwd("../A4.tophat.out")
bamFile <- 'accepted_hits.sorted.bam'
rep1 <- wrapKnown(bamFile, verbose=TRUE, genomeDB=hg19DB, readLength=76, rpkm=FALSE, priorq=2, mc.cores.int=mcr, mc.cores=mcr)
save(rep1,file = "expCasper.RData")

setwd("../A5.tophat.out")
bamFile <- 'accepted_hits.sorted.bam'
rep1 <- wrapKnown(bamFile, verbose=TRUE, genomeDB=hg19DB, readLength=76, rpkm=FALSE, priorq=2, mc.cores.int=mcr, mc.cores=mcr)
save(rep1,file = "expCasper.RData")

