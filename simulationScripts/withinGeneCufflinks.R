FinalBitSeqNames <- read.table("../../FinalBitSeqNames.txt",header = TRUE,fill = TRUE)


bitseq1 <- read.table("cuffA1.txt")[,3] 
nams <- read.table("cuffA1.txt")[,1]
bitseq2 <- read.table("../A1/transcriptNames.txt")[,2] # this will correspond to the true within gene expression
bitseq1 <- exp(bitseq1)

bitseq1 <- bitseq1/sum(bitseq1)
bitseq2 <- bitseq2/sum(bitseq2)

geneNames <- unique(FinalBitSeqNames[,1])



withinGenePointsBitseq <- array(data = NA,dim = dim(FinalBitSeqNames))

i <- 0
for(g in geneNames){
        transcripts <- which(FinalBitSeqNames[,1] == g)
        transcriptsID <- which(is.na(match(nams,FinalBitSeqNames[transcripts,2]))==FALSE)
        l <- length(transcriptsID)
        i1 <- i + 1
        i2 <- i + l
        i <- i2
	if(l>0){
        if(l==1){
                withinGenePointsBitseq[i,] <- c(1,1)
        }else{
                withinGenePointsBitseq[i1:i2,] <- cbind(bitseq1[transcriptsID]/sum(bitseq1[transcriptsID]),bitseq2[transcriptsID]/sum(bitseq2[transcriptsID]))
        }
	}
        #points(withinGenePointsBitseq[i1:i2,],pch = 16,cex = 0.3)
        if(i%%1000==0){print(i)}
}

#plot(withinGenePointsBitseq)
write.table(withinGenePointsBitseq,file = "withinGene.txt")



bitseq1 <- read.table("cuffA2.txt")[,3] 
nams <- read.table("cuffA2.txt")[,1]
bitseq2 <- read.table("../A2/transcriptNames.txt")[,2] # this will correspond to the true within gene expression
bitseq1 <- exp(bitseq1)

bitseq1 <- bitseq1/sum(bitseq1)
bitseq2 <- bitseq2/sum(bitseq2)

geneNames <- unique(FinalBitSeqNames[,1])



withinGenePointsBitseq <- array(data = NA,dim = dim(FinalBitSeqNames))

i <- 0
for(g in geneNames){
        transcripts <- which(FinalBitSeqNames[,1] == g)
        transcriptsID <- which(is.na(match(nams,FinalBitSeqNames[transcripts,2]))==FALSE)
        l <- length(transcriptsID)
        i1 <- i + 1
        i2 <- i + l
        i <- i2
	if(l>0){
        if(l==1){
                withinGenePointsBitseq[i,] <- c(1,1)
        }else{
                withinGenePointsBitseq[i1:i2,] <- cbind(bitseq1[transcriptsID]/sum(bitseq1[transcriptsID]),bitseq2[transcriptsID]/sum(bitseq2[transcriptsID]))
        }
	}
        #points(withinGenePointsBitseq[i1:i2,],pch = 16,cex = 0.3)
        if(i%%1000==0){print(i)}
}

#plot(withinGenePointsBitseq)
write.table(withinGenePointsBitseq,file = "withinGene2.txt")
####################
###################
bitseq1 <- read.table("cuffA3.txt")[,3] 
nams <- read.table("cuffA3.txt")[,1]
bitseq2 <- read.table("../A3/transcriptNames.txt")[,2] # this will correspond to the true within gene expression
bitseq1 <- exp(bitseq1)
bitseq1 <- bitseq1/sum(bitseq1)
bitseq2 <- bitseq2/sum(bitseq2)
geneNames <- unique(FinalBitSeqNames[,1])
withinGenePointsBitseq <- array(data = NA,dim = dim(FinalBitSeqNames))
i <- 0
for(g in geneNames){
        transcripts <- which(FinalBitSeqNames[,1] == g)
        transcriptsID <- which(is.na(match(nams,FinalBitSeqNames[transcripts,2]))==FALSE)
        l <- length(transcriptsID)
        i1 <- i + 1
        i2 <- i + l
        i <- i2
	if(l>0){
        if(l==1){
                withinGenePointsBitseq[i,] <- c(1,1)
        }else{
                withinGenePointsBitseq[i1:i2,] <- cbind(bitseq1[transcriptsID]/sum(bitseq1[transcriptsID]),bitseq2[transcriptsID]/sum(bitseq2[transcriptsID]))
        }
	}
        if(i%%1000==0){print(i)}
}
#plot(withinGenePointsBitseq)
write.table(withinGenePointsBitseq,file = "withinGene3.txt")
####################
###################
bitseq1 <- read.table("cuffA4.txt")[,3] 
nams <- read.table("cuffA4.txt")[,1]
bitseq2 <- read.table("../A4/transcriptNames.txt")[,2] # this will correspond to the true within gene expression
bitseq1 <- exp(bitseq1)
bitseq1 <- bitseq1/sum(bitseq1)
bitseq2 <- bitseq2/sum(bitseq2)
geneNames <- unique(FinalBitSeqNames[,1])
withinGenePointsBitseq <- array(data = NA,dim = dim(FinalBitSeqNames))
i <- 0
for(g in geneNames){
        transcripts <- which(FinalBitSeqNames[,1] == g)
        transcriptsID <- which(is.na(match(nams,FinalBitSeqNames[transcripts,2]))==FALSE)
        l <- length(transcriptsID)
        i1 <- i + 1
        i2 <- i + l
        i <- i2
	if(l>0){
        if(l==1){
                withinGenePointsBitseq[i,] <- c(1,1)
        }else{
                withinGenePointsBitseq[i1:i2,] <- cbind(bitseq1[transcriptsID]/sum(bitseq1[transcriptsID]),bitseq2[transcriptsID]/sum(bitseq2[transcriptsID]))
        }
	}
        if(i%%1000==0){print(i)}
}
#plot(withinGenePointsBitseq)
write.table(withinGenePointsBitseq,file = "withinGene4.txt")
####################
###################
bitseq1 <- read.table("cuffA5.txt")[,3] 
nams <- read.table("cuffA5.txt")[,1]
bitseq2 <- read.table("../A5/transcriptNames.txt")[,2] # this will correspond to the true within gene expression
bitseq1 <- exp(bitseq1)
bitseq1 <- bitseq1/sum(bitseq1)
bitseq2 <- bitseq2/sum(bitseq2)
geneNames <- unique(FinalBitSeqNames[,1])
withinGenePointsBitseq <- array(data = NA,dim = dim(FinalBitSeqNames))
i <- 0
for(g in geneNames){
        transcripts <- which(FinalBitSeqNames[,1] == g)
        transcriptsID <- which(is.na(match(nams,FinalBitSeqNames[transcripts,2]))==FALSE)
        l <- length(transcriptsID)
        i1 <- i + 1
        i2 <- i + l
        i <- i2
	if(l>0){
        if(l==1){
                withinGenePointsBitseq[i,] <- c(1,1)
        }else{
                withinGenePointsBitseq[i1:i2,] <- cbind(bitseq1[transcriptsID]/sum(bitseq1[transcriptsID]),bitseq2[transcriptsID]/sum(bitseq2[transcriptsID]))
        }
	}
        if(i%%1000==0){print(i)}
}
#plot(withinGenePointsBitseq)
write.table(withinGenePointsBitseq,file = "withinGene5.txt")


