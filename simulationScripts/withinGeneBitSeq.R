FinalBitSeqNames <- read.table("../../FinalBitSeqNames.txt",header = TRUE,fill = TRUE)
wtf <- which(FinalBitSeqNames[,2]== "")
FinalBitSeqNames[wtf,2] <- FinalBitSeqNames[wtf,1]


bitseq1 <- read.table("A1.thetaMeans")[,2] 

bitseq2 <- read.table("../A1/transcriptNames.txt")[,2] # this will correspond to the true within gene expression

#bitseq1 <- read.table("A1.thetaMeans")[,2] 
#bitseq2 <- read.table("A2.thetaMeans")[,2] 
#bitseq1 <- read.table("A1.thetaMeans")[,2] 

geneNames <- unique(FinalBitSeqNames[,1])

withinGenePointsBitseq <- array(data = NA,dim = dim(FinalBitSeqNames))

i <- 0
for(g in geneNames){
        transcripts <- which(FinalBitSeqNames[,1] == g)
        l <- length(transcripts)
        i1 <- i + 1
        i2 <- i + l
        i <- i2
        if(l==1){
                withinGenePointsBitseq[i,] <- c(1,1)
        }else{
                withinGenePointsBitseq[i1:i2,] <- cbind(bitseq1[transcripts]/sum(bitseq1[transcripts]),bitseq2[transcripts]/sum(bitseq2[transcripts]))
        }
        #points(withinGenePointsBitseq[i1:i2,],pch = 16,cex = 0.3)
        if(i%%1000==0){print(i)}
}

plot(withinGenePointsBitseq)
write.table(withinGenePointsBitseq,file = "withinGene.txt")
##


bitseq1 <- read.table("A2.thetaMeans")[,2] 

bitseq2 <- read.table("../A2/transcriptNames.txt")[,2] # this will correspond to the true within gene expression

#bitseq1 <- read.table("A1.thetaMeans")[,2] 
#bitseq2 <- read.table("A2.thetaMeans")[,2] 
#bitseq1 <- read.table("A1.thetaMeans")[,2] 

geneNames <- unique(FinalBitSeqNames[,1])

withinGenePointsBitseq <- array(data = NA,dim = dim(FinalBitSeqNames))

i <- 0
for(g in geneNames){
        transcripts <- which(FinalBitSeqNames[,1] == g)
        l <- length(transcripts)
        i1 <- i + 1
        i2 <- i + l
        i <- i2
        if(l==1){
                withinGenePointsBitseq[i,] <- c(1,1)
        }else{
                withinGenePointsBitseq[i1:i2,] <- cbind(bitseq1[transcripts]/sum(bitseq1[transcripts]),bitseq2[transcripts]/sum(bitseq2[transcripts]))
        }
        #points(withinGenePointsBitseq[i1:i2,],pch = 16,cex = 0.3)
        if(i%%1000==0){print(i)}
}

#plot(withinGenePointsBitseq)
write.table(withinGenePointsBitseq,file = "withinGene2.txt")

##

bitseq1 <- read.table("A3.thetaMeans")[,2] 
bitseq2 <- read.table("../A3/transcriptNames.txt")[,2] # this will correspond to the true within gene expression
geneNames <- unique(FinalBitSeqNames[,1])
withinGenePointsBitseq <- array(data = NA,dim = dim(FinalBitSeqNames))
i <- 0
for(g in geneNames){
        transcripts <- which(FinalBitSeqNames[,1] == g)
        l <- length(transcripts)
        i1 <- i + 1
        i2 <- i + l
        i <- i2
        if(l==1){
                withinGenePointsBitseq[i,] <- c(1,1)
        }else{
                withinGenePointsBitseq[i1:i2,] <- cbind(bitseq1[transcripts]/sum(bitseq1[transcripts]),bitseq2[transcripts]/sum(bitseq2[transcripts]))
        }
        if(i%%1000==0){print(i)}
}
#plot(withinGenePointsBitseq)
write.table(withinGenePointsBitseq,file = "withinGene3.txt")
##
bitseq1 <- read.table("A4.thetaMeans")[,2] 
bitseq2 <- read.table("../A4/transcriptNames.txt")[,2] # this will correspond to the true within gene expression
geneNames <- unique(FinalBitSeqNames[,1])
withinGenePointsBitseq <- array(data = NA,dim = dim(FinalBitSeqNames))
i <- 0
for(g in geneNames){
        transcripts <- which(FinalBitSeqNames[,1] == g)
        l <- length(transcripts)
        i1 <- i + 1
        i2 <- i + l
        i <- i2
        if(l==1){
                withinGenePointsBitseq[i,] <- c(1,1)
        }else{
                withinGenePointsBitseq[i1:i2,] <- cbind(bitseq1[transcripts]/sum(bitseq1[transcripts]),bitseq2[transcripts]/sum(bitseq2[transcripts]))
        }
        if(i%%1000==0){print(i)}
}

#plot(withinGenePointsBitseq)
write.table(withinGenePointsBitseq,file = "withinGene4.txt")
##

bitseq1 <- read.table("A5.thetaMeans")[,2] 
bitseq2 <- read.table("../A5/transcriptNames.txt")[,2] # this will correspond to the true within gene expression
geneNames <- unique(FinalBitSeqNames[,1])
withinGenePointsBitseq <- array(data = NA,dim = dim(FinalBitSeqNames))
i <- 0
for(g in geneNames){
        transcripts <- which(FinalBitSeqNames[,1] == g)
        l <- length(transcripts)
        i1 <- i + 1
        i2 <- i + l
        i <- i2
        if(l==1){
                withinGenePointsBitseq[i,] <- c(1,1)
        }else{
                withinGenePointsBitseq[i1:i2,] <- cbind(bitseq1[transcripts]/sum(bitseq1[transcripts]),bitseq2[transcripts]/sum(bitseq2[transcripts]))
        }
        if(i%%1000==0){print(i)}
}
#plot(withinGenePointsBitseq)
write.table(withinGenePointsBitseq,file = "withinGene5.txt")

