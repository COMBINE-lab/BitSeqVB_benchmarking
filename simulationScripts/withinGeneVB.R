FinalBitSeqNames <- read.table("../../FinalBitSeqNames.txt",header = TRUE,fill = TRUE)
geneNames <- unique(FinalBitSeqNames[,1])
vb1 <- read.table("A1.m_alphas")[,1] 

vb2 <- read.table("../A1/transcriptNames.txt")[,2] # this will correspond to the true within gene expression
vb1 <- vb1[-1]/(1-vb1[1])
#vb2 <- vb2/(1-vb2[2])

withinGenePointsVB <- array(data = NA,dim = dim(FinalBitSeqNames))

i <- 0
for(g in geneNames){
        transcripts <- which(FinalBitSeqNames[,1] == g)
        l <- length(transcripts)
        i1 <- i + 1
        i2 <- i + l
        i <- i2
        if(l==1){
                withinGenePointsVB[i,] <- c(1,1)
        }else{
                withinGenePointsVB[i1:i2,] <- cbind(vb1[transcripts]/sum(vb1[transcripts]),vb2[transcripts]/sum(vb2[transcripts]))
        }
        #points(withinGenePointsBitseq[i1:i2,],pch = 16,cex = 0.3)
        if(i%%1000==0){print(i)}
}
maeVB <- sum(abs(apply(withinGenePointsVB,1,diff)))/dim(withinGenePointsVB)[1]
#plot(withinGenePointsVB,col = 2,pch = 16,cex = 0.4,main = paste("VB: MAD = ",round(maeVB,3),sep = ""),xlab = "replicate 1",ylab = "true")
write.table(withinGenePointsVB,file = "withinGene.txt")



geneNames <- unique(FinalBitSeqNames[,1])
vb1 <- read.table("A2.m_alphas")[,1] 

vb2 <- read.table("../A2/transcriptNames.txt")[,2] # this will correspond to the true within gene expression
vb1 <- vb1[-1]/(1-vb1[1])
#vb2 <- vb2/(1-vb2[2])

withinGenePointsVB <- array(data = NA,dim = dim(FinalBitSeqNames))

i <- 0
for(g in geneNames){
        transcripts <- which(FinalBitSeqNames[,1] == g)
        l <- length(transcripts)
        i1 <- i + 1
        i2 <- i + l
        i <- i2
        if(l==1){
                withinGenePointsVB[i,] <- c(1,1)
        }else{
                withinGenePointsVB[i1:i2,] <- cbind(vb1[transcripts]/sum(vb1[transcripts]),vb2[transcripts]/sum(vb2[transcripts]))
        }
        #points(withinGenePointsBitseq[i1:i2,],pch = 16,cex = 0.3)
        if(i%%1000==0){print(i)}
}
maeVB <- sum(abs(apply(withinGenePointsVB,1,diff)))/dim(withinGenePointsVB)[1]
#plot(withinGenePointsVB,col = 2,pch = 16,cex = 0.4,main = paste("VB: MAD = ",round(maeVB,3),sep = ""),xlab = "replicate 1",ylab = "true")
write.table(withinGenePointsVB,file = "withinGene2.txt")
###################
###################
geneNames <- unique(FinalBitSeqNames[,1])
vb1 <- read.table("A3.m_alphas")[,1] 
vb2 <- read.table("../A3/transcriptNames.txt")[,2] # this will correspond to the true within gene expression
vb1 <- vb1[-1]/(1-vb1[1])
withinGenePointsVB <- array(data = NA,dim = dim(FinalBitSeqNames))
i <- 0
for(g in geneNames){
        transcripts <- which(FinalBitSeqNames[,1] == g)
        l <- length(transcripts)
        i1 <- i + 1
        i2 <- i + l
        i <- i2
        if(l==1){
                withinGenePointsVB[i,] <- c(1,1)
        }else{
                withinGenePointsVB[i1:i2,] <- cbind(vb1[transcripts]/sum(vb1[transcripts]),vb2[transcripts]/sum(vb2[transcripts]))
        }
        #points(withinGenePointsBitseq[i1:i2,],pch = 16,cex = 0.3)
        if(i%%1000==0){print(i)}
}
maeVB <- sum(abs(apply(withinGenePointsVB,1,diff)))/dim(withinGenePointsVB)[1]
#plot(withinGenePointsVB,col = 2,pch = 16,cex = 0.4,main = paste("VB: MAD = ",round(maeVB,3),sep = ""),xlab = "replicate 1",ylab = "true")
write.table(withinGenePointsVB,file = "withinGene3.txt")
###################
###################
geneNames <- unique(FinalBitSeqNames[,1])
vb1 <- read.table("A4.m_alphas")[,1] 
vb2 <- read.table("../A4/transcriptNames.txt")[,2] # this will correspond to the true within gene expression
vb1 <- vb1[-1]/(1-vb1[1])
withinGenePointsVB <- array(data = NA,dim = dim(FinalBitSeqNames))
i <- 0
for(g in geneNames){
        transcripts <- which(FinalBitSeqNames[,1] == g)
        l <- length(transcripts)
        i1 <- i + 1
        i2 <- i + l
        i <- i2
        if(l==1){
                withinGenePointsVB[i,] <- c(1,1)
        }else{
                withinGenePointsVB[i1:i2,] <- cbind(vb1[transcripts]/sum(vb1[transcripts]),vb2[transcripts]/sum(vb2[transcripts]))
        }
        #points(withinGenePointsBitseq[i1:i2,],pch = 16,cex = 0.3)
        if(i%%1000==0){print(i)}
}
maeVB <- sum(abs(apply(withinGenePointsVB,1,diff)))/dim(withinGenePointsVB)[1]
#plot(withinGenePointsVB,col = 2,pch = 16,cex = 0.4,main = paste("VB: MAD = ",round(maeVB,3),sep = ""),xlab = "replicate 1",ylab = "true")
write.table(withinGenePointsVB,file = "withinGene4.txt")
###################
###################
geneNames <- unique(FinalBitSeqNames[,1])
vb1 <- read.table("A5.m_alphas")[,1] 
vb2 <- read.table("../A5/transcriptNames.txt")[,2] # this will correspond to the true within gene expression
vb1 <- vb1[-1]/(1-vb1[1])
withinGenePointsVB <- array(data = NA,dim = dim(FinalBitSeqNames))
i <- 0
for(g in geneNames){
        transcripts <- which(FinalBitSeqNames[,1] == g)
        l <- length(transcripts)
        i1 <- i + 1
        i2 <- i + l
        i <- i2
        if(l==1){
                withinGenePointsVB[i,] <- c(1,1)
        }else{
                withinGenePointsVB[i1:i2,] <- cbind(vb1[transcripts]/sum(vb1[transcripts]),vb2[transcripts]/sum(vb2[transcripts]))
        }
        #points(withinGenePointsBitseq[i1:i2,],pch = 16,cex = 0.3)
        if(i%%1000==0){print(i)}
}
maeVB <- sum(abs(apply(withinGenePointsVB,1,diff)))/dim(withinGenePointsVB)[1]
#plot(withinGenePointsVB,col = 2,pch = 16,cex = 0.4,main = paste("VB: MAD = ",round(maeVB,3),sep = ""),xlab = "replicate 1",ylab = "true")
write.table(withinGenePointsVB,file = "withinGene5.txt")





