
withinGene <- read.table("bitseq/withinGene.txt")
aoua1 <- abs(apply(withinGene,1,diff))
Replicate1 <- withinGene[,1]
withinGene <- read.table("bitseq/withinGene2.txt")
aoua2 <- abs(apply(withinGene,1,diff))
Replicate2 <- withinGene[,1]
withinGene <- read.table("bitseq/withinGene3.txt")
aoua3 <- abs(apply(withinGene,1,diff))
Replicate3 <- withinGene[,1]
withinGene <- read.table("bitseq/withinGene4.txt")
aoua4 <- abs(apply(withinGene,1,diff))
Replicate4 <- withinGene[,1]
withinGene <- read.table("bitseq/withinGene5.txt")
aoua5 <- abs(apply(withinGene,1,diff))
Replicate5 <- withinGene[,1]


withinGene <- read.table("rsem/withinGene.txt")
aoua1 <- cbind(aoua1,abs(apply(withinGene,1,diff)))
Replicate1 <- cbind(Replicate1,withinGene[,1])
withinGene <- read.table("rsem/withinGene2.txt")
aoua2 <- cbind(aoua2,abs(apply(withinGene,1,diff)))
Replicate2 <- cbind(Replicate2,withinGene[,1])
withinGene <- read.table("rsem/withinGene3.txt")
aoua3 <- cbind(aoua3,abs(apply(withinGene,1,diff)))
Replicate3 <- cbind(Replicate3,withinGene[,1])
withinGene <- read.table("rsem/withinGene4.txt")
aoua4 <- cbind(aoua4,abs(apply(withinGene,1,diff)))
Replicate4 <- cbind(Replicate4,withinGene[,1])
withinGene <- read.table("rsem/withinGene5.txt")
aoua5 <- cbind(aoua5,abs(apply(withinGene,1,diff)))
Replicate5 <- cbind(Replicate5,withinGene[,1])

withinGene <- read.table("sailfish/withinGene.txt")
aoua1 <- cbind(aoua1,abs(apply(withinGene,1,diff)))
Replicate1 <- cbind(Replicate1,withinGene[,1])
withinGene <- read.table("sailfish/withinGene2.txt")
aoua2 <- cbind(aoua2,abs(apply(withinGene,1,diff)))
Replicate2 <- cbind(Replicate2,withinGene[,1])
withinGene <- read.table("sailfish/withinGene3.txt")
aoua3 <- cbind(aoua3,abs(apply(withinGene,1,diff)))
Replicate3 <- cbind(Replicate3,withinGene[,1])
withinGene <- read.table("sailfish/withinGene4.txt")
aoua4 <- cbind(aoua4,abs(apply(withinGene,1,diff)))
Replicate4 <- cbind(Replicate4,withinGene[,1])
withinGene <- read.table("sailfish/withinGene4.txt")
aoua5 <- cbind(aoua5,abs(apply(withinGene,1,diff)))
Replicate5 <- cbind(Replicate5,withinGene[,1])

withinGene <- read.table("tigar/withinGene.txt")
aoua1 <- cbind(aoua1,abs(apply(withinGene,1,diff)))
Replicate1 <- cbind(Replicate1,withinGene[,1])
withinGene <- read.table("tigar/withinGene2.txt")
aoua2 <- cbind(aoua2,abs(apply(withinGene,1,diff)))
Replicate2 <- cbind(Replicate2,withinGene[,1])
withinGene <- read.table("tigar/withinGene3.txt")
aoua3 <- cbind(aoua3,abs(apply(withinGene,1,diff)))
Replicate3 <- cbind(Replicate3,withinGene[,1])
withinGene <- read.table("tigar/withinGene4.txt")
aoua4 <- cbind(aoua4,abs(apply(withinGene,1,diff)))
Replicate4 <- cbind(Replicate4,withinGene[,1])
withinGene <- read.table("tigar/withinGene5.txt")
aoua5 <- cbind(aoua5,abs(apply(withinGene,1,diff)))
Replicate5 <- cbind(Replicate5,withinGene[,1])

withinGene <- read.table("vb/withinGene.txt")
aoua1 <- cbind(aoua1,abs(apply(withinGene,1,diff)))
Replicate1 <- cbind(Replicate1,withinGene[,1])
withinGene <- read.table("vb/withinGene2.txt")
aoua2 <- cbind(aoua2,abs(apply(withinGene,1,diff)))
Replicate2 <- cbind(Replicate2,withinGene[,1])
withinGene <- read.table("vb/withinGene3.txt")
aoua3 <- cbind(aoua3,abs(apply(withinGene,1,diff)))
Replicate3 <- cbind(Replicate3,withinGene[,1])
withinGene <- read.table("vb/withinGene4.txt")
aoua4 <- cbind(aoua4,abs(apply(withinGene,1,diff)))
Replicate4 <- cbind(Replicate4,withinGene[,1])
withinGene <- read.table("vb/withinGene5.txt")
aoua5 <- cbind(aoua5,abs(apply(withinGene,1,diff)))
Replicate5 <- cbind(Replicate5,withinGene[,1])

withinGene <- read.table("xpress/withinGene.txt")
aoua1 <- cbind(aoua1,abs(apply(withinGene,1,diff)))
Replicate1 <- cbind(Replicate1,withinGene[,1])
withinGene <- read.table("xpress/withinGene2.txt")
aoua2 <- cbind(aoua2,abs(apply(withinGene,1,diff)))
Replicate2 <- cbind(Replicate2,withinGene[,1])
withinGene <- read.table("xpress/withinGene3.txt")
aoua3 <- cbind(aoua3,abs(apply(withinGene,1,diff)))
Replicate3 <- cbind(Replicate3,withinGene[,1])
withinGene <- read.table("xpress/withinGene4.txt")
aoua4 <- cbind(aoua4,abs(apply(withinGene,1,diff)))
Replicate4 <- cbind(Replicate4,withinGene[,1])
withinGene <- read.table("xpress/withinGene5.txt")
aoua5 <- cbind(aoua5,abs(apply(withinGene,1,diff)))
Replicate5 <- cbind(Replicate5,withinGene[,1])

withinGene <- read.table("cufflinks/withinGene.txt")
aoua1 <- cbind(aoua1,abs(apply(withinGene,1,diff)))
Replicate1 <- cbind(Replicate1,withinGene[,1])
withinGene <- read.table("cufflinks/withinGene2.txt")
aoua2 <- cbind(aoua2,abs(apply(withinGene,1,diff)))
Replicate2 <- cbind(Replicate2,withinGene[,1])
withinGene <- read.table("cufflinks/withinGene3.txt")
aoua3 <- cbind(aoua3,abs(apply(withinGene,1,diff)))
Replicate3 <- cbind(Replicate3,withinGene[,1])
withinGene <- read.table("cufflinks/withinGene4.txt")
aoua4 <- cbind(aoua4,abs(apply(withinGene,1,diff)))
Replicate4 <- cbind(Replicate4,withinGene[,1])
withinGene <- read.table("cufflinks/withinGene5.txt")
aoua5 <- cbind(aoua5,abs(apply(withinGene,1,diff)))
Replicate5 <- cbind(Replicate5,withinGene[,1])

withinGene <- read.table("casper/withinGene.txt")
aoua1 <- cbind(aoua1,abs(apply(withinGene,1,diff)))
Replicate1 <- cbind(Replicate1,withinGene[,1])
withinGene <- read.table("casper/withinGene2.txt")
aoua2 <- cbind(aoua2,abs(apply(withinGene,1,diff)))
Replicate2 <- cbind(Replicate2,withinGene[,1])
withinGene <- read.table("casper/withinGene3.txt")
aoua3 <- cbind(aoua3,abs(apply(withinGene,1,diff)))
Replicate3 <- cbind(Replicate3,withinGene[,1])
withinGene <- read.table("casper/withinGene4.txt")
aoua4 <- cbind(aoua4,abs(apply(withinGene,1,diff)))
Replicate4 <- cbind(Replicate4,withinGene[,1])
withinGene <- read.table("casper/withinGene5.txt")
aoua5 <- cbind(aoua5,abs(apply(withinGene,1,diff)))
Replicate5 <- cbind(Replicate5,withinGene[,1])

colnames(aoua1) <- colnames(Replicate1) <- colnames(aoua2) <- colnames(Replicate2) <-colnames(aoua3) <- colnames(Replicate3) <-colnames(aoua4) <- colnames(Replicate4) <-colnames(aoua5) <- colnames(Replicate5) <-c("bitseqMCMC","rsem","sailfish","tigar","bitseqVB","eXpress","cufflinks","casper")
aoua <- (aoua1 + aoua2 + aoua3 + aoua4 + aoua5)/5

#boxplot(aoua)
nah <- apply(aoua,2,function(y)mean(y,na.rm = TRUE))
perm <- order(nah)
write.table(nah,"MeanWgeTrue.txt")
#######################################################################################################
#######################################################################################################

# interReplicate consistency plot
interReplicate <- abs(Replicate1 - Replicate2)
interReplicate <- interReplicate + abs(Replicate1 - Replicate3) + abs(Replicate1 - Replicate4) + abs(Replicate1 - Replicate5) + abs(Replicate2 - Replicate3) + abs(Replicate2 - Replicate4) + abs(Replicate2 - Replicate5) + abs(Replicate3 - Replicate4) + abs(Replicate3 - Replicate5) + abs(Replicate4 - Replicate5)
interReplicate <- interReplicate/10

nah <- apply(interReplicate,2,function(y)mean(y,na.rm = TRUE))
perm <- order(nah)
write.table(nah,"MeanWgeInterReplicate.txt")
##################

#perm <- c(1,5,2,8,4,3,7,6)
pdf(file = "spanki1vs2.pdf",width = 12,height = 6,pointsize = 15)
par(mfrow=c(2,4),mar = c(4,4,2,1))
for (k in 1:8){
#	plot(Replicate1[,k],Replicate2[,k],main = colnames(Replicate2)[k],xlab = "Replicate 1",ylab = "Replicate 2",pch = 16,cex = 0.5,col = 2)
	egkira <- length(which(is.na(abs(Replicate1[,perm[k]]-Replicate2[,perm[k]]))==FALSE))
	methodMAE = round(sum(abs(Replicate1[,perm[k]]-Replicate2[,perm[k]]),na.rm=TRUE)/egkira,4)
	smoothScatter(Replicate1[,perm[k]],Replicate2[,perm[k]],main = paste(colnames(Replicate2)[perm[k]],": MAE = ",methodMAE,sep=""),xlab = "Replicate 1",ylab = "Replicate 2",pch = 16,cex = 0.5,col = 2,nrpoints = 500)
}
dev.off()


pdf(file = "spanki1vsTrue.pdf",width = 12,height = 6,pointsize = 15)
par(mfrow=c(2,4),mar = c(4,4,2,1))
myDir <- c(
"bitseq/withinGene.txt",
"rsem/withinGene.txt",
"sailfish/withinGene.txt",
"tigar/withinGene.txt",
"vb/withinGene.txt",
"xpress/withinGene.txt",
"cufflinks/withinGene.txt",
"casper/withinGene.txt")
for (k in 1:8){
	withinGene <- read.table(myDir[perm[k]])
	egkira <- length(which(is.na(abs(apply(withinGene,1,diff)))==FALSE))
	methodMAE = round(sum(abs(apply(withinGene,1,diff)),na.rm=TRUE)/egkira,4)
	smoothScatter(withinGene,main = paste(colnames(Replicate2)[perm[k]],": MAE = ",methodMAE,sep=""),xlab = "Replicate 1",ylab = "true Value",pch = 16,cex = 0.5,col = 2,nrpoints = 500)
}
dev.off()
























