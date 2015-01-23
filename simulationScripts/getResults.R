myResults <- read.table("sailfish/sailA1.txt")
colnames(myResults)[3] <- "sailfish"
myResults <- cbind(myResults,read.table("xpress/xpressA1.txt")[,3])
colnames(myResults)[4] <- "eXpress"
myResults <- cbind(myResults,read.table("rsem/rsemA1.txt")[,3])
colnames(myResults)[5] <- "rsem"
myResults <- cbind(myResults,read.table("cufflinks/cuffA1.txt")[,3])
colnames(myResults)[6] <- "cufflinks"
myResults <- cbind(myResults,read.table("casper/sailA1.txt")[,3])
colnames(myResults)[7] <- "casper"
myResults <- cbind(myResults,log(read.table("bitseq/A1.thetaMeans")[,2]))
colnames(myResults)[8] <- "bitseq MCMC"
vb <- read.table("vb/A1.m_alphas")
vb[,1] <- vb[,1]/(1-vb[1,1])
vb <- log(vb[-1,1])
myResults <- cbind(myResults,vb)
colnames(myResults)[9] <- "bitseq VB"
myResults <- cbind(myResults,read.table("tigar/tigarA1.txt")[,3])
colnames(myResults)[10] <- "tigar"

ind <- which(is.infinite(myResults[,2])==FALSE)
apply(myResults[ind,-c(1,2)],2,function(y){sum(abs(y-myResults[ind,2]))})
for(k in 3:10){
ind1 <- which(is.infinite(myResults[,k])==FALSE)
ind <- intersect(ind,ind1)
print(length(ind))
}
apply(myResults[ind,-c(1,2)],2,function(y){sum(abs(y-myResults[ind,2]))})
aoua <- abs(myResults[ind,-c(1,2)] - myResults[ind,2])
boxplot(log(aoua))

aoua <- myResults[ind,-c(1,2)] - myResults[ind,2]
aoua1 <- exp(myResults[ind,-c(1,2)]) - exp(myResults[ind,2])
aoua2 <- log(myResults[ind,-c(1,2)]/myResults[ind,2])

#boxplot(aoua,ylim = c(-1,1))
#boxplot(abs(aoua),ylim = c(0,1))
#apply(abs(aoua),2,sum)

#aoua <- log(myResults[ind,-c(1,2)]/myResults[ind,2])
#boxplot(aoua)

mae1 <- apply(abs(aoua),2,sum)/dim(aoua)[1]
maExpe1 <- apply(abs(aoua1),2,sum)/dim(aoua1)[1]
MLR1 <- apply(abs(aoua2),2,sum)/dim(aoua1)[1]
perm <- order(mae1)
barplot(mae1[perm],ylab = "MAE")

####################################################################

myResults <- read.table("sailfish/sailA2.txt")
colnames(myResults)[3] <- "sailfish"
myResults <- cbind(myResults,read.table("xpress/xpressA2.txt")[,3])
colnames(myResults)[4] <- "eXpress"
myResults <- cbind(myResults,read.table("rsem/rsemA2.txt")[,3])
colnames(myResults)[5] <- "rsem"
myResults <- cbind(myResults,read.table("cufflinks/cuffA2.txt")[,3])
colnames(myResults)[6] <- "cufflinks"
myResults <- cbind(myResults,read.table("casper/sailA2.txt")[,3])
colnames(myResults)[7] <- "casper"
myResults <- cbind(myResults,log(read.table("bitseq/A2.thetaMeans")[,2]))
colnames(myResults)[8] <- "bitseq MCMC"
vb <- read.table("vb/A2.m_alphas")
vb[,1] <- vb[,1]/(1-vb[1,1])
vb <- log(vb[-1,1])
myResults <- cbind(myResults,vb)
colnames(myResults)[9] <- "bitseq VB"
myResults <- cbind(myResults,read.table("tigar/tigarA2.txt")[,3])
colnames(myResults)[10] <- "tigar"


ind <- which(is.infinite(myResults[,2])==FALSE)
apply(myResults[ind,-c(1,2)],2,function(y){sum(abs(y-myResults[ind,2]))})
for(k in 3:10){
ind1 <- which(is.infinite(myResults[,k])==FALSE)
ind <- intersect(ind,ind1)
print(length(ind))
}
apply(myResults[ind,-c(1,2)],2,function(y){sum(abs(y-myResults[ind,2]))})
aoua <- abs(myResults[ind,-c(1,2)] - myResults[ind,2])
boxplot(log(aoua))

aoua <- myResults[ind,-c(1,2)] - myResults[ind,2]
aoua1 <- exp(myResults[ind,-c(1,2)]) - exp(myResults[ind,2])
aoua2 <- log(myResults[ind,-c(1,2)]/myResults[ind,2])
#boxplot(aoua,ylim = c(-1,1))
#boxplot(abs(aoua),ylim = c(0,1))
#apply(abs(aoua),2,sum)

#aoua <- log(myResults[ind,-c(1,2)]/myResults[ind,2])
#boxplot(aoua)

mae2 <- apply(abs(aoua),2,sum)/dim(aoua)[1]
maExpe2 <- apply(abs(aoua1),2,sum)/dim(aoua1)[1]
MLR2 <- apply(abs(aoua2),2,sum)/dim(aoua1)[1]
perm <- order(mae2)
barplot(mae2[perm],ylab = "MAE")
####################################################################

myResults <- read.table("sailfish/sailA3.txt")
colnames(myResults)[3] <- "sailfish"
myResults <- cbind(myResults,read.table("xpress/xpressA3.txt")[,3])
colnames(myResults)[4] <- "eXpress"
myResults <- cbind(myResults,read.table("rsem/rsemA3.txt")[,3])
colnames(myResults)[5] <- "rsem"
myResults <- cbind(myResults,read.table("cufflinks/cuffA3.txt")[,3])
colnames(myResults)[6] <- "cufflinks"
myResults <- cbind(myResults,read.table("casper/sailA3.txt")[,3])
colnames(myResults)[7] <- "casper"
myResults <- cbind(myResults,log(read.table("bitseq/A3.thetaMeans")[,2]))
colnames(myResults)[8] <- "bitseq MCMC"
vb <- read.table("vb/A3.m_alphas")
vb[,1] <- vb[,1]/(1-vb[1,1])
vb <- log(vb[-1,1])
myResults <- cbind(myResults,vb)
colnames(myResults)[9] <- "bitseq VB"
myResults <- cbind(myResults,read.table("tigar/tigarA3.txt")[,3])
colnames(myResults)[10] <- "tigar"


ind <- which(is.infinite(myResults[,2])==FALSE)
apply(myResults[ind,-c(1,2)],2,function(y){sum(abs(y-myResults[ind,2]))})
for(k in 3:10){
ind1 <- which(is.infinite(myResults[,k])==FALSE)
ind <- intersect(ind,ind1)
print(length(ind))
}
apply(myResults[ind,-c(1,2)],2,function(y){sum(abs(y-myResults[ind,2]))})
aoua <- abs(myResults[ind,-c(1,2)] - myResults[ind,2])
boxplot(log(aoua))

aoua <- myResults[ind,-c(1,2)] - myResults[ind,2]
aoua1 <- exp(myResults[ind,-c(1,2)]) - exp(myResults[ind,2])
aoua2 <- log(myResults[ind,-c(1,2)]/myResults[ind,2])
#boxplot(aoua,ylim = c(-1,1))
#boxplot(abs(aoua),ylim = c(0,1))
#apply(abs(aoua),2,sum)

#aoua <- log(myResults[ind,-c(1,2)]/myResults[ind,2])
#boxplot(aoua)

mae3 <- apply(abs(aoua),2,sum)/dim(aoua)[1]
maExpe3 <- apply(abs(aoua1),2,sum)/dim(aoua1)[1]
MLR3 <- apply(abs(aoua2),2,sum)/dim(aoua1)[1]
perm <- order(mae3)
barplot(mae3[perm],ylab = "MAE")
####################################################################

myResults <- read.table("sailfish/sailA4.txt")
colnames(myResults)[3] <- "sailfish"
myResults <- cbind(myResults,read.table("xpress/xpressA4.txt")[,3])
colnames(myResults)[4] <- "eXpress"
myResults <- cbind(myResults,read.table("rsem/rsemA4.txt")[,3])
colnames(myResults)[5] <- "rsem"
myResults <- cbind(myResults,read.table("cufflinks/cuffA4.txt")[,3])
colnames(myResults)[6] <- "cufflinks"
myResults <- cbind(myResults,read.table("casper/sailA4.txt")[,3])
colnames(myResults)[7] <- "casper"
myResults <- cbind(myResults,log(read.table("bitseq/A4.thetaMeans")[,2]))
colnames(myResults)[8] <- "bitseq MCMC"
vb <- read.table("vb/A4.m_alphas")
vb[,1] <- vb[,1]/(1-vb[1,1])
vb <- log(vb[-1,1])
myResults <- cbind(myResults,vb)
colnames(myResults)[9] <- "bitseq VB"
myResults <- cbind(myResults,read.table("tigar/tigarA4.txt")[,3])
colnames(myResults)[10] <- "tigar"


ind <- which(is.infinite(myResults[,2])==FALSE)
apply(myResults[ind,-c(1,2)],2,function(y){sum(abs(y-myResults[ind,2]))})
for(k in 3:10){
ind1 <- which(is.infinite(myResults[,k])==FALSE)
ind <- intersect(ind,ind1)
print(length(ind))
}
apply(myResults[ind,-c(1,2)],2,function(y){sum(abs(y-myResults[ind,2]))})
aoua <- abs(myResults[ind,-c(1,2)] - myResults[ind,2])
boxplot(log(aoua))

aoua <- myResults[ind,-c(1,2)] - myResults[ind,2]
aoua1 <- exp(myResults[ind,-c(1,2)]) - exp(myResults[ind,2])
aoua2 <- log(myResults[ind,-c(1,2)]/myResults[ind,2])
#boxplot(aoua,ylim = c(-1,1))
#boxplot(abs(aoua),ylim = c(0,1))
apply(abs(aoua),2,sum)

#aoua <- log(myResults[ind,-c(1,2)]/myResults[ind,2])
#boxplot(aoua)

mae4 <- apply(abs(aoua),2,sum)/dim(aoua)[1]
maExpe4 <- apply(abs(aoua1),2,sum)/dim(aoua1)[1]
MLR4 <- apply(abs(aoua2),2,sum)/dim(aoua1)[1]
perm <- order(mae4)
barplot(mae4[perm],ylab = "MAE")
####################################################################

myResults <- read.table("sailfish/sailA5.txt")
colnames(myResults)[3] <- "sailfish"
myResults <- cbind(myResults,read.table("xpress/xpressA5.txt")[,3])
colnames(myResults)[4] <- "eXpress"
myResults <- cbind(myResults,read.table("rsem/rsemA5.txt")[,3])
colnames(myResults)[5] <- "rsem"
myResults <- cbind(myResults,read.table("cufflinks/cuffA5.txt")[,3])
colnames(myResults)[6] <- "cufflinks"
myResults <- cbind(myResults,read.table("casper/sailA5.txt")[,3])
colnames(myResults)[7] <- "casper"
myResults <- cbind(myResults,log(read.table("bitseq/A5.thetaMeans")[,2]))
colnames(myResults)[8] <- "bitseq MCMC"
vb <- read.table("vb/A5.m_alphas")
vb[,1] <- vb[,1]/(1-vb[1,1])
vb <- log(vb[-1,1])
myResults <- cbind(myResults,vb)
colnames(myResults)[9] <- "bitseq VB"
myResults <- cbind(myResults,read.table("tigar/tigarA5.txt")[,3])
colnames(myResults)[10] <- "tigar"

ind <- which(is.infinite(myResults[,2])==FALSE)
apply(myResults[ind,-c(1,2)],2,function(y){sum(abs(y-myResults[ind,2]))})
for(k in 3:10){
ind1 <- which(is.infinite(myResults[,k])==FALSE)
ind <- intersect(ind,ind1)
print(length(ind))
}
apply(myResults[ind,-c(1,2)],2,function(y){sum(abs(y-myResults[ind,2]))})
aoua <- abs(myResults[ind,-c(1,2)] - myResults[ind,2])
boxplot(log(aoua))

aoua <- myResults[ind,-c(1,2)] - myResults[ind,2]
aoua1 <- exp(myResults[ind,-c(1,2)]) - exp(myResults[ind,2])
aoua2 <- log(myResults[ind,-c(1,2)]/myResults[ind,2])
boxplot(abs(aoua),ylim = c(0,1))
apply(abs(aoua),2,sum)

#aoua <- log(myResults[ind,-c(1,2)]/myResults[ind,2])
#boxplot(aoua)

mae5 <- apply(abs(aoua),2,sum)/dim(aoua)[1]
maExpe5 <- apply(abs(aoua1),2,sum)/dim(aoua1)[1]
MLR5 <- apply(abs(aoua2),2,sum)/dim(aoua1)[1]
perm <- order(mae2)
barplot(mae5[perm],ylab = "MAE")


mae <- rbind(mae1,mae2, mae3, mae4, mae5)
write.table(mae,file = "mae.txt")

