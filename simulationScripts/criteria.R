mae <- read.table("mae.txt")
mae <- apply(mae,2,mean)
MeanWgeTrue <- read.table("MeanWgeTrue.txt")
rownames(MeanWgeTrue)[c(1,5)] <- c("bitseqMCMC","bitseqVB")
MeanWgeInterReplicate <- read.table("MeanWgeInterReplicate.txt")
names(mae)[c(6,7)] <- rownames(MeanWgeTrue)[c(1,5)]


criterion <- array(data = NA, dim =c(8,3))
colnames(criterion) <- c("Theta","WGE-True","WGE-Inter")
rownames(criterion) <- rownames(MeanWgeTrue)
criterion[,2] <- MeanWgeTrue[,1]
criterion[,3] <- MeanWgeInterReplicate[,1]
criterion[names(mae),1] <- mae


mm <- apply(criterion,1,mean)
perm <- order(mm,decreasing=TRUE)

sCrit <- apply(criterion,2,function(x)return(x/sum(x)))

mm <- apply(sCrit,1,mean)
pdf(file = "simCriteria.pdf",width = 10,height = 5)
mm <- apply(sCrit,1,mean)
perm <- order(mm,decreasing=TRUE)
par(mfrow = c(1,1),mar=c(4,6.0,1,1))
barplot(t(sCrit[perm,]),beside=TRUE,horiz = TRUE,legend.text=TRUE,xlab = "Mean Absolute Error",ylab = "",las=1)
dev.off()

#times <- read.table("runtime.txt",row.names = 1)
#perm <- order(times[,1])
#barplot(log(times[perm,1]/(5*60)), names.arg = rownames(times)[perm],ylab = "log(run-time)",)

#times <- read.table("runTime.txt",header=TRUE)
#times <- times/60
#meanTimes <- apply(times,2,mean)
#permTimes <- order(meanTimes)
#pdf(file = "spankiTime.pdf",width = 10,height = 5)
#par(mfrow = c(1,1),mar = c(3,4,1,1))
#boxplot(times[,permTimes],ylab = "run-time (minutes)",log="y")
#dev.off()

