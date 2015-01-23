system("grep \">\" ../../transcriptome_data/known.fa > names.txt")
myFile <- file("names.txt",open="r")
xpressNames <- array(data = NA, dim = c(48009,2))
i <- 0
while (length(oneLine <- readLines(myFile, n = 1, warn = FALSE)) > 0){
	i <- i + 1
	xpressNames[i,1] <- as.numeric(strsplit(strsplit(oneLine,split = " ")[[1]][1],split = ">")[[1]][2])
	xpressNames[i,2] <- strsplit(oneLine,split = " ")[[1]][2]
}
close(myFile)

A1 <- read.table("A1/results.xprs",header = TRUE)
ind <- A1[,2]
groundTruth <- read.table("../A1/transcriptNames.txt")
xpressA1 <- cbind(groundTruth,numeric(dim(groundTruth)[1]))
colnames(xpressA1) <- c("id","true","est")
for (i in 1:dim(groundTruth)[1]){
	j <- which(xpressNames[,2] == as.character(xpressA1[i,1]))
	k <- which(A1[,2] == xpressNames[j,1])
	xpressA1[i,3] <- A1[k,7]
	if(i%%1000 == 0){print(i)}
}
xpressA1[,2] <- log(xpressA1[,2]/sum(xpressA1[,2]))
xpressA1[,3] <- log(xpressA1[,3]/sum(xpressA1[,3]))
#smoothScatter(xpressA1[,-1]);points(c(-30,10),c(-30,10),type = "l",col = 2)
write.table(xpressA1,file = "xpressA1.txt")

A1 <- read.table("A2/results.xprs",header = TRUE)
ind <- A1[,2]
groundTruth <- read.table("../A2/transcriptNames.txt")
xpressA1 <- cbind(groundTruth,numeric(dim(groundTruth)[1]))
colnames(xpressA1) <- c("id","true","est")
for (i in 1:dim(groundTruth)[1]){
	j <- which(xpressNames[,2] == as.character(xpressA1[i,1]))
	k <- which(A1[,2] == xpressNames[j,1])
	xpressA1[i,3] <- A1[k,7]
	if(i%%1000 == 0){print(i)}
}
xpressA1[,2] <- log(xpressA1[,2]/sum(xpressA1[,2]))
xpressA1[,3] <- log(xpressA1[,3]/sum(xpressA1[,3]))
#smoothScatter(xpressA1[,-1]);points(c(-30,10),c(-30,10),type = "l",col = 2)
write.table(xpressA1,file = "xpressA2.txt")

A1 <- read.table("A3/results.xprs",header = TRUE)
ind <- A1[,2]
groundTruth <- read.table("../A3/transcriptNames.txt")
xpressA1 <- cbind(groundTruth,numeric(dim(groundTruth)[1]))
colnames(xpressA1) <- c("id","true","est")
for (i in 1:dim(groundTruth)[1]){
	j <- which(xpressNames[,2] == as.character(xpressA1[i,1]))
	k <- which(A1[,2] == xpressNames[j,1])
	xpressA1[i,3] <- A1[k,7]
	if(i%%1000 == 0){print(i)}
}
xpressA1[,2] <- log(xpressA1[,2]/sum(xpressA1[,2]))
xpressA1[,3] <- log(xpressA1[,3]/sum(xpressA1[,3]))
#smoothScatter(xpressA1[,-1]);points(c(-30,10),c(-30,10),type = "l",col = 2)
write.table(xpressA1,file = "xpressA3.txt")

A1 <- read.table("A4/results.xprs",header = TRUE)
ind <- A1[,2]
groundTruth <- read.table("../A4/transcriptNames.txt")
xpressA1 <- cbind(groundTruth,numeric(dim(groundTruth)[1]))
colnames(xpressA1) <- c("id","true","est")
for (i in 1:dim(groundTruth)[1]){
	j <- which(xpressNames[,2] == as.character(xpressA1[i,1]))
	k <- which(A1[,2] == xpressNames[j,1])
	xpressA1[i,3] <- A1[k,7]
	if(i%%1000 == 0){print(i)}
}
xpressA1[,2] <- log(xpressA1[,2]/sum(xpressA1[,2]))
xpressA1[,3] <- log(xpressA1[,3]/sum(xpressA1[,3]))
#smoothScatter(xpressA1[,-1]);points(c(-30,10),c(-30,10),type = "l",col = 2)
write.table(xpressA1,file = "xpressA4.txt")

A1 <- read.table("A5/results.xprs",header = TRUE)
ind <- A1[,2]
groundTruth <- read.table("../A5/transcriptNames.txt")
xpressA1 <- cbind(groundTruth,numeric(dim(groundTruth)[1]))
colnames(xpressA1) <- c("id","true","est")
for (i in 1:dim(groundTruth)[1]){
	j <- which(xpressNames[,2] == as.character(xpressA1[i,1]))
	k <- which(A1[,2] == xpressNames[j,1])
	xpressA1[i,3] <- A1[k,7]
	if(i%%1000 == 0){print(i)}
}
xpressA1[,2] <- log(xpressA1[,2]/sum(xpressA1[,2]))
xpressA1[,3] <- log(xpressA1[,3]/sum(xpressA1[,3]))
#smoothScatter(xpressA1[,-1]);points(c(-30,10),c(-30,10),type = "l",col = 2)
write.table(xpressA1,file = "xpressA5.txt")

