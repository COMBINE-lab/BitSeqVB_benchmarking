system("grep \">\" ../../transcriptome_data/known.fa > names.txt")
myFile <- file("names.txt",open="r")
sailNames <- array(data = NA, dim = c(48009,2))
i <- 0
while (length(oneLine <- readLines(myFile, n = 1, warn = FALSE)) > 0){
        i <- i + 1
        sailNames[i,1] <- as.numeric(strsplit(strsplit(oneLine,split = " ")[[1]][1],split = ">")[[1]][2])
        sailNames[i,2] <- strsplit(oneLine,split = " ")[[1]][2]
}
close(myFile)
par(mfrow = c(2,3))



load('../tophat/A1.tophat.out/expCasper.RData')
rep1 #load
head(fData(rep1$exp)) #shows the expected counts
head(exprs(rep1$exp)) #shows within gene est
A1 <- fData(rep1$exp)
wge <- exprs(rep1$exp)
groundTruth <- read.table("../A1/transcriptNames.txt")
sailA1 <- cbind(groundTruth,numeric(dim(groundTruth)[1]))
colnames(sailA1) <- c("id","true","est")
nas <- 0
for (i in 1:dim(groundTruth)[1]){
        j <- which(sailNames[,2] == as.character(sailA1[i,1]))
        k <- which(A1[,1] == sailNames[j,2])
	if(length(k)>0){
		sailA1[i,3] <- wge[k[1]]*A1[k[1],4] # (relative within gene expression)*(counts per gene) = counts per transcript
	}else{nas <- nas + 1}
        if(i%%1000 == 0){print(i);print(nas)}
}
sailA1[,2] <- log(sailA1[,2]/sum(sailA1[,2]))
sailA1[,3] <- log(sailA1[,3]/sum(sailA1[,3]))
#smoothScatter(sailA1[,-1]);points(c(-30,10),c(-30,10),type = "l",col = 2)
write.table(sailA1,file = "sailA1.txt")


load('../tophat/A2.tophat.out/expCasper.RData')
rep1 #load
head(fData(rep1$exp)) #shows the expected counts
head(exprs(rep1$exp)) #shows within gene est
A1 <- fData(rep1$exp)
wge <- exprs(rep1$exp)
groundTruth <- read.table("../A2/transcriptNames.txt")
sailA1 <- cbind(groundTruth,numeric(dim(groundTruth)[1]))
colnames(sailA1) <- c("id","true","est")
nas <- 0
for (i in 1:dim(groundTruth)[1]){
        j <- which(sailNames[,2] == as.character(sailA1[i,1]))
        k <- which(A1[,1] == sailNames[j,2])
	if(length(k)>0){
		sailA1[i,3] <- wge[k[1]]*A1[k[1],4] # (relative within gene expression)*(counts per gene) = counts per transcript
	}else{nas <- nas + 1}
        if(i%%1000 == 0){print(i);print(nas)}
}
sailA1[,2] <- log(sailA1[,2]/sum(sailA1[,2]))
sailA1[,3] <- log(sailA1[,3]/sum(sailA1[,3]))
#smoothScatter(sailA1[,-1]);points(c(-30,10),c(-30,10),type = "l",col = 2)
write.table(sailA1,file = "sailA2.txt")

load('../tophat/A3.tophat.out/expCasper.RData')
rep1 #load
head(fData(rep1$exp)) #shows the expected counts
head(exprs(rep1$exp)) #shows within gene est
A1 <- fData(rep1$exp)
wge <- exprs(rep1$exp)
groundTruth <- read.table("../A3/transcriptNames.txt")
sailA1 <- cbind(groundTruth,numeric(dim(groundTruth)[1]))
colnames(sailA1) <- c("id","true","est")
nas <- 0
for (i in 1:dim(groundTruth)[1]){
        j <- which(sailNames[,2] == as.character(sailA1[i,1]))
        k <- which(A1[,1] == sailNames[j,2])
	if(length(k)>0){
		sailA1[i,3] <- wge[k[1]]*A1[k[1],4] # (relative within gene expression)*(counts per gene) = counts per transcript
	}else{nas <- nas + 1}
        if(i%%1000 == 0){print(i);print(nas)}
}
sailA1[,2] <- log(sailA1[,2]/sum(sailA1[,2]))
sailA1[,3] <- log(sailA1[,3]/sum(sailA1[,3]))
#smoothScatter(sailA1[,-1]);points(c(-30,10),c(-30,10),type = "l",col = 2)
write.table(sailA1,file = "sailA3.txt")

load('../tophat/A4.tophat.out/expCasper.RData')
rep1 #load
head(fData(rep1$exp)) #shows the expected counts
head(exprs(rep1$exp)) #shows within gene est
A1 <- fData(rep1$exp)
wge <- exprs(rep1$exp)
groundTruth <- read.table("../A4/transcriptNames.txt")
sailA1 <- cbind(groundTruth,numeric(dim(groundTruth)[1]))
colnames(sailA1) <- c("id","true","est")
nas <- 0
for (i in 1:dim(groundTruth)[1]){
        j <- which(sailNames[,2] == as.character(sailA1[i,1]))
        k <- which(A1[,1] == sailNames[j,2])
	if(length(k)>0){
		sailA1[i,3] <- wge[k[1]]*A1[k[1],4] # (relative within gene expression)*(counts per gene) = counts per transcript
	}else{nas <- nas + 1}
        if(i%%1000 == 0){print(i);print(nas)}
}
sailA1[,2] <- log(sailA1[,2]/sum(sailA1[,2]))
sailA1[,3] <- log(sailA1[,3]/sum(sailA1[,3]))
#smoothScatter(sailA1[,-1]);points(c(-30,10),c(-30,10),type = "l",col = 2)
write.table(sailA1,file = "sailA4.txt")

load('../tophat/A5.tophat.out/expCasper.RData')
rep1 #load
head(fData(rep1$exp)) #shows the expected counts
head(exprs(rep1$exp)) #shows within gene est
A1 <- fData(rep1$exp)
wge <- exprs(rep1$exp)
groundTruth <- read.table("../A5/transcriptNames.txt")
sailA1 <- cbind(groundTruth,numeric(dim(groundTruth)[1]))
colnames(sailA1) <- c("id","true","est")
nas <- 0
for (i in 1:dim(groundTruth)[1]){
        j <- which(sailNames[,2] == as.character(sailA1[i,1]))
        k <- which(A1[,1] == sailNames[j,2])
	if(length(k)>0){
		sailA1[i,3] <- wge[k[1]]*A1[k[1],4] # (relative within gene expression)*(counts per gene) = counts per transcript
	}else{nas <- nas + 1}
        if(i%%1000 == 0){print(i);print(nas)}
}
sailA1[,2] <- log(sailA1[,2]/sum(sailA1[,2]))
sailA1[,3] <- log(sailA1[,3]/sum(sailA1[,3]))
#smoothScatter(sailA1[,-1]);points(c(-30,10),c(-30,10),type = "l",col = 2)
write.table(sailA1,file = "sailA5.txt")
















