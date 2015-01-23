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
A1 <- read.table("A1.cufflinks.out/isoforms.fpkm_tracking",header=TRUE)
groundTruth <- read.table("../A1/transcriptNames.txt")
sailA1 <- cbind(groundTruth,numeric(dim(groundTruth)[1]))
colnames(sailA1) <- c("id","true","est")
nas <- 0
for (i in 1:dim(groundTruth)[1]){
	j <- which(sailNames[,2] == as.character(sailA1[i,1]))
	k <- which(A1[,1] == sailNames[j,2])
	if(length(k)==0){
		sailA1[i,3] <- 0
		nas <- nas +1
		}else{
		sailA1[i,3] <- A1[k[1],10]*A1[k[1],8] #(fpkm*length)
	}
	if(i%%1000 == 0){print(i);print(nas)}
}
sail <- sailA1
sailA1[,2] <- log(sailA1[,2]/sum(sailA1[,2]))
sailA1[,3] <- log(sailA1[,3]/sum(sailA1[,3]))
#smoothScatter(sailA1[,-1],ylim = c(-15,-5));points(c(-30,10),c(-30,10),type = "l",col = 2)
write.table(sailA1,file = "cuffA1.txt")

A1 <- read.table("A2.cufflinks.out/isoforms.fpkm_tracking",header=TRUE)
groundTruth <- read.table("../A2/transcriptNames.txt")
sailA1 <- cbind(groundTruth,numeric(dim(groundTruth)[1]))
colnames(sailA1) <- c("id","true","est")
nas <- 0
for (i in 1:dim(groundTruth)[1]){
	j <- which(sailNames[,2] == as.character(sailA1[i,1]))
	k <- which(A1[,1] == sailNames[j,2])
	if(length(k)==0){
		sailA1[i,3] <- 0
		nas <- nas +1
		}else{
		sailA1[i,3] <- A1[k[1],10]*A1[k[1],8] #(fpkm*length)
	}
	if(i%%1000 == 0){print(i);print(nas)}
}
sail <- sailA1
sailA1[,2] <- log(sailA1[,2]/sum(sailA1[,2]))
sailA1[,3] <- log(sailA1[,3]/sum(sailA1[,3]))
#smoothScatter(sailA1[,-1],ylim = c(-15,-5));points(c(-30,10),c(-30,10),type = "l",col = 2)
write.table(sailA1,file = "cuffA2.txt")

A1 <- read.table("A3.cufflinks.out/isoforms.fpkm_tracking",header=TRUE)
groundTruth <- read.table("../A3/transcriptNames.txt")
sailA1 <- cbind(groundTruth,numeric(dim(groundTruth)[1]))
colnames(sailA1) <- c("id","true","est")
nas <- 0
for (i in 1:dim(groundTruth)[1]){
	j <- which(sailNames[,2] == as.character(sailA1[i,1]))
	k <- which(A1[,1] == sailNames[j,2])
	if(length(k)==0){
		sailA1[i,3] <- 0
		nas <- nas +1
		}else{
		sailA1[i,3] <- A1[k[1],10]*A1[k[1],8] #(fpkm*length)
	}
	if(i%%1000 == 0){print(i);print(nas)}
}
sail <- sailA1
sailA1[,2] <- log(sailA1[,2]/sum(sailA1[,2]))
sailA1[,3] <- log(sailA1[,3]/sum(sailA1[,3]))
#smoothScatter(sailA1[,-1],ylim = c(-15,-5));points(c(-30,10),c(-30,10),type = "l",col = 2)
write.table(sailA1,file = "cuffA3.txt")

A1 <- read.table("A4.cufflinks.out/isoforms.fpkm_tracking",header=TRUE)
groundTruth <- read.table("../A4/transcriptNames.txt")
sailA1 <- cbind(groundTruth,numeric(dim(groundTruth)[1]))
colnames(sailA1) <- c("id","true","est")
nas <- 0
for (i in 1:dim(groundTruth)[1]){
	j <- which(sailNames[,2] == as.character(sailA1[i,1]))
	k <- which(A1[,1] == sailNames[j,2])
	if(length(k)==0){
		sailA1[i,3] <- 0
		nas <- nas +1
		}else{
		sailA1[i,3] <- A1[k[1],10]*A1[k[1],8] #(fpkm*length)
	}
	if(i%%1000 == 0){print(i);print(nas)}
}
sail <- sailA1
sailA1[,2] <- log(sailA1[,2]/sum(sailA1[,2]))
sailA1[,3] <- log(sailA1[,3]/sum(sailA1[,3]))
#smoothScatter(sailA1[,-1],ylim = c(-15,-5));points(c(-30,10),c(-30,10),type = "l",col = 2)
write.table(sailA1,file = "cuffA4.txt")

A1 <- read.table("A5.cufflinks.out/isoforms.fpkm_tracking",header=TRUE)
groundTruth <- read.table("../A5/transcriptNames.txt")
sailA1 <- cbind(groundTruth,numeric(dim(groundTruth)[1]))
colnames(sailA1) <- c("id","true","est")
nas <- 0
for (i in 1:dim(groundTruth)[1]){
	j <- which(sailNames[,2] == as.character(sailA1[i,1]))
	k <- which(A1[,1] == sailNames[j,2])
	if(length(k)==0){
		sailA1[i,3] <- 0
		nas <- nas +1
		}else{
		sailA1[i,3] <- A1[k[1],10]*A1[k[1],8] #(fpkm*length)
	}
	if(i%%1000 == 0){print(i);print(nas)}
}
sail <- sailA1
sailA1[,2] <- log(sailA1[,2]/sum(sailA1[,2]))
sailA1[,3] <- log(sailA1[,3]/sum(sailA1[,3]))
#smoothScatter(sailA1[,-1],ylim = c(-15,-5));points(c(-30,10),c(-30,10),type = "l",col = 2)
write.table(sailA1,file = "cuffA5.txt")


