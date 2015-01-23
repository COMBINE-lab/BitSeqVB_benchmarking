# load names
txid <- read.table("../trNames.tr")[,1]
K <- length(txid)
#conditionA1

set.seed(1)
mus <- 10 + 200*runif(K)
mySize = 20
musA1 <- rnbinom(K,mu = mus,size = mySize) 
musA2 <- rnbinom(K,mu = mus,size = mySize) 
musA3 <- rnbinom(K,mu = mus,size = mySize) 
musA4 <- rnbinom(K,mu = mus,size = mySize) 
musA5 <- rnbinom(K,mu = mus,size = mySize) 
#plot(log(musA1),log(musA2),col = colInd,pch = 16,cex = 0.5)

rpk <- musA1
tr_File_1 <- data.frame(txid,rpk)
write.table(tr_File_1,file = "tr_File_A1.tr",quote = FALSE,row.names = FALSE,sep = "\t")
#conditionA2
rpk <- musA2
tr_File_1 <- data.frame(txid,rpk)
write.table(tr_File_1,file = "tr_File_A2.tr",quote = FALSE,row.names = FALSE,sep = "\t")
rpk <- musA3
tr_File_1 <- data.frame(txid,rpk)
write.table(tr_File_1,file = "tr_File_A3.tr",quote = FALSE,row.names = FALSE,sep = "\t")
#conditionA2
rpk <- musA4
tr_File_1 <- data.frame(txid,rpk)
write.table(tr_File_1,file = "tr_File_A4.tr",quote = FALSE,row.names = FALSE,sep = "\t")
rpk <- musA5
tr_File_1 <- data.frame(txid,rpk)
write.table(tr_File_1,file = "tr_File_A5.tr",quote = FALSE,row.names = FALSE,sep = "\t")


