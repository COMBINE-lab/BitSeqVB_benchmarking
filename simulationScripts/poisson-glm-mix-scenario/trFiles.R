
# load names
txid <- read.table("../trNames.tr")[,1]
KTR <- length(txid)
#conditionA1

set.seed(10)
mus <- 30*runif(KTR)


set.seed(100)
# number of simulated datasets
iterations <- 1

runs<-vector("list",length=iterations)

# ktrue: number of mixture components
ktrue <- 20

alpha.true<-vector("list",length=1)
beta.true<-vector("list",length=1)
gamma.true<-vector("list",length=1)
weights.true<-vector("list",length=1)


# extracting the reference from real data
n<-KTR
#
p <- 5;
data <- array(data=NA,dim=c(KTR,p+1))
data[,1]<- log(mus+1)
tau<-1
q<-p
simdata<-array(data=NA, dim =c(iterations,n,p+1))

K<-ktrue
# simulating the weights form a Dirichlet distribution with mean rep(1/K,K)
weights<-rgamma(K,scale=1,shape=1)
weights<-weights/sum(weights)

#       this is the q times K matrix containing the constant terms.
alphatrue <- array(data = NA,dim = c(q,K)) #- 1.5
betatrue <- array(data = NA, dim = c(q,K,tau))
if(2*floor(K/2)<K){
alphatrue[1,(K+1)/2]<-0
betatrue[1,(K+1)/2,1]<-1
for(k in 1:((K-1)/2)){
alphatrue[1,k]<- (((K+1)/2)-k)/2
betatrue[1,k,1]<- (10 - alphatrue[1,k])/10
}
for(k in (1+(K+1)/2):K){
alphatrue[1,k]<- (((K+1)/2) - k)/2
betatrue[1,k,1]<- (10 - alphatrue[1,k])/10
}
alphatrue[2,]<- alphatrue[1,] #+ rnorm(length(alphatrue[1,]))
alphatrue[3,]<- alphatrue[1,] #+ rnorm(length(alphatrue[1,]))
alphatrue[4,]<- alphatrue[1,] #+ rnorm(length(alphatrue[1,]))
alphatrue[5,]<- alphatrue[1,] #+ rnorm(length(alphatrue[1,]))
betatrue[2,,1]<- betatrue[1,,1] + 0.05*rnorm(length(betatrue[1,,1]))
betatrue[3,,1]<- betatrue[1,,1] + 0.05*rnorm(length(betatrue[1,,1]))
betatrue[4,,1]<- betatrue[1,,1] + 0.05*rnorm(length(betatrue[1,,1]))
betatrue[5,,1]<- betatrue[1,,1] + 0.05*rnorm(length(betatrue[1,,1]))
}else{
K<-K+1

alphatrue <- array(data = NA,dim = c(q,K)) #- 1.5
betatrue <- array(data = NA, dim = c(q,K,tau))

alphatrue[1,(K+1)/2]<-0
betatrue[1,(K+1)/2,1]<-1
for(k in 1:((K-1)/2)){
alphatrue[1,k]<- (((K+1)/2)-k)/2
betatrue[1,k,1]<- (10 - alphatrue[1,k])/10
}
for(k in (1+(K+1)/2):K){
alphatrue[1,k]<- (((K+1)/2) - k)/2
betatrue[1,k,1]<- (10 - alphatrue[1,k])/10
}
alphatrue[2,]<- alphatrue[1,] #+ 0.1*rnorm(length(alphatrue[1,]))
alphatrue[3,]<- alphatrue[1,] #+ 0.1*rnorm(length(alphatrue[1,]))
alphatrue[4,]<- alphatrue[1,] #+ 0.1*rnorm(length(alphatrue[1,]))
alphatrue[5,]<- alphatrue[1,] #+ 0.1*rnorm(length(alphatrue[1,]))
betatrue[2,,1]<- betatrue[1,,1] + 0.001*rnorm(length(betatrue[1,,1]))
betatrue[3,,1]<- betatrue[1,,1] + 0.001*rnorm(length(betatrue[1,,1]))
betatrue[4,,1]<- betatrue[1,,1] + 0.001*rnorm(length(betatrue[1,,1]))
betatrue[5,,1]<- betatrue[1,,1] + 0.001*rnorm(length(betatrue[1,,1]))

K<-K-1
alphatrue <- array(alphatrue[,1:K],dim=c(q,K))
betatrue<-array(betatrue[,1:K,],dim=c(q,K,tau))

}
L<-c(1,1,1,1,1)
qq<- sum(L)
cumL<- cumsum(L)
index <- c(1,2,3,4,5)
ml<-max(L)
gammatrue <- array(data = NA,dim = c(q,max(L)))
for(j in 1:q){gammatrue[j,1:L[j]]<-runif(L[j]);gammatrue[j,1:L[j]]<-gammatrue[j,1:L[j]]-mean(gammatrue[j,1:L[j]])}
indcum<-c(0,cumsum(L))

x<-array(data = NA, dim = c(n,tau))
for(t in 1:tau)x[,t] <- data[,t]
#       simulating the true allocations
sam<-rmultinom(n, 1, weights)

weights<-table(apply(sam,2,which.max))/n
weights<-as.vector(weights)

truez<-numeric(n)
for (i in 1:n){
truez[i]<-order(sam[,i])[K]
}



# simulating y's
diskoles <- c(0,0)
zmean <- matrix(data = 0, nrow = n, ncol = K, byrow = T)
#for (fora in 1:iterations){
fora <- 1
print("##############################")
print(fora)
print("##############################")
ss <- 1:n
y <- array(data = NA, dim = c(n,sum(L)))
for(k in 1:K){
ind <- ss[sam[k,]==1]
for (j in 1:q){
u<-numeric(length(ind))
for(t in 1:tau){u <- u + betatrue[j,k,t]*x[ind,t]}
for (l in 1:L[j]){
        y[ind,index[j]-1+l] <- rpois(length(ind),exp(alphatrue[j,k] + gammatrue[j,l] + u))
}
}
}
simdata[fora,,]<-cbind(x,y)
#}
#hist(simdata[1,,2:6])
print(length(which(simdata[1,,2]==0)))
#plot(log(y[,1:2]))


K<-KTR
mySize <- 50
musA1 <- rnbinom(K,mu = y[,1],size = mySize) 
musA2 <- rnbinom(K,mu = y[,2],size = mySize) 
musA3 <- rnbinom(K,mu = y[,3],size = mySize) 
musA4 <- rnbinom(K,mu = y[,4],size = mySize) 
musA5 <- rnbinom(K,mu = y[,5],size = mySize) 
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


