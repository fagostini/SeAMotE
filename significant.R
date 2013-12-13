Args <- commandArgs()

th <- as.double(Args[6])
out <- as.integer(Args[7])

DFall <- as.data.frame(read.table(paste("tmp/motifs_", Args[5],"nt.dat", sep="")))
# DFdis <- as.data.frame(read.table(paste("tmp/motifs_", Args[5],"nt_Pdistribution.dat", sep="")))

DF <- as.data.frame(DFall[which(DFall$V2/DFall$V3 >= th),])
DF$V3 <- DF$V3-DF$V2;
DF$V5 <- DF$V5-DF$V4;

#Pval <- vector(length=dim(DF)[1], mode="numeric")
#Nval <- vector(length=dim(DF)[1])
#for( i in seq(1,dim(DF)[1],1) ){
#	m <- matrix(c(DF[i,]$V2, DF[i,]$V3, DF[i,]$V4, DF[i,]$V5), 2, 2)
#	Pval[i] <- fisher.test(m, alternative="greater")$p.value
#	Nval[i] <- fisher.test(m, alternative="less")$p.value
#}
#DF$Pval <- Pval
# DF$Nval <- Nval

library(multicore)
myFun <- function(x){ m <- matrix(c(DF[x,]$V2, DF[x,]$V3, DF[x,]$V4, DF[x,]$V5), 2, 2); return(fisher.test(m, alternative="greater")$p.value) }
Pval <- mclapply(seq(1,dim(DF)[1],1), myFun, mc.cores=2)
DF$Pval <- unlist(Pval)

# if( length(Pval) < 1000 ){
# 	if( length(Pval) <= 50 ){
# 		Pquant <- max(Pval)
# 		Nquant <- max(Nval)
# 	}else{
# 		Pquant <- Pval[order(Pval, decreasing=FALSE)[50]]
# 		Nquant <- Nval[order(Nval, decreasing=FALSE)[50]]
# 	}
# }else{
# 	Pquant <- quantile(Pval, 0.1)
# 	Nquant <- quantile(Nval, 0.1)	
# }

# if( length(Pval) < 150 ){
# 	select <- order(Pval)
# }else{
# 	select <- order(Pval)[1:150]
# }

if( out == 1 ){
#	write.table(DF[which( (DF$V2/(DF$V2+DF$V3) >= th+0.1 & DF$Pval <= min(Pval) ) ),],"tmp/best_motifs.dat", row.names=FALSE, col.names=FALSE, quote = FALSE, append=TRUE)
	write.table(DF[which( (DF$V2/(DF$V2+DF$V3) >= th) ),], paste("tmp/motifs_",Args[5],"nt.dat", sep=""), row.names=FALSE, col.names=FALSE, quote = FALSE)
}else{
	write.table(DF[which( DF$Pval <= min(unlist(Pval)) ),],"tmp/best_motifs.dat", row.names=FALSE, col.names=FALSE, quote = FALSE, append=TRUE)
#	write.table(DF[which(DF$Pval <= Pquant ),], paste("tmp/motifs_",Args[5],"nt.dat", sep=""), row.names=FALSE, col.names=FALSE, quote = FALSE)
#	write.table(DF[which(row.names(DF) %in% select),], paste("tmp/motifs_",Args[5],"nt.dat", sep=""), row.names=FALSE, col.names=FALSE, quote = FALSE)
	write.table(DF[which( (DF$V2/(DF$V2+DF$V3) >= th) ),], paste("tmp/motifs_",Args[5],"nt.dat", sep=""), row.names=FALSE, col.names=FALSE, quote = FALSE)
}

