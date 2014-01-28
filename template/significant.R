Args <- commandArgs()

th <- as.double(Args[6])
out <- as.integer(Args[7])

DFall <- as.data.frame(read.table(paste("tmp/motifs_", Args[5],"nt.dat", sep="")))
# DFdis <- as.data.frame(read.table(paste("tmp/motifs_", Args[5],"nt_Pdistribution.dat", sep="")))

DF <- as.data.frame(DFall[which(DFall$V2/DFall$V3 >= th),])
DF$V6 <- DF$V2/DF$V3-DF$V4/DF$V5
DF$V3 <- DF$V3-DF$V2;
DF$V5 <- DF$V5-DF$V4;
DFbest <- DF[order(DF$V6, decreasing=T)[1:10],]

myFun <- function(x){ m <- matrix(c(DFbest[x,]$V2, DFbest[x,]$V3, DFbest[x,]$V4, DFbest[x,]$V5), 2, 2); return(fisher.test(m, alternative="greater")$p.value) }
Pval <- lapply(seq(1,10,1), myFun)
DFbest$Pval <- unlist(Pval)

# library(multicore)
# myFun <- function(x){ m <- matrix(c(DF[x,]$V2, DF[x,]$V3, DF[x,]$V4, DF[x,]$V5), 2, 2); return(fisher.test(m, alternative="greater")$p.value) }
# Pval <- mclapply(seq(1,dim(DF)[1],1), myFun, mc.cores=4)
# DF$Pval <- unlist(Pval)

write.table(DFbest,"tmp/best_motifs.dat", row.names=FALSE, col.names=FALSE, quote = FALSE, append=TRUE)
# write.table(DF[order(unlist(Pval))[1:10],],"tmp/best_motifs.dat", row.names=FALSE, col.names=FALSE, quote = FALSE, append=TRUE)
DF <- DF[which( (DF$V2/(DF$V2+DF$V3) >= th) & (DF$V6 >= 0) ),]
if( dim(DF)[1] > 500 ){
	DF <- DF[order(DF$V6, decreasing=T)[1:500],]	
}
write.table(DF, paste("tmp/motifs_",Args[5],"nt.dat", sep=""), row.names=FALSE, col.names=FALSE, quote = FALSE)
