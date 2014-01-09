Args <- commandArgs()

th <- as.double(Args[6])

DFall <- as.data.frame(read.table(paste("tmp/motifs_", Args[5],"nt_exp.dat", sep="")))

DF <- as.data.frame(DFall[which(DFall$V2/DFall$V3 >= th),])
DF$V3 <- DF$V3-DF$V2;
DF$V5 <- DF$V5-DF$V4;

library(multicore)
myFun <- function(x){ m <- matrix(c(DF[x,]$V2, DF[x,]$V3, DF[x,]$V4, DF[x,]$V5), 2, 2); return(fisher.test(m, alternative="greater")$p.value) }
Pval <- mclapply(seq(1,dim(DF)[1],1), myFun, mc.cores=4)
DF$Pval <- unlist(Pval)

write.table(DF[order(unlist(Pval))[1:10],],"tmp/best_motifs.dat", row.names=FALSE, col.names=FALSE, quote = FALSE, append=TRUE)
write.table(DF, paste("tmp/motifs_",Args[5],"nt_exp.dat", sep=""), row.names=FALSE, col.names=FALSE, quote = FALSE)
