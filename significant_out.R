Args <- commandArgs()

th <- as.double(Args[6])

DFall <- as.data.frame(read.table(paste("tmp/motifs_", Args[5],"nt_exp.dat", sep="")))

DF <- as.data.frame(DFall[which(DFall$V2/DFall$V3 >= th),])
DF$V3 <- DF$V3-DF$V2;
DF$V5 <- DF$V5-DF$V4;

Pval <- vector(length=dim(DF)[1])
for( i in seq(1,dim(DF)[1],1) ){
	m <- matrix(c(DF[i,]$V2, DF[i,]$V3, DF[i,]$V4, DF[i,]$V5), 2, 2)
	Pval[i] <- fisher.test(m, alternative="greater")$p.value
}
DF$Pval <- Pval

write.table(DF[order(Pval)[1:10],],"tmp/best_motifs.dat", row.names=FALSE, col.names=FALSE, quote = FALSE, append=TRUE)
write.table(DF, paste("tmp/motifs_",Args[5],"nt_exp.dat", sep=""), row.names=FALSE, col.names=FALSE, quote = FALSE)
