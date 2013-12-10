# library(qvalue)
#  
# DF <- as.data.frame(read.table("a"))
#  
# qvalue(DF$V2, fdr.level=0.01)

library(Exact)
# library(Barnard)
DF <- as.data.frame(read.table("tmp/motifs_5nt.dat"))
DF$V3 <- DF$V3-DF$V2;
DF$V5 <- DF$V5-DF$V4;

Pval <- vector(length=dim(DF)[1])
for( i in seq(1,dim(DF)[1],1) ){
	m <- matrix(c(DF[i,]$V2, DF[i,]$V3, DF[i,]$V4, DF[i,]$V5), 2, 2)
#	Pval[i] <- fisher.test(m, alternative="greater")$p.value
	Pval[i] <- exact.test(m, alternative="greater",  method="Boschloo", to.plot = FALSE)$p.value
#	Pval[i] <- barnardw.test(DF[i,]$V2, DF[i,]$V3, DF[i,]$V4, DF[i,]$V5)$p.value[1]
}
DF$Pval <- Pval

print.data.frame(DF, row.names=FALSE, col.names=FALSE)

