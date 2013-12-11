# library(qvalue)
#  
# DF <- as.data.frame(read.table("a"))
#  
# qvalue(DF$V2, fdr.level=0.01)

Args <- commandArgs()

library(Exact)
# library(Barnard)
DF <- as.data.frame(read.table(Args[5]))
DF$V3 <- DF$V3-DF$V2+1;
DF$V5 <- DF$V5-DF$V4+1;
DF$V2 <- DF$V2+1
DF$V4 <- DF$V4+1

Pval <- vector(length=dim(DF)[1])
for( i in seq(1,dim(DF)[1],1) ){
	m <- matrix(c(DF[i,]$V2, DF[i,]$V3, DF[i,]$V4, DF[i,]$V5), 2, 2)
#	Pval[i] <- fisher.test(m, alternative="greater")$p.value
	Pval[i] <- exact.test(m, alternative="greater",  method="Boschloo", to.plot = FALSE)$p.value
#	Pval[i] <- barnardw.test(DF[i,]$V2, DF[i,]$V3, DF[i,]$V4, DF[i,]$V5)$p.value[1]
}
DF$Pval <- Pval

print.data.frame(DF, row.names=FALSE, col.names=FALSE)
write.table(DF, Args[5], row.names=FALSE, col.names=FALSE)
