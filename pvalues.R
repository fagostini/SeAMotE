# library(qvalue)
#  
# DF <- as.data.frame(read.table("a"))
#  
# qvalue(DF$V2, fdr.level=0.01)

Args <- commandArgs()

library(Exact)
library(qvalue)
# library(Barnard)

th <- as.double(Args[6])
DFall <- as.data.frame(read.table(Args[5]))
# DF <- as.data.frame(DFall[which((DFall$V2/DFall$V3 >= th) | (DFall$V4/DFall$V5 >= th)),])
DF <- as.data.frame(DFall[which(DFall$V2/DFall$V3 >= th),])
DF$V3 <- DF$V3-DF$V2;
DF$V5 <- DF$V5-DF$V4;
# DF$V2 <- DF$V2
# DF$V4 <- DF$V4


Pval <- vector(length=dim(DF)[1])
for( i in seq(1,dim(DF)[1],1) ){
	m <- matrix(c(DF[i,]$V2, DF[i,]$V3, DF[i,]$V4, DF[i,]$V5), 2, 2)
# 	if( (DF[i,]$V3 == 0) & (DF[i,]$V5 == 0) ){
# 		Pval[i] <- 0.0000000000000001
# 	}else{
		Pval[i] <- fisher.test(m, alternative="greater")$p.value
# 	}
# 	if( (DF[i,]$V3 == 0) & (DF[i,]$V5 == 0) ){
#  		Pval[i] <- 0.0000000000000001
# 	}else{
# 		Pval[i] <- exact.test(m, alternative="greater",  method="Boschloo", to.plot = FALSE)$p.value
# 	}

#	Pval[i] <- barnardw.test(DF[i,]$V2, DF[i,]$V3, DF[i,]$V4, DF[i,]$V5)$p.value[1]
}
DF$Pval <- Pval

DF <- DF[which(DF$Pval <= quantile(qvalue(p.adjust(Pval, method = "holm", n = length(Pval)), fdr.level=0.01)$qvalue, 0.05)),]

# print.data.frame(DF, row.names=FALSE, col.names=FALSE)
# write.table(DF[which(DF$Pval < 0.05),], Args[5], row.names=FALSE, col.names=FALSE, quote = FALSE)
# write.table(DF[which(DF$Pval < 0.05),], "tmp/last.txt", row.names=FALSE, col.names=FALSE, quote = FALSE)
write.table(DF, Args[5], row.names=FALSE, col.names=FALSE, quote = FALSE)


