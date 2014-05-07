# cat motif_filter.R | R --slave --vanilla --args <seamote output>
library(stringdist)

Args <- commandArgs()

myTab <- data.frame(read.table(Args[5]), row.names=1)

for( l in seq(3,7)){ 
	myVec <- myTab[which(nchar(rownames(myTab))==l),]
	myOrd <- myVec[order(myVec$V7),]
	
	myFirst <- rownames(myOrd)[1]
	mat <- stringdistmatrix(rownames(myOrd)[1], rownames(myOrd), method="hamming")
	m <- max(mat)
	sel <- unique(which(mat / m == 1 , arr.ind = TRUE)[,1])
	myVec <- myOrd[sel,]
	myOrd <- myVec[order(myVec$V7),]
	
	mySecond <- rownames(myOrd)[1]
	mat <- stringdistmatrix(rownames(myOrd)[1], rownames(myOrd), method="hamming")
	m <- max(mat)
	sel <- unique(which(mat / m == 1 , arr.ind = TRUE)[,1])
	myVec <- myOrd[sel,]
	myOrd <- myVec[order(myVec$V7),]

	myThird <- rownames(myOrd)[1]	
	
	top3 <- myTab[which( rownames(myTab) %in% c(myFirst, mySecond, myThird)),]
	
	if( l == 3 ){
		final <- top3
	}else{
		final <- rbind(final,top3)
	}
}
	
motSelection <- final[order(final$V7),]
	
write.table(motSelection, "selected_motifs.dat", col.names=FALSE, quote = FALSE)
