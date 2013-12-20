b <- read.table("b", row.names=1)
c <- read.table("c", row.names=1)
d <- which( b$V1002-c$V1002 > 0 )

Pvec <- rep(0,length=dim(b)[2]-1)
Nvec <- rep(0,length=dim(b)[2]-1)
selected <- vector(length=dim(b)[2]-1)
excluded <- vector(length=dim(b)[2]-1)
s=1
e=1
for (m in d){
	if(b[m,1001] > c[m,1001]){
#		oM <- matrix(c(sum(Pvec), sum(Nvec), 1000-sum(Pvec), 1000-sum(Nvec)), 2, 2)
#		oF <- fisher.test(oM, alternative="greater")$p.value
		Ptmp <- Pvec
		Ptmp <- Ptmp + b[m,1:1000]
		Ptmp[which(Ptmp == 2)] = 1
		Ntmp <- Nvec
		Ntmp <- Ntmp + c[m,1:1000]
		Ntmp[which(Ntmp == 2)] = 1
#		tM <- matrix(c(sum(Ptmp), sum(Ntmp), 1000-sum(Ptmp), 1000-sum(Ntmp)), 2, 2)
#		tF <- fisher.test(tM, alternative="greater")$p.value
#		if( tF < oF ){
		if( (sum(Ptmp)-sum(Ntmp)) > (sum(Pvec)-sum(Nvec)) ){
			selected[s] <- m
			s <- s+1
			Pvec <- Ptmp
			Nvec <- Ntmp
		}else{
			excluded[e] <- m
			e <- e+1
		}
		print((s+e-2)/length(d)*100)
	}
}

write.table(rownames(b[selected,]), "Rmotifs.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)