library(boot)
library(multicore)

DF <- as.data.frame(read.table("/Users/fagostini/algorithm-Motif_enrichment_testing/tmp/motifs_17nt.dat", row.names=1))
p <- read.table("/Users/fagostini/algorithm-Motif_enrichment_testing/tmp/motifs_17nt_Pdistribution.dat", row.names=1)
n <- read.table("/Users/fagostini/algorithm-Motif_enrichment_testing/tmp/motifs_17nt_Ndistribution.dat", row.names=1)

plen <- length(p[1,])
nlen <- length(n[1,])

enrichment <- function(d){	
	pool <- sample(d, plen, replace=FALSE)
	ssum <- sum(pool[1:plen])
}

bootstrap <- function(x){
	vec <- c(as.vector(unlist(p[x,])), as.vector(unlist(n[x,])))
	psum <- sum(vec[1:plen])
	bres <- boot(vec, enrichment, R=1000, sim="parametric")
	1-sum(psum >= bres$t)/1000
}

# DF$boo <- sapply(c(seq(1,dim(p)[1],1)),bootstrap)
DF$boo <- unlist(mclapply(c(seq(1,dim(p)[1],1)),bootstrap, mc.cores=4))
