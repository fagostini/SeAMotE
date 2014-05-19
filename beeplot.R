# Rscript beeplot.R --vanilla --slave
library(beeswarm)

d <- read.table("dreme_res")
s <- read.table("seamote_res")
e <- read.table("decod_res")
x <- read.table("xxmotif_res")

da <- data.frame("DREME", d)
sa <- data.frame("SeAMotE", s)
e$V1[which(e$V1<0)] = 0
ep <- vector(mode="numeric",length=length(e$V1))
for (i in seq(1,length(e$V1))){
	ep[i] <- fisher.test(matrix(as.numeric(as.character(e[i,1:4])),2,2))$p.val
	if( ep[i] == 0 ){
		ep[i] = 5e-324
	}

}
ea <- data.frame("DECOD", e$V5, ep)
xp <- vector(mode="numeric",length=length(x$V1))
for (i in seq(1,length(x$V1))){
	xp[i] <- fisher.test(matrix(as.numeric(as.character(x[i,1:4])),2,2))$p.val
	if( xp[i] == 0 ){
		xp[i] = 5e-324
	}
}
xa <- data.frame("XXmotif", x$V5, xp)

colnames(da) <- c("ID", "V1", "V2")
colnames(sa) <- c("ID", "V1", "V2")
colnames(ea) <- c("ID", "V1", "V2")
colnames(xa) <- c("ID", "V1", "V2")

ds <- rbind(da,sa)
ds1 <- rbind(da,sa,ea,xa)
# da$col <- "red"
# sa$col <- "darkcyan"
# ds1 <- rbind(da,sa)
# ds2 <- rbind(da,sa)
# ds1$ID <- "disc"
# ds1$V1 <- ds1$V1*100
# ds1$V2 <- 0
# ds2$ID <- "pval"
# ds2$V1 <- -log10(ds2$V2)
# ds2$V2 <- 0

pdf("plot_DiscSign.pdf", width=12, height=6, pointsize=12)
par(mfrow=c(1,2), bg="white")
beeswarm(    V1*100 ~ ID, data=ds1, method="center", pch = 16, xlab="", ylab="Discrimination (%)", ylim=c(0,100), col=c(2,"darkcyan",1,3), main="Motif discrimination ability", cex=0.75);
beeswarm(-log10(V2) ~ ID, data=ds1, method="center", pch=16, xlab="", ylab="-log10(p-value) (Fisher's exact test)", ylim=c(0,325), col=c(2,"darkcyan",1,3), main="Motif significance", cex=0.75)
dev.off()

# wilcox.test(s$V1, d$V1, alternative="greater")$p.val
# wilcox.test(s$V1, e$V5, alternative="greater")$p.val
# wilcox.test(s$V1, x$V5, alternative="greater")$p.val