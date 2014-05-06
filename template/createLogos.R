# Rscript createLogo.R --vanilla --slave

library(Biostrings)
library(seqLogo)

motifTable <- read.table("outputs/best_motifs.txt", row.names=1)

cmd = with(motifTable, paste("grep -o -e ",V2," outputs/positive.seq > matches.txt"))

for( l in seq(1,length(cmd),1)){

	system(cmd[l])
	myMatches <- read.table("matches.txt")
	mySet <- DNAStringSet(myMatches$V1)
	m <- consensusMatrix(mySet, as.prob=TRUE, baseOnly=TRUE)[1:4,]
	
	len=mySet[[1]]@length*80+240
	nameF <- paste("logos/",rownames(motifTable)[l],"_logo.png", sep="")
	png(nameF, width=len, height=240, res=72, pointsize=10, type=c("cairo-png"))
	seqLogo(m, xfontsize=14, yfontsize=14)
	dev.off()

	nameF <- paste("logos/",rownames(motifTable)[l],"_pwm.txt", sep="")
	write(paste("IUPAC motif:", rownames(motifTable)[l], "\n"), file=nameF, append = FALSE)
	write(paste("Consensus sequence:", consensusString(m), "\n"), file=nameF, append = TRUE)
	write("Position weighted matrix (pwm):\n#     A        C        G        T", file=nameF, append = TRUE)
	write.table(format(t(m), digits=6), file=nameF, append = TRUE, col.names = FALSE, quote = FALSE)
}