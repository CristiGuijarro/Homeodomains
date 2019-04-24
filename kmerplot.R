library(ggplot2)
kmers <- read.table('/home/cristig/Desktop/Homeodomains/kmerProfileHbx5.tsv', header=FALSE, sep = "\t")
kmerp <- ggplot(kmers, aes(x=kmers$V4, y=kmers$V5)) + geom_point(stat = "identity", position=position_dodge(width=1), aes(colour = factor(kmers$V1)), size=0.75, shape = 20, alpha = 0.75) + theme(panel.spacing = unit(1, "lines")) + theme_minimal() + theme(legend.position="bottom", panel.grid = element_blank(), axis.text.y = element_blank()) + theme(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_text(size = 8)) + facet_grid(kmers$V3 ~ ., scales="free_x") + xlab("Amino acid pair") + ylab("Median occurance by Phylum") + labs(color="") + scale_colour_manual(values = c("Consensus" = "black", "Annelida" = "violetred1", "Arthropoda" = "darkorange", "Brachiopoda" = "violetred4", "Chordata" = "darkgreen", "Cnidaria" = "rosybrown", "Ctenophora" = "rosybrown1", "Echinodermata" = "darkolivegreen3", "Hemichordata" = "darkolivegreen4", "Mollusca" = "violet", "Nematoda" = "chocolate", "Orthonectida" = "steelblue", "Placozoa" = "rosybrown2", "Platyhelminthes" = "orchid4", "Porifera" = "rosybrown4", "Rotifera" = "slateblue2", "Tardigrada" = "tan1"))
kmerp

library(ggplot2)
library(ggforce)
basevar <- read.table('/home/cristig/Desktop/Homeodomains/baseVariable.tsv', header=FALSE, sep = "\t")
bv <- ggplot(basevar, aes(x=basevar$V3, y=basevar$V5)) + geom_bar(stat="identity", position="stack", aes(fill=factor(basevar$V4)))
#bv <- bv + facet_grid(basevar$V1 ~ basevar$V2, scales="free_x")
bv <- bv + facet_grid_paginate(basevar$V1 ~ basevar$V2, ncol = 2, nrow = 16, page = 1)
n_pages(bv)
bv

pdf("baseVar.pdf")
for(i in 1:13){
	print(ggplot(basevar, aes(x=basevar$V3, y=basevar$V5)) + geom_bar(stat="identity", position="stack", aes(fill=factor(basevar$V4))) +
	theme(panel.spacing = unit(1, "lines")) + theme_minimal() + theme(panel.grid = element_blank(), axis.text.y = element_blank()) +
	theme(strip.background = element_blank(), strip.text.y = element_text(size = 8, angle = 360)) +
	xlab("Position") + ylab("Ratio of positions (%)") + labs(fill="Amino acid") +
	facet_grid_paginate(basevar$V1 ~ basevar$V2, ncol = 1, nrow = 16, page = i))
}
dev.off()

bv2 <- ggplot(basevar, aes(x=basevar$V3, y=basevar$V5) +      
  geom_logo(aes(label=basevar$V4, fill=basevar$V4))
bv2
