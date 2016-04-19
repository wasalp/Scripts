library(igraph)
library(robust)

#g <- barabasi.game(10000)
d <- degree(g, mode="in")
hd <- hist(d)
pdf("gene_name.pdf", width=6, height=4)
plot(log10(hd$counts) ~ log10(seq(1,length(hd$counts))))
defv.off()
lm_pl <- lmrob( log10(hd$counts+.1) ~ log10(seq(1,length(hd$counts))) )
summary(lm_pl)

