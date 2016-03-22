###################################################################
##Plot network
#####################################################################
#Koel et al.(2013)
#################################################################################################################################A&A&

setwd("/Users/louis/Desktop/dsDNA/hpv4/YP_001129462.1")
pairAMat_adj_sortedLPOL3 <- read.csv("Edge_support1.csv")

#lkl_res_adj_sor <-  lkl_res_adj_sorted[1:27,] 
pairAMat_adj_sor <- pairAMat_adj_sortedLPOL3[1:100,]


#distance <- c(36.919, 49.486, 19.499, 21.234, 24.974, 32.731, 13.411, 27.945, 35.015, 9.925, 24.456, 30.895, 31.242, 16.831, 4.624, 52.618, 17.176, 27.828, 5.496, 12.111, 30.678, 21.977, 27.456, 18.106, 21.057, 71.646, 30.682)

#lkl_res_adj_sordetected <- cbind(lkl_res_adj_sor, distance )

#log10pvalue <- -log10(lkl_res_adj_sordetected$pvalue_fdr)


#log10pvalue <- round(log10pvalue, 3)
#lkl_res_adj_sordetected <- cbind(lkl_res_adj_sordetected, log10pvalue)

#corrected_site_2 <- as.numeric(lkl_res_adj_sordetected$site_2)-16
#corrected_site_1 <- as.numeric(lkl_res_adj_sordetected$site_1)-16
corrected_site_1 <- pairAMat_adj_sor$Site1
corrected_site_2 <- pairAMat_adj_sor$Site2



#combined_sites <- paste(corrected_site_1,corrected_site_2, sep="vs")
combined_sites <- paste(corrected_site_1, corrected_site_2, sep="vs")

#lkl_res_adj_sordetected <- cbind(lkl_res_adj_sordetected, combined_sites)
pairAMat_adj_sor <- cbind(pairAMat_adj_sor, combined_sites)


# generate random data
npairs <- 25

site1 <-  corrected_site_1

site2 <- corrected_site_2



#probs <- lkl_res_adj_sordetected$pvalue_fdr
probs <- pairAMat_adj_sor$Site1...Site2


rand_pairs <- cbind(site1, site2, probs)
write.csv(rand_pairs, "rand_pairs.csv", quote=F, row.names=F)

# read network from file
net_dat <- read.csv("rand_pairs.csv")
net_dat
n_pairs <- length(net_dat$site1)

list_of_uniq_sites <- unique(union(net_dat$site1, net_dat$site2))
list_len <- length(list_of_uniq_sites)

mat_form <- matrix(0, list_len, list_len) 
colnames(mat_form) <- rownames(mat_form) <- list_of_uniq_sites

# fill up matrix -- non-sym form
for(i in 1:n_pairs){
	ii <- which(colnames(mat_form) == net_dat$site1[i])
	jj <- which(colnames(mat_form) == net_dat$site2[i])
	mat_form[ii,jj] <- net_dat$probs[i]
	#mat_form[ii,jj] <- mat_form[jj,ii] <- net_dat$probs[i]
}
mat_form

# now draw network: solution 1
#library(diagram)
#plotmat(mat_form, curve=0, lwd=-2*log10(mat_form), box.type="multi", arr.type="curved", arr.width=0, cex=0.0001)

# now draw network: solution 2
#library("qgraph")
#qgraph(mat_form,esize=5,gray=TRUE)

# now draw network: solution 3
library(igraph)
g <- graph.adjacency(mat_form, mode="undirected", weighted=TRUE)
log10_weights <- -log10(E(g)$weight)
E(g)$width <- log10_weights + min(log10_weights) + 1
E(g)$color <- "black" 
V(g)$color <-"white"
V(g)$label.color="black"
plot(g)

# now draw network: solution 3'
plot(g, layout=layout.fruchterman.reingold(g))







