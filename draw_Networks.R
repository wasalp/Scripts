###################################################################
##Plot network
#####################################################################
#Koel et al.(2013)
#################################################################################################################################A&A&
library(igraph)
library(gridExtra)
library(ggplot2)
setwd("/Users/louis/Desktop/spidermeme2")

filelist <-list.files(recursive = TRUE,no.. = TRUE,full.names = TRUE,pattern = "Edge_support.csv")

color <- 1
palette(rainbow(11))
sumframe <- data.frame(virtype= c(),average=c(),transitivity=c(),between=c(),central=c(),meandegree=c(),SpiderTreshold=c())
for(files in filelist){
  
  pairAMat_adj_sortedLPOL3 <- read.csv(files)
  
  #lkl_res_adj_sor <-  lkl_res_adj_sorted[1:27,] 
  orderedpairMat <- pairAMat_adj_sortedLPOL3[order(-pairAMat_adj_sortedLPOL3$Site1...Site2),]
  
  top <- max(orderedpairMat$Site1...Site2)
  thres_list <- seq(0.05,top,0.05)

    for(f in thres_list){
    virtype1<-average1 <- transitivity1 <- between1 <- alpha1 <- central1 <- meandegree1 <- powerlaw_slope1 <- powerlaw_alpha1 <- c()
    tempsub <- subset(orderedpairMat,orderedpairMat$Site1...Site2 >= f)
    pairAMat_adj_sor <- tempsub
    
    
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
    
    site1 <- corrected_site_1
    
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

    g <- graph.adjacency(mat_form, mode="undirected", weighted=TRUE)
    log10_weights <- -log10(E(g)$weight)
    E(g)$width <- log10_weights + min(log10_weights) + 1
    E(g)$color <- "black" 
    V(g)$color <-"white"
    V(g)$label.color="black"
    #plot(g,main=c("treshold = ", tf))
    # now draw network: solution 3'
    #plot(g, layout=layout.fruchterman.reingold(g))
    meandegree1 <- c(meandegree1, mean(degree(g)))
    alpha1 <- c(alpha1,mean(alpha_centrality(g, nodes = V(g), loops = FALSE,exo = 1, weights = NULL, tol = 1e-07, sparse = TRUE)))	
    central1 <- c(central1,mean(closeness(g)))
    between1<- c(between1,mean(betweenness(g)))
    average1 <- c(average1,average.path.length(g, directed=TRUE, unconnected=TRUE))
    transitivity1 <- c(transitivity1,mean(transitivity(g)))
    virtype1 <- c(virtype1,strsplit(files,"/")[[1]][2])
    tempfram <- data.frame(virtype= virtype1,average=average1,transitivity=transitivity1,between=between1,central=central1,meandegree=meandegree1,alpha=alpha1,SpiderTreshold=f)
    sumframe <- rbind(sumframe,tempfram)
    }
  

  
}
dsDNAsum <- subset(sumframe,sumframe$virtype == "dsDNA")
ssRNAsum <- subset(sumframe,sumframe$virtype == "ssRNA-")
dsDNAmean <- data.frame(average=c(),transitivity=c(),between=c(),central=c(),meandegree=c(),SpiderTreshold=c())
ssRNAmean <- data.frame(average=c(),transitivity=c(),between=c(),central=c(),meandegree=c(),SpiderTreshold=c())
for(i in seq(0.05,0.95,0.05)){
  tempds <- subset(dsDNAsum,dsDNAsum$SpiderTreshold == i)
  tempss <- subset(ssRNAsum,dsDNAsum$SpiderTreshold == i)
  #DNA
  meandegree1 <- mean(tempds$meandegree)
  alpha1 <- mean(tempds$alpha)	
  central1 <- mean(tempds$central)
  between1<- mean(tempds$between)
  average1 <- mean(tempds$average)
  transitivity1 <- mean(tempds$transitivity)
  tempfram <- data.frame(average=average1,transitivity=transitivity1,between=between1,central=central1,meandegree=meandegree1,SpiderTreshold=i)
  dsDNAmean <- rbind(dsDNAmean,tempfram)
  #RNA
  meandegree1 <- mean(tempss$meandegree)
  alpha1 <- mean(tempss$alpha)	
  central1 <- mean(tempss$central)
  between1<- mean(tempss$between)
  average1 <- mean(tempss$average)
  transitivity1 <- mean(tempss$transitivity)
  tempfram <- data.frame(average=average1,transitivity=transitivity1,between=between1,central=central1,meandegree=meandegree1,SpiderTreshold=i)
  ssRNAmean <- rbind(ssRNAmean,tempfram)
}
pdf("network.pdf", width=24, height=16)
panel <- c()
for(i in c(1:5)){
q <- ggplot()
q <- q + geom_smooth(data=dsDNAmean, aes(x=SpiderTreshold, y=dsDNAmean[i]), fill=1,
          colour=1, size=1) 
q <- q + geom_smooth(data=ssRNAmean, aes(x=SpiderTreshold, y=ssRNAmean[i]), fill=7,
                     colour=7, size=1) 
q <- q + labs(title=colnames(dsDNAmean[i]),y=c(colnames(dsDNAmean[i])))

}
print(panel)
grid.arrange(panel)
dev.off()