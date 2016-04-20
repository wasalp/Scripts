###################################################################
##Plot network
#####################################################################
#Koel et al.(2013)
#################################################################################################################################A&A&
library(igraph)
library(gridExtra)
library(ggplot2)
library(stringr)
setwd("/Users/louis/Desktop/spidermeme2")

filelist <-list.files(recursive = TRUE,no.. = TRUE,full.names = TRUE,pattern = "Edge_support.csv")

color <- 1
palette(rainbow(11))
sumframe <- data.frame(virtype= c(),average=c(),transitivity=c(),between=c(),central=c(),meandegree=c(),SpiderTreshold=c(),gentype=c())
pdf("graph.pdf", width=12, height=8)
for(files in filelist){
  
  pairAMat_adj_sortedLPOL3 <- read.csv(files)
  
  #lkl_res_adj_sor <-  lkl_res_adj_sorted[1:27,] 
  orderedpairMat <- pairAMat_adj_sortedLPOL3[order(-pairAMat_adj_sortedLPOL3$Site1...Site2),]
  
  top <- max(orderedpairMat$Site1...Site2)
  thres_list <- seq(0.05,top,0.025)

    for(f in thres_list){
    assortativity_degree1<- eccentricity1 <- strength1 <- diameter1 <- gentype1<-virtype1<-average1 <- transitivity1 <- between1 <- alpha1 <- central1 <- meandegree1 <- powerlaw_slope1 <- powerlaw_alpha1 <- c()
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
    #V(g)$label.main=c("treshold = ", f,strsplit(files,"/")[[1]][2],sttkrsplit(files,"/")[[1]][3],strsplit(files,"/")[[1]][5])
    V(g)$label.color="black"

    if(isTRUE(all.equal(f,0.85))|isTRUE(all.equal(f,0.1))){
      #plot(g,main=c("treshold = ", f,strsplit(files,"/")[[1]][2],sttkrsplit(files,"/")[[1]][3],strsplit(files,"/")[[1]][5]))
      # now draw network: solution 3'
      title <- paste("treshold =", f,strsplit(files,"/")[[1]][2],strsplit(files,"/")[[1]][3],strsplit(files,"/")[[1]][5])
      plot(g, layout=layout.fruchterman.reingold(g),vertex.label=NA, main=title)
    }
    meandegree1 <- c(meandegree1, mean(degree(g)))
    alpha1 <- c(alpha1,mean(alpha_centrality(g, nodes = V(g), loops = FALSE,exo = 1, weights = NULL, tol = 1e-07, sparse = TRUE)))	
    central1 <- c(central1,mean(closeness(g)))
    between1<- c(between1,mean(betweenness(g)))
    average1 <- c(average1,average.path.length(g, directed=TRUE, unconnected=TRUE))
    transitivity1 <- c(transitivity1,mean(transitivity(g)))
    diameter1 <- c(diameter1,mean(diameter(g)))
    strength1 <- c(strength1,mean(strength(g)))
    eccentricity1 <- c(eccentricity1,mean(eccentricity(g)))
    assortativity_degree1 <- c(assortativity_degree1,mean(assortativity_degree(g)))
    virtype1 <- c(virtype1,strsplit(files,"/")[[1]][2])
    gentype1 <- c(gentype1,strsplit(files,"/")[[1]][5])
    tempfram <- data.frame(virtype= virtype1,average=average1,transitivity=transitivity1,between=between1,
                           central=central1,meandegree=meandegree1,diameter=diameter1,strength=strength1,eccentricity=eccentricity1,assortativity_degree=assortativity_degree1, alpha=alpha1,SpiderTreshold=f,
                           gentype=gentype1)
    sumframe <- rbind(sumframe,tempfram)
    }
  

  
}
dev.off()
dsDNAsum <- subset(sumframe,sumframe$virtype == "dsDNA")
ssRNAsum <- subset(sumframe,sumframe$virtype == "ssRNA-")

dsDNAother <- subset(dsDNAsum,dsDNAsum$gentype == "other")
dsDNAGP <- subset(dsDNAsum,dsDNAsum$gentype == "GP")
ssRNAother <- subset(ssRNAsum,ssRNAsum$gentype == "other")
ssRNAGP <- subset(ssRNAsum,ssRNAsum$gentype == "GP")

dsDNAothermean <- data.frame(average=c(),transitivity=c(),between=c(),central=c(),meandegree=c(),diameter=c(),strength=c(),eccentricity=c(),assortativity_degree=c(),SpiderTreshold=c())
dsDNAGPmean <- data.frame(average=c(),transitivity=c(),between=c(),central=c(),meandegree=c(),diameter=c(),strength=c(),eccentricity=c(),assortativity_degree=c(),SpiderTreshold=c())
ssRNAothermean <- data.frame(average=c(),transitivity=c(),between=c(),central=c(),meandegree=c(),diameter=c(),strength=c(),eccentricity=c(),assortativity_degree=c(),SpiderTreshold=c())
ssRNAGPmean <- data.frame(average=c(),transitivity=c(),between=c(),central=c(),meandegree=c(),diameter=c(),strength=c(),eccentricity=c(),assortativity_degree=c(),SpiderTreshold=c())

for(i in seq(0.05,0.95,0.05)){
  tempdsother <- subset(dsDNAother,dsDNAother$SpiderTreshold == i)
  tempdsGP <- subset(dsDNAGP,dsDNAGP$SpiderTreshold == i)
  
  tempssother <- subset(ssRNAother,ssRNAother$SpiderTreshold == i)
  tempssGP <- subset(ssRNAGP,ssRNAGP$SpiderTreshold == i)
  
  #DNA
    #other
  meandegree1 <- mean(tempdsother$meandegree)
  alpha1 <- mean(tempdsother$alpha)	
  central1 <- mean(tempdsother$central)
  between1<- mean(tempdsother$between)
  average1 <- mean(tempdsother$average)
  transitivity1 <- mean(tempdsother$transitivity)
  diameter1 <- mean(tempdsother$diameter)
  strength1 <- mean(tempdsother$strength)
  eccentricity1 <- mean(tempdsother$eccentricity)
  assortativity_degree1 <- mean(tempdsother$assortativity_degree)
  tempfram <- data.frame(average=average1,transitivity=transitivity1,between=between1,central=central1,meandegree=meandegree1,diameter=diameter1,strength=strength1,eccentricity=eccentricity1,assortativity_degree=assortativity_degree1,SpiderTreshold=i)
  dsDNAothermean <- rbind(dsDNAothermean,tempfram)
    #GP
  meandegree1 <- mean(tempdsGP$meandegree)
  alpha1 <- mean(tempdsGP$alpha)	
  central1 <- mean(tempdsGP$central)
  between1<- mean(tempdsGP$between)
  average1 <- mean(tempdsGP$average)
  transitivity1 <- mean(tempdsGP$transitivity)
  diameter1 <- mean(tempdsGP$diameter)
  strength1 <- mean(tempdsGP$strength)
  eccentricity1 <- mean(tempdsGP$eccentricity)
  assortativity_degree1 <- mean(tempdsGP$assortativity_degree)
  tempfram <- data.frame(average=average1,transitivity=transitivity1,between=between1,central=central1,meandegree=meandegree1,diameter=diameter1,strength=strength1,eccentricity=eccentricity1,assortativity_degree=assortativity_degree1,SpiderTreshold=i)
  dsDNAGPmean <- rbind(dsDNAGPmean,tempfram)
  
  #RNA
    #other
  meandegree1 <- mean(tempssother$meandegree)
  alpha1 <- mean(tempssother$alpha)	
  central1 <- mean(tempssother$central)
  between1<- mean(tempssother$between)
  average1 <- mean(tempssother$average)
  transitivity1 <- mean(tempssother$transitivity)
  diameter1 <- mean(tempssother$diameter)
  strength1 <- mean(tempssother$strength)
  eccentricity1 <- mean(tempssother$eccentricity)
  assortativity_degree1 <- mean(tempssother$assortativity_degree)
  tempfram <- data.frame(average=average1,transitivity=transitivity1,between=between1,central=central1,meandegree=meandegree1,diameter=diameter1,strength=strength1,eccentricity=eccentricity1,assortativity_degree=assortativity_degree1,SpiderTreshold=i)
  ssRNAothermean <- rbind(ssRNAothermean,tempfram)
    #GP
  meandegree1 <- mean(tempssGP$meandegree)
  alpha1 <- mean(tempssGP$alpha)	
  central1 <- mean(tempssGP$central)
  between1<- mean(tempssGP$between)
  average1 <- mean(tempssGP$average)
  transitivity1 <- mean(tempssGP$transitivity)
  diameter1 <- mean(tempssGP$diameter)
  strength1 <- mean(tempssGP$strength)
  eccentricity1 <- mean(tempssGP$eccentricity)
  assortativity_degree1 <- mean(tempssGP$assortativity_degree)
  tempfram <- data.frame(average=average1,transitivity=transitivity1,between=between1,central=central1,meandegree=meandegree1,diameter=diameter1,strength=strength1,eccentricity=eccentricity1,assortativity_degree=assortativity_degree1,SpiderTreshold=i)
  ssRNAGPmean <- rbind(ssRNAGPmean,tempfram)
}
pdf("network.pdf", width=24, height=16)

    q1 <- ggplot()+
    geom_smooth(data=dsDNAothermean, aes(x=SpiderTreshold, y=dsDNAothermean[1],xmin=0,xmax=1),
                                   fill = 1,colour=1,size=1) +
    geom_smooth(data=dsDNAGPmean, aes(x=SpiderTreshold, y=dsDNAGPmean[1],xmin=0,xmax=1),
                                   fill =2,colour=2, size=1) + 
    geom_smooth(data=ssRNAothermean, aes(x=SpiderTreshold, y=ssRNAothermean[1],xmin=0,xmax=1),
                           fill = 4,colour=4,size=1) +
    geom_smooth(data=ssRNAGPmean, aes(x=SpiderTreshold, y=ssRNAGPmean[1],xmin=0,xmax=1),
                           fill =6,colour=6, size=1) +
      labs(title=colnames(dsDNAothermean[1]),y=colnames(dsDNAothermean[1]),x="Posterior probability threshold")+
      theme_bw()

    q2 <- ggplot()+
      geom_smooth(data=dsDNAothermean, aes(x=SpiderTreshold, y=dsDNAothermean[2],xmin=0,xmax=1),
                  fill = 1,colour=1,size=1) +
      geom_smooth(data=dsDNAGPmean, aes(x=SpiderTreshold, y=dsDNAGPmean[2],xmin=0,xmax=1),
                  fill =2,colour=2, size=1) + 
      geom_smooth(data=ssRNAothermean, aes(x=SpiderTreshold, y=ssRNAothermean[2],xmin=0,xmax=1),
                  fill = 4,colour=4,size=1) +
      geom_smooth(data=ssRNAGPmean, aes(x=SpiderTreshold, y=ssRNAGPmean[2],xmin=0,xmax=1),
                  fill =6,colour=6, size=1) +
      labs(title=colnames(dsDNAothermean[2]),y=colnames(dsDNAothermean[2]),x="Posterior probability threshold")+
      theme_bw()

    q3 <- ggplot()+
      geom_smooth(data=dsDNAothermean, aes(x=SpiderTreshold, y=dsDNAothermean[3],xmin=0,xmax=1),
                  fill = 1,colour=1,size=1) +
      geom_smooth(data=dsDNAGPmean, aes(x=SpiderTreshold, y=dsDNAGPmean[3],xmin=0,xmax=1),
                  fill =2,colour=2, size=1) + 
      geom_smooth(data=ssRNAothermean, aes(x=SpiderTreshold, y=ssRNAothermean[3],xmin=0,xmax=1),
                  fill = 4,colour=4,size=1) +
      geom_smooth(data=ssRNAGPmean, aes(x=SpiderTreshold, y=ssRNAGPmean[3],xmin=0,xmax=1),
                  fill =6,colour=6, size=1) +
      labs(title=colnames(dsDNAothermean[3]),y=colnames(dsDNAothermean[3]),x="Posterior probability threshold")+
      theme_bw()

    q4 <- ggplot()+
      geom_smooth(data=dsDNAothermean, aes(x=SpiderTreshold, y=dsDNAothermean[4],xmin=0,xmax=1),
                  fill = 1,colour=1,size=1) +
      geom_smooth(data=dsDNAGPmean, aes(x=SpiderTreshold, y=dsDNAGPmean[4],xmin=0,xmax=1),
                  fill =2,colour=2, size=1) + 
      geom_smooth(data=ssRNAothermean, aes(x=SpiderTreshold, y=ssRNAothermean[4],xmin=0,xmax=1),
                  fill = 4,colour=4,size=1) +
      geom_smooth(data=ssRNAGPmean, aes(x=SpiderTreshold, y=ssRNAGPmean[4],xmin=0,xmax=1),
                  fill =6,colour=6, size=1) +
      labs(title=colnames(dsDNAothermean[4]),y=colnames(dsDNAothermean[4]),x="Posterior probability threshold")+
      theme_bw()
  

    q5 <- ggplot()+
      geom_smooth(data=dsDNAothermean, aes(x=SpiderTreshold, y=dsDNAothermean[5],xmin=0,xmax=1),
                  fill = 1,colour=1,size=1) +
      geom_smooth(data=dsDNAGPmean, aes(x=SpiderTreshold, y=dsDNAGPmean[5],xmin=0,xmax=1),
                  fill =2,colour=2, size=1) + 
      geom_smooth(data=ssRNAothermean, aes(x=SpiderTreshold, y=ssRNAothermean[5],xmin=0,xmax=1),
                  fill = 4,colour=4,size=1) +
      geom_smooth(data=ssRNAGPmean, aes(x=SpiderTreshold, y=ssRNAGPmean[5],xmin=0,xmax=1),
                  fill =6,colour=6, size=1) +
      labs(title=colnames(dsDNAothermean[5]),y=colnames(dsDNAothermean[5]),x="Posterior probability threshold")+
      theme_bw()
  
  
  q6 <- ggplot()+
    geom_smooth(data=dsDNAothermean, aes(x=SpiderTreshold, y=dsDNAothermean[6],xmin=0,xmax=1),
                fill = 1,colour=1,size=1) +
    geom_smooth(data=dsDNAGPmean, aes(x=SpiderTreshold, y=dsDNAGPmean[6],xmin=0,xmax=1),
                fill =2,colour=2, size=1) + 
    geom_smooth(data=ssRNAothermean, aes(x=SpiderTreshold, y=ssRNAothermean[6],xmin=0,xmax=1),
                fill = 4,colour=4,size=1) +
    geom_smooth(data=ssRNAGPmean, aes(x=SpiderTreshold, y=ssRNAGPmean[6],xmin=0,xmax=1),
                fill =6,colour=6, size=1) +
    labs(title=colnames(dsDNAothermean[6]),y=colnames(dsDNAothermean[6]),x="Posterior probability threshold")+
    theme_bw()
  
  q7 <- ggplot()+
    geom_smooth(data=dsDNAothermean, aes(x=SpiderTreshold, y=dsDNAothermean[7],xmin=0,xmax=1),
                fill = 1,colour=1,size=1) +
    geom_smooth(data=dsDNAGPmean, aes(x=SpiderTreshold, y=dsDNAGPmean[7],xmin=0,xmax=1),
                fill =2,colour=2, size=1) + 
    geom_smooth(data=ssRNAothermean, aes(x=SpiderTreshold, y=ssRNAothermean[7],xmin=0,xmax=1),
                fill = 4,colour=4,size=1) +
    geom_smooth(data=ssRNAGPmean, aes(x=SpiderTreshold, y=ssRNAGPmean[7],xmin=0,xmax=1),
                fill =6,colour=6, size=1) +
    labs(title=colnames(dsDNAothermean[7]),y=colnames(dsDNAothermean[7]),x="Posterior probability threshold")+
    theme_bw()
  
  q8 <- ggplot()+
    geom_smooth(data=dsDNAothermean, aes(x=SpiderTreshold, y=dsDNAothermean[8],xmin=0,xmax=1),
                fill = 1,colour=1,size=1) +
    geom_smooth(data=dsDNAGPmean, aes(x=SpiderTreshold, y=dsDNAGPmean[8],xmin=0,xmax=1),
                fill =2,colour=2, size=1) + 
    geom_smooth(data=ssRNAothermean, aes(x=SpiderTreshold, y=ssRNAothermean[8],xmin=0,xmax=1),
                fill = 4,colour=4,size=1) +
    geom_smooth(data=ssRNAGPmean, aes(x=SpiderTreshold, y=ssRNAGPmean[8],xmin=0,xmax=1),
                fill =6,colour=6, size=1) +
    labs(title=colnames(dsDNAothermean[8]),y=colnames(dsDNAothermean[8]),x="Posterior probability threshold")+
    theme_bw()
  
  q9 <- ggplot()+
    geom_smooth(data=dsDNAothermean, aes(x=SpiderTreshold, y=dsDNAothermean[9],xmin=0,xmax=1),
                fill = 1,colour=1,size=1) +
    geom_smooth(data=dsDNAGPmean, aes(x=SpiderTreshold, y=dsDNAGPmean[9],xmin=0,xmax=1),
                fill =2,colour=2, size=1) + 
    geom_smooth(data=ssRNAothermean, aes(x=SpiderTreshold, y=ssRNAothermean[9],xmin=0,xmax=1),
                fill = 4,colour=4,size=1) +
    geom_smooth(data=ssRNAGPmean, aes(x=SpiderTreshold, y=ssRNAGPmean[9],xmin=0,xmax=1),
                fill =6,colour=6, size=1) +
    labs(title=colnames(dsDNAothermean[9]),y=colnames(dsDNAothermean[9]),x="Posterior probability threshold")+
    theme_bw()
  
          


grid.arrange(q1,q2,q3,q4,q5,q6,q7,q8,q9)
dev.off()

