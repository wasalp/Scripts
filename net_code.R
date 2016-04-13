library(igraph)
library(ggplot2)
library(gridExtra)


###########################################
# maybe not the fastest way to threshold matrix b/c I/O overhead, but it works
csv2mat_form <- function(filename, threshold = .95){
	# reads network from file and format data for downstream analyses
	net_dat <- read.csv(filename)
	n_pairs <- length(net_dat[,1])
	
	list_of_uniq_sites <- unique(union(net_dat[,1], net_dat[,2]))
	list_len <- length(list_of_uniq_sites)
	
	mat_form <- matrix(0, list_len, list_len) 
	colnames(mat_form) <- rownames(mat_form) <- list_of_uniq_sites
	
	# fills up matrix -- non-sym form
	for(i in 1:n_pairs){
		ii <- which(colnames(mat_form) == net_dat[i,1])
		jj <- which(colnames(mat_form) == net_dat[i,2])
		if(net_dat[i,5] > threshold){
			mat_form[ii,jj] <- net_dat[i,5]
		}else{
			mat_form[ii,jj] <- 0
		}
	}
	return(mat_form)
}


											#######################
											####### HA data #######
											#######################

###########################################
# posterior probability threshold analysis
thres_list <- seq(.01,.99,.01)
#thres_list <- seq(.01,.99,.2)
# pandemic clade
meanStrength1 <- meanEccentric1 <- diameter1 <- dyad1 <- authorScore1 <- assort_deg1 <- avepathlen1 <- transitivity1 <- between1 <- alpha1 <- central1 <- meandegree1 <- powerlaw_alpha1 <- c()
for(i in 1:length(thres_list)){
	print(paste0("Now doing threshold ", thres_list[i]))
	my_mat_form <- csv2mat_form("ha_pandemic/LAST_FILE_PATH.csv", thres_list[i])
	g <- graph.adjacency(my_mat_form, mode="undirected", weighted=T)
	d <- degree(g)
	meandegree1[i] <- mean(d)
	alpha1[i] <- mean(alpha_centrality(g, nodes = V(g), loops = F,exo = 1, weights = NULL, tol = 1e-07, sparse = T))	
	central1[i] <- mean(closeness(g))
	between1[i]<- mean(betweenness(g))
	avepathlen1[i] <- average.path.length(g, directed=F, unconnected=T)
	transitivity1[i] <- mean(transitivity(g))
	powerlaw_alpha1[i] <- fit_power_law(d+1, 10)$alpha
	assort_deg1[i] <- assortativity_degree(g)
	authorScore1[i] <- authority_score(g)$value # same as hub_score(g)$value, Kleinberg’s hub centrality scores
	dyad1[i] <- dyad.census(g)$mut
	diameter1[i] <- diameter(g)
	meanEccentric1[i] <- mean(eccentricity(g))
	meanStrength1[i] <- mean(strength(g))
}
# seasonal clade
meanStrength2 <- meanEccentric2 <- diameter2 <- dyad2 <- authorScore2 <- assort_deg2 <- avepathlen2 <-transitivity2 <- between2 <- alpha2 <- central2 <- meandegree2 <- powerlaw_alpha2 <- c()
for(i in 1:length(thres_list)){
	print(paste0("Now doing threshold ", thres_list[i]))
	my_mat_form <- csv2mat_form("ha_seasonal/LAST_FILE_PATH.csv", thres_list[i])
	g <- graph.adjacency(my_mat_form, mode="undirected", weighted=T)
	d <- degree(g)
	meandegree2[i] <- mean(d)
	alpha2[i] <- mean(alpha_centrality(g, nodes = V(g), loops = F, exo = 1, weights = NULL, tol = 1e-07, sparse = T))
	central2[i] <- mean(closeness(g))
	between2[i]<- mean(betweenness(g))
	avepathlen2[i] <- average.path.length(g, directed=F, unconnected=T)
	transitivity2[i] <- mean(transitivity(g))
	powerlaw_alpha2[i] <- fit_power_law(d+1, 10)$alpha
	assort_deg2[i] <- assortativity_degree(g)
	authorScore2[i] <- authority_score(g)$value
	dyad2[i] <- dyad.census(g)$mut
	diameter2[i] <- diameter(g)
	meanEccentric2[i] <- mean(eccentricity(g))
	meanStrength2[i] <- mean(strength(g))
}

save.image("ha.RData")


dat_ha_pandemic <- data.frame(thres_list, alpha1, meanStrength1, meanEccentric1, diameter1, dyad1, authorScore1, assort_deg1, avepathlen1, transitivity1, between1, central1, meandegree1, powerlaw_alpha1)
dat_ha_seasonal <- data.frame(thres_list, alpha2, meanStrength2, meanEccentric2, diameter2, dyad2, authorScore2, assort_deg2, avepathlen2, transitivity2, between2, central2, meandegree2, powerlaw_alpha2)

###########################################	
# plots
##########
panel <- list()

# Alpha centrality 
panel[[1]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Alpha centrality") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_ha_pandemic, aes(x=thres_list, y=alpha1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_ha_seasonal, aes(x=thres_list, y=alpha2), fill="blue",
      colour="darkblue", size=1)

# meandegree 
panel[[2]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Average degree") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_ha_pandemic, aes(x=thres_list, y= meandegree1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_ha_seasonal, aes(x=thres_list, y= meandegree2), fill="blue",
      colour="darkblue", size=1)

# meanStrength 
panel[[3]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Average strength") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_ha_pandemic, aes(x=thres_list, y= meanStrength1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_ha_seasonal, aes(x=thres_list, y= meanStrength2), fill="blue",
      colour="darkblue", size=1)

# meanEccentric 
panel[[4]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Average eccentricity") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_ha_pandemic, aes(x=thres_list, y= meanEccentric1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_ha_seasonal, aes(x=thres_list, y= meanEccentric2), fill="blue",
      colour="darkblue", size=1)

# diameter 
panel[[5]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Average diameter") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_ha_pandemic, aes(x=thres_list, y= diameter1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_ha_seasonal, aes(x=thres_list, y= diameter2), fill="blue",
      colour="darkblue", size=1)

# dyad 
panel[[6]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Number of dyads with mutual connections") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_ha_pandemic, aes(x=thres_list, y= dyad1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_ha_seasonal, aes(x=thres_list, y= dyad2), fill="blue",
      colour="darkblue", size=1)

# authorScore 
panel[[7]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Authority score") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_ha_pandemic, aes(x=thres_list, y= authorScore1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_ha_seasonal, aes(x=thres_list, y= authorScore2), fill="blue",
      colour="darkblue", size=1)

# assort_deg 
panel[[8]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Assortative degree") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_ha_pandemic, aes(x=thres_list, y= assort_deg1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_ha_seasonal, aes(x=thres_list, y= assort_deg2), fill="blue",
      colour="darkblue", size=1)

# avepathlen 
panel[[9]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Average path length") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_ha_pandemic, aes(x=thres_list, y= avepathlen1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_ha_seasonal, aes(x=thres_list, y= avepathlen2), fill="blue",
      colour="darkblue", size=1)

# transitivity2 
panel[[10]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Transitivity") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_ha_pandemic, aes(x=thres_list, y= transitivity1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_ha_seasonal, aes(x=thres_list, y= transitivity2), fill="blue",
      colour="darkblue", size=1)

# between 
panel[[11]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Average betweenness") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_ha_pandemic, aes(x=thres_list, y= between1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_ha_seasonal, aes(x=thres_list, y= between2), fill="blue",
      colour="darkblue", size=1)

# central 
panel[[12]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Average centrality") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_ha_pandemic, aes(x=thres_list, y= central1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_ha_seasonal, aes(x=thres_list, y= central2), fill="blue",
      colour="darkblue", size=1)

# powerlaw_alpha 
panel[[13]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Power law (alpha)") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_ha_pandemic, aes(x=thres_list, y= powerlaw_alpha1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_ha_seasonal, aes(x=thres_list, y= powerlaw_alpha2), fill="blue",
      colour="darkblue", size=1)


# print all figures
pdf("fig_influenza_ha.pdf" ,width=24, height=16)
do.call("grid.arrange", panel)  
dev.off()



											#######################
											####### NA data #######
											#######################

###########################################
# posterior probability threshold analysis
thres_list <- seq(.01,.99,.01)
#thres_list <- seq(.01,.99,.2)
# pandemic clade
meanStrength1 <- meanEccentric1 <- diameter1 <- dyad1 <- authorScore1 <- assort_deg1 <- avepathlen1 <- transitivity1 <- between1 <- alpha1 <- central1 <- meandegree1 <- powerlaw_alpha1 <- c()
for(i in 1:length(thres_list)){
	print(paste0("Now doing threshold ", thres_list[i]))
	my_mat_form <- csv2mat_form("na_pandemic/LAST_FILE_PATH.csv", thres_list[i])
	g <- graph.adjacency(my_mat_form, mode="undirected", weighted=T)
	d <- degree(g)
	meandegree1[i] <- mean(d)
	alpha1[i] <- mean(alpna_centrality(g, nodes = V(g), loops = F,exo = 1, weights = NULL, tol = 1e-07, sparse = T))	
	central1[i] <- mean(closeness(g))
	between1[i]<- mean(betweenness(g))
	avepathlen1[i] <- average.path.length(g, directed=F, unconnected=T)
	transitivity1[i] <- mean(transitivity(g))
	powerlaw_alpha1[i] <- fit_power_law(d+1, 10)$alpha
	assort_deg1[i] <- assortativity_degree(g)
	authorScore1[i] <- authority_score(g)$value # same as hub_score(g)$value, Kleinberg’s hub centrality scores
	dyad1[i] <- dyad.census(g)$mut
	diameter1[i] <- diameter(g)
	meanEccentric1[i] <- mean(eccentricity(g))
	meanStrength1[i] <- mean(strength(g))
}
# seasonal clade
meanStrength2 <- meanEccentric2 <- diameter2 <- dyad2 <- authorScore2 <- assort_deg2 <- avepathlen2 <-transitivity2 <- between2 <- alpha2 <- central2 <- meandegree2 <- powerlaw_alpha2 <- c()
for(i in 1:length(thres_list)){
	print(paste0("Now doing threshold ", thres_list[i]))
	my_mat_form <- csv2mat_form("na_seasonal/LAST_FILE_PATH.csv", thres_list[i])
	g <- graph.adjacency(my_mat_form, mode="undirected", weighted=T)
	d <- degree(g)
	meandegree2[i] <- mean(d)
	alpha2[i] <- mean(alpna_centrality(g, nodes = V(g), loops = F, exo = 1, weights = NULL, tol = 1e-07, sparse = T))
	central2[i] <- mean(closeness(g))
	between2[i]<- mean(betweenness(g))
	avepathlen2[i] <- average.path.length(g, directed=F, unconnected=T)
	transitivity2[i] <- mean(transitivity(g))
	powerlaw_alpha2[i] <- fit_power_law(d+1, 10)$alpha
	assort_deg2[i] <- assortativity_degree(g)
	authorScore2[i] <- authority_score(g)$value
	dyad2[i] <- dyad.census(g)$mut
	diameter2[i] <- diameter(g)
	meanEccentric2[i] <- mean(eccentricity(g))
	meanStrength2[i] <- mean(strength(g))
}

save.image("na.RData")


dat_na_pandemic <- data.frame(thres_list, alpha1, meanStrength1, meanEccentric1, diameter1, dyad1, authorScore1, assort_deg1, avepathlen1, transitivity1, between1, central1, meandegree1, powerlaw_alpha1)
dat_na_seasonal <- data.frame(thres_list, alpha2, meanStrength2, meanEccentric2, diameter2, dyad2, authorScore2, assort_deg2, avepathlen2, transitivity2, between2, central2, meandegree2, powerlaw_alpha2)

###########################################	
# plots
##########
panel <- list()

# Alpha centrality 
panel[[1]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Alpha centrality") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_na_pandemic, aes(x=thres_list, y=alpha1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_na_seasonal, aes(x=thres_list, y=alpha2), fill="blue",
      colour="darkblue", size=1)

# meandegree 
panel[[2]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Average degree") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_na_pandemic, aes(x=thres_list, y= meandegree1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_na_seasonal, aes(x=thres_list, y= meandegree2), fill="blue",
      colour="darkblue", size=1)

# meanStrength 
panel[[3]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Average strength") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_na_pandemic, aes(x=thres_list, y= meanStrength1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_na_seasonal, aes(x=thres_list, y= meanStrength2), fill="blue",
      colour="darkblue", size=1)

# meanEccentric 
panel[[4]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Average eccentricity") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_na_pandemic, aes(x=thres_list, y= meanEccentric1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_na_seasonal, aes(x=thres_list, y= meanEccentric2), fill="blue",
      colour="darkblue", size=1)

# diameter 
panel[[5]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Average diameter") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_na_pandemic, aes(x=thres_list, y= diameter1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_na_seasonal, aes(x=thres_list, y= diameter2), fill="blue",
      colour="darkblue", size=1)

# dyad 
panel[[6]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Number of dyads with mutual connections") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_na_pandemic, aes(x=thres_list, y= dyad1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_na_seasonal, aes(x=thres_list, y= dyad2), fill="blue",
      colour="darkblue", size=1)

# authorScore 
panel[[7]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Authority score") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_na_pandemic, aes(x=thres_list, y= authorScore1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_na_seasonal, aes(x=thres_list, y= authorScore2), fill="blue",
      colour="darkblue", size=1)

# assort_deg 
panel[[8]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Assortative degree") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_na_pandemic, aes(x=thres_list, y= assort_deg1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_na_seasonal, aes(x=thres_list, y= assort_deg2), fill="blue",
      colour="darkblue", size=1)

# avepathlen 
panel[[9]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Average path length") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_na_pandemic, aes(x=thres_list, y= avepathlen1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_na_seasonal, aes(x=thres_list, y= avepathlen2), fill="blue",
      colour="darkblue", size=1)

# transitivity2 
panel[[10]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Transitivity") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_na_pandemic, aes(x=thres_list, y= transitivity1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_na_seasonal, aes(x=thres_list, y= transitivity2), fill="blue",
      colour="darkblue", size=1)

# between 
panel[[11]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Average betweenness") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_na_pandemic, aes(x=thres_list, y= between1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_na_seasonal, aes(x=thres_list, y= between2), fill="blue",
      colour="darkblue", size=1)

# central 
panel[[12]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Average centrality") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_na_pandemic, aes(x=thres_list, y= central1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_na_seasonal, aes(x=thres_list, y= central2), fill="blue",
      colour="darkblue", size=1)

# powerlaw_alpha 
panel[[13]] <- ggplot() + theme_bw() + xlab("Posterior probability threshold (strength of coevolution)") + ylab("Power law (alpha)") + scale_x_continuous(limits = c(0, 1)) +
      # for pandemic clade
      geom_smooth(data=dat_na_pandemic, aes(x=thres_list, y= powerlaw_alpha1), fill="red",
      colour="darkred", size=1) +
      # for seasonal clade
      geom_smooth(data=dat_na_seasonal, aes(x=thres_list, y= powerlaw_alpha2), fill="blue",
      colour="darkblue", size=1)


# print all figures
pdf("fig_influenza_ha.pdf" ,width=24, height=16)
do.call("grid.arrange", panel)  
dev.off()


save.image("all.RData")
