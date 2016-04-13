library(igraph)
library(robust)
###########################################
# generates random data
npairs <- 2500
seqlen <- 20
site1 <- floor(runif(npairs, 1, seqlen))
site2 <- floor(runif(npairs, 1, seqlen))
probs <- rbeta(npairs, .5, .2)
rand_pairs <- cbind(site1, site2, probs)
write.csv(rand_pairs, "rand_pairs1.csv", quote=F, row.names=F)

npairs <- 2500
seqlen <- 20
site1 <- floor(runif(npairs, 1, seqlen))
site2 <- floor(runif(npairs, 1, seqlen))
probs <- rbeta(npairs, .2, .5)
rand_pairs <- cbind(site1, site2, probs)
write.csv(rand_pairs, "rand_pairs2.csv", quote=F, row.names=F)
###########################################
# maybe not the fastest way to threshold matrix b/c I/O overhead, but it works
csv2mat_form <- function(filename, threshold = .95){
	# reads network from file and format data for downstream analyses
	net_dat <- read.csv(filename)
	n_pairs <- length(net_dat$site1)
	
	list_of_uniq_sites <- unique(union(net_dat$site1, net_dat$site2))
	list_len <- length(list_of_uniq_sites)
	
	mat_form <- matrix(0, list_len, list_len) 
	colnames(mat_form) <- rownames(mat_form) <- list_of_uniq_sites
	
	# fills up matrix -- non-sym form
	for(i in 1:n_pairs){
		ii <- which(colnames(mat_form) == net_dat$site1[i])
		jj <- which(colnames(mat_form) == net_dat$site2[i])
		if(net_dat$probs[i] > threshold){
			mat_form[ii,jj] <- net_dat$probs[i]
		}else{
			mat_form[ii,jj] <- 0
		}
	}
	return(mat_form)
}

my_mat_form <- csv2mat_form("Clade_Rose_Spidermonkey_BGM_report.csv", .5)
my_mat_form

###########################################
# draw network
g <- graph.adjacency(mat_form, mode="undirected", weighted=TRUE)
log10_weights <- -2*log10(E(g)$weight)
E(g)$width <- log10_weights + min(log10_weights) + 1
igraph_options(plot.layout=layout_as_tree, label.cex=.2)
plot(g, vertex.color="orange", vertex.size=10) # see ?igraph.plotting for options


###########################################
# posterior probability threshold analysis
thres_list <- seq(.05,.95,.05)
# clade 1
average1 <- transitivity1 <- between1 <- alpha1 <- central1 <- meandegree1 <- powerlaw_slope1 <- powerlaw_alpha1 <- c()
for(i in 1:length(thres_list)){
	my_mat_form <- csv2mat_form("Clade_Rose_Spidermonkey_BGM_report.csv", thres_list[i])
	g <- graph.adjacency(my_mat_form, mode="undirected", weighted=TRUE)
	d <- degree(g)
	meandegree1[i] <- mean(d)
	alpha1[i] <- mean(alpha_centrality(g, nodes = V(g), loops = FALSE,exo = 1, weights = NULL, tol = 1e-07, sparse = TRUE))	
	central1[i] <- mean(closeness(g))
	between1[i]<- mean(betweenness(g))
	average1[i] <- average.path.length(g, directed=TRUE, unconnected=TRUE)
	transitivity1[i] <- mean(transitivity(g))
	hd <- hist(d, plot=F)
	#lm_pl <- lmrob(log10(hd$counts+.1) ~ log10(seq(1,length(hd$counts))) )
	#ummary(lm_pl)
    # slope of fitted power law
	#powerlaw_slope1[i] <- summary(lm_pl)$coefficients[2,1]
	powerlaw_alpha1[i] <- fit_power_law(d+1, 10)$alpha
}
# clade 2
 <- average2 <-transitivity2 <- between2 <- alpha2 <- central2 <- meandegree2 <- powerlaw_slope2 <- powerlaw_alpha2 <- c()
for(i in 1:length(thres_list)){
	my_mat_form <- csv2mat_form("Clade_Orange_Spidermonkey_BGM_report.csv", thres_list[i])
	g <- graph.adjacency(my_mat_form, mode="undirected", weighted=TRUE)
	d <- degree(g)
	meandegree2[i] <- mean(d)
	alpha2[i] <- mean(alpha_centrality(g, nodes = V(g), loops = FALSE, exo = 1, weights = NULL, tol = 1e-07, sparse = TRUE))
	central2[i] <- mean(closeness(g))
	between2[i]<- mean(betweenness(g))
	average2[i] <- average.path.length(g, directed=TRUE, unconnected=TRUE)
	transitivity2[i] <- mean(transitivity(g))
		hd <- hist(d, plot=F)
	#lm_pl <- lmrob( log10(hd$counts+.1) ~ log10(seq(1,length(hd$counts))) )
    #summary(lm_pl)
	# slope of fitted power law
	#powerlaw_slope2[i] <- summary(lm_pl)$coefficients[2,1]
	powerlaw_alpha2[i] <- fit_power_law(d+1, 10)$alpha
}


average.path.length(graph, directed=TRUE, unconnected=TRUE)

###########################################	
###GGplot#
##########
###dataframe_alpha
##1
zika_data.frame_Brazil <- data.frame(alpha1, thres_list)
c <- ggplot(zika_data.frame_Brazil, aes(thres_list, alpha1))
c + stat_smooth() + geom_smooth(fill="blue", colour="blue", size=1)
##2
zika_data.frame_senegal <- data.frame(alpha2, thres_list)
c <- ggplot(zika_data.frame_senegal, aes(thres_list, alpha2))
c + stat_smooth() + geom_smooth(fill="red", colour="red", size=1)

###method 1
p <- ggplot() +
      # blue plot
        geom_smooth(data=zika_data.frame_Brazil, aes(x=thres_list, y=alpha1), fill="blue",
        colour="darkblue", size=1) +
      # red plot
        geom_smooth(data=zika_data.frame_senegal, aes(x=thres_list, y=alpha2), fill="red",
        colour="red", size=1)
###########################################
###dataframe_mean degree
##1
df_degree_b <- data.frame(meandegree1, thres_list)
c <- ggplot(df_degree_b, aes(thres_list, meandegree1))
c + stat_smooth() + geom_smooth(fill="blue", colour="blue", size=1)
##2
df_degree_s <- data.frame(meandegree2, thres_list)
c <- ggplot(df_degree_s, aes(thres_list, meandegree2))
c + stat_smooth() + geom_smooth(fill="red", colour="red", size=1)

###method 1
p <- ggplot() +
      # blue plot
        geom_smooth(data=df_degree_b, aes(x=thres_list, y=meandegree1), fill="blue",
        colour="darkblue", size=1) +
      # red plot
        geom_smooth(data=df_degree_s, aes(x=thres_list, y=meandegree2), fill="red",
        colour="red", size=1)
###########################################
###dataframe_central_closeness
##1
df_central_b <- data.frame(central1, thres_list)
c <- ggplot(df_central_b, aes(thres_list, central1))
c + stat_smooth() + geom_smooth(fill="blue", colour="blue", size=1)
##2
df_central_s <- data.frame(central2, thres_list)
c <- ggplot(df_central_s, aes(thres_list, central2))
c + stat_smooth() + geom_smooth(fill="red", colour="red", size=1)

###method 1
p <- ggplot() +
      # blue plot
        geom_smooth(data=df_central_b, aes(x=thres_list, y=central1), fill="blue",
        colour="darkblue", size=1) +
      # red plot
        geom_smooth(data=df_central_s, aes(x=thres_list, y=central2), fill="red",
        colour="red", size=1)
###########################################
###dataframe_betwenesss
##1
df_between_b <- data.frame(between1, thres_list)
c <- ggplot(df_between_b, aes(thres_list, between1))
c + stat_smooth() + geom_smooth(fill="blue", colour="blue", size=1)
##2
df_between_s <- data.frame(between2, thres_list)
c <- ggplot(df_between_s, aes(thres_list, between2))
c + stat_smooth() + geom_smooth(fill="red", colour="red", size=1)

###method 1
p <- ggplot() +
      # blue plot
        geom_smooth(data=df_between_b, aes(x=thres_list, y=between1), fill="blue",
        colour="darkblue", size=1) +
      # red plot
        geom_smooth(data=df_between_s, aes(x=thres_list, y=between2), fill="red",
        colour="red", size=1)
		
###########################################
###dataframe_transitivity
##1
df_transitivity_b <- data.frame(transitivity1, thres_list)
c <- ggplot(df_transitivity_b, aes(thres_list, transitivity1))
c + stat_smooth() + geom_smooth(fill="blue", colour="blue", size=1)
##2
df_transitivity_s <- data.frame(transitivity2, thres_list)
c <- ggplot(df_transitivity_s, aes(thres_list, transitivity2))
c + stat_smooth() + geom_smooth(fill="red", colour="red", size=1)

###method 1
p <- ggplot() +
      # blue plot
        geom_smooth(data=df_transitivity_b, aes(x=thres_list, y=transitivity1), fill="blue",
        colour="darkblue", size=1) +
      # red plot
        geom_smooth(data=df_transitivity_s, aes(x=thres_list, y=transitivity2), fill="red",
        colour="red", size=1)
###########################################
###dataframe_average path lenght
##1
df_average_b <- data.frame(average1, thres_list)
c <- ggplot(df_average_b, aes(thres_list, average1))
c + stat_smooth() + geom_smooth(fill="blue", colour="blue", size=1)
##2
df_average_s <- data.frame(average2, thres_list)
c <- ggplot(df_average_s, aes(thres_list, average2))
c + stat_smooth() + geom_smooth(fill="red", colour="red", size=1)

###method 1
p <- ggplot() +
      # blue plot
        geom_smooth(data=df_average_b, aes(x=thres_list, y=average1), fill="blue",
        colour="darkblue", size=1) +
      # red plot
        geom_smooth(data=df_average_s, aes(x=thres_list, y=average2), fill="red",
        colour="red", size=1)


###########################################
####BETTER RESOLUTION
pdf("meandegreelog.pdf", width=6, height=4)
fudge <- .001
plot(log10(meandegree1+fudge) ~ thres_list, pch=20, xlab="Posterior probability threshold", ylab="Log 10 Mean node degree", col="red", ylim=c(-3,3))
points(log10(meandegree2+fudge) ~ thres_list, col="blue")
legend("topright", c("Clade Brazil","Clade Senegal"), pch=20, col=c("red","blue"), bty="y", cex=.7)
dev.off()

pdf("meandegree.pdf", width=6, height=4)
plot(meandegree1 ~ thres_list, pch=20, xlab="Posterior probability threshold", ylab="Mean node degree", col="red", ylim=c(min(meandegree1, meandegree2), max(meandegree1, meandegree2)))
points(meandegree2 ~ thres_list, pch=20, col="blue")
legend("topright", c("Clade Brazil","Clade Senegal"), pch=20, col=c("red","blue"), bty="y", cex=.7)
dev.off()

pdf("powerlaw_alpha.pdf", width=6, height=4)
plot(log10(powerlaw_alpha1) ~ thres_list, pch=20, xlab="Posterior probability threshold", ylab="Mean node degree", col="red", ylim=c(min(log10(powerlaw_alpha1), log10(powerlaw_alpha2)), max(log10(powerlaw_alpha1), log10(powerlaw_alpha2))))
points(log10(powerlaw_alpha2) ~ thres_list, pch=20, col="blue")
legend("topright", c("Clade Brazil","Clade Senegal"), pch=20, col=c("red","blue"), bty="y", cex=.7)
dev.off()

###Central closeness
pdf("Central_Closeness.pdf", width=6, height=4)
plot(central1 ~ thres_list, pch=20, xlab="Posterior probability threshold", ylab="Mean node centrality", col="red", ylim=c(min(central1, central2), max(central1, central2)))
points(central2 ~ thres_list, pch=20, col="blue")
legend("topright", c("Clade Brazil","Clade Senegal"), pch=20, col=c("red","blue"), bty="y", cex=.7)
dev.off()

###Better resolution
pdf("Central_Closeness_Log.pdf", width=6, height=4)
plot(log10(central1) ~ thres_list, pch=20, xlab="Posterior probability threshold", ylab="Mean node centrality", col="red", ylim=c(min(log10(central1), log10(central2)), max(log10(central1), log10(central2))))
points(log10(central2) ~ thres_list, pch=20, col="blue")
legend("topright", c("Clade Brazil","Clade Senegal"), pch=20, col=c("red","blue"), bty="y", cex=.7)
dev.off()
##Alpha_centrality
pdf("Alpha_centrality.pdf", width=6, height=4)
plot(alpha1 ~ thres_list , pch=20, xlab="Posterior probability threshold", ylab="Mean Alpha Centrality", col="red", ylim=c(min(alpha1, alpha2), max(alpha1, alpha2)))
points(alpha2 ~ thres_list, pch=20, col="blue")
legend("bottomright", c("Clade Brazil","Clade Senegal"), pch=20, col=c("red","blue"), bty="y", cex=.7)
dev.off()

##Betweenness
pdf("mean Betweeness.pdf", width=6, height=4)
plot(between1 ~ thres_list, pch=20, xlab="Posterior probability threshold", ylab="Mean node degree", col="red", ylim=c(min(between1, between2), max(between1, between2)))
points(between2 ~ thres_list, pch=20, col="blue")
legend("topright", c("Clade Brazil","Clade Senegal"), pch=20, col=c("red","blue"), bty="y", cex=.7)
dev.off()

##transitivity
pdf("transitivity.pdf", width=6, height=4)
plot(transitivity1 ~ thres_list, pch=20, xlab="Posterior probability threshold", ylab=" transitivity", col="red", ylim=c(-2,2))
points(transitivity2 ~ thres_list, pch=20, col="blue")
legend("topright", c("Clade Brazil","Clade Senegal"), pch=20, col=c("red","blue"), bty="y", cex=.7)
dev.off()

##average
pdf("average.pdf", width=6, height=4)
plot(average1 ~ thres_list, pch=20, xlab="Posterior probability threshold", ylab="average", col="red", ylim=c(0,6))
points(average2 ~ thres_list, pch=20, col="blue")
legend("topright", c("Clade Brazil","Clade Senegal"), pch=20, col=c("red","blue"), bty="y", cex=.7)
dev.off()


###########################################


plot(powerlaw_slope1 ~ thres_list, pch=20, xlab="Posterior probability threshold", ylab="MM slope of power law")
plot(powerlaw_alpha1 ~ thres_list, pch=20, xlab="Posterior probability threshold", ylab=expression(paste(alpha," of power law")))

###########################################
######Reseau

	pdf()
	i <- 2
	thres_list[i]
	my_mat_form <- csv2mat_form("Clade_Orange_Spidermonkey_BGM_report.csv", thres_list[i])
	g <- graph.adjacency(my_mat_form, mode="undirected", weighted=TRUE)
	plot(g, layout=layout_with_kk, vertex.color="green")
	i <- 4
	thres_list[i]
	my_mat_form <- csv2mat_form("Clade_Orange_Spidermonkey_BGM_report.csv", thres_list[i])
	g <- graph.adjacency(my_mat_form, mode="undirected", weighted=TRUE)
	plot(g, layout=layout_with_kk, vertex.color="orange")
	i <- 8
	thres_list[i]
	my_mat_form <- csv2mat_form("Clade_Orange_Spidermonkey_BGM_report.csv", thres_list[i])
	g <- graph.adjacency(my_mat_form, mode="undirected", weighted=TRUE)
	plot(g, layout=layout_with_kk, vertex.color="orange")
	dev.off()
	
	
	
	i <- 8
	thres_list[i]
	my_mat_form <- csv2mat_form("Clade_Rose_Spidermonkey_BGM_report.csv", thres_list[i])
	g <- graph.adjacency(my_mat_form, mode="undirected", weighted=TRUE)
	plot(g, layout=layout_with_kk, vertex.color="pink")
	
####Find clique
i=5
my_mat_form <- csv2mat_form("Clade_Orange_Spidermonkey_BGM_report.csv", thres_list[i])
	g <- graph.adjacency(my_mat_form, mode="undirected", weighted=TRUE)
mc <- maximal.cliques(g)
length(mc)
sapply(mc, length)
col <- rep("blue", length(V(g)))
col[mc[[length(mc)]]] <- "red"
lay.kk <- layout.kamada.kawai(g, kkconst=5)
plot.igraph(g, layout=lay.kk, vertex.color=col,vertex.label=NA)

lay.kk <- layout.kamada.kawai(g, kkconst=5)
plot.igraph(g, layout=lay.kk, vertex.label=NA, vertex.size=10, vertex.color="black")


###########################################	
##### LTT analysis
	b <- read.nexus(file = "combined_out_medians.trees", tree.names = NULL)
ltt.plot(b)
###########################################





Ref:
http://stackoverflow.com/questions/21192002/how-to-combine-2-plots-ggplot-into-one-plot
http://download.springer.com/static/pdf/678/art%253A10.1007%252Fs00227-016-2826-x.pdf?originUrl=http%3A%2F%2Flink.springer.com%2Farticle%2F10.1007%2Fs00227-016-2826-x&token2=exp=1459198945~acl=%2Fstatic%2Fpdf%2F678%2Fart%25253A10.1007%25252Fs00227-016-2826-x.pdf%3ForiginUrl%3Dhttp%253A%252F%252Flink.springer.com%252Farticle%252F10.1007%252Fs00227-016-2826-x*~hmac=371f409ae54d8b1763a12e9636e49e873ae71ced9e0e65404115b10f8b3ba47f
