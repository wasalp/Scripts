setwd("/Users/louis/Desktop/project_temp")
library(plyr)
treshEv <- read.csv("coEv.csv")
palette("default")

plot(0,0,xlim = c(0.5,11),ylim = c(0,22) ,
     main="Number of hits", xlab="Gene",
     ylab="Number of sites found",xaxt="n")

with(treshEv,points(gene[order(gene)],hitnum,col=func))
legend('topright', legend=levels(treshEv$func),col=1:6,cex=0.8,pch=1)
with(treshEv,text(hitnum~gene[order(gene)],labels = gene, cex=0.6, pos=3, col="black"))