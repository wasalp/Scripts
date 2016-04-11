setwd("/Users/louis/Desktop/project_temp")

coEv <- read.csv("AllCoEvSel.csv")
palette(rainbow(22))


plot(NULL,xlim=c(0,0.95),ylim=c(0,100),ylab = "hit mean", xlab = "Spider treshold")
spidertreshlevel <- unique(as.factor(coEv$tresholdSPIDER))
print(spidertreshlevel)
memelevel <- unique(as.factor(coEv$tresholdMEME))
print(levels(memelevel))
colvar <- 1
for(b in memelevel){
  
  memesub <- subset(coEv, coEv$tresholdMEME == b)
  vectormcgee <- data.frame(a=c(NA),b=c(NA))
  for(i in spidertreshlevel){
    
    vectormcgee <- rbind(vectormcgee,c(as.numeric(i),mean(memesub$intersection[memesub$tresholdSPIDER == i])))
    #print(as.numeric(i))
    #print(mean(intersection[tresholdSPIDER == i]))
    
    
  }
  lines(vectormcgee$a,vectormcgee$b,col=colvar)
  colvar = colvar +1
}
legend('topright', legend=levels(as.factor(memelevel)),col=1:colvar,cex=0.8,pch=1)



