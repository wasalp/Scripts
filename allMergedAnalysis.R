library(ggplot2)
library(gridExtra)

setwd("F:/bioinformatics/research_project")

allMerged <- read.csv("info.csv")

dsDNA <- subset(allMerged, Type == "dsDNA")

ssRNA <- subset(allMerged, Type == "ssRNA-")

dsDNA <- subset(dsDNA, seqnum > 99)

ssRNA <- subset(ssRNA, seqnum > 99)

ssRNA <- ssRNA[order(ssRNA$Hits),]
dsDNA <- dsDNA[order(dsDNA$Hits),]

a <-ggplot()+ 
  geom_point(aes(y=ssRNA$Hits, x=ssRNA[order(ssRNA$Protein_id),]$Protein_id),fill=2,colour=2,size=1)+
  theme_bw()+labs(title="ssRNA blast hits",x="Genes ordered from least to most hits",y="Blast hits")+ theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    axis.ticks = element_blank(),
                    axis.text.x = element_blank())
a <- a + ylim(0,6000)

b <- ggplot() 
b <- b + geom_point(aes(y=dsDNA$Hits, x=dsDNA[order(dsDNA$Protein_id),]$Protein_id),fill=3,colour=3,size=1)
b <- b + theme_bw()+labs(title="dsDNA blast hits",x="Genes ordered from least to most hits",y="Blast hits")

b <- b + theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_blank()) 
b <- b+ylim(0,6000)


ssRNA1 <- ssRNA[order(ssRNA$seqnum),]
dsDNA1 <- dsDNA[order(dsDNA$seqnum),]
c <-ggplot()+ 
  geom_point(aes(y=ssRNA1$seqnum, x=ssRNA1[order(ssRNA1$Protein_id),]$Protein_id),fill=2,colour=2,size=1)+
  theme_bw()+labs(title="ssRNA retrieved sequences",x="Genes ordered from least to most sequences",y="number of sequences")+ theme(axis.line = element_line(colour = "black"),
                                                                                                            panel.grid.major = element_blank(),
                                                                                                            panel.grid.minor = element_blank(),
                                                                                                            panel.border = element_blank(),
                                                                                                            panel.background = element_blank(),
                                                                                                            axis.ticks = element_blank(),
                                                                                                            axis.text.x = element_blank())
c <- c + ylim(0,4000)

d <- ggplot() 
d <- d + geom_point(aes(y=dsDNA1$seqnum, x=dsDNA1[order(dsDNA1$Protein_id),]$Protein_id),fill=3,colour=3,size=1)
d <- d + theme_bw()+labs(title="dsDNA retrieved sequences",x="Genes ordered from least to most sequences",y="number of sequences")

d <- d + theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_blank()) 
d <- d+ylim(0,4000)

e <- ggplot() 
e <- e + geom_point(aes(y=mean(dsDNA1$seqnum), x=dsDNA1$Type,fill=3,colour=3,size=1))
e <- e + geom_point(aes(y=mean(ssRNA1$seqnum), x=ssRNA1$Type,fill=3,colour=3,size=1))
e <- e + theme_bw()+labs(title="dsDNA retrieved sequences",x="Genes ordered from least to most sequences",y="number of sequences")

e <- e + theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.position="none")

f <- ggplot() 
f <- f + geom_point(aes(y=length(dsDNA1$Protein_id), x=dsDNA1$Type,fill=3,colour=3,size=1))
f <- f + geom_point(aes(y=length(ssRNA1$Protein_id), x=ssRNA1$Type,fill=3,colour=3,size=1))
f <- f + theme_bw()+labs(title="dsDNA retrieved sequences",x="Genes ordered from least to most sequences",y="number of sequences")

f <- f + theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.position="none")
grid.arrange(a,b,c,d,e,f)