library(data.table)

setwd("F:/bioinformatics/research_project")

allMerged <- read.csv("allMerged.csv")

dsDNA.dt <- data.table(subset(allMerged, Type == "dsDNA"))

ssRNA.dt <- data.table(subset(allMerged, Type == "ssRNA+" | Type == "ssRNA-"))


MeanRNA <- aggregate(ssRNA.dt$Hits, list(ssRNA.dt$Virus),mean)
MeanDNA <- aggregate(dsDNA.dt$Hits, list(dsDNA.dt$Virus),mean)
MeanAll <- aggregate(allMerged$Hits, list(allMerged$Virus),mean)
