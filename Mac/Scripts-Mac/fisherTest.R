setwd("F:/bioinformatics/research_project")
library(stringr)
bob <- list()
GOTerm <- read.csv("CountAllGO.csv")
BP <- subset(GOTerm,str_detect(GOTerm$Term,":biological_process"))
MF <- subset(GOTerm,str_detect(GOTerm$Term,":molecular_function"))
CC <- subset(GOTerm,str_detect(GOTerm$Term,":cellular_component"))

TotalCountBP <- c(sum(BP$HighCount),sum(BP$allCount))
for( i in 1:nrow(BP)){
  row <- BP[i,]
  HT <- c(row$HighCount, row$allCount)
  fish <- fisher.test(matrix(c(HT,TotalCountBP-HT),ncol = 2, nrow = 2), alternative = "greater")
  bob <- rbind(bob, list(as.character(row$Term),as.double(fish[1])))
}

write.csv(bob, file = "fisherBP.csv")
bob <- list()
TotalCountMF <- c(sum(MF$HighCount),sum(MF$allCount))
for( i in 1:nrow(MF)){
  row <- MF[i,]
  HT <- c(row$HighCount, row$allCount)
  fish <- fisher.test(matrix(c(HT,TotalCountMF-HT),ncol = 2, nrow = 2), alternative = "greater")
  bob <- rbind(bob, list(as.character(row$Term),as.double(fish[1])))
}

write.csv(bob, file = "fisherMF.csv")
bob <- list()
TotalCountCC <- c(sum(CC$HighCount),sum(CC$allCount))
for( i in 1:nrow(CC)){
  row <- CC[i,]
  HT <- c(row$HighCount, row$allCount)
  fish <- fisher.test(matrix(c(HT,TotalCountCC-HT),ncol = 2, nrow = 2), alternative = "greater")
  bob <- rbind(bob, list(as.character(row$Term),as.double(fish[1])))
}

write.csv(bob, file = "fisherCC.csv")