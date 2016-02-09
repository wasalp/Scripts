setwd("D:/bioinformatics/research project")
bob <- list()
data <- read.csv("CountAllGO.csv")
TotalCount <- c(sum(data$HighCount),sum(data$allCount))
for( i in 1:nrow(data)){
  row <- data[i,]
  HT <- c(row$HighCount, row$allCount)
  fish <- fisher.test(matrix(c(HT,TotalCount-HT),ncol = 2, nrow = 2), alternative = "greater")
  bob <- rbind(bob, list(as.character(row$Term),as.double(fish[1])))
}
write.csv(bob, file = "fisherTest.csv")