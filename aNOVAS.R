setwd("/Users/louis/Desktop/project_temp")

coEv <- read.csv("AllCo.csv")

fit <- aov(hitnum ~ seqlen*seqnum, data = coEv)
fitcov <- aov(hitnum ~ seqlen+seqnum, data = coEv)

summary(fit)
summary(fitcov)