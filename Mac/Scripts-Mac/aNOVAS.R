setwd("/Users/louis/Desktop/project_temp")

coEv <- read.csv("AllCoEvSel.csv")

fitcov <- aov(hitnum ~ seqlen+seqnum, data = coEv)

summary(fitcov)

summary(fit)
summary(fitcov)