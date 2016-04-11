setwd("/Users/louis/Desktop/project_temp")
treshEv <- read.csv("AllCoEvSel.csv")
setwd("/Users/louis/Desktop/plots_stephane")

pdf(paste("all.pdf",sep=""))
treshs = 0.95
treshm = 0.95
treshEv <- subset(treshEv, treshEv$tresholdSPIDER == treshs)
treshEv <- subset(treshEv, treshEv$tresholdMEME == treshm)
GP <- subset(treshEv,treshEv$function. == "GP")
other <- subset(treshEv,treshEv$function. == "Other")

par(mfrow = c(1,2))
boxplot(GP$intersection,x=GP$type,xlab="",main=c("GP spider=",treshs, " meme=", treshm))
text(x= GP$type, labels = c("dsDNA","ssRNA"),pos=1)
boxplot(other$intersection, x= other$type, xlab="", main=c("Other spider=",treshs, " meme=", treshm))
text(x= GP$type, labels = c("dsDNA","ssRNA"),pos=1)


setwd("/Users/louis/Desktop/project_temp")
treshEv <- read.csv("AllCoEvSel.csv")
setwd("/Users/louis/Desktop/plots_stephane")

treshs = 0.95
treshm = 0.50
treshEv <- subset(treshEv, treshEv$tresholdSPIDER == treshs)
treshEv <- subset(treshEv, treshEv$tresholdMEME == treshm)
GP <- subset(treshEv,treshEv$function. == "GP")
other <- subset(treshEv,treshEv$function. == "Other")

par(mfrow = c(1,2))
boxplot(GP$intersection,x=GP$type,xlab="",main=c("GP spider=",treshs, " meme=", treshm))
text(x= GP$type, labels = c("dsDNA","ssRNA"),pos=1)
boxplot(other$intersection, x= other$type, xlab="", main=c("Other spider=",treshs, " meme=", treshm))
text(x= GP$type, labels = c("dsDNA","ssRNA"),pos=1)


setwd("/Users/louis/Desktop/project_temp")
treshEv <- read.csv("AllCoEvSel.csv")
setwd("/Users/louis/Desktop/plots_stephane")

treshs = 0.50
treshm = 0.95
treshEv <- subset(treshEv, treshEv$tresholdSPIDER == treshs)
treshEv <- subset(treshEv, treshEv$tresholdMEME == treshm)
GP <- subset(treshEv,treshEv$function. == "GP")
other <- subset(treshEv,treshEv$function. == "Other")

par(mfrow = c(1,2))
boxplot(GP$intersection,x=GP$type,xlab="",main=c("GP spider=",treshs, " meme=", treshm))
text(x= GP$type, labels = c("dsDNA","ssRNA"),pos=1)
boxplot(other$intersection, x= other$type, xlab="", main=c("Other spider=",treshs, " meme=", treshm))
text(x= GP$type, labels = c("dsDNA","ssRNA"),pos=1)


setwd("/Users/louis/Desktop/project_temp")
treshEv <- read.csv("AllCoEvSel.csv")
setwd("/Users/louis/Desktop/plots_stephane")

treshs = 0.50
treshm = 0.50
treshEv <- subset(treshEv, treshEv$tresholdSPIDER == treshs)
treshEv <- subset(treshEv, treshEv$tresholdMEME == treshm)
GP <- subset(treshEv,treshEv$function. == "GP")
other <- subset(treshEv,treshEv$function. == "Other")

par(mfrow = c(1,2))
boxplot(GP$intersection,x=GP$type,xlab="",main=c("GP spider=",treshs, " meme=", treshm))
text(x= GP$type, labels = c("dsDNA","ssRNA"),pos=1)
boxplot(other$intersection, x= other$type, xlab="", main=c("Other spider=",treshs, " meme=", treshm))
text(x= GP$type, labels = c("dsDNA","ssRNA"),pos=1)


setwd("/Users/louis/Desktop/project_temp")
treshEv <- read.csv("AllCoEvSel.csv")
setwd("/Users/louis/Desktop/plots_stephane")

treshs = 0.70
treshm = 0.95
treshEv <- subset(treshEv, treshEv$tresholdSPIDER == treshs)
treshEv <- subset(treshEv, treshEv$tresholdMEME == treshm)
GP <- subset(treshEv,treshEv$function. == "GP")
other <- subset(treshEv,treshEv$function. == "Other")

par(mfrow = c(1,2))
boxplot(GP$intersection,x=GP$type,xlab="",main=c("GP spider=",treshs, " meme=", treshm))
text(x= GP$type, labels = c("dsDNA","ssRNA"),pos=1)
boxplot(other$intersection, x= other$type, xlab="", main=c("Other spider=",treshs, " meme=", treshm))
text(x= GP$type, labels = c("dsDNA","ssRNA"),pos=1)
dev.off()



