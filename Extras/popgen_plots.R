library(ggplot2)
library(plotly)
setwd("/home/wbooker/Dropbox/Documents/Research/ChrysoscelisVersicolor/SumStats/Obs")
data_EC_AE <- read.table("EC_AE_Allstat.txt", header = T)
data_EC_SW <- read.table("EC_SW_Allstat.txt", header = T)
data_EC_MW <- read.table("EC_MW_Allstat.txt", header = T)
data_CMW_MW <- read.table("CMW_MW_Allstat.txt", header = T)
data_CMW_SW <- read.table("CMW_SW_Allstat.txt", header = T)
data_CMW_AE <- read.table("CMW_AE_Allstat.txt", header = T)
data_CSW_SW <- read.table("CSW_SW_Allstat.txt", header = T)
data_CSW_MW <- read.table("CSW_MW_Allstat.txt", header = T)
data_CSW_AE <- read.table("CSW_AE_Allstat.txt", header = T)


##################### Pi #################
piMat <- c(data_EC_AE$piA, data_CMW_MW$piA, data_CSW_SW$piA, data_EC_AE$piB, data_CMW_MW$piB, data_CSW_SW$piB)
piMat <- cbind(piMat, c(rep("EC", nrow(data_EC_AE)), rep("CMW", nrow(data_CMW_MW)), rep("CSW", nrow(data_CSW_SW)), rep("AE", nrow(data_EC_AE)), rep("MW", nrow(data_CMW_MW)), rep("SW", nrow(data_CSW_SW))))


a <- as.data.frame(piMat)

colnames(a) <- c("Pi","Lineage")
a$Lineage <- as.factor(a$Lineage)
a$Pi <- as.numeric(as.character(a$Pi))
Pi <- ggplot(a, aes(x=Lineage, y = Pi)) + geom_violin()  + theme_classic() + #scale_y_continuous(breaks=seq(0.0, 0.6, 0.1), limits=c(0, 0.6)) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.2) + stat_summary(fun.y=mean, geom="point", color = "red", size = 2) + scale_x_discrete(limits=c("EC","CMW","CSW", "AE", "MW", "SW"))

ggplotly(Pi)
ggsave('Pi.png', width = 10, height = 7, dpi = 500)

##################### Dtaj ################
DtajMat <- c(data_EC_AE$DtajA, data_CMW_MW$DtajA, data_CSW_SW$DtajA, data_EC_AE$DtajB, data_CMW_MW$DtajB, data_CSW_SW$DtajB)
DtajMat <- cbind(DtajMat, c(rep("EC", nrow(data_EC_AE)), rep("CMW", nrow(data_CMW_MW)), rep("CSW", nrow(data_CSW_SW)), rep("AE", nrow(data_EC_AE)), rep("MW", nrow(data_CMW_MW)), rep("SW", nrow(data_CSW_SW))))


a <- as.data.frame(DtajMat)

colnames(a) <- c("Dtaj","Lineage")
a$Lineage <- as.factor(a$Lineage)
a$Dtaj <- as.numeric(as.character(a$Dtaj))
Dtaj <- ggplot(a, aes(x=Lineage, y = Dtaj)) + geom_violin()  + theme_classic() + #scale_y_continuous(breaks=seq(0.0, 0.6, 0.1), limits=c(0, 0.6)) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.2) + stat_summary(fun.y=mean, geom="point", color = "red", size = 2) + scale_x_discrete(limits=c("EC","CMW","CSW", "AE", "MW", "SW"))

ggplotly(Dtaj)
ggsave('Dtaj.png', width = 10, height = 7, dpi = 500)

############## Theta
thetaMat <- c(data_EC_AE$thetaA, data_CMW_MW$thetaA, data_CSW_SW$thetaA, data_EC_AE$thetaB, data_CMW_MW$thetaB, data_CSW_SW$thetaB)
thetaMat <- cbind(thetaMat, c(rep("EC", nrow(data_EC_AE)), rep("CMW", nrow(data_CMW_MW)), rep("CSW", nrow(data_CSW_SW)), rep("AE", nrow(data_EC_AE)), rep("MW", nrow(data_CMW_MW)), rep("SW", nrow(data_CSW_SW))))


a <- as.data.frame(thetaMat)

colnames(a) <- c("theta","Lineage")
a$Lineage <- as.factor(a$Lineage)
a$theta <- as.numeric(as.character(a$theta))
theta <- ggplot(a, aes(x=Lineage, y = theta)) + geom_violin()  + theme_classic() + #scale_y_continuous(breaks=seq(0.0, 0.6, 0.1), limits=c(0, 0.6)) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.2) + stat_summary(fun.y=mean, geom="point", color = "red", size = 2) + scale_x_discrete(limits=c("EC","CMW","CSW", "AE", "MW", "SW"))

ggplotly(theta)
ggsave('theta.png', width = 10, height = 7, dpi = 500)

################# D3a
D3aMat <- c(data_EC_AE$D3a, data_CMW_MW$D3a, data_CSW_SW$D3a)
D3aMat <- cbind(D3aMat, c(rep("EC_AE", nrow(data_EC_AE)), rep("CMW_MW", nrow(data_CMW_MW)), rep("CSW_SW", nrow(data_CSW_SW))))


a <- as.data.frame(D3aMat)

colnames(a) <- c("D3a","Lineage")
a$Lineage <- as.factor(a$Lineage)
a$D3a <- as.numeric(as.character(a$D3a))
D3a <- ggplot(a, aes(x=Lineage, y = D3a)) + geom_violin()  + theme_classic() + #scale_y_continuous(breaks=seq(0.0, 0.6, 0.1), limits=c(0, 0.6)) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.2) + stat_summary(fun.y=mean, geom="point", color = "red", size = 2) + scale_x_discrete(limits=c("EC_AE", "CMW_MW", "CSW_SW"))

ggplotly(D3a)
ggsave('D3a.png', width = 10, height = 7, dpi = 500)

################# divAB
divABMat <- c(data_EC_AE$divAB, data_CMW_MW$divAB, data_CSW_SW$divAB)
divABMat <- cbind(divABMat, c(rep("EC_AE", nrow(data_EC_AE)), rep("CMW_MW", nrow(data_CMW_MW)), rep("CSW_SW", nrow(data_CSW_SW))))


a <- as.data.frame(divABMat)

colnames(a) <- c("divAB","Lineage")
a$Lineage <- as.factor(a$Lineage)
a$divAB <- as.numeric(as.character(a$divAB))
divAB <- ggplot(a, aes(x=Lineage, y = divAB)) + geom_violin()  + theme_classic() + #scale_y_continuous(breaks=seq(0.0, 0.6, 0.1), limits=c(0, 0.6)) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.2) + stat_summary(fun.y=mean, geom="point", color = "red", size = 2) + scale_x_discrete(limits=c("EC_AE", "CMW_MW", "CSW_SW"))

ggplotly(divAB)
ggsave('divAB.png', width = 10, height = 7, dpi = 500)

################# netdivAB
netdivABMat <- c(data_EC_AE$netdivAB, data_CMW_MW$netdivAB, data_CSW_SW$netdivAB)
netdivABMat <- cbind(netdivABMat, c(rep("EC_AE", nrow(data_EC_AE)), rep("CMW_MW", nrow(data_CMW_MW)), rep("CSW_SW", nrow(data_CSW_SW))))


a <- as.data.frame(netdivABMat)

colnames(a) <- c("netdivAB","Lineage")
a$Lineage <- as.factor(a$Lineage)
a$netdivAB <- as.numeric(as.character(a$netdivAB))
netdivAB <- ggplot(a, aes(x=Lineage, y = netdivAB)) + geom_violin()  + theme_classic() + #scale_y_continuous(breaks=seq(0.0, 0.6, 0.1), limits=c(0, 0.6)) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.2) + stat_summary(fun.y=mean, geom="point", color = "red", size = 2) + scale_x_discrete(limits=c("EC_AE", "CMW_MW", "CSW_SW"))

ggplotly(netdivAB)
ggsave('netdivAB.png', width = 10, height = 7, dpi = 500)

################# FST
FSTMat <- c(data_EC_AE$FST, data_CMW_MW$FST, data_CSW_SW$FST)
FSTMat <- cbind(FSTMat, c(rep("EC_AE", nrow(data_EC_AE)), rep("CMW_MW", nrow(data_CMW_MW)), rep("CSW_SW", nrow(data_CSW_SW))))


a <- as.data.frame(FSTMat)

colnames(a) <- c("FST","Lineage")
a$Lineage <- as.factor(a$Lineage)
a$FST <- as.numeric(as.character(a$FST))
FST <- ggplot(a, aes(x=Lineage, y = FST)) + geom_violin()  + theme_classic() + #scale_y_continuous(breaks=seq(0.0, 0.6, 0.1), limits=c(0, 0.6)) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.2) + stat_summary(fun.y=mean, geom="point", color = "red", size = 2) + scale_x_discrete(limits=c("EC_AE", "CMW_MW", "CSW_SW"))

ggplotly(FST)
ggsave('FST.png', width = 10, height = 7, dpi = 500)

################# FST AE
FSTMat <- c(data_EC_AE$FST, data_CMW_AE$FST, data_CSW_AE$FST)
FSTMat <- cbind(FSTMat, c(rep("EC_AE", nrow(data_EC_AE)), rep("CMW_AE", nrow(data_CMW_MW)), rep("CSW_AE", nrow(data_CSW_SW))))


a <- as.data.frame(FSTMat)

colnames(a) <- c("FST","Lineage")
a$Lineage <- as.factor(a$Lineage)
a$FST <- as.numeric(as.character(a$FST))
FST <- ggplot(a, aes(x=Lineage, y = FST)) + geom_violin()  + theme_classic() + #scale_y_continuous(breaks=seq(0.0, 0.6, 0.1), limits=c(0, 0.6)) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.2) + stat_summary(fun.y=mean, geom="point", color = "red", size = 2) + scale_x_discrete(limits=c("EC_AE", "CMW_AE", "CSW_AE"))

ggplotly(FST)
ggsave('FST_AE.png', width = 10, height = 7, dpi = 500)

################# FST MW
FSTMat <- c(data_EC_MW$FST, data_CMW_MW$FST, data_CSW_MW$FST)
FSTMat <- cbind(FSTMat, c(rep("EC_MW", nrow(data_EC_MW)), rep("CMW_MW", nrow(data_CMW_MW)), rep("CSW_MW", nrow(data_CSW_SW))))


a <- as.data.frame(FSTMat)

colnames(a) <- c("FST","Lineage")
a$Lineage <- as.factor(a$Lineage)
a$FST <- as.numeric(as.character(a$FST))
FST <- ggplot(a, aes(x=Lineage, y = FST)) + geom_violin()  + theme_classic() + #scale_y_continuous(breaks=seq(0.0, 0.6, 0.1), limits=c(0, 0.6)) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.2) + stat_summary(fun.y=mean, geom="point", color = "red", size = 2) + scale_x_discrete(limits=c("EC_MW", "CMW_MW", "CSW_MW"))

ggplotly(FST)
ggsave('FST_MW.png', width = 10, height = 7, dpi = 500)

################# FST SW
FSTMat <- c(data_EC_SW$FST, data_CSW_SW$FST, data_CSW_SW$FST)
FSTMat <- cbind(FSTMat, c(rep("EC_SW", nrow(data_EC_SW)), rep("CMW_SW", nrow(data_CSW_SW)), rep("CSW_SW", nrow(data_CSW_SW))))


a <- as.data.frame(FSTMat)

colnames(a) <- c("FST","Lineage")
a$Lineage <- as.factor(a$Lineage)
a$FST <- as.numeric(as.character(a$FST))
FST <- ggplot(a, aes(x=Lineage, y = FST)) + geom_violin()  + theme_classic() + #scale_y_continuous(breaks=seq(0.0, 0.6, 0.1), limits=c(0, 0.6)) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.2) + stat_summary(fun.y=mean, geom="point", color = "red", size = 2) + scale_x_discrete(limits=c("EC_SW", "CMW_SW", "CSW_SW"))

ggplotly(FST)
ggsave('FST_SW.png', width = 10, height = 7, dpi = 500)

################# divAB AE
divABMat <- c(data_EC_AE$divAB, data_CMW_AE$divAB, data_CSW_AE$divAB)
divABMat <- cbind(divABMat, c(rep("EC_AE", nrow(data_EC_AE)), rep("CMW_AE", nrow(data_CMW_MW)), rep("CSW_AE", nrow(data_CSW_SW))))


a <- as.data.frame(divABMat)

colnames(a) <- c("divAB","Lineage")
a$Lineage <- as.factor(a$Lineage)
a$divAB <- as.numeric(as.character(a$divAB))
divAB <- ggplot(a, aes(x=Lineage, y = divAB)) + geom_violin()  + theme_classic() + #scale_y_continuous(breaks=seq(0.0, 0.6, 0.1), limits=c(0, 0.6)) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.2) + stat_summary(fun.y=mean, geom="point", color = "red", size = 2) + scale_x_discrete(limits=c("EC_AE", "CMW_AE", "CSW_AE"))

ggplotly(divAB)
ggsave('divAB_AE.png', width = 10, height = 7, dpi = 500)

################# divAB MW
divABMat <- c(data_EC_MW$divAB, data_CMW_MW$divAB, data_CSW_MW$divAB)
divABMat <- cbind(divABMat, c(rep("EC_MW", nrow(data_EC_MW)), rep("CMW_MW", nrow(data_CMW_MW)), rep("CSW_MW", nrow(data_CSW_SW))))


a <- as.data.frame(divABMat)

colnames(a) <- c("divAB","Lineage")
a$Lineage <- as.factor(a$Lineage)
a$divAB <- as.numeric(as.character(a$divAB))
divAB <- ggplot(a, aes(x=Lineage, y = divAB)) + geom_violin()  + theme_classic() + #scale_y_continuous(breaks=seq(0.0, 0.6, 0.1), limits=c(0, 0.6)) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.2) + stat_summary(fun.y=mean, geom="point", color = "red", size = 2) + scale_x_discrete(limits=c("EC_MW", "CMW_MW", "CSW_MW"))

ggplotly(divAB)
ggsave('divAB_MW.png', width = 10, height = 7, dpi = 500)

################# divAB SW
divABMat <- c(data_EC_SW$divAB, data_CSW_SW$divAB, data_CSW_SW$divAB)
divABMat <- cbind(divABMat, c(rep("EC_SW", nrow(data_EC_SW)), rep("CMW_SW", nrow(data_CSW_SW)), rep("CSW_SW", nrow(data_CSW_SW))))


a <- as.data.frame(divABMat)

colnames(a) <- c("divAB","Lineage")
a$Lineage <- as.factor(a$Lineage)
a$divAB <- as.numeric(as.character(a$divAB))
divAB <- ggplot(a, aes(x=Lineage, y = divAB)) + geom_violin()  + theme_classic() + #scale_y_continuous(breaks=seq(0.0, 0.6, 0.1), limits=c(0, 0.6)) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.2) + stat_summary(fun.y=mean, geom="point", color = "red", size = 2) + scale_x_discrete(limits=c("EC_SW", "CMW_SW", "CSW_SW"))

ggplotly(divAB)
ggsave('divAB_SW.png', width = 10, height = 7, dpi = 500)

################# netdivAB AE
netdivABMat <- c(data_EC_AE$netdivAB, data_CMW_AE$netdivAB, data_CSW_AE$netdivAB)
netdivABMat <- cbind(netdivABMat, c(rep("EC_AE", nrow(data_EC_AE)), rep("CMW_AE", nrow(data_CMW_MW)), rep("CSW_AE", nrow(data_CSW_SW))))


a <- as.data.frame(netdivABMat)

colnames(a) <- c("netdivAB","Lineage")
a$Lineage <- as.factor(a$Lineage)
a$netdivAB <- as.numeric(as.character(a$netdivAB))
netdivAB <- ggplot(a, aes(x=Lineage, y = netdivAB)) + geom_violin()  + theme_classic() + #scale_y_continuous(breaks=seq(0.0, 0.6, 0.1), limits=c(0, 0.6)) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.2) + stat_summary(fun.y=mean, geom="point", color = "red", size = 2) + scale_x_discrete(limits=c("EC_AE", "CMW_AE", "CSW_AE"))

ggplotly(netdivAB)
ggsave('netdivAB_AE.png', width = 10, height = 7, dpi = 500)

################# netdivAB MW
netdivABMat <- c(data_EC_MW$netdivAB, data_CMW_MW$netdivAB, data_CSW_MW$netdivAB)
netdivABMat <- cbind(netdivABMat, c(rep("EC_MW", nrow(data_EC_MW)), rep("CMW_MW", nrow(data_CMW_MW)), rep("CSW_MW", nrow(data_CSW_SW))))


a <- as.data.frame(netdivABMat)

colnames(a) <- c("netdivAB","Lineage")
a$Lineage <- as.factor(a$Lineage)
a$netdivAB <- as.numeric(as.character(a$netdivAB))
netdivAB <- ggplot(a, aes(x=Lineage, y = netdivAB)) + geom_violin()  + theme_classic() + #scale_y_continuous(breaks=seq(0.0, 0.6, 0.1), limits=c(0, 0.6)) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.2) + stat_summary(fun.y=mean, geom="point", color = "red", size = 2) + scale_x_discrete(limits=c("EC_MW", "CMW_MW", "CSW_MW"))

ggplotly(netdivAB)
ggsave('netdivAB_MW.png', width = 10, height = 7, dpi = 500)

################# netdivAB SW
netdivABMat <- c(data_EC_SW$netdivAB, data_CSW_SW$netdivAB, data_CSW_SW$netdivAB)
netdivABMat <- cbind(netdivABMat, c(rep("EC_SW", nrow(data_EC_SW)), rep("CMW_SW", nrow(data_CSW_SW)), rep("CSW_SW", nrow(data_CSW_SW))))


a <- as.data.frame(netdivABMat)

colnames(a) <- c("netdivAB","Lineage")
a$Lineage <- as.factor(a$Lineage)
a$netdivAB <- as.numeric(as.character(a$netdivAB))
netdivAB <- ggplot(a, aes(x=Lineage, y = netdivAB)) + geom_violin()  + theme_classic() + #scale_y_continuous(breaks=seq(0.0, 0.6, 0.1), limits=c(0, 0.6)) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.2) + stat_summary(fun.y=mean, geom="point", color = "red", size = 2) + scale_x_discrete(limits=c("EC_SW", "CMW_SW", "CSW_SW"))

ggplotly(netdivAB)
ggsave('netdivAB_SW.png', width = 10, height = 7, dpi = 500)

################# D3a AE
D3aMat <- c(data_EC_AE$D3a, data_CMW_AE$D3a, data_CSW_AE$D3a)
D3aMat <- cbind(D3aMat, c(rep("EC_AE", nrow(data_EC_AE)), rep("CMW_AE", nrow(data_CMW_MW)), rep("CSW_AE", nrow(data_CSW_SW))))


a <- as.data.frame(D3aMat)

colnames(a) <- c("D3a","Lineage")
a$Lineage <- as.factor(a$Lineage)
a$D3a <- as.numeric(as.character(a$D3a))
D3a <- ggplot(a, aes(x=Lineage, y = D3a)) + geom_violin()  + theme_classic() + #scale_y_continuous(breaks=seq(0.0, 0.6, 0.1), limits=c(0, 0.6)) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.2) + stat_summary(fun.y=mean, geom="point", color = "red", size = 2) + scale_x_discrete(limits=c("EC_AE", "CMW_AE", "CSW_AE"))

ggplotly(D3a)
ggsave('D3a_AE.png', width = 10, height = 7, dpi = 500)

################# D3a MW
D3aMat <- c(data_EC_MW$D3a, data_CMW_MW$D3a, data_CSW_MW$D3a)
D3aMat <- cbind(D3aMat, c(rep("EC_MW", nrow(data_EC_MW)), rep("CMW_MW", nrow(data_CMW_MW)), rep("CSW_MW", nrow(data_CSW_SW))))


a <- as.data.frame(D3aMat)

colnames(a) <- c("D3a","Lineage")
a$Lineage <- as.factor(a$Lineage)
a$D3a <- as.numeric(as.character(a$D3a))
D3a <- ggplot(a, aes(x=Lineage, y = D3a)) + geom_violin()  + theme_classic() + #scale_y_continuous(breaks=seq(0.0, 0.6, 0.1), limits=c(0, 0.6)) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.2) + stat_summary(fun.y=mean, geom="point", color = "red", size = 2) + scale_x_discrete(limits=c("EC_MW", "CMW_MW", "CSW_MW"))

ggplotly(D3a)
ggsave('D3a_MW.png', width = 10, height = 7, dpi = 500)

################# D3a SW
D3aMat <- c(data_EC_SW$D3a, data_CSW_SW$D3a, data_CSW_SW$D3a)
D3aMat <- cbind(D3aMat, c(rep("EC_SW", nrow(data_EC_SW)), rep("CMW_SW", nrow(data_CSW_SW)), rep("CSW_SW", nrow(data_CSW_SW))))


a <- as.data.frame(D3aMat)

colnames(a) <- c("D3a","Lineage")
a$Lineage <- as.factor(a$Lineage)
a$D3a <- as.numeric(as.character(a$D3a))
D3a <- ggplot(a, aes(x=Lineage, y = D3a)) + geom_violin()  + theme_classic() + #scale_y_continuous(breaks=seq(0.0, 0.6, 0.1), limits=c(0, 0.6)) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.2) + stat_summary(fun.y=mean, geom="point", color = "red", size = 2) + scale_x_discrete(limits=c("EC_SW", "CMW_SW", "CSW_SW"))

ggplotly(D3a)
ggsave('D3a_SW.png', width = 10, height = 7, dpi = 500)
