library(abcrf)
setwd("/home/wbooker/Dropbox/Documents/Research/ChrysoscelisVersicolor/PolyploidABC_versicolor/PolyploidABC/EC_AE/")
ss=c(2:46) ## These are the rows to pull out of the sim data for ABC analysis/comparison
statObs=read.table("EC_AE_ObsData.txt", h = T)[ss]
nRep = 40
allo_disomic_noMig=allo_heterosomic_noMig=allo_tetrasomic_noMig=auto_disomic_noMig=auto_heterosomic_noMig=auto_tetrasomic_noMig=polyP_disomic_noMig=polyP_heterosomic_noMig=polyP_tetrasomic_noMig=allo_disomic_migAB=allo_heterosomic_migAB=allo_tetrasomic_migAB=auto_disomic_migAB=auto_heterosomic_migAB=auto_tetrasomic_migAB=polyP_disomic_migAB=polyP_heterosomic_migAB=polyP_tetrasomic_migAB=NULL
for(i in 1:nRep){
  #allo_disomic_noMig = rbind(allo_disomic_noMig, read.table(paste("allo_disomic_noMig/allo_disomic_noMig_",i,"/ABCstat.txt", sep = ""), skip=0, h = T)[,ss])
  #allo_heterosomic_noMig = rbind(allo_heterosomic_noMig, read.table(paste("allo_heterosomic_noMig/allo_heterosomic_noMig_",i,"/ABCstat.txt", sep = ""), skip=0, h = T)[,ss])
  #allo_tetrasomic_noMig = rbind(allo_tetrasomic_noMig, read.table(paste("allo_tetrasomic_noMig/allo_tetrasomic_noMig_",i,"/ABCstat.txt", sep = ""), skip=0, h = T)[,ss])
  #auto_disomic_noMig = rbind(auto_disomic_noMig, read.table(paste("auto_disomic_noMig/auto_disomic_noMig_",i,"/ABCstat.txt", sep = ""), skip=0, h = T)[,ss])
  #auto_heterosomic_noMig = rbind(auto_heterosomic_noMig, read.table(paste("auto_heterosomic_noMig/auto_heterosomic_noMig_",i,"/ABCstat.txt", sep = ""), skip=0, h = T)[,ss])
  #auto_tetrasomic_noMig = rbind(auto_tetrasomic_noMig, read.table(paste("auto_tetrasomic_noMig/auto_tetrasomic_noMig_",i,"/ABCstat.txt", sep = ""), skip=0, h = T)[,ss])
  #polyP_disomic_noMig = rbind(polyP_disomic_noMig, read.table(paste("polyP_disomic_noMig/polyP_disomic_noMig_",i,"/ABCstat.txt", sep = ""), skip=0, h = T)[,ss])
  #olyP_heterosomic_noMig = rbind(polyP_heterosomic_noMig, read.table(paste("polyP_heterosomic_noMig/polyP_heterosomic_noMig_",i,"/ABCstat.txt", sep = ""), skip=0, h = T)[,ss])
  #polyP_tetrasomic_noMig = rbind(polyP_tetrasomic_noMig, read.table(paste("polyP_tetrasomic_noMig/polyP_tetrasomic_noMig_",i,"/ABCstat.txt", sep = ""), skip=0, h = T)[,ss])
  allo_disomic_migAB = rbind(allo_disomic_migAB, read.table(paste("allo_disomic_migAB/allo_disomic_migAB_",i,"/ABCstat.txt", sep = ""), skip=0, h = T)[,ss])
  allo_heterosomic_migAB = rbind(allo_heterosomic_migAB, read.table(paste("allo_heterosomic_migAB/allo_heterosomic_migAB_",i,"/ABCstat.txt", sep = ""), skip=0, h = T)[,ss])
  allo_tetrasomic_migAB = rbind(allo_tetrasomic_migAB, read.table(paste("allo_tetrasomic_migAB/allo_tetrasomic_migAB_",i,"/ABCstat.txt", sep = ""), skip=0, h = T)[,ss])
  auto_disomic_migAB = rbind(auto_disomic_migAB, read.table(paste("auto_disomic_migAB/auto_disomic_migAB_",i,"/ABCstat.txt", sep = ""), skip=0, h = T)[,ss])
  auto_heterosomic_migAB = rbind(auto_heterosomic_migAB, read.table(paste("auto_heterosomic_migAB/auto_heterosomic_migAB_",i,"/ABCstat.txt", sep = ""), skip=0, h = T)[,ss])
  auto_tetrasomic_migAB = rbind(auto_tetrasomic_migAB, read.table(paste("auto_tetrasomic_migAB/auto_tetrasomic_migAB_",i,"/ABCstat.txt", sep = ""), skip=0, h = T)[,ss])
  polyP_disomic_migAB = rbind(polyP_disomic_migAB, read.table(paste("polyP_disomic_migAB/polyP_disomic_migAB_",i,"/ABCstat.txt", sep = ""), skip=0, h = T)[,ss])
  polyP_heterosomic_migAB = rbind(polyP_heterosomic_migAB, read.table(paste("polyP_heterosomic_migAB/polyP_heterosomic_migAB_",i,"/ABCstat.txt", sep = ""), skip=0, h = T)[,ss])
  polyP_tetrasomic_migAB = rbind(polyP_tetrasomic_migAB, read.table(paste("polyP_tetrasomic_migAB/polyP_tetrasomic_migAB_",i,"/ABCstat.txt", sep = ""), skip=0, h = T)[,ss])
  
}

for(i in 1:ncol(allo_disomic_migAB)){
  #allo_disomic_noMig[which(allo_disomic_noMig[,i]=="NaN"),i]=mean(allo_disomic_noMig[,i], na.rm=T)
  #allo_heterosomic_noMig[which(allo_heterosomic_noMig[,i]=="NaN"),i]=mean(allo_heterosomic_noMig[,i], na.rm=T)
  #allo_tetrasomic_noMig[which(allo_tetrasomic_noMig[,i]=="NaN"),i]=mean(allo_tetrasomic_noMig[,i], na.rm=T)
  #auto_disomic_noMig[which(auto_disomic_noMig[,i]=="NaN"),i]=mean(auto_disomic_noMig[,i], na.rm=T)
  #auto_heterosomic_noMig[which(auto_heterosomic_noMig[,i]=="NaN"),i]=mean(auto_heterosomic_noMig[,i], na.rm=T)
  #auto_tetrasomic_noMig[which(auto_tetrasomic_noMig[,i]=="NaN"),i]=mean(auto_tetrasomic_noMig[,i], na.rm=T)
  #polyP_disomic_noMig[which(polyP_disomic_noMig[,i]=="NaN"),i]=mean(polyP_disomic_noMig[,i], na.rm=T)
  #polyP_heterosomic_noMig[which(polyP_heterosomic_noMig[,i]=="NaN"),i]=mean(polyP_heterosomic_noMig[,i], na.rm=T)
  #polyP_tetrasomic_noMig[which(polyP_tetrasomic_noMig[,i]=="NaN"),i]=mean(polyP_tetrasomic_noMig[,i], na.rm=T)
  allo_disomic_migAB[which(allo_disomic_migAB[,i]=="NaN"),i]=mean(allo_disomic_migAB[,i], na.rm=T)
  allo_heterosomic_migAB[which(allo_heterosomic_migAB[,i]=="NaN"),i]=mean(allo_heterosomic_migAB[,i], na.rm=T)
  allo_tetrasomic_migAB[which(allo_tetrasomic_migAB[,i]=="NaN"),i]=mean(allo_tetrasomic_migAB[,i], na.rm=T)
  auto_disomic_migAB[which(auto_disomic_migAB[,i]=="NaN"),i]=mean(auto_disomic_migAB[,i], na.rm=T)
  auto_heterosomic_migAB[which(auto_heterosomic_migAB[,i]=="NaN"),i]=mean(auto_heterosomic_migAB[,i], na.rm=T)
  auto_tetrasomic_migAB[which(auto_tetrasomic_migAB[,i]=="NaN"),i]=mean(auto_tetrasomic_migAB[,i], na.rm=T)
  polyP_disomic_migAB[which(polyP_disomic_migAB[,i]=="NaN"),i]=mean(polyP_disomic_migAB[,i], na.rm=T)
  polyP_heterosomic_migAB[which(polyP_heterosomic_migAB[,i]=="NaN"),i]=mean(polyP_heterosomic_migAB[,i], na.rm=T)
  polyP_tetrasomic_migAB[which(polyP_tetrasomic_migAB[,i]=="NaN"),i]=mean(polyP_tetrasomic_migAB[,i], na.rm=T)
}


#################Reduce metrics here

nScenarios=9
nRef= nrow(allo_disomic_migAB)
nTrees_in_forest=750
nSS=ncol(statObs)
modIndex <- NULL
for(i in 1:nScenarios){
  modIndex <- c(modIndex, rep(i,nRef))
}
modIndex <- as.factor(modIndex)
sumstat=rbind(allo_disomic_noMig,allo_heterosomic_noMig,allo_tetrasomic_noMig,auto_disomic_noMig,auto_heterosomic_noMig,auto_tetrasomic_noMig, polyP_disomic_noMig, polyP_heterosomic_noMig, polyP_tetrasomic_noMig,allo_disomic_migAB,allo_heterosomic_migAB,allo_tetrasomic_migAB,auto_disomic_migAB,auto_heterosomic_migAB,auto_tetrasomic_migAB, polyP_disomic_migAB, polyP_heterosomic_migAB, polyP_tetrasomic_migAB)
data1 <- data.frame(modIndex, sumstat)
model.rf1 <- abcrf(modIndex~., data = data1, ntree=nTrees_in_forest, paral = T, ncores =12, lda = T)

###################Plot model and lda axes
source("../plotabcrf_v2.R")
plot.abcrf_all(model.rf1, data1, statObs, n.var=30, pdf = T)

########### Asess model fit and posterior probability
predRF_n <- predict(model.rf1, statObs, data1, ntree=nTrees_in_forest, paral = T, ncores = 12)

##############Prior Info
allo_disomic_noMig_prior=allo_heterosomic_noMig_prior=allo_tetrasomic_noMig_prior=auto_disomic_noMig_prior=auto_heterosomic_noMig_prior=auto_tetrasomic_noMig_prior=polyP_disomic_noMig_prior=polyP_heterosomic_noMig_prior=polyP_tetrasomic_noMig_prior=allo_disomic_migAB_prior=allo_heterosomic_migAB_prior=allo_tetrasomic_migAB_prior=auto_disomic_migAB_prior=auto_heterosomic_migAB_prior=auto_tetrasomic_migAB_prior=polyP_disomic_migAB_prior=polyP_heterosomic_migAB_prior=polyP_tetrasomic_migAB_prior=NULL
priorCol = c(1:6)
for(i in 1:nRep){
  allo_disomic_noMig_prior = rbind(allo_disomic_noMig_prior, read.table(paste("allo_disomic_noMig/allo_disomic_noMig_",i,"/output.txt", sep = ""), skip=0, h = T)[,priorCol])
  allo_heterosomic_noMig_prior = rbind(allo_heterosomic_noMig_prior, read.table(paste("allo_heterosomic_noMig/allo_heterosomic_noMig_",i,"/output.txt", sep = ""), skip=0, h = T)[,priorCol])
  allo_tetrasomic_noMig_prior = rbind(allo_tetrasomic_noMig_prior, read.table(paste("allo_tetrasomic_noMig/allo_tetrasomic_noMig_",i,"/output.txt", sep = ""), skip=0, h = T)[,priorCol])
  auto_disomic_noMig_prior = rbind(auto_disomic_noMig_prior, read.table(paste("auto_disomic_noMig/auto_disomic_noMig_",i,"/output.txt", sep = ""), skip=0, h = T)[,priorCol])
  auto_heterosomic_noMig_prior = rbind(auto_heterosomic_noMig_prior, read.table(paste("auto_heterosomic_noMig/auto_heterosomic_noMig_",i,"/output.txt", sep = ""), skip=0, h = T)[,priorCol])
  auto_tetrasomic_noMig_prior = rbind(auto_tetrasomic_noMig_prior, read.table(paste("auto_tetrasomic_noMig/auto_tetrasomic_noMig_",i,"/output.txt", sep = ""), skip=0, h = T)[,priorCol])
  polyP_disomic_noMig_prior = rbind(polyP_disomic_noMig_prior, read.table(paste("polyP_disomic_noMig/polyP_disomic_noMig_",i,"/output.txt", sep = ""), skip=0, h = T)[,priorCol])
  polyP_heterosomic_noMig_prior = rbind(polyP_heterosomic_noMig_prior, read.table(paste("polyP_heterosomic_noMig/polyP_heterosomic_noMig_",i,"/output.txt", sep = ""), skip=0, h = T)[,priorCol])
  polyP_tetrasomic_noMig_prior = rbind(polyP_tetrasomic_noMig_prior, read.table(paste("polyP_tetrasomic_noMig/polyP_tetrasomic_noMig_",i,"/output.txt", sep = ""), skip=0, h = T)[,priorCol])
  
}
selectedModel <- "15"
priors=rbind(allo_disomic_noMig_prior,allo_heterosomic_noMig_prior,allo_tetrasomic_noMig_prior,auto_disomic_noMig_prior,auto_heterosomic_noMig_prior,auto_tetrasomic_noMig_prior, polyP_disomic_noMig_prior, polyP_heterosomic_noMig_prior, polyP_tetrasomic_noMig_prior,allo_disomic_migAB_prior,allo_heterosomic_migAB_prior,allo_tetrasomic_migAB_prior,auto_disomic_migAB_prior,auto_heterosomic_migAB_prior,auto_tetrasomic_migAB_prior, polyP_disomic_migAB_prior, polyP_heterosomic_migAB_prior, polyP_tetrasomic_migAB_prior)
sumstatModel <-sumstat[modIndex==selectedModel,]

priors = cbind(modIndex,priors)

Tsplit <- priors$Tsplit[modIndex==selectedModel]
data2 <- cbind(Tsplit, sumstatModel)
model.rf.Tsplit <- regAbcrf(Tsplit~., data2, ntree=500, paral = T, ncores = 12)
predTsplit <- predict(model.rf.Tsplit, statObs, data2, paral = T, ncores = 12)

TWGD <- priors$Twgd[modIndex==selectedModel]
data2 <- cbind(TWGD, sumstatModel)
model.rf.TWGD <- regAbcrf(TWGD~., data2, ntree=500, paral = T, ncores = 12)
predTWGD<- predict(model.rf.TWGD, statObs, data2, paral = T, ncores = 12)

NAsize <- priors$NA.[modIndex==selectedModel]
data2 <- cbind(NAsize, sumstatModel)
model.rf.NAsize <- regAbcrf(NAsize~., data2, ntree=500, paral = T, ncores = 12)
predNAsize<- predict(model.rf.NAsize, statObs, data2, paral = T, ncores = 12)

N2Asize <- priors$N2A[modIndex==selectedModel]
data2 <- cbind(N2Asize, sumstatModel)
model.rf.N2Asize <- regAbcrf(N2Asize~., data2, ntree=500, paral = T, ncores = 12)
predN2Asize<- predict(model.rf.N2Asize, statObs, data2, paral = T, ncores = 12)

N1size <- priors$N1[modIndex==selectedModel]
data2 <- cbind(N1size, sumstatModel)
model.rf.N1size <- regAbcrf(N1size~., data2, ntree=500, paral = T, ncores = 12)
predN1size<- predict(model.rf.N1size, statObs, data2, paral = T, ncores = 12)

N2size <- priors$N2[modIndex==selectedModel]
data2 <- cbind(N2size, sumstatModel)
model.rf.N2size <- regAbcrf(N2size~., data2, ntree=500, paral = T, ncores = 12)
predN2size<- predict(model.rf.N2size, statObs, data2, paral = T, ncores = 12)


plot(x = model.rf.TWGD, n.var = 22)
