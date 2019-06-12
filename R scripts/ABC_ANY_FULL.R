library(nnet)
dataset <- "EC_AE"
setwd(paste(c("/gpfs/research/scratch/wwb15/ABC_WGD/50Loci/Exponential/",dataset,"/"), collapse = ""))
source("../cv4abc.R")
ss=c(4:50) ## These are the rows to pull out of the sim data for ABC analysis/comparison, removed biallelic sites column since we removed some snps converting to ms format for observed data
EpReplicates <- 4500
nModels <- 9 ########## Number of models analyzed for each individual analysis (eg noMig)
nModelsFinal <- 15 ######### Number of models analyzed in the final analysis of 3 per individual. Will require some tweaking in that section as well, however 
nnetReplicates = 20
nFolders = 40
target=read.table(paste(c(dataset,"_ObsData.txt"),collapse = ""), h = T)[ss]
allo_disomic_noMig=allo_heterosomic_noMig=allo_tetrasomic_noMig=auto_disomic_noMig=auto_heterosomic_noMig=auto_tetrasomic_noMig=polyP_disomic_noMig=polyP_heterosomic_noMig=polyP_tetrasomic_noMig=allo_disomic_migAB_1W=allo_heterosomic_migAB_1W=allo_tetrasomic_migAB_1W=auto_disomic_migAB_1W=auto_heterosomic_migAB_1W=auto_tetrasomic_migAB_1W=polyP_disomic_migAB_1W=polyP_heterosomic_migAB_1W=polyP_tetrasomic_migAB_1W=allo_disomic_migA_1W=allo_heterosomic_migA_1W=allo_tetrasomic_migA_1W=auto_disomic_migA_1W=auto_heterosomic_migA_1W=auto_tetrasomic_migA_1W=polyP_disomic_migA_1W=polyP_heterosomic_migA_1W=polyP_tetrasomic_migA_1W=allo_disomic_migAB=allo_heterosomic_migAB=allo_tetrasomic_migAB=auto_disomic_migAB=auto_heterosomic_migAB=auto_tetrasomic_migAB=polyP_disomic_migAB=polyP_heterosomic_migAB=polyP_tetrasomic_migAB=allo_disomic_migA=allo_heterosomic_migA=allo_tetrasomic_migA=auto_disomic_migA=auto_heterosomic_migA=auto_tetrasomic_migA=polyP_disomic_migA=polyP_heterosomic_migA=polyP_tetrasomic_migA=NULL
for(i in 1:nFolders){
  allo_disomic_noMig = rbind(allo_disomic_noMig, read.table(paste("allo_disomic_noMig/allo_disomic_noMig_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  allo_heterosomic_noMig = rbind(allo_heterosomic_noMig, read.table(paste("allo_heterosomic_noMig/allo_heterosomic_noMig_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  allo_tetrasomic_noMig = rbind(allo_tetrasomic_noMig, read.table(paste("allo_tetrasomic_noMig/allo_tetrasomic_noMig_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  auto_disomic_noMig = rbind(auto_disomic_noMig, read.table(paste("auto_disomic_noMig/auto_disomic_noMig_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  auto_heterosomic_noMig = rbind(auto_heterosomic_noMig, read.table(paste("auto_heterosomic_noMig/auto_heterosomic_noMig_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  auto_tetrasomic_noMig = rbind(auto_tetrasomic_noMig, read.table(paste("auto_tetrasomic_noMig/auto_tetrasomic_noMig_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  polyP_disomic_noMig = rbind(polyP_disomic_noMig, read.table(paste("polyP_disomic_noMig/polyP_disomic_noMig_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  polyP_heterosomic_noMig = rbind(polyP_heterosomic_noMig, read.table(paste("polyP_heterosomic_noMig/polyP_heterosomic_noMig_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  polyP_tetrasomic_noMig = rbind(polyP_tetrasomic_noMig, read.table(paste("polyP_tetrasomic_noMig/polyP_tetrasomic_noMig_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
}

for(i in 1:ncol(allo_disomic_noMig)){
  allo_disomic_noMig[which(allo_disomic_noMig[,i]=="NaN"),i]=mean(allo_disomic_noMig[,i], na.rm=T)
  allo_heterosomic_noMig[which(allo_heterosomic_noMig[,i]=="NaN"),i]=mean(allo_heterosomic_noMig[,i], na.rm=T)
  allo_tetrasomic_noMig[which(allo_tetrasomic_noMig[,i]=="NaN"),i]=mean(allo_tetrasomic_noMig[,i], na.rm=T)
  auto_disomic_noMig[which(auto_disomic_noMig[,i]=="NaN"),i]=mean(auto_disomic_noMig[,i], na.rm=T)
  auto_heterosomic_noMig[which(auto_heterosomic_noMig[,i]=="NaN"),i]=mean(auto_heterosomic_noMig[,i], na.rm=T)
  auto_tetrasomic_noMig[which(auto_tetrasomic_noMig[,i]=="NaN"),i]=mean(auto_tetrasomic_noMig[,i], na.rm=T)
  polyP_disomic_noMig[which(polyP_disomic_noMig[,i]=="NaN"),i]=mean(polyP_disomic_noMig[,i], na.rm=T)
  polyP_heterosomic_noMig[which(polyP_heterosomic_noMig[,i]=="NaN"),i]=mean(polyP_heterosomic_noMig[,i], na.rm=T)
  polyP_tetrasomic_noMig[which(polyP_tetrasomic_noMig[,i]=="NaN"),i]=mean(polyP_tetrasomic_noMig[,i], na.rm=T)
}

x=c(rep(1:nModels, each=nrow(allo_disomic_noMig)))
obs = NULL
for(i in 1:nnetReplicates){
  obs <- rbind(obs,target)
}

res=model_selection_abc_nnet(target=obs, x=x, sumstat=rbind(allo_disomic_noMig,allo_heterosomic_noMig,allo_tetrasomic_noMig,auto_disomic_noMig,auto_heterosomic_noMig,auto_tetrasomic_noMig, polyP_disomic_noMig, polyP_heterosomic_noMig, polyP_tetrasomic_noMig
                                                            
), tol=EpReplicates/(nModels*nrow(allo_disomic_noMig))#250/(6*nrow(allo_disomic_noMig))
, noweight=F, rejmethod=F, nb.nnet=35, size.nnet=10, output=paste(c(dataset,"_NoMig_0.05.txt"),collapse = ""))
#################### migAB_1W
for(i in 1:nFolders){
  allo_disomic_migAB_1W = rbind(allo_disomic_migAB_1W, read.table(paste("allo_disomic_migAB_1W/allo_disomic_migAB_1W_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  allo_heterosomic_migAB_1W = rbind(allo_heterosomic_migAB_1W, read.table(paste("allo_heterosomic_migAB_1W/allo_heterosomic_migAB_1W_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  allo_tetrasomic_migAB_1W = rbind(allo_tetrasomic_migAB_1W, read.table(paste("allo_tetrasomic_migAB_1W/allo_tetrasomic_migAB_1W_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  auto_disomic_migAB_1W = rbind(auto_disomic_migAB_1W, read.table(paste("auto_disomic_migAB_1W/auto_disomic_migAB_1W_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  auto_heterosomic_migAB_1W = rbind(auto_heterosomic_migAB_1W, read.table(paste("auto_heterosomic_migAB_1W/auto_heterosomic_migAB_1W_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  auto_tetrasomic_migAB_1W = rbind(auto_tetrasomic_migAB_1W, read.table(paste("auto_tetrasomic_migAB_1W/auto_tetrasomic_migAB_1W_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  polyP_disomic_migAB_1W = rbind(polyP_disomic_migAB_1W, read.table(paste("polyP_disomic_migAB_1W/polyP_disomic_migAB_1W_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  polyP_heterosomic_migAB_1W = rbind(polyP_heterosomic_migAB_1W, read.table(paste("polyP_heterosomic_migAB_1W/polyP_heterosomic_migAB_1W_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  polyP_tetrasomic_migAB_1W = rbind(polyP_tetrasomic_migAB_1W, read.table(paste("polyP_tetrasomic_migAB_1W/polyP_tetrasomic_migAB_1W_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
}#

for(i in 1:ncol(allo_disomic_migAB_1W)){
  allo_disomic_migAB_1W[which(allo_disomic_migAB_1W[,i]=="NaN"),i]=mean(allo_disomic_migAB_1W[,i], na.rm=T)
  allo_heterosomic_migAB_1W[which(allo_heterosomic_migAB_1W[,i]=="NaN"),i]=mean(allo_heterosomic_migAB_1W[,i], na.rm=T)
  allo_tetrasomic_migAB_1W[which(allo_tetrasomic_migAB_1W[,i]=="NaN"),i]=mean(allo_tetrasomic_migAB_1W[,i], na.rm=T)
  auto_disomic_migAB_1W[which(auto_disomic_migAB_1W[,i]=="NaN"),i]=mean(auto_disomic_migAB_1W[,i], na.rm=T)
  auto_heterosomic_migAB_1W[which(auto_heterosomic_migAB_1W[,i]=="NaN"),i]=mean(auto_heterosomic_migAB_1W[,i], na.rm=T)
  auto_tetrasomic_migAB_1W[which(auto_tetrasomic_migAB_1W[,i]=="NaN"),i]=mean(auto_tetrasomic_migAB_1W[,i], na.rm=T)
  polyP_disomic_migAB_1W[which(polyP_disomic_migAB_1W[,i]=="NaN"),i]=mean(polyP_disomic_migAB_1W[,i], na.rm=T)
  polyP_heterosomic_migAB_1W[which(polyP_heterosomic_migAB_1W[,i]=="NaN"),i]=mean(polyP_heterosomic_migAB_1W[,i], na.rm=T)
  polyP_tetrasomic_migAB_1W[which(polyP_tetrasomic_migAB_1W[,i]=="NaN"),i]=mean(polyP_tetrasomic_migAB_1W[,i], na.rm=T)
}

x=c(rep(1:nModels, each=nrow(allo_disomic_migAB_1W)))
obs = NULL
for(i in 1:nnetReplicates){
  obs <- rbind(obs,target)
}

res=model_selection_abc_nnet(target=obs, x=x, sumstat=rbind(allo_disomic_migAB_1W,allo_heterosomic_migAB_1W,allo_tetrasomic_migAB_1W,
                                                            auto_disomic_migAB_1W,auto_heterosomic_migAB_1W,auto_tetrasomic_migAB_1W,
                                                            polyP_disomic_migAB_1W,polyP_heterosomic_migAB_1W,polyP_tetrasomic_migAB_1W
                                                            
), tol=EpReplicates/(nModels*nrow(allo_disomic_migAB_1W))#250/(6*nrow(allo_disomic_noMig))
, noweight=F, rejmethod=F, nb.nnet=35, size.nnet=10, output=paste(c(dataset,"_migAB_1W_0.05.txt"),collapse = ""))



##################  migA_1W
for(i in 1:nFolders){
  allo_disomic_migA_1W = rbind(allo_disomic_migA_1W, read.table(paste("allo_disomic_migA_1W/allo_disomic_migA_1W_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  allo_heterosomic_migA_1W = rbind(allo_heterosomic_migA_1W, read.table(paste("allo_heterosomic_migA_1W/allo_heterosomic_migA_1W_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  allo_tetrasomic_migA_1W = rbind(allo_tetrasomic_migA_1W, read.table(paste("allo_tetrasomic_migA_1W/allo_tetrasomic_migA_1W_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  auto_disomic_migA_1W = rbind(auto_disomic_migA_1W, read.table(paste("auto_disomic_migA_1W/auto_disomic_migA_1W_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  auto_heterosomic_migA_1W = rbind(auto_heterosomic_migA_1W, read.table(paste("auto_heterosomic_migA_1W/auto_heterosomic_migA_1W_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  auto_tetrasomic_migA_1W = rbind(auto_tetrasomic_migA_1W, read.table(paste("auto_tetrasomic_migA_1W/auto_tetrasomic_migA_1W_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  polyP_disomic_migA_1W = rbind(polyP_disomic_migA_1W, read.table(paste("polyP_disomic_migA_1W/polyP_disomic_migA_1W_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  polyP_heterosomic_migA_1W = rbind(polyP_heterosomic_migA_1W, read.table(paste("polyP_heterosomic_migA_1W/polyP_heterosomic_migA_1W_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  polyP_tetrasomic_migA_1W = rbind(polyP_tetrasomic_migA_1W, read.table(paste("polyP_tetrasomic_migA_1W/polyP_tetrasomic_migA_1W_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
}

for(i in 1:ncol(allo_disomic_migA_1W)){
  allo_disomic_migA_1W[which(allo_disomic_migA_1W[,i]=="NaN"),i]=mean(allo_disomic_migA_1W[,i], na.rm=T)
  allo_heterosomic_migA_1W[which(allo_heterosomic_migA_1W[,i]=="NaN"),i]=mean(allo_heterosomic_migA_1W[,i], na.rm=T)
  allo_tetrasomic_migA_1W[which(allo_tetrasomic_migA_1W[,i]=="NaN"),i]=mean(allo_tetrasomic_migA_1W[,i], na.rm=T)
  auto_disomic_migA_1W[which(auto_disomic_migA_1W[,i]=="NaN"),i]=mean(auto_disomic_migA_1W[,i], na.rm=T)
  auto_heterosomic_migA_1W[which(auto_heterosomic_migA_1W[,i]=="NaN"),i]=mean(auto_heterosomic_migA_1W[,i], na.rm=T)
  auto_tetrasomic_migA_1W[which(auto_tetrasomic_migA_1W[,i]=="NaN"),i]=mean(auto_tetrasomic_migA_1W[,i], na.rm=T)
  polyP_disomic_migA_1W[which(polyP_disomic_migA_1W[,i]=="NaN"),i]=mean(polyP_disomic_migA_1W[,i], na.rm=T)
  polyP_heterosomic_migA_1W[which(polyP_heterosomic_migA_1W[,i]=="NaN"),i]=mean(polyP_heterosomic_migA_1W[,i], na.rm=T)
  polyP_tetrasomic_migA_1W[which(polyP_tetrasomic_migA_1W[,i]=="NaN"),i]=mean(polyP_tetrasomic_migA_1W[,i], na.rm=T)
}

x=c(rep(1:nModels, each=nrow(allo_disomic_migA_1W)))
obs = NULL
for(i in 1:nnetReplicates){
  obs <- rbind(obs,target)
}

res=model_selection_abc_nnet(target=obs, x=x, sumstat=rbind(allo_disomic_migA_1W,allo_heterosomic_migA_1W,allo_tetrasomic_migA_1W,
                                                            auto_disomic_migA_1W,auto_heterosomic_migA_1W,auto_tetrasomic_migA_1W,
                                                            polyP_disomic_migA_1W,polyP_heterosomic_migA_1W,polyP_tetrasomic_migA_1W
                                                            
), tol=EpReplicates/(nModels*nrow(allo_disomic_migA_1W))#250/(6*nrow(allo_disomic_noMig))
, noweight=F, rejmethod=F, nb.nnet=35, size.nnet=10, output=paste(c(dataset,"_migA_1W_0.05.txt"),collapse = ""))



#################### migAB
for(i in 1:nFolders){
  allo_disomic_migAB = rbind(allo_disomic_migAB, read.table(paste("allo_disomic_migAB/allo_disomic_migAB_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  allo_heterosomic_migAB = rbind(allo_heterosomic_migAB, read.table(paste("allo_heterosomic_migAB/allo_heterosomic_migAB_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  allo_tetrasomic_migAB = rbind(allo_tetrasomic_migAB, read.table(paste("allo_tetrasomic_migAB/allo_tetrasomic_migAB_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  auto_disomic_migAB = rbind(auto_disomic_migAB, read.table(paste("auto_disomic_migAB/auto_disomic_migAB_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  auto_heterosomic_migAB = rbind(auto_heterosomic_migAB, read.table(paste("auto_heterosomic_migAB/auto_heterosomic_migAB_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  auto_tetrasomic_migAB = rbind(auto_tetrasomic_migAB, read.table(paste("auto_tetrasomic_migAB/auto_tetrasomic_migAB_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  polyP_disomic_migAB = rbind(polyP_disomic_migAB, read.table(paste("polyP_disomic_migAB/polyP_disomic_migAB_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  polyP_heterosomic_migAB = rbind(polyP_heterosomic_migAB, read.table(paste("polyP_heterosomic_migAB/polyP_heterosomic_migAB_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  polyP_tetrasomic_migAB = rbind(polyP_tetrasomic_migAB, read.table(paste("polyP_tetrasomic_migAB/polyP_tetrasomic_migAB_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
}#

for(i in 1:ncol(allo_disomic_migAB)){
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

x=c(rep(1:nModels, each=nrow(allo_disomic_migAB)))
obs = NULL
for(i in 1:nnetReplicates){
  obs <- rbind(obs,target)
}

res=model_selection_abc_nnet(target=obs, x=x, sumstat=rbind(allo_disomic_migAB,allo_heterosomic_migAB,allo_tetrasomic_migAB,
                                                            auto_disomic_migAB,auto_heterosomic_migAB,auto_tetrasomic_migAB,
                                                            polyP_disomic_migAB,polyP_heterosomic_migAB,polyP_tetrasomic_migAB
                                                            
), tol=EpReplicates/(nModels*nrow(allo_disomic_migAB))#250/(6*nrow(allo_disomic_noMig))
, noweight=F, rejmethod=F, nb.nnet=35, size.nnet=10, output=paste(c(dataset,"_migAB_0.05.txt"),collapse = ""))



##################  migA
for(i in 1:nFolders){
  allo_disomic_migA = rbind(allo_disomic_migA, read.table(paste("allo_disomic_migA/allo_disomic_migA_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  allo_heterosomic_migA = rbind(allo_heterosomic_migA, read.table(paste("allo_heterosomic_migA/allo_heterosomic_migA_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  allo_tetrasomic_migA = rbind(allo_tetrasomic_migA, read.table(paste("allo_tetrasomic_migA/allo_tetrasomic_migA_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  auto_disomic_migA = rbind(auto_disomic_migA, read.table(paste("auto_disomic_migA/auto_disomic_migA_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  auto_heterosomic_migA = rbind(auto_heterosomic_migA, read.table(paste("auto_heterosomic_migA/auto_heterosomic_migA_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  auto_tetrasomic_migA = rbind(auto_tetrasomic_migA, read.table(paste("auto_tetrasomic_migA/auto_tetrasomic_migA_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  polyP_disomic_migA = rbind(polyP_disomic_migA, read.table(paste("polyP_disomic_migA/polyP_disomic_migA_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  polyP_heterosomic_migA = rbind(polyP_heterosomic_migA, read.table(paste("polyP_heterosomic_migA/polyP_heterosomic_migA_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  polyP_tetrasomic_migA = rbind(polyP_tetrasomic_migA, read.table(paste("polyP_tetrasomic_migA/polyP_tetrasomic_migA_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
}

for(i in 1:ncol(allo_disomic_migA)){
  allo_disomic_migA[which(allo_disomic_migA[,i]=="NaN"),i]=mean(allo_disomic_migA[,i], na.rm=T)
  allo_heterosomic_migA[which(allo_heterosomic_migA[,i]=="NaN"),i]=mean(allo_heterosomic_migA[,i], na.rm=T)
  allo_tetrasomic_migA[which(allo_tetrasomic_migA[,i]=="NaN"),i]=mean(allo_tetrasomic_migA[,i], na.rm=T)
  auto_disomic_migA[which(auto_disomic_migA[,i]=="NaN"),i]=mean(auto_disomic_migA[,i], na.rm=T)
  auto_heterosomic_migA[which(auto_heterosomic_migA[,i]=="NaN"),i]=mean(auto_heterosomic_migA[,i], na.rm=T)
  auto_tetrasomic_migA[which(auto_tetrasomic_migA[,i]=="NaN"),i]=mean(auto_tetrasomic_migA[,i], na.rm=T)
  polyP_disomic_migA[which(polyP_disomic_migA[,i]=="NaN"),i]=mean(polyP_disomic_migA[,i], na.rm=T)
  polyP_heterosomic_migA[which(polyP_heterosomic_migA[,i]=="NaN"),i]=mean(polyP_heterosomic_migA[,i], na.rm=T)
  polyP_tetrasomic_migA[which(polyP_tetrasomic_migA[,i]=="NaN"),i]=mean(polyP_tetrasomic_migA[,i], na.rm=T)
}

x=c(rep(1:nModels, each=nrow(allo_disomic_migA)))
obs = NULL
for(i in 1:nnetReplicates){
  obs <- rbind(obs,target)
}

res=model_selection_abc_nnet(target=obs, x=x, sumstat=rbind(allo_disomic_migA,allo_heterosomic_migA,allo_tetrasomic_migA,
                                                            auto_disomic_migA,auto_heterosomic_migA,auto_tetrasomic_migA,
                                                            polyP_disomic_migA,polyP_heterosomic_migA,polyP_tetrasomic_migA
                                                            
), tol=EpReplicates/(nModels*nrow(allo_disomic_migA))#250/(6*nrow(allo_disomic_noMig))
, noweight=F, rejmethod=F, nb.nnet=35, size.nnet=10, output=paste(c(dataset,"_migA_0.05.txt"),collapse = ""))



########### Test of top 3 of each model, currently set up for 15 models so will need manual tweaking if that number is different 
model1 = model2 = model3 = model4 = model5 = model6 = model7 = model8 = model9 = model10 = model11 = model12 = model13 = model14 = model15 = NULL
modelOrder <- c("allo_disomic","allo_heterosomic","allo_tetrasomic",
                "auto_disomic","auto_heterosomic","auto_tetrasomic",
                "polyP_disomic","polyP_heterosomic","polyP_tetrasomic")
datN <- order(colMeans(read.table(paste(c(dataset,"_NoMig_0.05.txt"), collapse = ""))), decreasing = T)
datAB1w <- order(colMeans(read.table(paste(c(dataset,"_migAB_1W_0.05.txt"), collapse = ""))), decreasing = T)
datA1w <- order(colMeans(read.table(paste(c(dataset,"_migA_1W_0.05.txt"), collapse = ""))), decreasing = T)
datAB <- order(colMeans(read.table(paste(c(dataset,"_migAB_0.05.txt"), collapse = ""))), decreasing = T)
datA <- order(colMeans(read.table(paste(c(dataset,"_migA_0.05.txt"), collapse = ""))), decreasing = T)
models <- c(paste(c(modelOrder[datN[1]], "_noMig"), collapse = ""),paste(c(modelOrder[datN[2]], "_noMig"), collapse = ""),paste(c(modelOrder[datN[3]], "_noMig"), collapse = ""),
            paste(c(modelOrder[datAB1w[1]], "_migAB_1W"), collapse = ""),paste(c(modelOrder[datAB1w[2]], "_migAB_1W"), collapse = ""),paste(c(modelOrder[datAB1w[3]], "_migAB_1W"), collapse = ""),
            paste(c(modelOrder[datA1w[1]], "_migA_1W"), collapse = ""),paste(c(modelOrder[datA1w[2]], "_migA_1W"), collapse = ""),paste(c(modelOrder[datA1w[3]], "_migA_1W"), collapse = ""),
            paste(c(modelOrder[datAB[1]], "_migAB"), collapse = ""),paste(c(modelOrder[datAB[2]], "_migAB"), collapse = ""),paste(c(modelOrder[datAB[3]], "_migAB"), collapse = ""),
            paste(c(modelOrder[datA[1]], "_migA"), collapse = ""),paste(c(modelOrder[datA[2]], "_migA"), collapse = ""),paste(c(modelOrder[datA[3]], "_migA"), collapse = ""))
for(i in 1:nFolders){
  model1 = rbind(model1, read.table(paste(models[1],"/",models[1],"_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  model2 = rbind(model2, read.table(paste(models[2],"/",models[2],"_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  model3 = rbind(model3, read.table(paste(models[3],"/",models[3],"_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  model4 = rbind(model4, read.table(paste(models[4],"/",models[4],"_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  model5 = rbind(model5, read.table(paste(models[5],"/",models[5],"_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  model6 = rbind(model6, read.table(paste(models[6],"/",models[6],"_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  model7 = rbind(model7, read.table(paste(models[7],"/",models[7],"_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  model8 = rbind(model8, read.table(paste(models[8],"/",models[8],"_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  model9 = rbind(model9, read.table(paste(models[9],"/",models[9],"_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  model10 = rbind(model10, read.table(paste(models[10],"/",models[10],"_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  model11 = rbind(model11, read.table(paste(models[11],"/",models[11],"_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  model12 = rbind(model12, read.table(paste(models[12],"/",models[12],"_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  model13 = rbind(model13, read.table(paste(models[13],"/",models[13],"_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  model14 = rbind(model14, read.table(paste(models[14],"/",models[14],"_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
  model15 = rbind(model15, read.table(paste(models[15],"/",models[15],"_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss])
}
for(i in 1:ncol(model1)){
  model1[which(model1[,i]=="NaN"),i]=mean(model1[,i], na.rm=T)
  model2[which(model2[,i]=="NaN"),i]=mean(model2[,i], na.rm=T)
  model3[which(model3[,i]=="NaN"),i]=mean(model3[,i], na.rm=T)
  model4[which(model4[,i]=="NaN"),i]=mean(model4[,i], na.rm=T)
  model5[which(model5[,i]=="NaN"),i]=mean(model5[,i], na.rm=T)
  model6[which(model6[,i]=="NaN"),i]=mean(model6[,i], na.rm=T)
  model7[which(model7[,i]=="NaN"),i]=mean(model7[,i], na.rm=T)
  model8[which(model8[,i]=="NaN"),i]=mean(model8[,i], na.rm=T)
  model9[which(model9[,i]=="NaN"),i]=mean(model9[,i], na.rm=T)
  model10[which(model10[,i]=="NaN"),i]=mean(model10[,i], na.rm=T)
  model11[which(model11[,i]=="NaN"),i]=mean(model11[,i], na.rm=T)
  model12[which(model12[,i]=="NaN"),i]=mean(model12[,i], na.rm=T)
  model13[which(model13[,i]=="NaN"),i]=mean(model13[,i], na.rm=T)
  model14[which(model14[,i]=="NaN"),i]=mean(model14[,i], na.rm=T)
  model15[which(model15[,i]=="NaN"),i]=mean(model15[,i], na.rm=T)
  
}
########### Test of top 3 of each model
x=c(rep(1:nModels, each=nrow(model1)))
obs = NULL
for(i in 1:nnetReplicates){
  obs <- rbind(obs,target)
}

res=model_selection_abc_nnet(target=obs, x=x, sumstat= rbind(model1,model2, model3, model4, model5, model6, model7, model8, model9, model10, model11, model12, model13, model14, model15)
                             , tol=(EpReplicates*nModelsFinal/nModels)/(nModels*nrow(model1))
                             , noweight=F, rejmethod=F, nb.nnet=35, size.nnet=10, output=paste(c(dataset,"_top3_0.05.txt"),collapse = ""))
write.table(models, file = paste(c("Models_",dataset,".txt"),collapse = ""), col.names = F, row.names = F, quote = F, sep = "\t")

################### generating posterior distribution of top 3
datFinal <- order(colMeans(read.table(paste(c(dataset,"_top3_0.05.txt"), collapse = ""))), decreasing = T)
source("../cv4estimations.R")
target = NULL
tmp=read.table(paste(c(dataset,"_ObsData.txt"),collapse = ""), skip = 1, h = F)[,ss]
priorCol = 7
for(i in 1:10){target=rbind(target, tmp)}
for(i in 1:3){
  prior = NULL
  post = NULL
  chosenModel <- models[datFinal[i]]
  print(chosenModel)
  if("migAB" %in% strsplit(chosenModel, "_")[[1]] | "migA" %in% strsplit(chosenModel, "_")[[1]]){
    priorCol <- 23
  }
  if("noMig" %in% strsplit(chosenModel, "_")[[1]]){
    priorCol <- 7
  }
  priorHead <-read.table(paste(chosenModel,"/",chosenModel,"_",i,"/output.txt", sep = ""), skip = 0, h = F)[1,1:priorCol]
  for(i in 1:nFolders){
    tmpPost <- read.table(paste(chosenModel,"/",chosenModel,"_",i,"/ABCstat.txt", sep = ""), skip=1, h = F)[,ss]
    tmpPrior <-read.table(paste(chosenModel,"/",chosenModel,"_",i,"/output.txt", sep = ""), skip = 1, h = F)[,1:priorCol]
    if(nrow(tmpPost)==nrow(tmpPrior)){
      post <- rbind(post, tmpPost)
      prior <- rbind(prior,tmpPrior)
    }
  }
  bb=rbind(apply(prior, MARGIN=2, FUN="min"), apply(prior, MARGIN=2, FUN="max"))
  remCol = NULL
  for(i in 1:ncol(bb)){
    if(bb[1,i] >= bb[2,i]){
      remCol <- c(remCol,i)
    }
  }
  prior <- prior[,-remCol]
  priorHead <- priorHead[,-remCol]
  
  mod <- abc_nnet_multivar(target=target, x=prior, sumstat=post, tol=10000/nrow(prior), rejmethod=F, noweight=F, transf=rep("logit", ncol(prior)), bb=rbind(apply(prior, MARGIN=2, FUN="min"), 
                                                                                                                                                            apply(prior, MARGIN=2, FUN="max")), nb.nnet=25, size.nnet=10, trace=T, output=paste(c("Posterior_",dataset,chosenModel), collapse = ""),apply(prior, MARGIN=2, FUN="max"))
  write.table(priorHead, file = paste(c("Posterior_",dataset,chosenModel,"_Header"), collapse = ""), col.names = F, row.names = F, quote = F, sep = "\t")
}
