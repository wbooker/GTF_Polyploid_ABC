abc_nnet_multivar <- function(target,x,sumstat,tol,gwt,rejmethod=F,noweight=F,transf="none",bb=c(0,0),nb.nnet=10,size.nnet=5,trace=T, output="Posterior_pseudo")
  
{
  
  require(nnet)
  
  
  
  # target is the set of target summary stats
  
  # x is the parameter vector (long vector of numbers from the simulations) and is the dependent variable for the regression
  
  # x can also be a matrix for multi-dimensional models. Each column corresponds to a parameter and each row to a simulation.
  
  # sumstat is an array of simulated summary stats (i.e. independent variables).
  
  # tol is the required proportion of points nearest the target values
  
  # gwt is a vector with T/F weights, weighting out any 'bad' values (determined by the simulation program - i.e. nan's etc)
  
  # if noweight=T, no Epanechnikov weights are calculated
  
  # if rejmethod=T it doesn't bother with the regression, and just does rejection.
  
  # nb.nnet>1 is the number of trained neural networks, the more neural nets the more robust is the inference
  
  # size.nnet is the number of hidden network in the regression. Typically >= the number of parameters in the model
  
  # transf the vector of transformation for the parameter, ex:transf=c("none","logit","log")
  
  # bb the vector of bounds for the logit transformation ex:bb=cbind(c(0,0),c(0,2),c(0,0)) (The second column is the only one to be taken into account)
  
  # If trace=T print messages during the algorithm
  
  
  
  
  
  
  
  # If rejmethod=F it returns a list with the following components:-
  
  
  
  # $x regression adjusted values
  
  # $vals - unadjusted values in rejection region (i.e. normal rejection)
  
  # $wt - the regression weight (i.e. the Epanechnikov weight)
  
  # $ss - the sumstats corresponding to these points
  
  # $predmean - estimate of the posterior mean
  
  
  
  if(class(x)=="numeric")
    
  {
    
    bb<-cbind(bb)
    
    x<-cbind(x)
    
  }
  
  if(rejmethod)
    
    transf<-rep("none", dim(x)[2])
  
  normalise <- function(x,y){
    
    
    
    if(mad(y) == 0)
      
      return (x)
    
    else
      
      return (x/mad(y))
    
  }
  
  
  
  ####Define the weight-decay paramaeter
  
  repet<-floor(nb.nnet/3)+1
  
  the_decay<-rep(c(10^(-4),10^(-3),10^(-2)),repet)[1:nb.nnet]
  
  
  
  lt<-dim(x)[2]
  
  nb_simu<-dim(sumstat)[1]
  
  for (i in 1:lt)
    
  {
    
    if(sum(transf[i] == c("none","log","logit")) == 0){
      
      stop("transf must be none, log, or logit")
      
      
      
    }
    
    if(transf[i]=="logit"){
      
      if(bb[1,i] >= bb[2,i]){
        
        stop("bounds wrong for logit")
        
        
        
      }
      
    }
    
  }
  
  if(missing(gwt))gwt <- rep(T,length(sumstat[,1]))
  
  
  
  nss <- length(sumstat[1,])
  
  
  
  
  
  # scale everything
  
  
  
  scaled.sumstat <- sumstat
  
  
  
  for(j in 1:nss){
    
    
    
    scaled.sumstat[,j] <- normalise(sumstat[,j],sumstat[,j][gwt])
    
  }
  
  target.s.tmp <- target
  
  
  
  for(j in 1:nss){
    
    
    
    target.s.tmp[,j] <- normalise(target[,j],sumstat[,j][gwt])
    
  }
  
  
  #déplacé (origine : après calcul de abstol)
  for (i in 1:lt)
    
  {
    
    if(transf[i] == "log"){
      
      if(min(x[,i]) <= 0){
        
        print("log transform: val out of bounds - correcting")
        
        x.tmp <- ifelse(x[,i] <= 0,max(x[,i]),x[,i])
        
        x.tmp.min <- min(x.tmp)
        
        xx[,i] <- ifelse(x[,i] <= 0, x.tmp.min,x[,i])
        
      }
      
      x[,i] <- log(x[,i])
      
    }
    
    else if(transf[i] == "logit"){
      
      if(min(x[,i]) <= bb[1,i]){
        
        x.tmp <- ifelse(x[,i] <= bb[1,i],max(x[,i]),x[,i])
        
        x.tmp.min <- min(x.tmp)
        
        x[,i] <- ifelse(x[,i] <= bb[1,i], x.tmp.min,x[,i])
        
        
        
      }
      
      if(max(x[,i]) >= bb[2,i]){
        
        x.tmp <- ifelse(x[,i] >= bb[2,i],min(x[,i]),x[,i])
        
        x.tmp.max <- max(x.tmp)
        
        x[,i] <- ifelse(x[,i] >= bb[2,i], x.tmp.max,x[,i])
        
        
        
      }
      
      x[,i] <- (x[,i]-bb[1,i])/(bb[2,i]-bb[1,i])
      
      x[,i] <- log(x[,i]/(1-x[,i]))
      
    }
    
  }
  
  
  
  
  #boucle le long des targets
  for(L in 1:nrow(target.s.tmp)){
    target.s=as.numeric(target.s.tmp[L,])
    sum1=dst=abstol=wt1=regwt=l1=fit1=ll=predmean=array_pred=mean_averaged=my_residuals=predvar=var_averaged=the_sd=predsd=x.tmp=x.tmp.min=xx=x.tmp.max=NULL
    
    
    
    # calc euclidean distance
    
    sum1 <- 0
    
    for(j in 1:nss){
      
      sum1 <- sum1 + (scaled.sumstat[,j]-target.s[j])^2
      
    }
    
    dst <- sqrt(sum1)
    
    # includes the effect of gwt in the tolerance
    
    dst[!gwt] <- floor(max(dst[gwt])+10)
    
    
    
    
    
    # wt1 defines the region we're interested in
    
    abstol <- quantile(dst,tol)
    
    wt1 <- dst <= abstol
    
    
    
    
    
    if(rejmethod){
      
      regwt <- 1-dst[wt1]^2/abstol^2
      
      l1 <- list(x=cbind(x[wt1,]),wt=regwt,ind=wt1,dst=dst)
      
    }
    
    else{
      
      regwt <- 1-dst[wt1]^2/abstol^2
      
      if(noweight)
        
        regwt <- rep(1,length(regwt))
      
      
      
      ll<-NULL
      
      if(trace==TRUE)
        
        cat("Regression of the mean ")
      
      for (i in 1:nb.nnet)
        
      {
        
        if(trace==TRUE)
          
          cat(i," ")
        
        fit1 <- nnet(scaled.sumstat[wt1,],x[wt1,],weights=regwt,decay=the_decay[i],size=size.nnet,linout=T,maxit=500,trace=F)
        
        ll<-c(ll,list(fit1))
        
      }
      
      #Compute the residuals
      
      predmean<-NULL
      
      array_pred<-array(dim=c(nb.nnet,sum(wt1),lt))
      
      for (i in 1:nb.nnet)
        
      {
        
        array_pred[i,,]<-ll[[i]]$fitted.values
        
        predmean<-cbind(predmean,as.numeric(predict(ll[[i]],data.frame(rbind(target.s)))))	
        
      }
      
      mean_averaged<-NULL
      
      for (j in 1:lt)
        
        mean_averaged<-cbind(mean_averaged,apply(array_pred[,,j],FUN=median,MARGIN=2))
      
      predmean<-apply(predmean,FUN=median,MARGIN=1)
      
      
      
      
      
      my_residuals<-(x[wt1,]-mean_averaged)
      
      
      
      #Fit a neural network for predicting the conditional variance
      
      if(trace==TRUE)
        
        cat("\nRegression of the variance ")
      
      ll<-NULL
      
      for (i in 1:nb.nnet)
        
      {
        
        if(trace==TRUE)
          
          cat(i," ")
        
        fit2 <- nnet(scaled.sumstat[wt1,],log(my_residuals^2),weights=regwt,decay=the_decay[i],size=size.nnet,linout=T,maxit=500,trace=F)
        
        ll<-c(ll,list(fit2))
        
      }
      
      if(trace==TRUE)
        
        cat("\n")
      
      predvar<-NULL
      
      array_pred<-array(dim=c(nb.nnet,sum(wt1),lt))
      
      for (i in 1:nb.nnet)
        
      {
        
        array_pred[i,,]<-ll[[i]]$fitted.values
        
        predvar<-cbind(predvar,as.numeric(predict(ll[[i]],data.frame(rbind(target.s)))))	
        
      }
      
      var_averaged<-NULL
      
      for (j in 1:lt)
        
        var_averaged<-cbind(var_averaged,apply(array_pred[,,j],FUN=median,MARGIN=2))
      
      the_sd<-sqrt(exp(var_averaged))
      
      predsd<-sqrt(exp(apply(predvar,FUN=median,MARGIN=1)))
      
      
      
      
      
      res_correc<-sapply(1:lt,FUN=function(i){predmean[i]+ ((predsd[i]*my_residuals[,i])/the_sd[,i]) })
      
      
      
      
      
      l1 <- list(x=cbind(res_correc),vals=cbind(x[wt1,]),wt=regwt,ss=sumstat[wt1,],predmean=predmean)
      
    }
    
    for (i in 1:lt)
      
    {
      
      if(transf[i] == "log"){
        
        l1$x[,i] <- exp(l1$x[,i])
        
        l1$vals[,i] <- exp(l1$vals[,i])
        
      }
      
      if(transf[i] == "logit"){
        
        l1$x[,i] <- exp(l1$x[,i])/(1+exp(l1$x[,i]))
        
        l1$x[,i] <- l1$x[,i]*(bb[2,i]-bb[1,i])+bb[1,i]
        
        l1$vals[,i] <- exp(l1$vals[,i])/(1+exp(l1$vals[,i]))
        
        l1$vals[,i] <- l1$vals[,i]*(bb[2,i]-bb[1,i])+bb[1,i]
        
      }
      
    }
    
    #l1
    write.table(l1$x, col.names=F, row.names=F, file=paste(output,L,sep=""))
    
  }
  
}
