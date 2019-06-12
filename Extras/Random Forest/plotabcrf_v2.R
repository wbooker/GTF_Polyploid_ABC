plot.abcrf_v2 <- function(x, training, obs=NULL, n.var=20, pdf=FALSE, xlim=NULL, lda1=1, lda2=2, xlda=NULL, ...)
{
  if (!inherits(obs, "data.frame") && !is.null(obs) ) 
    stop("obs needs to be a data.frame object or NULL")
  if (!inherits(training, "data.frame"))
    stop("training needs to be a data.frame object")
  if (!is.null(xlim) && !is.numeric(xlim))
    stop("xlim needs to be a numeric object or NULL")

  if (length(x$group)!=0)
  {
    ngroup <- length(x$group)
    varn <- x$formula[[2]]
    training[[as.character(varn)]] <- as.vector(training[[as.character(varn)]])
    allmod <- unique(training[[as.character(varn)]])
    for (k in 1:ngroup) for (l in 1:length(x$group[[k]])) 
      training[[as.character(varn)]][which(training[[as.character(varn)]]==x$group[[k]][l])] <- paste("g",k,sep="")
    if (!setequal(allmod,unlist(x$group)))
    {
      diffe <- setdiff(allmod,unlist(x$group))
      for (l in 1:length(diffe)) training <- training[-which(training[[as.character(varn)]]==diffe[l]),]
    }
    training[[as.character(varn)]] <- as.factor(training[[as.character(varn)]])
  }
  
  old.par <- par(no.readonly = TRUE)
  
  if (length(x$model.rf$variable.importance)<20) 
    n.var <- length(x$model.rf$variable.importance)
  
  mf <- match.call(expand.dots=FALSE)
  mf <- mf[1]
  mf$formula <- x$formula
  mf$data <- training
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame() )
  mt <- attr(mf, "terms")
  modindex <- model.response(mf)
  
  if (x$lda) {
    if (pdf) { 
      pdf("graph_varImpPlot.pdf")
      variableImpPlot(x, n.var=n.var, xlim=xlim)
      dev.off()
    }
    variableImpPlot(x, n.var=n.var, xlim=xlim)
    nmod <- length(x$model.rf$forest$levels)
    nstat <- x$model.rf$num.independent.variables
    projections <- predict(x$model.lda, training)$x
    if  (!is.null(obs)) projobs <- predict(x$model.lda, obs)$x
    coloris <- rainbow(nmod)
    colo <- coloris[modindex]
    if (nmod > 2) {
      print(projobs[lda1])
      print(projobs[lda2])
      if (pdf)
      {
        pdf(paste(c("graph_",lda1,"_",lda2,".pdf"), collapse = ""))
        par(mar=par()$mar+c(0,0,0,5), xpd=TRUE)
        plot(projections[,c(lda1,lda2)], col=colo, pch=3)
        legend("topright", legend = as.character(x$model.rf$forest$levels), col = coloris, 
               pch = 15, bty = "o", pt.cex = 2, cex = .8, horiz = FALSE, 
               inset = c(-.16, 0), ncol = 2, title = "Models", bg = "white")
        if  (!is.null(obs)) points(projobs[lda1],projobs[lda2],pch="*",cex=5.3)
        
        dev.off()
      }
      par(mar=par()$mar+c(0,0,0,5), xpd=TRUE)
      plot(projections[,c(lda1,lda2)], col=colo, pch=3)
      legend("topright", legend = as.character(x$model.rf$forest$levels), col = coloris, 
             pch = 15, bty = "o", pt.cex = 2, cex = .8, horiz = FALSE, 
             inset = c(-.16, 0), ncol = 2, title = "Models", bg = "white")
      if  (!is.null(obs)) points(projobs[lda1],projobs[lda2],pch="*",cex=5.3)
    } else {
      
      l1 <- x$model.rf$forest$levels[lda1]
      l2 <- x$model.rf$forest$levels[lda2]
      d1 <- density(projections[modindex == paste(c("l",lda1),collapse = "")])
      d2 <- density(projections[modindex == paste(c("l",lda2),collapse = "")])
      coloris <- c("blue", "orange")
      if(is.null(xlda)){
        xlda <- range(c(d1$x, d2$x))
        ylda <- c(0, 1.2*max(c(d1$y, d2$y)))
      }
      else{
        ylda <- xlda
      }
      xrange <- range(c(d1$x, d2$x))
      yrange <- c(0, 1.2*max(c(d1$y, d2$y)))
      if (pdf)
      {
        pdf(paste(c("graph_lda",lda1,"_lda",lda2,".pdf"), collapse = ""))
        plot(d1, xlim = xlda, ylim = ylda,
             col=coloris[1], main="", xlab="")
        lines(d2, col=coloris[2])
        legend("topleft", legend = as.character(x$model.rf$forest$levels), col = coloris, 
               cex = .8, horiz = TRUE, lty=1, bty="o",
               inset = c(.01, .01), title = "Models", bg = "white")
        if  (!is.null(obs)) abline(v=projobs)
        dev.off()
      }
      plot(d1, xlim = xlda, ylim = ylda,
           col=coloris[1], main="", xlab="")
      lines(d2, col=coloris[2])
      legend("topleft", legend = as.character(x$model.rf$forest$levels), col = coloris, 
             cex = .8, horiz = TRUE, lty=1, bty="o",
             inset = c(.01, .01), title = "Models", bg = "white")
      if  (!is.null(obs)) abline(v=projobs)
    }
  } else {
    if (pdf)
    {
      pdf("graph_varImpPlot.pdf")
      variableImpPlot(x , n.var=n.var, xlim=xlim)
      dev.off()
    }
    variableImpPlot(x , n.var=n.var, xlim=xlim)
  }
  par(old.par)
}

plot.abcrf_all <- function(x, training, obs=NULL, n.var=20, pdf=FALSE, xlim=NULL, xlda=NULL, ...)
{
  lda1_list <- c(1,1,1,2,2,3)
  lda2_list <- c(2,3,4,3,4,4)
  if (!inherits(obs, "data.frame") && !is.null(obs) ) 
    stop("obs needs to be a data.frame object or NULL")
  if (!inherits(training, "data.frame"))
    stop("training needs to be a data.frame object")
  if (!is.null(xlim) && !is.numeric(xlim))
    stop("xlim needs to be a numeric object or NULL")
  for(i in 1:length(lda1_list)){
    lda1 <- lda1_list[i]
    lda2 <- lda2_list[i]
    if (length(x$group)!=0)
    {
      ngroup <- length(x$group)
      varn <- x$formula[[2]]
      training[[as.character(varn)]] <- as.vector(training[[as.character(varn)]])
      allmod <- unique(training[[as.character(varn)]])
      for (k in 1:ngroup) for (l in 1:length(x$group[[k]])) 
        training[[as.character(varn)]][which(training[[as.character(varn)]]==x$group[[k]][l])] <- paste("g",k,sep="")
      if (!setequal(allmod,unlist(x$group)))
      {
        diffe <- setdiff(allmod,unlist(x$group))
        for (l in 1:length(diffe)) training <- training[-which(training[[as.character(varn)]]==diffe[l]),]
      }
      training[[as.character(varn)]] <- as.factor(training[[as.character(varn)]])
    }
    
    old.par <- par(no.readonly = TRUE)
    
    if (length(x$model.rf$variable.importance)<20) 
      n.var <- length(x$model.rf$variable.importance)
    
    mf <- match.call(expand.dots=FALSE)
    mf <- mf[1]
    mf$formula <- x$formula
    mf$data <- training
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame() )
    mt <- attr(mf, "terms")
    modindex <- model.response(mf)
    
    if (x$lda) {
      if(i == 1){
        if (pdf) { 
          pdf("graph_varImpPlot.pdf")
          variableImpPlot(x, n.var=n.var, xlim=xlim)
          dev.off()
        }
        variableImpPlot(x, n.var=n.var, xlim=xlim)
      }
      nmod <- length(x$model.rf$forest$levels)
      nstat <- x$model.rf$num.independent.variables
      projections <- predict(x$model.lda, training)$x
      if  (!is.null(obs)) projobs <- predict(x$model.lda, obs)$x
      #coloris <- c("#F41802", "#F44B02","#F48A02","#E0FF02", "#44F402", "#02F474", "#0212FF", "#A702FF", "#FF02C2")
      coloris <- rainbow(nmod)
      colo <- coloris[modindex]
      if (nmod > 2) {
        print(projobs[lda1])
        print(projobs[lda2])
        if (pdf)
        {
          pdf(paste(c("graph_",lda1,"_",lda2,".pdf"), collapse = ""))
          par(mar=par()$mar+c(0,0,0,5), xpd=TRUE)
          plot(projections[,c(lda1,lda2)], col=colo, pch=3)
          legend("topright", legend = as.character(x$model.rf$forest$levels), col = coloris, 
                 pch = 15, bty = "o", pt.cex = 2, cex = .8, horiz = FALSE, 
                 inset = c(-.16, 0), ncol = 2, title = "Models", bg = "white")
          if  (!is.null(obs)) points(projobs[lda1],projobs[lda2],pch="*",cex=5.3)
          
          dev.off()
        }
        par(mar=par()$mar+c(0,0,0,5), xpd=TRUE)
        plot(projections[,c(lda1,lda2)], col=colo, pch=3)
        legend("topright", legend = as.character(x$model.rf$forest$levels), col = coloris, 
               pch = 15, bty = "o", pt.cex = 2, cex = .8, horiz = FALSE, 
               inset = c(-.16, 0), ncol = 2, title = "Models", bg = "white")
        if  (!is.null(obs)) points(projobs[lda1],projobs[lda2],pch="*",cex=5.3)
      } else {
        
        l1 <- x$model.rf$forest$levels[lda1]
        l2 <- x$model.rf$forest$levels[lda2]
        d1 <- density(projections[modindex == paste(c("l",lda1),collapse = "")])
        d2 <- density(projections[modindex == paste(c("l",lda2),collapse = "")])
        coloris <- c("blue", "orange")
        if(is.null(xlda)){
          xlda <- range(c(d1$x, d2$x))
          ylda <- c(0, 1.2*max(c(d1$y, d2$y)))
        }
        else{
          ylda <- xlda
        }
        xrange <- range(c(d1$x, d2$x))
        yrange <- c(0, 1.2*max(c(d1$y, d2$y)))
        if (pdf)
        {
          pdf(paste(c("graph_lda",lda1,"_lda",lda2,".pdf"), collapse = ""))
          plot(d1, xlim = xlda, ylim = ylda,
               col=coloris[1], main="", xlab="")
          lines(d2, col=coloris[2])
          legend("topleft", legend = as.character(x$model.rf$forest$levels), col = coloris, 
                 cex = .8, horiz = TRUE, lty=1, bty="o",
                 inset = c(.01, .01), title = "Models", bg = "white")
          if  (!is.null(obs)) abline(v=projobs)
          dev.off()
        }
        plot(d1, xlim = xlda, ylim = ylda,
             col=coloris[1], main="", xlab="")
        lines(d2, col=coloris[2])
        legend("topleft", legend = as.character(x$model.rf$forest$levels), col = coloris, 
               cex = .8, horiz = TRUE, lty=1, bty="o",
               inset = c(.01, .01), title = "Models", bg = "white")
        if  (!is.null(obs)) abline(v=projobs)
      }
    } else {
      if (pdf)
      {
        pdf("graph_varImpPlot.pdf")
        variableImpPlot(x , n.var=n.var, xlim=xlim)
        dev.off()
      }
      variableImpPlot(x , n.var=n.var, xlim=xlim)
    }
    par(old.par)
  }  
  
}