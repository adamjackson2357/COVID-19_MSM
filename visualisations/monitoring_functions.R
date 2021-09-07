MonitorParams=function(models, param, legCol=NULL, name=expression(theta), init=1, legend=FALSE, factor=1, log_scale=FALSE, sigma=TRUE){
  # Prepare the seed colours
  if (is.null(legCol)){
    legCol <- c('darkgreen', 'orange', 'darkred', 'skyblue', "goldenrod2")
    legCol=colorRampPalette(legCol)(length(models))
  }
  
  # Sequence of data points to consider (skip points for efficiency with factor>1)
  Myseq <- seq(init,NIter,by=factor)
  
  # Log-transformation
  if (log_scale){
    for (k in 1:length(models)){
      x=eval(parse(text=models[k]))
      x[[param]]=log(x[[param]])
      x[[param]][is.infinite(x[[param]])]=NA
      assign(models[k], x)
    }
  }
  
  # Prepare limits of the Y-axis
  Ymin=min(eval(parse(text=paste0(models[1],"$", param)))[Myseq], na.rm=TRUE)
  for (k in 1:length(models)){
    Ymin=min(Ymin, min(eval(parse(text=paste0(models[k],"$", param)))[Myseq], na.rm=TRUE))
  }
  Ymax=min(eval(parse(text=paste0(models[1],"$", param)))[Myseq], na.rm=TRUE)
  for (k in 1:length(models)){
    Ymax=max(Ymax, max(eval(parse(text=paste0(models[k],"$", param)))[Myseq], na.rm=TRUE))
  }
  
  # Prepare the legend
  if (legend){
    Ymax=Ymax+(Ymax-Ymin)/5
  }
  
  # Making plot
  par(mar=c(5,5,1,5))
  # if (log_scale){
  #   plot(eval(parse(text=paste0(models[1],"$Iter")))[Myseq]/1000,log(eval(parse(text=paste0(models[1], "$", param)))[Myseq]),
  #        ylim=c(Ymin,Ymax), type="l", 
  #        ylab=name,xlab="Iteration (x1000)",main="",cex.lab=1.5,cex.axis=1.5,cex.main=1,col=legCol[1])  
  # } else {
  plot(eval(parse(text=paste0(models[1],"$Iter")))[Myseq]/1000,eval(parse(text=paste0(models[1], "$", param)))[Myseq],
       ylim=c(Ymin,Ymax), type="l", 
       ylab=eval(parse(text=paste0("expression(",name,")"))),
       xlab="Iteration (x1000)",main="",cex.lab=1.5,cex.axis=1.5,cex.main=1,col=legCol[1])
  # }
  abline(v=axTicks(1),lty=3,col='grey',lwd=1.5)
  abline(h=axTicks(2),lty=3,col='grey',lwd=1.5)
  for (k in 1:length(models)){
    # if (log_scale){
    #   lines(eval(parse(text=paste0(models[k],"$Iter")))[Myseq]/1000,log(eval(parse(text=paste0(models[k], "$", param)))[Myseq]), col=legCol[k])
    # } else {
    lines(eval(parse(text=paste0(models[k],"$Iter")))[Myseq]/1000,eval(parse(text=paste0(models[k], "$", param)))[Myseq], col=legCol[k])
    # }
  }
  if (sigma){
    par(new=TRUE)
    plot(eval(parse(text=paste0(models[1],"$Iter")))[Myseq]/1000,eval(parse(text=paste0(models[1], "$sigma_", tolower(param), "_", tolower(param))))[Myseq],
         type="l", xaxt="n", yaxt="n",
         ylab="",xlab="",main="",cex.lab=1.5,cex.axis=1.5,cex.main=1,col="red")
    axis(side=4, at=axTicks(2), cex.axis=1.5, col="red", col.axis="red")
    mtext(text=eval(parse(text=paste0("expression(sigma[",name,"])"))),side=4,line=3,cex=1.5,col="red")
  }
  
  if (legend){
    legTxt <- paste('Seed', c(1:length(models)))
    legend("top",col=legCol,legend=legTxt,lty=1,lwd=3,border=F,cex=1.5,bg='white',horiz=TRUE)
  }
}


GetPosterior=function(models, param, legCol=NULL, name=expression(theta), init=1, legend=FALSE, factor=1){
  # Prepare the seed colours
  if (is.null(legCol)){
    legCol <- c('darkgreen', 'orange', 'darkred', 'skyblue', "goldenrod2")
    legCol=colorRampPalette(legCol)(length(models))
  }
  
  # Sequence of data points to consider (skip points for efficiency with factor>1)
  Myseq <- seq(init,NIter,by=factor)
  
  # Read the MCMC output
  for (k in 1:length(models)){
    tmp=density(eval(parse(text=paste0(models[k],"$", param)))[Myseq])
    assign(paste0("Density",k),tmp)
  }
  
  # Prepare limits of the X-axis
  Xmin=min(eval(parse(text=paste0(models[1],"$", param)))[Myseq])
  for (k in 1:length(models)){
    Xmin=min(Xmin, min(eval(parse(text=paste0(models[k],"$", param)))[Myseq]))
  }
  Xmax=min(eval(parse(text=paste0(models[1],"$", param)))[Myseq])
  for (k in 1:length(models)){
    Xmax=max(Xmax, max(eval(parse(text=paste0(models[k],"$", param)))[Myseq]))
  }
  
  # Prepare limits of the Y-axis
  Ymax=max(Density1$y)
  for (k in 1:length(models)){
    Ymax=max(Ymax, max(eval(parse(text=paste0("Density",k,"$y")))))
  }
  
  # Prepare the legend
  if (legend){
    Ymax=Ymax+Ymax/5
    legTxt=NULL
    for (k in 1:length(models)){
      legTxt=c(legTxt, paste('Mean=',round(mean(eval(parse(text=paste0(models[k], "$", param)))[Myseq]),digits=2),
                             ' [',round(quantile(eval(parse(text=paste0(models[k], "$", param)))[Myseq],probs = 0.025),digits=2),
                             ';',round(quantile(eval(parse(text=paste0(models[k], "$", param)))[Myseq],probs = 0.975),digits=2),']',sep=''))
    }
  }
  
  # Making plot
  plot(Density1,xlim=c(Xmin,Xmax), ylim=c(0,Ymax),type="l",lwd=2,
       xlab=name,main="",cex.lab=1.5,cex.axis=1.5,cex.main=1,col=legCol[1])
  abline(v=axTicks(1),lty=3,col='grey',lwd=1.5)
  abline(h=axTicks(2),lty=3,col='grey',lwd=1.5)
  for (k in 1:length(models)){
    lines(eval(parse(text=paste0("Density",k))),col=legCol[k],lwd=2)
  }
  
  if(legend){
    legend('top',col=legCol,legend=legTxt,lty=1,lwd=2,bg="white",cex=1.5,ncol=ifelse(length(models)>2,yes=2,no=1))
  }
}


MonitorLikelihood=function(models, param, legCol=NULL, name=expression(theta), init=1, legend=FALSE, factor=1){
  # Prepare the seed colours
  if (is.null(legCol)){
    legCol <- c('darkgreen', 'orange', 'darkred', 'skyblue', "goldenrod2")
    legCol=colorRampPalette(legCol)(length(models))
  }
  
  # Sequence of data points to consider (skip points for efficiency with factor>1)
  Myseq <- seq(init,NIter,by=factor)
  
  # Read the MCMC output
  Xmin=Xmax=Ymin=Ymax=NULL
  for (k in 1:length(models)){
    TmpX=eval(parse(text=paste0(models[k],"$", param)))[Myseq]
    Xmin=c(Xmin,min(TmpX))
    Xmax=c(Xmax,max(TmpX))
    assign(paste0("TmpX",k),TmpX)
  }
  Xmin=min(Xmin)
  Xmax=max(Xmax)
  for (k in 1:length(models)){
    TmpY=eval(parse(text=paste0(models[k],"$L")))[Myseq] # likelihood at the end of serial updates
    Ymin=c(Ymin, min(TmpY))
    Ymax=c(Ymax, max(TmpY))
    assign(paste0("TmpY",k),TmpY)
  }
  Ymin=min(Ymin)
  Ymax=max(Ymax)
  
  # Making plot
  plot(TmpX1,TmpY1,xlim=c(Xmin,Xmax),ylim=c(Ymin,Ymax), type="p", 
       xlab=name,ylab="Likelihood",main="",cex.lab=1.5,cex.axis=1.5,cex.main=1,
       col=adjustcolor(legCol[1],alpha.f = 0.6),pch=19)
  abline(v=axTicks(1),lty=3,col='grey',lwd=1.5)
  abline(h=axTicks(2),lty=3,col='grey',lwd=1.5)
  for (k in 1:length(models)){
    points(eval(parse(text=paste0("TmpX",k))),eval(parse(text=paste0("TmpY",k))),type="p",
           col=adjustcolor(legCol[k],alpha.f = 0.6),pch=19)
  }
  # Adding posterior mean and 5%-confidence regions
  for (k in 1:length(models)){
    abline(v=mean(eval(parse(text=paste0("TmpX",k)))),col=legCol[k],lty=1,lwd=1)
    abline(v=quantile(eval(parse(text=paste0("TmpX",k))),probs = 0.025),col=legCol[k],lty=2)
    abline(v=quantile(eval(parse(text=paste0("TmpX",k))),probs = 0.975),col=legCol[k],lty=2)
  }
}

