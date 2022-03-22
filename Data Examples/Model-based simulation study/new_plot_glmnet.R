
myplot.glmnet=function(x, xvar=c("norm","lambda","dev"),label=FALSE,...){
  xvar=match.arg(xvar)
  myplotCoef(x$beta,lambda=x$lambda,df=x$df,dev=x$dev.ratio,label=label,xvar=xvar,...)
}

myplotCoef=function(beta,norm,lambda,df,dev,mycols=1,mylty=1,label=FALSE,xvar=c("norm","lambda","dev"),xlab=iname,ylab="Coefficients",...){
  ##beta should be in "dgCMatrix" format
  #which=nonzeroCoef(beta)
  which = 1:nrow(beta)
  mycols <- mycols[which]
  nwhich=length(which)
  switch(nwhich+1,#we add one to make switch work
         "0"={
           warning("No plot produced since all coefficients zero")
           return()
         },
         "1"=warning("1 or less nonzero coefficients; glmnet plot is not meaningful")
  )
  beta=as.matrix(beta[which,,drop=FALSE])
  xvar=match.arg(xvar)
  switch(xvar,
         "norm"={
           index=if(missing(norm))apply(abs(beta),2,sum)else norm
           iname="L1 Norm"
           approx.f=1
         },
         "lambda"={
           index=log(lambda)
           iname="Log Lambda"
           approx.f=0
         },
         "dev"= {
           index=dev
           iname="Fraction Deviance Explained"
           approx.f=1
         }
  )
  dotlist=list(...)
  type=dotlist$type
  if(is.null(type))
    matplot(index,t(beta),lty=mylty,xlab=xlab,ylab=ylab,type="l",col=mycols,...)
  else matplot(index,t(beta),lty=mylty,xlab=xlab,ylab=ylab, col=mycols, ...)
  atdf=pretty(index)
  ### compute df by interpolating to df at next smaller lambda
  ### thanks to Yunyang Qian
  prettydf=approx(x=index,y=df,xout=atdf,rule=2,method="constant",f=approx.f)$y
  # prettydf=ceiling(approx(x=index,y=df,xout=atdf,rule=2)$y)
 # axis(3,at=atdf,labels=prettydf,tcl=NA)
  if(label){
    nnz=length(which)
    xpos=max(index)
    pos=4
    if(xvar=="lambda"){
      xpos=min(index)
      pos=2
    }
    xpos=rep(xpos,nnz)
    ypos=beta[,ncol(beta)]
    text(xpos,ypos,paste(which),cex=.5,pos=pos)
  }
  
}


environment(myplot.glmnet) <- asNamespace('glmnet')
environment(myplotCoef) <- asNamespace('glmnet')