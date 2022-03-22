###Simulation study variable selection
#NOTE 1: in Rstudio; press ALT+o for overview of code section headings
#This file may be used to replicate the analysis for Supplementary Section 3.1.4

#make sure working directory is set to where this file is
setwd("") 
pathResults <- "" #if desired, separate folder can be specified to save results..
pathFigures <- "" #..and figures

#Load libraries----
if(!requireNamespace("glmnet")) install.packages("glmnet")
library(glmnet)
if(!requireNamespace("ecpc")) install.packages("ecpc")
library(ecpc)
if(!requireNamespace("squeezy")) install.packages("squeezy")
library(squeezy)
if(!requireNamespace("ggpubr")) install.packages("ggpubr")
library(ggpubr)
if(!requireNamespace("gridGraphics")) install.packages("gridGraphics")
library(gridGraphics)
if(!requireNamespace("mvtnorm")) install.packages("mvtnorm")
library(mvtnorm)
source("new_plot_glmnet.R")

#Scenario 1----
set.seed(1234)
#Set simulation variables####
n<-100 #sample size of each training data set
p<-400 #number of covariates
pstar <- 24 #number of covariates in the model
p0 <- p - pstar
s2 <- 1
rho <- 0.8
beta <- c(rep(c(-1,1),pstar/2), rep(0,p0))

#Simulate data and co-data####
#Generate X 
ppair <- pstar/2 #number of variables which have a strongly correlated friend, not part of the model
indtruecor <- 1:ppair
indtruenoncor <- (1:pstar)[-indtruecor]
indfalsecor <- ((pstar + 1):(pstar+ppair))
indfalsenoncor <- ((pstar+ppair+1):p)

#generate bi-correlated 
sigma <- matrix(c(s2,rho,rho,s2),nrow=2)
Xcor <- rmvnorm(n=n,sigma=sigma)
Xmat <- matrix(nrow=n,ncol=p)
Xuncor <- matrix(rnorm(n*(p-ppair*2)),nrow=n)
dim(Xuncor)
Xmat[,c(indtruenoncor,indfalsenoncor)] <- matrix(rnorm(n*(p-ppair*2)),nrow=n)
for(j in 1:ppair) {
  Xcor <- rmvnorm(n=n,sigma=sigma)
  Xmat[,c(indtruecor[j],indfalsecor[j])] <- Xcor
}


#Generate response 
y <- Xmat %*% matrix(beta,ncol=1) + rnorm(n)

save(y,Xmat,n,p,pstar,rho,s2,beta,file=paste(pathResults,"simulp400n100.Rdata",sep=""))


#Co-data indices 
Gtrue <- 1:pstar
Gnottrue <- pstar + ppair + 1:pstar
Co1 <- c(Gtrue,Gnottrue)
Nonco1 <- (1:p)[-Co1]

#Weaker co-data 
Co2 <- sort(c(sample(1:pstar, pstar/2),Gnottrue))
Nonco2 <- (1:p)[-Co2]

#Random co-data 
Co3 <- sort(c(sample(1:pstar, pstar/16),Gnottrue))
Nonco3 <- (1:p)[-Co3]


#Fit models ########
#glmnet (ordinary elastic net)
gnfit <- glmnet(Xmat,y,family="gaussian",intercept=FALSE, alpha=1)

#squeezy strong co-data
grp <- list(Co1,Nonco1)
squeezyfit <- squeezy(y,Xmat,groupset=grp,intrcpt=FALSE,reCV=FALSE, alpha=1)
pfgr1 <- squeezyfit$lambdaApprox
pfgr1
pf <- rep(1,p)
pf[Co1] <- pfgr1[1];pf[Nonco1]<-pfgr1[2]
sqfit1 <- glmnet(Xmat,y,family="gaussian",intercept=FALSE,penalty.factor = pf, alpha=1)

#squeezy weaker co-data
grp <- list(Co2,Nonco2)
squeezyfit <- squeezy(y,Xmat,groupset=grp,intrcpt=FALSE,reCV=FALSE, alpha=1)
pfgr2 <- squeezyfit$lambdaApprox
pfgr2
pf <- rep(1,p)
pf[Co2] <- pfgr2[1];pf[Nonco2]<-pfgr2[2]
sqfit2 <- glmnet(Xmat,y,family="gaussian",intercept=FALSE,penalty.factor = pf, alpha=1)

#squeezy non-informative co-data
grp <- list(Co3,Nonco3)
squeezyfit <- squeezy(y,Xmat,groupset=grp,intrcpt=FALSE,reCV=FALSE, alpha=1)
pfgr3 <- squeezyfit$lambdaApprox
pfgr3
pf <- rep(1,p)
pf[Co3] <- pfgr3[1];pf[Nonco3]<-pfgr3[2]
sqfit3 <- glmnet(Xmat,y,family="gaussian",intercept=FALSE,penalty.factor = pf, alpha=1)

mycolors <- rep("red",p)
mycolors[1:pstar] <- "green"
myltype <- rep(1,p); myltype[mycolors=="red"] <- 2
lw <- rep(1,p)
lw[1:pstar] <- 2
save(sqfit1,sqfit2,sqfit3,gnfit,mycolors,lw,myltype,
     file=paste(pathResults,"inputplot_glmnet.Rdata",sep=""))

#Obtain roc curve information####
load(paste(pathResults,"inputplot_glmnet.Rdata",sep=""))

#coefficient matrices; rows are features, columns are penalties
coefs1 <- coefficients(sqfit1)[-1,] #remove intercept
coefs2 <- coefficients(sqfit2)[-1,] #remove intercept
coefs3 <- coefficients(sqfit3)[-1,] #remove intercept
coefs4 <- coefficients(gnfit)[-1,] #remove intercept

#computes FPR and TPR
TPRf <- function(vec,nnonzero){
  nbeta <- length(vec)
  betanz <- vec[1:nnonzero]
  betaz <- vec[-(1:nnonzero)]
  TPR <- length(which(abs(betanz) >= 10^(-10)))/nnonzero
  FPR <- length(which(abs(betaz) >= 10^(-10)))/(nbeta-nnonzero)
  return(c(FPR,TPR))
} 


#ROC
ROC1 <- rbind(t(apply(coefs1,2,TPRf, nnonzero=24)),c(1,1))
ROC2 <- rbind(t(apply(coefs2,2,TPRf, nnonzero=24)),c(1,1))
ROC3 <- rbind(t(apply(coefs3,2,TPRf, nnonzero=24)),c(1,1))
ROC4 <- rbind(t(apply(coefs4,2,TPRf, nnonzero=24)),c(1,1))


#Plots: general parameters----
wdth<-600
hght<-wdth*5/8
wdthpdf <- wdth/75
hghtpdf <- hght/75
ts <- 16 #basis text size in figures
ls <- 1.5 #basis line size in figures
ps <- 2 #basis point size in figures
sz <- 2 #point size
strk <- 1.5 #stroke size
palette <- "Dark2"

#Scenario 1 Load data plots----
load(paste(pathResults,"inputplot_glmnet.Rdata",sep=""))

#Scenario 1 plot traceplot----
figname<-paste(pathFigures,"FigSimsTraceplots_sc1.pdf",sep="")
pdf(width = wdthpdf, height = hghtpdf*1.4,
    file = figname)
layout(matrix(c(1,3,2,4),ncol=2,byrow=F))
par(mar=c(5, 4, 4, 2) + 0.1)
myplot.glmnet(sqfit1, xvar="lambda", label=FALSE, mycol=mycolors, mylty=myltype,
              lwd=2,ylab="Regression coefficient estimates",xlab=expression(log(lambda)),
              main="Scenario 1, case A")
legend("bottomright", title="Effect",legend=c("strong", "none"),
       col=unique(mycolors),lty=c(1,3),lwd=2, bg="white" )
myplot.glmnet(sqfit2, xvar="lambda", label=FALSE, mycol=mycolors, mylty=myltype,
              lwd=2,ylab="Regression coefficient estimates",xlab=expression(log(lambda)),
              main="Scenario 1, case B")
legend("bottomright", title="Effect",legend=c("strong", "none"),
       col=unique(mycolors),lty=c(1,3),lwd=2, bg="white")
myplot.glmnet(sqfit3, xvar="lambda", label=FALSE, mycol=mycolors, mylty=myltype,
              lwd=2,ylab="Regression coefficient estimates",xlab=expression(log(lambda)),
              main="Scenario 1, case C")
legend("bottomright", title="Effect",legend=c("strong", "none"),
       col=unique(mycolors),lty=c(1,3),lwd=2, bg="white")
myplot.glmnet(gnfit, xvar="lambda", label=FALSE, mycol=mycolors, mylty=myltype,
              lwd=2,ylab="Regression coefficient estimates",xlab=expression(log(lambda)),
              main="Scenario 1, case D")
legend("bottomright", title="Effect",legend=c("strong", "none"),
       col=unique(mycolors),lty=c(1,3),lwd=2, bg="white")
dev.off()



#Scenario 1 plot traceplot pairs glmnet----
colors <- rep("grey",p)
lw_temp <- rep(1,p)

pdf(paste(pathFigures,"traceplots12pairs_glmnet.pdf",sep=""),width=14,height=7)
par(mfrow=c(3,4), mar=c(2,2,1,1))
for(j in 1:ppair){
  #j <-7
  mycolors <- rep("grey",p)
  lw_temp <- rep(1,p)
  lw_temp[c(indtruecor[j],indfalsecor[j])] <- 3
  mycolors[indtruecor[j]] <- "green"
  mycolors[indfalsecor[j]] <- "red"
  myplot.glmnet(gnfit, xvar="lambda", label=FALSE, mycols=mycolors,
                lwd=lw_temp,xlab="",ylab="",ylim=c(-1.35,1.35))
  true <- beta[j]
  abline(h=true,lwd=2)
}
dev.off()

#Scenario 1 plot traceplot pairs squeezy----
pdf(paste(pathFigures,"traceplots12pairs_sqfit1.pdf",sep=""),width=14,height=7)
par(mfrow=c(3,4), mar=c(2,2,2,1))
for(j in 1:ppair){
  #j <-7
  mycolors <- rep("grey",p)
  lw_temp <- rep(1,p)
  lw_temp[c(indtruecor[j],indfalsecor[j])] <- 3
  mycolors[indtruecor[j]] <- "green"
  mycolors[indfalsecor[j]] <- "red"
  myplot.glmnet(sqfit1, xvar="lambda", label=FALSE, mycol=mycolors,
                lwd=lw_temp ,xlab="",ylab="", ylim=c(-1.35,1.35))
  true <- beta[j]
  abline(h=true,lwd=2)
}
dev.off()


#Scenario 1 plot ROC-curves (FPR vs TPR) ----
par(mar=c(4,4,1,1),mfrow=c(1,1))
plot(ROC1[,1],ROC1[,2],type="l",lwd=2,col=1,xlab="FPR", ylab="TPR")
points(ROC2[,1],ROC2[,2],type="l",lwd=2,col=2)
points(ROC3[,1],ROC3[,2],type="l",lwd=2,col=3)
points(ROC4[,1],ROC4[,2],type="l",lwd=2,col=4, lty=2)
legend(x=0.65,y=0.2, legend=c("strong", "moderate", "non-informative","glmnet"),col=1:4,lty=c(1,1,1,2),lwd=2 )

#Scenario 2----
set.seed(12345)

#Set simulation variables----
n<-100 #sample size of each training data set
p<-400 #number of covariates
pstar <- 24 #number of covariates in the model
p0 <- p - pstar
s2 <- 1
beta <- c(rep(c(-1,1),pstar/4),rep(c(-1/3,1/3),pstar/4), rep(0,p0))


#Simulate data and co-data####
#Generate X 
indtruestrong <- 1:(pstar/2)
indtrueweak <- (pstar/2 + 1):pstar
indfalse <- (pstar+1):p
Xmat <- matrix(rnorm(n*p),nrow=n)


#Generate response 
y <- Xmat %*% matrix(beta,ncol=1) + rnorm(n)

save(y,Xmat,n,p,pstar,beta,file=paste(pathResults,"simulp400n100_sc2.Rdata",sep=""))


#Co-data indices 
Gtrue <- 1:pstar
Gnottrue <- pstar + 1:pstar
#Strong co-data
Co1 <- c(Gtrue,Gnottrue)
Nonco1 <- (1:p)[-Co1]

#Weaker co-data 
Co2 <- sort(c(sample(1:(pstar/2), pstar/4),sample((pstar/2+1):pstar, pstar/4),Gnottrue))
Nonco2 <- (1:p)[-Co2]

#Random co-data 
Co3 <- sort(c(sample(1:pstar, pstar/16),Gnottrue))
Nonco3 <- (1:p)[-Co3]


#Fit models ########
#glmnet
gnfit_2 <- glmnet(Xmat,y,family="gaussian",intercept=FALSE, alpha=1)

#squeezy with strong co-data
grp <- list(Co1,Nonco1)
squeezyfit <- squeezy(y,Xmat,groupset=grp,intrcpt=FALSE,reCV=FALSE, alpha=1)
pfgr1 <- squeezyfit$lambdaApprox
pfgr1
pf <- rep(1,p)
pf[Co1] <- pfgr1[1];pf[Nonco1]<-pfgr1[2]
sqfit1_2 <- glmnet(Xmat,y,family="gaussian",intercept=FALSE,penalty.factor = pf, alpha=1)

#squeezy with moderate co-data
grp <- list(Co2,Nonco2)
squeezyfit <- squeezy(y,Xmat,groupset=grp,intrcpt=FALSE,reCV=FALSE, alpha=1)
pfgr2 <- squeezyfit$lambdaApprox
pfgr2
pf <- rep(1,p)
pf[Co2] <- pfgr2[1];pf[Nonco2]<-pfgr2[2]
sqfit2_2 <- glmnet(Xmat,y,family="gaussian",intercept=FALSE,penalty.factor = pf, alpha=1)

#squeezy with non-informative co-data
grp <- list(Co3,Nonco3)
squeezyfit <- squeezy(y,Xmat,groupset=grp,intrcpt=FALSE,reCV=FALSE, alpha=1)
pfgr3 <- squeezyfit$lambdaApprox
pfgr3
pf <- rep(1,p)
pf[Co3] <- pfgr3[1];pf[Nonco3]<-pfgr3[2]
sqfit3_2 <- glmnet(Xmat,y,family="gaussian",intercept=FALSE,penalty.factor = pf, alpha=1)

mycolors_2 <- rep("red",p)
mycolors_2[1:(pstar/2)] <- "lightgreen"
mycolors_2[(pstar/2+1):pstar] <- "darkgreen"
myltype_2 <- rep(2,p); myltype_2[mycolors_2=="red"] <- 3; myltype_2[mycolors_2=="lightgreen"] <- 1
lw_2 <- lt_2 <- rep(1,p)
lw_2[1:pstar] <- 2
lt_2[1:pstar] <- 1
save(sqfit1_2,sqfit2_2,sqfit3_2,gnfit_2,mycolors_2,lw_2,myltype_2,
     file=paste(pathResults,"inputplot_glmnet_sc2.Rdata",sep=""))


#Obtain roc curve information----
#load("inputplot_glmnet_sc2.Rdata") 

#coefficient matrices; rows are features, columns are penalties
coefs1 <- coefficients(sqfit1_2)[-1,] #remove intercept
coefs2 <- coefficients(sqfit2_2)[-1,] #remove intercept
coefs3 <- coefficients(sqfit3_2)[-1,] #remove intercept
coefs4 <- coefficients(gnfit_2)[-1,] #remove intercept

#computes FPR and TPR
TPRf <- function(vec,nnonzero){
  nbeta <- length(vec)
  betanz <- vec[1:nnonzero]
  betaz <- vec[-(1:nnonzero)]
  TPR <- length(which(abs(betanz) >= 10^(-10)))/nnonzero
  FPR <- length(which(abs(betaz) >= 10^(-10)))/(nbeta-nnonzero)
  return(c(FPR,TPR))
} 

#ROC
ROC1_2 <- rbind(t(apply(coefs1,2,TPRf, nnonzero=24)),c(1,1))
ROC2_2 <- rbind(t(apply(coefs2,2,TPRf, nnonzero=24)),c(1,1))
ROC3_2 <- rbind(t(apply(coefs3,2,TPRf, nnonzero=24)),c(1,1))
ROC4_2 <- rbind(t(apply(coefs4,2,TPRf, nnonzero=24)),c(1,1))

#Scenario 2 Load data plots----
load(paste(pathResults,"inputplot_glmnet_sc2.Rdata",sep=""))

#Scenario 2 plot traceplot----
#pdf(paste(pathFigures,"traceplots_sc2.pdf",sep=""),width=14,height=7)
par(mfrow=c(2,2),mar=c(3,3,2,1))
myplot.glmnet(sqfit1_2, xvar="lambda", label=FALSE, mycol=mycolors_2,
              lwd=lw_2, ylab="",xlab="",main="(a)")
#text(x=4.35,y=1,"(A)")
myplot.glmnet(sqfit2_2, xvar="lambda", label=FALSE, mycol=mycolors_2,
              lwd=lw_2,ylab="",xlab="",main="(b)")
#text(x=2,y=1.3,"(B)")
myplot.glmnet(sqfit3_2, xvar="lambda", label=FALSE, mycol=mycolors_2,
              lwd=lw_2,ylab="",xlab="",main="(c)")
#text(x=1.1,y=0.88,"(C)")
myplot.glmnet(gnfit_2, xvar="lambda", label=FALSE, mycol=mycolors_2,
              lwd=lw_2,ylab="",xlab="",main="(d)")
#text(x=0.7,y=0.92,"(D)")
#dev.off()


#Scenario 2 ROC-curves (FPR vs TPR) ----
#ROC curves
par(mar=c(4,4,1,1), mfrow=c(1,1))
plot(ROC1_2[,1],ROC1_2[,2],type="l",lwd=2,col=1,xlab="FPR", ylab="TPR")
points(ROC2_2[,1],ROC2_2[,2],type="l",lwd=2,col=2)
points(ROC3_2[,1],ROC3_2[,2],type="l",lwd=2,col=3)
points(ROC4_2[,1],ROC4_2[,2],type="l",lwd=2,col=4, lty=2)
legend(x=0.65,y=0.2, legend=c("strong", "moderate", "non-informative","glmnet"),col=1:4,lty=c(1,1,1,2),lwd=2 )

#Scenario 1&2 roc + scenario 2 traceplots----
figname<-paste(pathFigures,"FigSimsVarSel.pdf",sep="")
pdf(width = wdthpdf, height = hghtpdf,
    file = figname)
layout(matrix(c(1,2,3,5,4,6),ncol=3,byrow=F))
par(mar=c(4, 4, 3, 2) + 0.1)
#ROC curve scenario 1
plot(ROC1[,1],ROC1[,2],type="l",lwd=2,col=1,xlab="FPR", ylab="TPR", main="Scenario 1",lty=1)
points(ROC2[,1],ROC2[,2],type="l",lwd=2,col=2, lty=2)
points(ROC3[,1],ROC3[,2],type="l",lwd=2,col=3, lty=3)
points(ROC4[,1],ROC4[,2],type="l",lwd=2,col=4, lty=4)
legend("bottomright", title= "Co-data", legend=c("strong", "moderate", "non-informative","none (glmnet)"),
       col=1:4,lty=c(1,2,3,4),lwd=2, seg.len=3,cex=1)
#ROC curve scenario 2
plot(ROC1_2[,1],ROC1_2[,2],type="l",lwd=2,col=1,xlab="FPR", ylab="TPR", main="Scenario 2", lty=1)
points(ROC2_2[,1],ROC2_2[,2],type="l",lwd=2,col=2, lty=2)
points(ROC3_2[,1],ROC3_2[,2],type="l",lwd=2,col=3, lty=3)
points(ROC4_2[,1],ROC4_2[,2],type="l",lwd=2,col=4, lty=4)
legend("bottomright", title= "Co-data", legend=c("strong", "moderate", "non-informative","none (glmnet)"),
       col=1:4,lty=c(1,2,3,4),lwd=2, seg.len=3, cex=1)

myplot.glmnet(sqfit1_2, xvar="lambda", label=FALSE, mycol=mycolors_2, mylty=myltype_2,
              lwd=2,ylab="Regression coefficient estimates",xlab=expression(log(lambda)),
              main="Scenario 2, case A")
legend("bottomright", title="Effect",legend=c("strong", "weak", "none"),
       col=unique(mycolors_2),lty=c(1,2,3),lwd=2, bg="white", cex=1)
myplot.glmnet(sqfit2_2, xvar="lambda", label=FALSE, mycol=mycolors_2, mylty=myltype_2,
              lwd=2,ylab="Regression coefficient estimates",xlab=expression(log(lambda)),main="Scenario 2, case B")
legend("bottomright", title="Effect",legend=c("strong", "weak", "none"),
       col=unique(mycolors_2),lty=c(1,2,3),lwd=2, bg="white", cex=1)
myplot.glmnet(sqfit3_2, xvar="lambda", label=FALSE, mycol=mycolors_2, mylty=myltype_2,
              lwd=2,ylab="Regression coefficient estimates",xlab=expression(log(lambda)),main="Scenario 2, case C")
legend("bottomright", title="Effect",legend=c("strong", "weak", "none"),
       col=unique(mycolors_2),lty=c(1,2,3),lwd=2, bg="white", cex=1)
myplot.glmnet(gnfit_2, xvar="lambda", label=FALSE, mycol=mycolors_2, mylty=myltype_2,
              lwd=2,ylab="Regression coefficient estimates",xlab=expression(log(lambda)),main="Scenario 2, case D")
legend("bottomright", title="Effect",legend=c("strong", "weak", "none"),
       col=unique(mycolors_2),lty=c(1,2,3),lwd=2, bg="white", cex=1)
dev.off()
