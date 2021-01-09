###Simulation study
#NOTE 1: in Rstudio; press ALT+o for overview of code section headings

#Load libraries----
setwd("C:\\Users\\VNOB-0728\\Documents\\Server\\Analysis squeezy")
pathResults <- "./SimulationResults/" 
pathFigures <- "./SimulationFigures/" 
#setwd("C:\\Synchr\\Rscripts\\EB\\Squeezy")
library(ecpc)
library(MASS)
library(penalized)
library(glmnet)
library(mvtnorm)
library(mgcv)
library(CVXR)
library(foreach)
library(doParallel)

library(gren)
library(ipflasso)
library(fwelnet)

library(multiridge)
source("multiridge.R")
source("squeezy.R")

#library(Matrix)
#library(GRridge)

#load libraries needed for storing and plotting results
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(tidyr)

#optional: setup parallel backend to use many processors
# cores=detectCores()
if(0){
  cl <- makeCluster(7) #not to overload your computer
  registerDoParallel(cl)
}

#generate seeds for reproducibility
#set.seed(1)
seeds <- round(runif(10000)*10^6)
#save(seeds,file=paste(pathResults,"seeds.R",sep=""))
load(paste(pathResults,"seeds.R",sep=""))

#Set simulation variables----
nSim <- 100
n<-150
p<-1200
n2<-300 #sample size of each test data set
q<-0 #set q*100% of the betas to 0
alphaEN <- 0.3

#Generate correlated data----
#set.seed(101) #used for initial values comparison
Nblock <- 10     #correlation blocks in X
CorX <- 0.2
setting <- "3"
for(setting in as.character(1)){
  switch(setting,
         "1"={ #setting 1: no signal (single lambda setting)
           taugrp <- rep(0.1,5)  #single lambda setting
           sigmasq <- 2
           n<-150
           p<-1200
           nSim <- 100
           G<-5 #number of groups for which co-data is generated
           
           set.seed(seeds[1])
         },
         "2"={ #setting 2: medium signal (multigroup setting with small differences)
           taugrp <- c(0.01,0.05,0.1,0.4,0.8)
           taugrp <- mean(taugrp) + (taugrp-mean(taugrp))/2
           sigmasq <- 2
           n<-150
           p<-1200
           nSim <- 100
           G<-5 #number of groups for which co-data is generated
           
           set.seed(seeds[2])
         },
         "3"={ #setting 3: high signal (multigroup setting with large differences)
           taugrp <- c(0.01,0.05,0.1,0.4,0.8)
           sigmasq <- 2
           n<-150
           p<-1200
           nSim <- 100
           G<-5 #number of groups for which co-data is generated
           
           set.seed(seeds[3])
         },
         "4"={ #setting 4: time dependence n, fixed p, G, small n
           taugrp <- c(0.01,0.05,0.1,0.4,0.8)
           sigmasq <- 2
           n<-75
           p<-1200
           nSim <- 5
           G<-5 #number of groups for which co-data is generated
           
           set.seed(seeds[3])
         },
         "5"={ #setting 5: time dependence n, fixed p, G, large n
           taugrp <- c(0.01,0.05,0.1,0.4,0.8)
           sigmasq <- 2
           n<-300
           p<-1200
           nSim <- 5
           G<-5 #number of groups for which co-data is generated
           
           set.seed(seeds[3])
         },
         "6"={ #setting 6: time dependence p, fixed n, G, small p
           taugrp <- c(0.01,0.05,0.1,0.4,0.8)*2
           sigmasq <- 2
           n<-150
           p<-600
           nSim <- 5
           G<-5 #number of groups for which co-data is generated
           
           set.seed(seeds[3])
         },
         "7"={ #setting 7: time dependence p, fixed n, G, large p
           taugrp <- c(0.01,0.05,0.1,0.4,0.8)/2
           sigmasq <- 2
           n<-150
           p<-2400
           nSim <- 5
           G<-5 #number of groups for which co-data is generated
           
           set.seed(seeds[3])
         },
         "8"={ #setting 8: time dependence G, fixed n, p, small G
           taugrp <- c(0.01,0.05,0.1,0.4,0.8)
           sigmasq <- 2
           n<-150
           p<-1200
           nSim <- 5
           G<-2 #number of groups for which co-data is generated
           
           set.seed(seeds[3])
         },
         "9"={ #setting 9: time dependence G, fixed n, p, large G
           taugrp <- c(0.01,0.05,0.1,0.4,0.8)
           sigmasq <- 2
           n<-150
           p<-1200
           nSim <- 5
           G<-10 #number of groups for which co-data is generated
           
           set.seed(seeds[3])
         },
         "10"={ #setting 10: time dependence central n,p,G
           taugrp <- c(0.01,0.05,0.1,0.4,0.8)
           sigmasq <- 2
           n<-150
           p<-1200
           nSim <- 5
           G<-5 #number of groups for which co-data is generated
           
           set.seed(seeds[3])
         },
         "11"={ #setting 7: time dependence p, fixed n, G, largest p
           taugrp <- c(0.01,0.05,0.1,0.4,0.8)/4
           sigmasq <- 2
           n<-150
           p<-4800
           nSim <- 5
           G<-5 #number of groups for which co-data is generated
           
           set.seed(seeds[3])
         }
  )
  
  #Make random grouping (used below to make informative grouping) ----
  #set.seed(74743)
  rankRandom <- sample(1:p,p,replace = F) #ranking used for random groups
  ind <- list() #first make list with for each number of groups, a px1 vector with group number
  size <- floor(p/G) #group size
  rest <- p- size*G #first rest number of groups obtain 1 covariate more than size
  if(rest>0){
    ind<-c(rep(1:rest,each=(size+1)),rep((rest+1):G,each=size))
  }else{
    ind<-rep(1:G,each=size)
  }
  #for each number of groups G[g], create list element with list of G[g] random groups
  #e.g. indRandom[[3]][[1]] is the vector with the covariate indices belonging to group 1 in the grouping with G[3] random groups
  #indRandom <- lapply(1:max(ind),function(x){rankRandom[which(ind==x)]})
  
  tauglobal2 <- rep(taugrp,each=p/G)
  if(G!=5){
    tauglobal2 <- rep(taugrp,each=p/5)
  }
  lambda1 <- sqrt(tauglobal2/2)
  dist <- "laplace"
  model <- "linear"
  
  betas <- list()
  X <- list()
  X2 <- list()
  Y <- list()
  Y2 <- list()
  rankBeta <- list(); indRank <- list()
  for(i in 1:nSim){
    #simulate betas
    if(dist=="gaussian"){
      betas[[i]] <- rnorm(p,0,sqrt(tauglobal2))
      rankBeta[[i]] <- order(abs(betas[[i]]))
    } 
    if(dist=="laplace"){
      lambda1 <- sqrt(tauglobal2/2)
      # print(paste("Simulating betas from laplace(0,",signif(lambda1,3),")",sep=""))
      #A Laplace(0,b) variate can also be generated as the difference of two i.i.d. 
      #Exponential(1/b) random variables
      betas[[i]] <-   rexp(p, 1/lambda1) -  rexp(p, 1/lambda1)
      rankBeta[[i]] <- order(abs(betas[[i]]))
      if(q > 0){
        set0 <- sample(1:p,size=floor(q*p),replace=F)
        betas[[i]][set0] <- 0
      }
    }
    #simulate co-data: list with G rank-based groups (using rank from betas before randomly setting some to 0)
    #indRank[[i]] <- lapply(1:max(ind),function(x){rankBeta[[i]][which(ind==x)]}) 
    indGrp <- rep(1:G,each=p/G)
    indRank[[i]] <- lapply(1:G,function(g) which(indGrp==g))
    
    #simulate observed data
    pblock <- p/Nblock
    X[[i]] <- Reduce(cbind,lapply(1:Nblock, function(z) matrix(rep(rnorm(n,sd=sqrt(CorX/(1-CorX))),times=pblock),n,pblock))) + 
      matrix(rnorm(n*p),n,p)
    X2[[i]] <- Reduce(cbind,lapply(1:Nblock, function(z) matrix(rep(rnorm(n2,sd=sqrt(CorX/(1-CorX))),times=pblock),n2,pblock))) + 
      matrix(rnorm(n2*p),n2,p)
    
    #X <- t((t(X) - apply(t(X),1,mean))/apply(t(X),1,sd))
    #Beta <- rnorm(G*p,mean=0,sd=sqrt(vbetas))
    #simulate response data
    if(model=="linear"){
      Y[[i]] <- as.numeric(X[[i]]%*%betas[[i]] + rnorm(n,sd=sqrt(sigmasq)))
      Y2[[i]] <- as.numeric(X2[[i]]%*%betas[[i]] + rnorm(n2,sd=sqrt(sigmasq)))
      
      fml <- "gaussian"
    }else if(model=="logistic"){
      probs <- 1/(1+exp(-as.numeric(X[[i]]%*%betas[[i]])))
      probs2 <- 1/(1+exp(-as.numeric(X2[[i]]%*%betas[[i]])))
      Y[[i]] <- sapply(probs,function(x) rbinom(n=1,size=1,prob=x))
      Y2[[i]] <- sapply(probs2,function(x) rbinom(n=1,size=1,prob=x))
      
      fml <- "binomial"
    }
  }
  
#Simulations all----
logname<-paste("logSimulations.txt",sep="")
write(paste(Sys.time(),logname),file=logname,append=T)

allMSEsim <- allctsim <- models <- c()
simpara <- list(n=n,ntest=n2,p=p,G=G,CorX=CorX,Nblock=Nblock,
                alphaEN=alphaEN, model=model, dist= dist,taugrp=taugrp,sigmasq=sigmasq)
#fname <- paste(pathResults,paste("simres",alphaEN,setting,sep="_"),".Rdata",sep="")
fname <- paste(pathResults,paste("simreshat",alphaEN,setting,sep="_"),".Rdata",sep="") #alphat

nSim2 <- nSim
#nSim2 <- 5

df <- data.frame()

#for(i in 1:nSim2){
finalMatrix <- foreach(i=1:nSim2, .combine=rbind,
                       .packages = c("glmnet","penalized","mvtnorm","gglasso","mgcv",
                                     "CVXR","GRridge","expm","Rsolnp","ecpc",
                                     "ipflasso","fwelnet","gren")) %dopar% {
  write(paste(Sys.time(),i,logname),file=logname,append=T)

  #Fit elastic net with glmnet----
  alphat <- alphaEN
  sd_y <- sqrt(var(Y[[i]])*(n-1)/n)[1]
  alphat <- 1/(1+sd_y*(1-alphaEN)/alphaEN)
  
  pmt <- proc.time()[[3]]
  glmGR <-  cv.glmnet(X[[i]],Y[[i]],alpha=alphat,family=fml, intercept = T, standardize = F, thresh=10^-10)
  lambda1approx <- glmGR$lambda.min
  betaEN <- coef(glmGR,s=lambda1approx,exact=T)[-1]
  a0EN <- coef(glmGR,s=lambda1approx)[1]
  MSEEN <- mean((Y2[[i]]-X2[[i]]%*%betaEN-a0EN)^2)
  ctglm <- proc.time()[[3]]-pmt

  MSEEN

  #Fit gren (logistic only)----
  #i <- 1 #simulation data set number
  # ind.vector <- rep(1,p)
  # for(j in 1:length(indRank[[i]])){
  #   ind.vector[indRank[[i]][[j]]] <- j
  # }
  # resgren <- NaN
  # if(model=="logistic"){
  #   resgren <- gren(x=X[[i]],y=Y[[i]],partitions=ind.vector,compare=F,alpha=alphaEN)
  #   betagren <- resgren$freq.model$groupreg$beta
  #   a0gren <- resgren$freq.model$groupreg$a0
  # }


  #Fit ipflasso----
  #i <- 1 #simulation data set number
  ind.vector <- rep(1,p)
  for(j in 1:length(indRank[[i]])){
   ind.vector[indRank[[i]][[j]]] <- j
  }
  #make list with proposal lambdas in the form of ratios 1:1:..:1, 0.5:..:1, etc.
  #weights vary from 2^-x to 2^x
  ratios <- 2^(0:1)
  A <- matrix(ratios,c(length(ratios),1))
  for(k in 2:G){
   l<-sapply(ratios,function(x){
     cbind(A,rep(x,length(ratios)))
   },simplify=F)
   A<-do.call(rbind,l)
  }
  pflist <- lapply(1:dim(A)[1],function(x)c(A[x,]))

  pmt <- proc.time()[[3]]
  if(model=="linear"){
   resipf <- cvr2.ipflasso(X=X[[i]],Y=Y[[i]],blocks=indRank[[i]],alpha=alphat,standardize = F,
                           family=fml,type.measure="mse",nfolds=10,ncv=1,pflist=pflist)
  }else{
   resipf <- cvr2.ipflasso(X=X[[i]],Y=Y[[i]],blocks=indRank[[i]],alpha=alphat,standardize = F,
                           family=fml,type.measure="auc",nfolds=10,ncv=1,pflist=pflist)
  }
  a0.ipf <- resipf$coeff[1,resipf$ind.bestlambda]
  beta.ipf <- c(resipf$coeff[-1,resipf$ind.bestlambda])
  predipf <- X2[[i]]%*%beta.ipf+a0.ipf
  MSEipf <- mean((Y2[[i]]-predipf)^2)
  ctipf <- proc.time()[[3]]-pmt

  # #ipflasso with more weights
  # if(i %in% 1:5){
  #   #make list with proposal lambdas in the form of ratios 1:1:..:1, 0.5:..:1, etc.
  #   #weights vary from 2^-x to 2^x
  #   ratios <- 2^(0:2)
  #   A <- matrix(ratios,c(length(ratios),1))
  #   for(k in 2:G){
  #     l<-sapply(ratios,function(x){
  #       cbind(A,rep(x,length(ratios)))
  #     },simplify=F)
  #     A<-do.call(rbind,l)
  #   }
  #   pflist <- lapply(1:dim(A)[1],function(x)c(A[x,]))
  #
  #   pmt <- proc.time()[[3]]
  #   if(model=="linear"){
  #     resipf <- cvr2.ipflasso(X=X[[i]],Y=Y[[i]],blocks=indRank[[i]],alpha=alphaEN,standardize = F,
  #                             family=fml,type.measure="mse",nfolds=10,ncv=1,pflist=pflist)
  #   }else{
  #     resipf <- cvr2.ipflasso(X=X[[i]],Y=Y[[i]],blocks=indRank[[i]],alpha=alphaEN,standardize = F,
  #                             family=fml,type.measure="auc",nfolds=10,ncv=1,pflist=pflist)
  #   }
  #   a0.ipf <- resipf$coeff[1,resipf$ind.bestlambda]
  #   beta.ipf <- c(resipf$coeff[-1,resipf$ind.bestlambda])
  #   predipf <- X2[[i]]%*%beta.ipf+a0.ipf
  #   MSEipf2 <- mean((Y2[[i]]-predipf)^2)
  #   ctipf2 <- proc.time()[[3]]-pmt
  # }


  #Fit fwelnet----
  #i <- 1 #simulation data set number
  Z <- matrix(rep(0,p*G),c(p,G))
  for(j in 1:G){
   Z[indRank[[i]][[j]],j] <- 1
  }

  pmt <- proc.time()[[3]]
  #inclusion of intercept not optional and included by default
  fit.fwEN <-  cv.fwelnet(x=X[[i]],y=Y[[i]],z=Z,
                         alpha=alphat,family=fml,
                         standardize = F, thresh=10^-10)
  ind.minlam <- which(fit.fwEN$lambda==fit.fwEN$lambda.min)

  betafwEN <- fit.fwEN$glmfit$beta[,ind.minlam]
  a0fwEN <- fit.fwEN$glmfit$a0[ind.minlam]
  MSEfwEN <- mean((Y2[[i]]-X2[[i]]%*%betafwEN-a0fwEN)^2)
  ctfwEN <- proc.time()[[3]]-pmt
  MSEfwEN

  #summarise results----

  df2 <- data.frame("time"=c(ctglm, ctipf,ctfwEN),
                   "MSE"=c(MSEEN, MSEipf,MSEfwEN),
                   "Method"=c("EN", "ipf","fwEN"))
  df2$Dataset <- i
  df2$n <- n
  df2$G <- G
  df2$p <- p
  #df <- rbind(df,df2)
  # if(i %in% 1:5){
  #   df3 <- data.frame("time"= ctipf2,
  #                     "MSE"= MSEipf2,
  #                     "Method"= "ipf2",
  #                     "Dataset"=i,
  #                     "n"=n,
  #                     "G"=G,
  #                     "p"=p)
  #   df2 <- rbind(df2,df3)
  # }

  # allct <- c(i, ctecpc, ctsqueezy, ctglm, ctipf,ctfwEN)
  # print(allct)
  # allctsim <- rbind(allctsim,allct)
  # allMSE <- c(i, mses_ecpcsq, MSEEN, MSEipf,MSEfwEN)
  # allMSEsim <- rbind(allMSEsim,allMSE)
  # colnames(allMSEsim) <- c("i", "ecpc", "sqAIC_MML_b", "sqAIC_MML_m", "sqAIC_MML_s", "EN", "ipf","fwEN")
  # colnames(allctsim) <- c("i", "ecpc", "sqAIC_MML", "EN", "ipf","fwEN")
  # print(allMSEsim)
  # print(models)


  list("df"=df2)
  }
df2 <- lapply(1:nSim2,function(i) finalMatrix[i,1][[1]])
df <- df2[[1]]; for(i in 2:nSim2) df <- rbind(df,df2[[i]])

save(df,simpara,file=fname)
}

#Simulations squeezy----
logname<-paste("logSimulations.txt",sep="")
write(paste(Sys.time(),logname),file=logname,append=T)

allMSEsim <- allctsim <- models <- c()
simpara <- list(n=n,ntest=n2,p=p,G=G,CorX=CorX,Nblock=Nblock,
                alphaEN=alphaEN, model=model, dist= dist,taugrp=taugrp,sigmasq=sigmasq)
fname <- paste(pathResults,paste("simresSqueezy",alphaEN,setting,sep="_"),".Rdata",sep="")

nSim2 <- nSim
#nSim2 <- 5

df <- data.frame()

#for(i in 1:nSim2){
finalMatrix <- foreach(i=1:nSim2, .combine=rbind,
                       .packages = c("glmnet","penalized","mvtnorm","gglasso","mgcv",
                                     "CVXR","GRridge","expm","Rsolnp","ecpc",
                                     "ipflasso","fwelnet","gren")) %dopar% {
       write(paste(Sys.time(),i,logname),file=logname,append=T)
       
       #Fit ecpc+squeezy defaults----
       #i <- 1 #simulation data set number
       sd_y <- sqrt(var(Y[[i]])*(n-1)/n)[1]
       #indRank[[i]] <- list(1:p)
       #indRank[[i]] <- indRandom
       
       #first fit ecpc to find group-ridge penalty estimates
       pmt <- proc.time()[[3]]
       res.ecpc<-ecpc(Y[[i]],X[[i]], #observed data and response to train model
                      groupings=list(indRank[[i]]), #informative grouping
                      Y2=Y2[[i]],X2=X2[[i]], #test set
                      model=model,
                      hypershrinkage="none",postselection = F)
       #then use multivariate normal approximation to fit model for group-regularised elastic net
       ctecpc <- proc.time()[[3]] - pmt
       ctecpc
       
       #multigroup, no reCV
       pmt <- proc.time()[[3]]
       res.ecpcsqueezy <- squeezy(Y[[i]],X[[i]],grouping=indRank[[i]],Y2=Y2[[i]],X2=X2[[i]],
                                  model=model,intrcpt=T,
                                  alpha=alphaEN, #elastic net parameter for which group-regularised EN has to be fit
                                  method="MML", #one of c("ecpcEN","MML","MML.noDeriv","CV")
                                  fit.ecpc=res.ecpc) #fit of ecpc function
       ctsqueezy <- proc.time()[[3]]-pmt + ctecpc;
       ctsqueezy
       
       #multigroup, with reCV
       pmt <- proc.time()[[3]]
       res.ecpcsqueezy2 <- squeezy(Y[[i]],X[[i]],grouping=indRank[[i]],Y2=Y2[[i]],X2=X2[[i]],
                                  model=model,intrcpt=T,
                                  alpha=alphaEN, #elastic net parameter for which group-regularised EN has to be fit
                                  method="MML", #one of c("ecpcEN","MML","MML.noDeriv","CV")
                                  fit.ecpc=res.ecpc,reCV=T) #fit of ecpc function
       ctsqueezy2 <- proc.time()[[3]]-pmt + ctecpc;
       ctsqueezy2
       
       #first fit ecpc to find group-ridge penalty estimates
       pmt <- proc.time()[[3]]
       res.ecpc2<-ecpc(Y[[i]],X[[i]], #observed data and response to train model
                      groupings=list(list(1:p)), #informative grouping
                      Y2=Y2[[i]],X2=X2[[i]], #test set
                      model=model,
                      hypershrinkage="none",postselection = F)
       #then use multivariate normal approximation to fit model for group-regularised elastic net
       ctecpc2 <- proc.time()[[3]] - pmt
       ctecpc2
       
       #one group, no reCV
       pmt <- proc.time()[[3]]
       res.ecpcsqueezy3 <- squeezy(Y[[i]],X[[i]],grouping=list(1:p),Y2=Y2[[i]],X2=X2[[i]],
                                  model=model,intrcpt=T,
                                  alpha=alphaEN, #elastic net parameter for which group-regularised EN has to be fit
                                  method="MML", #one of c("ecpcEN","MML","MML.noDeriv","CV")
                                  fit.ecpc=res.ecpc2) #fit of ecpc function
       ctsqueezy3 <- proc.time()[[3]]-pmt +ctecpc2;
       ctsqueezy3
       
       #one group, with reCV
       pmt <- proc.time()[[3]]
       res.ecpcsqueezy4 <- squeezy(Y[[i]],X[[i]],grouping=list(1:p),Y2=Y2[[i]],X2=X2[[i]],
                                   model=model,intrcpt=T,
                                   alpha=alphaEN, #elastic net parameter for which group-regularised EN has to be fit
                                   method="MML", #one of c("ecpcEN","MML","MML.noDeriv","CV")
                                   fit.ecpc=res.ecpc2,reCV=T) #fit of ecpc function
       ctsqueezy4 <- proc.time()[[3]]-pmt +ctecpc2;
       ctsqueezy4
       
       #default ecpcEN
       pmt <- proc.time()[[3]]
       res.ecpcsqueezy5 <- squeezy(Y[[i]],X[[i]],grouping=indRank[[i]],Y2=Y2[[i]],X2=X2[[i]],
                                   model=model,intrcpt=T,
                                   alpha=alphaEN, #elastic net parameter for which group-regularised EN has to be fit
                                   method="ecpcEN", #one of c("ecpcEN","MML","MML.noDeriv","CV")
                                   fit.ecpc=res.ecpc) #fit of ecpc function
       ctsqueezy5 <- proc.time()[[3]]-pmt +ctecpc;
       ctsqueezy5
       
       MSEmulti <- res.ecpcsqueezy$MSEApprox #e.g. MSE on independent test set of EN-approximation
       MSEmultireCV <- res.ecpcsqueezy2$MSEApprox #e.g. MSE on independent test set of EN-approximation
       MSEsingle <- res.ecpcsqueezy3$MSEApprox #e.g. MSE on independent test set of EN-approximation
       MSEsinglereCV <- res.ecpcsqueezy4$MSEApprox #e.g. MSE on independent test set of EN-approximation
       MSEecpcEN <- res.ecpcsqueezy5$MSEApprox #e.g. MSE on independent test set of EN-approximation
       
       mses_ecpcsq <- c(MSEmulti,MSEmultireCV, MSEsingle, MSEsinglereCV,MSEecpcEN)
       
       #summarise results----
       
       df2 <- data.frame("time"=c(ctsqueezy, ctsqueezy2,ctsqueezy3,ctsqueezy4,ctsqueezy5),
                         "MSE"=mses_ecpcsq,
                         "Method"=c("squeezy.m","squeezy.m.reCV","squeezy.s","squeezy.s.reCV","ecpcEN"))
       df2$Dataset <- i
       df2$n <- n
       df2$G <- G
       df2$p <- p
      
       
       list("df"=df2)
     }
df2 <- lapply(1:nSim2,function(i) finalMatrix[i,1][[1]])
df <- df2[[1]]; for(i in 2:nSim2) df <- rbind(df,df2[[i]])

save(df,simpara,file=fname)
}

#Simulations ipflasso larger grid----
fname <- paste(pathResults,paste("simresIPF2",alphaEN,setting,sep="_"),".Rdata",sep="")

#nSim2 <- nSim
nSim2 <- 5

df <- data.frame()

#for(i in 1:nSim2){
finalMatrix <- foreach(i=1:nSim2, .combine=rbind,
                       .packages = c("glmnet","penalized","mvtnorm","gglasso","mgcv",
                                     "CVXR","GRridge","expm","Rsolnp","ecpc",
                                     "ipflasso","fwelnet","gren")) %dopar% {
             write(paste(Sys.time(),i,logname),file=logname,append=T)
             
             
             #ipflasso with more weights
             if(i %in% 1:5){
               #make list with proposal lambdas in the form of ratios 1:1:..:1, 0.5:..:1, etc.
               #weights vary from 2^-x to 2^x
               ratios <- 2^(0:2)
               A <- matrix(ratios,c(length(ratios),1))
               for(k in 2:G){
                 l<-sapply(ratios,function(x){
                   cbind(A,rep(x,length(ratios)))
                 },simplify=F)
                 A<-do.call(rbind,l)
               }
               pflist <- lapply(1:dim(A)[1],function(x)c(A[x,]))

               pmt <- proc.time()[[3]]
               if(model=="linear"){
                 resipf <- cvr2.ipflasso(X=X[[i]],Y=Y[[i]],blocks=indRank[[i]],alpha=alphaEN,standardize = F,
                                         family=fml,type.measure="mse",nfolds=10,ncv=1,pflist=pflist)
               }else{
                 resipf <- cvr2.ipflasso(X=X[[i]],Y=Y[[i]],blocks=indRank[[i]],alpha=alphaEN,standardize = F,
                                         family=fml,type.measure="auc",nfolds=10,ncv=1,pflist=pflist)
               }
               a0.ipf <- resipf$coeff[1,resipf$ind.bestlambda]
               beta.ipf <- c(resipf$coeff[-1,resipf$ind.bestlambda])
               predipf <- X2[[i]]%*%beta.ipf+a0.ipf
               MSEipf2 <- mean((Y2[[i]]-predipf)^2)
               ctipf2 <- proc.time()[[3]]-pmt
             }
             
             
            
             if(i %in% 1:5){
               df2 <- data.frame("time"= ctipf2,
                                 "MSE"= MSEipf2,
                                 "Method"= "ipf2",
                                 "Dataset"=i,
                                 "n"=n,
                                 "G"=G,
                                 "p"=p)
             }
                                       
            write(paste(Sys.time(),setting,i,"ipflasso2 done"),file=logname)
            
             list("df"=df2)
  }

  if(nSim2==1){
    df <- finalMatrix[[1]]
  }else{
    df2 <- lapply(1:nSim2,function(i) finalMatrix[i,1][[1]])
    df <- df2[[1]]; for(i in 2:nSim2) df <- rbind(df,df2[[i]])
  }

save(df,simpara,file=fname)
}

#Simulations multiple alpha----
logname<-paste("logSimulations.txt",sep="")
write(paste(Sys.time(),logname),file=logname,append=T)

reCV <- T
allMSEsim <- allctsim <- models <- c()
simpara <- list(n=n,ntest=n2,p=p,G=G,CorX=CorX,Nblock=Nblock,
                alphaEN=alphaEN, model=model, dist= dist,taugrp=taugrp,sigmasq=sigmasq)
fname <- paste(pathResults,paste("simressqueezytemp",alphaEN,setting,sep="_"),".Rdata",sep="")
load(fname)

nSim2 <- nSim
#nSim2 <- 5

df <- data.frame()
res.ecpcsqueezy <- list() 
#for(i in 1:nSim2){
finalMatrix <- foreach(i=1:nSim2, .combine=rbind,
                       .packages = c("glmnet","penalized","mvtnorm","gglasso","mgcv",
                                     "CVXR","GRridge","expm","Rsolnp","ecpc",
                                     "ipflasso","fwelnet","gren")) %dopar% {
  write(paste(Sys.time(),i,logname),file=logname,append=T)
  #Fit ecpc+squeezy defaults----
  #i <- 1 #simulation data set number
  sd_y <- sqrt(var(Y[[i]])*(n-1)/n)[1]
  #indRank[[i]] <- list(1:p)
  indGrp <- rep(1:5,each=p/G)
  indRank[[i]] <- lapply(1:G,function(g) which(indGrp==g))
  #indRank[[i]] <- indRandom

  #first fit ecpc to find group-ridge penalty estimates
  pmt <- proc.time()[[3]]
  res.ecpc<-ecpc(Y[[i]],X[[i]], #observed data and response to train model
                groupings=list(indRank[[i]]), #informative grouping
                Y2=Y2[[i]],X2=X2[[i]], #test set
                model=model,
                hypershrinkage="none",postselection = F)
  #then use multivariate normal approximation to fit model for group-regularised elastic net
  ctecpc <- proc.time()[[3]] - pmt
  ctecpc
  pmt <- proc.time()[[3]]
  res.ecpcsqueezy <- squeezy(Y[[i]],X[[i]],grouping=indRank[[i]],Y2=Y2[[i]],X2=X2[[i]],
                            model=model,intrcpt=T,
                            alpha=alphaEN, #elastic net parameter for which group-regularised EN has to be fit
                            method="MML", #one of c("ecpcEN","MML","MML.noDeriv","CV")
                            fit.ecpc=res.ecpc) 
  res.ecpcsqueezy$Dataset <- i
  ctsqueezy <- proc.time()[[3]]-pmt;
  ctsqueezy
 

  list("res.ecpcsqueezy"=res.ecpcsqueezy)
  }
res.ecpcsqueezy <- lapply(1:nSim2,function(i) finalMatrix[[i]])

#save(res.ecpcsqueezy,simpara,file=fname)

load(fname)

alp <- seq(0,1,length.out=100)
df <- data.frame()
dfGrps <- data.frame()
#for(i in 1:nSim2){
finalMatrix <- foreach(i=1:nSim2, .combine=rbind,
                       .packages = c("glmnet","penalized","mvtnorm","gglasso","mgcv",
                                     "CVXR","GRridge","expm","Rsolnp","ecpc",
                                     "ipflasso","fwelnet","gren")) %dopar% {
     #write(paste(Sys.time(),i,logname),file=logname,append=T)
       #print(i)
       df <- data.frame()
       dfGrps <- data.frame()
       for(j in 1:length(alp)){
         #print(c(i,j))
         #without reCV
         pmt <- proc.time()[[3]]
         res.ecpcsqueezyAlp <- squeezy(Y[[i]],X[[i]],grouping=indRank[[i]],Y2=Y2[[i]],X2=X2[[i]],
                                       model=model,intrcpt=T,
                                       alpha=alp[j], #elastic net parameter for which group-regularised EN has to be fit
                                       method="MML", #one of c("ecpcEN","MML","MML.noDeriv","CV")
                                       lambdas=res.ecpcsqueezy[[i]]$lambdaMR,
                                       lambdaglobal=res.ecpcsqueezy[[i]]$lambdaglobal,
                                       sigmasq = res.ecpcsqueezy[[i]]$sigmahat,
                                       reCV=F) #T/F to return results of squeezy fit of both models
         res.ecpcsqueezy[[i]]$Dataset <- i
         ctsqueezy <- proc.time()[[3]]-pmt;
         ctsqueezy
         
         #with reCV
         pmt <- proc.time()[[3]]
         res.ecpcsqueezyAlp2 <- squeezy(Y[[i]],X[[i]],grouping=indRank[[i]],Y2=Y2[[i]],X2=X2[[i]],
                                       model=model,intrcpt=T,
                                       alpha=alp[j], #elastic net parameter for which group-regularised EN has to be fit
                                       method="MML", #one of c("ecpcEN","MML","MML.noDeriv","CV")
                                       lambdas=res.ecpcsqueezy[[i]]$lambdaMR,
                                       lambdaglobal=res.ecpcsqueezy[[i]]$lambdaglobal,
                                       sigmasq = res.ecpcsqueezy[[i]]$sigmahat,
                                       reCV=T) #T/F to return results of squeezy fit of both models
         
         ctsqueezy <- proc.time()[[3]]-pmt;
         ctsqueezy
         
         df2 <- data.frame(Alpha=rep(alp[j],2),
                          Dataset=rep(i,2),
                          MSEApprox=c(res.ecpcsqueezyAlp$MSEApprox,res.ecpcsqueezyAlp2$MSEApprox),
                          reCV =c(F,T))
         df <- rbind(df,df2)
         
         dfGrps2 <- data.frame(lambdaApprox = res.ecpcsqueezyAlp$lambdaApprox,
                           Group = 1:G)
         dfGrps2$Alpha <- alp[j]
         dfGrps2$Dataset <- i
         dfGrps <- rbind(dfGrps,dfGrps2)
       }
       
       
       
       list("df"=df,"dfGrps"=dfGrps)
}
df2 <- lapply(1:nSim2,function(i) finalMatrix[i,1][[1]])
df <- df2[[1]]; for(i in 2:nSim2) df <- rbind(df,df2[[i]])
dfGrps2 <- lapply(1:nSim2,function(i) finalMatrix[i,2][[1]])
dfGrps <- dfGrps2[[1]]; for(i in 2:nSim2) dfGrps <- rbind(dfGrps,dfGrps2[[i]])

fname <- paste(pathResults,paste("simressqueezyalpha",setting,sep="_"),".Rdata",sep="")
save(df,dfGrps,file=fname)

#intercept only model
dfIntercept <- data.frame()
for(i in 1:nSim2){
  a0 <- mean(Y[[i]])
  MSEApprox <- mean((Y2[[i]]-a0)^2)
  dfIntercept <- rbind(dfIntercept,
                       data.frame(Alpha=NaN,Dataset=i,MSEApprox=MSEApprox))
}
save(df,dfGrps,dfIntercept,file=fname)

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
#display.brewer.all(3,colorblindFriendly=T)
colpal <- "Dark2"
colsfill <- brewer.pal(3,"Dark2")[1:2]
colsAUC <- brewer.pal(8,"Dark2")

#load data----
fname <- paste(pathResults,paste("simres",alphaEN,p,n,taugrp[1],taugrp[5],sigmasq,sep="_"),".Rdata",sep="")
load(fname)

#plot performance
df2 <- data.frame()
for(setting in c("1","2","3")){
  fname <- paste(pathResults,paste("simres",alphaEN,setting,sep="_"),".Rdata",sep="")
  load(fname)
  df$Setting <- setting
  df2 <- rbind(df2,df)
}
df2$Setting2 <- as.numeric(df2$Setting)
df2$Setting <- factor(df2$Setting,levels=c("1","2","3"),
                      labels=c("No groups","Weakly informative groups","Highly informative groups"))

summBox <- df2 %>% group_by(Setting2) %>% summarise("low"=boxplot.stats(MSE)$stats[1],
                                                    "high"=boxplot.stats(MSE)$stats[5]) %>% ungroup()

df2$low <- summBox$low[df2$Setting2]
df2$high <- summBox$high[df2$Setting2]
df2NoOutliers <- df2[df2$MSE>=df2$low & df2$MSE<=df2$high,]


#Plot: prediction all with outliers----
# figname<-paste(pathFigures,"FigSimsPredictionAllOutliers.pdf",sep="")
# pdf(width = wdthpdf*1.5, height = hghtpdf*1.5,
#     file = figname)
ggplot(df2)+aes(x=Method,y=MSE)+
  geom_boxplot(aes(fill=Method))+
  facet_wrap(.~Setting,scales="free")+
  labs(y="MSE",x="Method")+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
# dev.off()

#Plot: prediction performance all no outliers----
figname<-paste(pathFigures,"FigSimsPredictionAllNoOutliers.pdf",sep="")
pdf(width = wdthpdf*1.7, height = hghtpdf,
    file = figname)
ggplot(df2NoOutliers)+aes(x=Method,y=MSE)+
  geom_boxplot(aes(fill=Method),outlier.shape = NA)+
  facet_wrap(.~Setting,scales="free")+
  labs(y="MSE",x="Method")+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
dev.off()

#Plot: time----
# figname<-paste("FigSimsPrediction.pdf",sep="")
# pdf(width = wdthpdf*1.5, height = hghtpdf*1.5,
#     file = figname)
ggplot(df)+aes(x=Method,y=time)+
  geom_boxplot(aes(fill=Method))+
  #facet_wrap(.~Method,scales="free_x",nrow=2)+
  labs(y="Time",x="Method")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
# dev.off()
