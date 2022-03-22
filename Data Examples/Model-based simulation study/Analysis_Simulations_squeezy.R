###Simulation study - linear regression
#NOTE: in Rstudio; press ALT+o for overview of code section headings
#This file may be used to replicate the analysis for Section 3.1.1-3.1.3

#make sure working directory is set to where this file is
setwd("") 
pathResults <- "" #if desired, separate folder can be specified to save results..
pathFigures <- "" #..and figures

#data simulation settings to run (one/multiple of "1"-"11")
Settings <- as.character(1:11)

#set T/F to run the following methods
runEN <- T #glmnet, fwelnet, ipflasso (small grid)
run_squeezy <- T #ecpc+squeezy, squeezy
run_ipflarge <- T #ipflasso (large grid)
run_rangeAlpha <- T #run squeezy for a range of alpha
alphaEN <- 0.3
#Or for specific sections, use the following simulation variables (uncomment relevant lines)
#Section 3.1.1:
# alphaEN <- 0.3 #or alphaEN <- 0.8 
# Settings <- as.character(1:3)
# runEN <- run_squeezy <- run_rangeAlpha <- T
# run_ipflarge <- F

#Section 3.1.2:
# alphaEN <- 1
# Settings <- as.character(1:3)
# run_squeezy <- T
# runEN <- run_ipflarge <- run_rangeAlpha <- F

#Section 3.1.3:
# alphaEN <- 0.3
# Settings <- as.character(4:11)
# runEN <- run_squeezy <- run_ipflarge <- T
# run_rangeAlpha <- F

#optional: set up parallel back-end to use multiple processors
runParallel <- T

#Install and load libraries----
#for parallel computing
#check whether each package has been installed, if not install package
if(!requireNamespace("foreach")) install.packages("foreach") 
if(!requireNamespace("doParallel")) install.packages("doParallel")
library(foreach)
library(doParallel)

#for installing development versions of packages on github
if(!requireNamespace("devtools")) install.packages("devtools")
library(devtools)

#methods to compare
if(!requireNamespace("glmnet")) install.packages("glmnet")
library(glmnet)
if(!requireNamespace("gren")) install.packages("gren")
library(gren)
if(!requireNamespace("ipflasso")) install.packages("ipflasso")
library(ipflasso)
if(!requireNamespace("fwelnet")) install_github("kjytay/fwelnet")
library(fwelnet)
if(!requireNamespace("multiridge")) install.packages("multiridge")
library(multiridge)
if(!requireNamespace("ecpc")) install.packages("ecpc")
library(ecpc)
if(!requireNamespace("squeezy")) install.packages("squeezy")
library(squeezy)

#for storing and plotting results
if(!requireNamespace("dplyr")) install.packages("dplyr")
library(dplyr)
if(!requireNamespace("ggplot2")) install.packages("ggplot2")
library(ggplot2)
if(!requireNamespace("RColorBrewer")) install.packages("RColorBrewer")
library(RColorBrewer)
if(!requireNamespace("scales")) install.packages("scales")
library(scales)
if(!requireNamespace("tidyr")) install.packages("tidyr")
library(tidyr)
if(!requireNamespace("ggpubr")) install.packages("ggpubr")
library(ggpubr)

#optional: setup parallel backend to use many processors
if(runParallel){
  cores <- detectCores() #check how many cores may be used
  cl <- makeCluster(cores-1) #not to overload your computer
  registerDoParallel(cl)
}

#generate seeds for reproducibility
#set.seed(1)
seeds <- round(runif(10000)*10^6)
#save(seeds,file=paste(pathResults,"seeds.Rdata",sep=""))
load(paste(pathResults,"seeds.Rdata",sep=""))

#Settings simulation----
#set some general simulation variables
n2<-300 #sample size of each test data set
Nblock <- 10     #correlation blocks in X
CorX <- 0.2 #correlation strength
q<-0 #set q*100% of the betas to 0
nSim <- 100 #number of simulated training&test data sets (reset inside for-loop)
n<-150 #number of training samples (reset inside for-loop)
p<-1200 #number of covariates (reset inside for-loop)


#Do simulations----
# setting <- "1"
for(setting in Settings){
  #Set specific simulation setting variables----
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
  
  #Make random groupset (used below to make informative groupset) ----
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
  #e.g. indRandom[[3]][[1]] is the vector with the covariate indices belonging to group 1 in the groupset with G[3] random groups
  #indRandom <- lapply(1:max(ind),function(x){rankRandom[which(ind==x)]})
  
  tauglobal2 <- rep(taugrp,each=p/G)
  if(G!=5){
    tauglobal2 <- rep(taugrp,each=p/5)
  }
  lambda1 <- sqrt(tauglobal2/2)
  dist <- "laplace"
  model <- "linear"
  
  #Generate correlated data----
  betas <- list()
  X <- list()
  X2 <- list()
  Y <- list()
  Y2 <- list()
  sigma2_emp <- c() #empirical residual variance in simulated linear regression
  lambda_emp <- list() #empirical residual variance in simulated linear regression with lasso prior
  tau2_emp <- list() #empirical residual variance in simulated linear regression with lasso prior
  dfGrps <- list()
  rankBeta <- list(); indRank <- list()
  for(i in 1:nSim){
    #simulate betas
    if(dist=="gaussian"){
      betas[[i]] <- rnorm(p,0,sqrt(tauglobal2))
      rankBeta[[i]] <- order(abs(betas[[i]]))
      lambdaApproxTrue <- sigmasq/taugrp
    } 
    if(dist=="laplace"){
      lambda1 <- sqrt(tauglobal2/2)
      # print(paste("Simulating betas from laplace(0,",signif(lambda1,3),")",sep=""))
      #A Laplace(0,b) variate can also be generated as the difference of two i.i.d. 
      #Exponential(1/b) random variables
      betas[[i]] <-   rexp(p, 1/lambda1) -  rexp(p, 1/lambda1)
      lambdaApproxTrue <- 2*sigmasq/sqrt(taugrp/2)
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
    
    #empirical 'true' estimates of the simulations for lasso prior
    sigma2_emp[i] <- n/(n-1)*mean((Y[[i]] - as.numeric(X[[i]]%*%betas[[i]]))^2) #residual variance for mean fixed at 0
    tau2_emp[[i]] <- sapply(indRank[[i]], function(x) 1/2/length(x)*sum(abs(betas[[i]][x]))) #0.5*average |\beta_j| in each group
    lambda_emp[[i]] <- sigma2_emp[i]/tau2_emp[[i]]
    dfGrps2 <- data.frame("Group"=1:G,"Lambda"=lambda_emp[[i]],"Sigma2"=rep(sigma2_emp[i],G),
                        "Method"=rep("True_emp",G),"Tau2"=tau2_emp[[i]])
    dfGrps2$Dataset <- i
    dfGrps2$n <- n
    dfGrps2$G <- G
    dfGrps2$p <- p
    dfGrps2$Setting <- setting
    dfGrps <- rbind(dfGrps,dfGrps2)
  }
  save(dfGrps,file=paste(pathResults,paste("simresSqueezy_GroupsEmpTrue",alphaEN,setting,sep="_"),".Rdata",sep=""))
  
  
  #Fit glmnet, fwelnet and ipflasso (set runEN=TRUE to run)----
  if(runEN){
    logname<-paste("logSimulations.txt",sep="")
    write(paste(Sys.time(),logname),file=logname,append=T)
    
    simpara <- list(n=n,ntest=n2,p=p,G=G,CorX=CorX,Nblock=Nblock,
                    alphaEN=alphaEN, model=model, dist= dist,taugrp=taugrp,sigmasq=sigmasq)
    fname <- paste(pathResults,paste("simresEN",alphaEN,setting,sep="_"),".Rdata",sep="")
    
    nSim2 <- nSim
    #nSim2 <- 5
    
    df <- data.frame()
    
    #for(i in 1:nSim2){
    finalMatrix <- foreach(i=1:nSim2, .combine=rbind,
                           .packages = c("glmnet","ipflasso","fwelnet")) %dopar% {
                                           write(paste(Sys.time(),i,logname),file=logname,append=T)
                                           
           #Fit elastic net with glmnet----
           alphat <- alphaEN
           # sd_y <- sqrt(var(Y[[i]])*(n-1)/n)[1]
           # alphat <- 1/(1+sd_y*(1-alphaEN)/alphaEN)
           
           pmt <- proc.time()[[3]]
           glmGR <-  cv.glmnet(X[[i]],Y[[i]],alpha=alphat,family=fml, intercept = T, standardize = F, thresh=10^-10)
           lambda1approx <- glmGR$lambda.min
           betaEN <- coef(glmGR,s=lambda1approx,exact=T)[-1]
           a0EN <- coef(glmGR,s=lambda1approx,exact=T)[1]
           MSEEN <- mean((Y2[[i]]-X2[[i]]%*%betaEN-a0EN)^2)
           ctglm <- proc.time()[[3]]-pmt
           
           MSEEN
           
           
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
           
           list("df"=df2)
         }
    df2 <- lapply(1:nSim2,function(i) finalMatrix[i,1][[1]])
    df <- df2[[1]]; for(i in 2:nSim2) df <- rbind(df,df2[[i]])
    
    save(df,simpara,file=fname)
  }


  #Fit squeezy (set run_squeezy=TRUE to run)----
  if(run_squeezy){
    logname<-paste("logSimulations.txt",sep="")
    write(paste(Sys.time(),logname),file=logname,append=T)
    
    simpara <- list(n=n,ntest=n2,p=p,G=G,CorX=CorX,Nblock=Nblock,
                    alphaEN=alphaEN, model=model, dist= dist,taugrp=taugrp,sigmasq=sigmasq)
    fname <- paste(pathResults,paste("simresSqueezy",alphaEN,setting,sep="_"),".Rdata",sep="")
    
    nSim2 <- nSim
    #nSim2 <- 5
    
    df <- data.frame()
    dfGrps <- data.frame()
    
    #for(i in 1:nSim2){
    finalMatrix <- foreach(i=1:nSim2, .combine=rbind,
                           .packages = c("ecpc","squeezy")) %dopar% {
             write(paste(Sys.time(),i,logname),file=logname,append=T)
             
             #Fit ecpc+squeezy defaults----
             
             #first fit ecpc to find group-ridge penalty estimates
             pmt <- proc.time()[[3]]
             res.ecpc<-ecpc(Y[[i]],X[[i]], #observed data and response to train model
                            groupsets=list(indRank[[i]]), #informative groupset
                            Y2=Y2[[i]],X2=X2[[i]], #test set
                            model=model,
                            hypershrinkage="none",postselection = F)
             #then use multivariate normal approximation to fit model for group-regularised elastic net
             ctecpc <- proc.time()[[3]] - pmt
             ctecpc
             
             #multigroup, no reCV
             pmt <- proc.time()[[3]]
             res.ecpcsqueezy <- squeezy(Y[[i]],X[[i]],groupset=indRank[[i]],Y2=Y2[[i]],X2=X2[[i]],
                                        model=model,intrcpt=T,
                                        alpha=alphaEN, #elastic net parameter for which group-regularised EN has to be fit
                                        method="MML", #one of c("ecpcEN","MML","MML.noDeriv","CV")
                                        fit.ecpc=res.ecpc, reCV=F) #fit of ecpc function
             ctsqueezy <- proc.time()[[3]]-pmt + ctecpc;
             ctsqueezy
             
             #multigroup, with reCV
             pmt <- proc.time()[[3]]
             res.ecpcsqueezy2 <- squeezy(Y[[i]],X[[i]],groupset=indRank[[i]],Y2=Y2[[i]],X2=X2[[i]],
                                         model=model,intrcpt=T,
                                         alpha=alphaEN, #elastic net parameter for which group-regularised EN has to be fit
                                         method="MML", #one of c("ecpcEN","MML","MML.noDeriv","CV")
                                         fit.ecpc=res.ecpc,reCV=T) #fit of ecpc function
             ctsqueezy2 <- proc.time()[[3]]-pmt + ctecpc;
             ctsqueezy2
             
             #first fit ecpc to find group-ridge penalty estimates
             pmt <- proc.time()[[3]]
             res.ecpc2<-ecpc(Y[[i]],X[[i]], #observed data and response to train model
                             groupsets=list(list(1:p)), #informative groupset
                             Y2=Y2[[i]],X2=X2[[i]], #test set
                             model=model,
                             hypershrinkage="none",postselection = F)
             #then use multivariate normal approximation to fit model for group-regularised elastic net
             ctecpc2 <- proc.time()[[3]] - pmt
             ctecpc2
             
             #one group, no reCV
             pmt <- proc.time()[[3]]
             res.ecpcsqueezy3 <- squeezy(Y[[i]],X[[i]],groupset=list(1:p),Y2=Y2[[i]],X2=X2[[i]],
                                         model=model,intrcpt=T,
                                         alpha=alphaEN, #elastic net parameter for which group-regularised EN has to be fit
                                         method="MML", #one of c("ecpcEN","MML","MML.noDeriv","CV")
                                         fit.ecpc=res.ecpc2, reCV=F) #fit of ecpc function
             ctsqueezy3 <- proc.time()[[3]]-pmt +ctecpc2;
             ctsqueezy3
             
             #one group, with reCV
             pmt <- proc.time()[[3]]
             res.ecpcsqueezy4 <- squeezy(Y[[i]],X[[i]],groupset=list(1:p),Y2=Y2[[i]],X2=X2[[i]],
                                         model=model,intrcpt=T,
                                         alpha=alphaEN, #elastic net parameter for which group-regularised EN has to be fit
                                         method="MML", #one of c("ecpcEN","MML","MML.noDeriv","CV")
                                         fit.ecpc=res.ecpc2,reCV=T) #fit of ecpc function
             ctsqueezy4 <- proc.time()[[3]]-pmt +ctecpc2;
             ctsqueezy4
             
             #directly transform moment ecpc-estimates with method="ecpcEN"
             pmt <- proc.time()[[3]]
             res.ecpcsqueezy5 <- squeezy(Y[[i]],X[[i]],groupset=indRank[[i]],Y2=Y2[[i]],X2=X2[[i]],
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
             
             
             
             mses_ecpcsq <- c(MSEmulti,MSEmultireCV, MSEsingle, MSEsinglereCV, MSEecpcEN)
             
             #summarise results----
             
             df2 <- data.frame("time"=c(ctsqueezy, ctsqueezy2,ctsqueezy3,ctsqueezy4,ctsqueezy5),
                               "MSE"=mses_ecpcsq,
                               "Method"=c("squeezy.m","squeezy.m.reCV","squeezy.s","squeezy.s.reCV","ecpcEN"))
             df2$Dataset <- i
             df2$n <- n
             df2$G <- G
             df2$p <- p
             
             dfGrps <- data.frame("Group"=rep(1:G,3),
                                  "Lambda"=c(res.ecpcsqueezy$lambdaApprox,res.ecpcsqueezy5$lambdaApprox,lambdaApproxTrue),
                                  "Sigma2"=rep(c(res.ecpcsqueezy$sigmahat,res.ecpcsqueezy5$sigmahat,sigmasq),each=G),
                                  "Method"=rep(c("squeezy.m","ecpcEN","Truth"),each=G))
             dfGrps$Tau2 <- dfGrps$Sigma2/dfGrps$Lambda
             dfGrps$Dataset <- i
             dfGrps$n <- n
             dfGrps$G <- G
             dfGrps$p <- p
             dfGrps$Setting <- setting
             
             list("df"=df2, "dfGrps"=dfGrps)
           }
    df2 <- lapply(1:nSim2,function(i) finalMatrix[i,1][[1]])
    dfGrps2 <- lapply(1:nSim2,function(i) finalMatrix[i,2][[1]])
    df <- df2[[1]]; for(i in 2:nSim2) df <- rbind(df,df2[[i]])
    dfGrps <- dfGrps2[[1]]; for(i in 2:nSim2) dfGrps <- rbind(dfGrps,dfGrps2[[i]])
    
    save(df,dfGrps,simpara,file=fname)
  }
  


  #Fit ipflasso on larger grid (set run_ipflarge=TRUE to run)----
  if(run_ipflarge & !(setting%in%as.character(c(1,2,3)))){
    fname <- paste(pathResults,paste("simresIPF2",alphaEN,setting,sep="_"),".Rdata",sep="")

    #nSim2 <- nSim
    #nSim2 <- 5

    df <- data.frame()

    #for(i in 1:nSim2){
    finalMatrix <- foreach(i=1:nSim2, .combine=rbind,
                           .packages = c("glmnet","ipflasso")) %dopar% {
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


  #Fit squeezy for range of alpha----
  if(run_rangeAlpha & setting=="3"){
    simpara <- list(n=n,ntest=n2,p=p,G=G,CorX=CorX,Nblock=Nblock,
                    alphaEN=alphaEN, model=model, dist= dist,taugrp=taugrp,sigmasq=sigmasq)
    fname <- paste(pathResults,paste("simressqueezytemp2",alphaEN,setting,sep="_"),".Rdata",sep="")
    
    nSim2 <- nSim
    #nSim2 <- 5
    
    df <- data.frame()
    res.ecpcsqueezy <- list() 
    #for(i in 1:nSim2){
    finalMatrix <- foreach(i=1:nSim2, .combine=rbind,
                           .packages = c("glmnet","penalized","mvtnorm","gglasso","mgcv",
                                         "CVXR","GRridge","expm","Rsolnp","ecpc",
                                         "ipflasso","fwelnet","gren","squeezy")) %dopar% {
                   write(paste(Sys.time(),i,logname),file=logname,append=T)
                   #Fit ecpc+squeezy defaults----
                   
                   #first fit ecpc to find group-ridge penalty estimates
                   pmt <- proc.time()[[3]]
                   res.ecpc<-ecpc(Y[[i]],X[[i]], #observed data and response to train model
                                  groupsets=list(indRank[[i]]), #informative groupset
                                  Y2=Y2[[i]],X2=X2[[i]], #test set
                                  model=model,
                                  hypershrinkage="none",postselection = F)
                   #then use multivariate normal approximation to fit model for group-regularised elastic net
                   ctecpc <- proc.time()[[3]] - pmt
                   ctecpc
                   pmt <- proc.time()[[3]]
                   res.ecpcsqueezy <- squeezy(Y[[i]],X[[i]],groupset=indRank[[i]],Y2=Y2[[i]],X2=X2[[i]],
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
    
    save(res.ecpcsqueezy,simpara,file=fname)
    
    #load(fname)
    
    alp <- seq(0,1,length.out=100)
    df <- data.frame()
    dfGrps <- data.frame()
    #for(i in 1:nSim2){
    finalMatrix <- foreach(i=1:nSim2, .combine=rbind,
                           .packages = c("glmnet","penalized","mvtnorm","gglasso","mgcv",
                                         "CVXR","GRridge","expm","Rsolnp","ecpc",
                                         "ipflasso","fwelnet","gren","squeezy")) %dopar% {
                     write(paste(Sys.time(),i,logname),file=logname,append=T)
                     #print(i)
                     df <- data.frame()
                     dfGrps <- data.frame()
                     for(j in 1:length(alp)){
                       #print(c(i,j))
                       #without reCV
                       pmt <- proc.time()[[3]]
                       res.ecpcsqueezyAlp <- squeezy(Y[[i]],X[[i]],groupset=indRank[[i]],Y2=Y2[[i]],X2=X2[[i]],
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
                       res.ecpcsqueezyAlp2 <- squeezy(Y[[i]],X[[i]],groupset=indRank[[i]],Y2=Y2[[i]],X2=X2[[i]],
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
    
  }
}






