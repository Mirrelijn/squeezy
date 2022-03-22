##Farkas data (n=37, p=40000)
#NOTE: in Rstudio; press ALT+o for overview of code section headings
#This file may be used to replicate the analysis for Section 3.3

#make sure working directory is set to where this file is
setwd("") 
pathResults <- "" #if desired, separate folder can be specified to save results..
pathFigures <- "" #..and figures
pathData <- "" #path to data and co-data

#set T/F to run the following analyses/methods
#(all should be rerun to reproduce the figures/tables below)
run_normalitycheck <- T #run squeezy for normality check
run_squeezy <- T #ecpc+squeezy, squeezy
runEN <- T #glmnet
run_fwelnet <- T #fwelnet on groups
run_ipf <- T #ipflasso (small grid)
run_ipflarge <- T #ipflasso (large grid)
run_gren <- T #gren

#optional: set up parallel back-end to use multiple processors
runParallel <- F

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
if(!requireNamespace("gren")) install_github("magnusmunch/gren/rpackage")
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

#for data (contained in GRridge package)
if(!requireNamespace("GRridge")){
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("GRridge")
}
library(GRridge)

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
if(!requireNamespace("ggrepel")) install.packages("ggrepel")
library(ggrepel)
if(!requireNamespace("ggraph")) install.packages("ggraph")
library(ggraph)
if(!requireNamespace("igraph")) install.packages("igraph")
library(igraph)

#optional: setup parallel backend to use many processors
if(runParallel){
  cores <- detectCores() #check how many cores may be used
  cl <- makeCluster(cores-1) #not to overload your computer
  registerDoParallel(cl)
}

#Load data----
data(dataFarkas) #contained in GRridge-package

means <- apply(as.matrix(datcenFarkas),1,mean) #abundance: average gene expression
Xcen <- as.matrix(datcenFarkas)-means
dim(Xcen) #centered data
sds <- apply(Xcen,1,sd) #standard deviations
Xstd<- t(Xcen/sds) #standardized gene expression

n <- dim(Xstd)[1] #number of samples
p <- dim(Xstd)[2] #number of covariates

Y <- respFarkas
Y<-as.numeric(Y)-1 #response: numeric values of 0/1

#Define co-data groupsets and types of hypershrinkage----
GroupsetCpG <- createGroupset(CpGannFarkas) #groups for squeezy

#variables for ecpc
GroupsetsAll <- list(GroupsetCpG)
hypershrinkage <- "none"

G <- length(GroupsetCpG)

#variables for fwelnet
Z <- matrix(rep(0,p*G),c(p,G))
for(j in 1:G){
  Z[GroupsetCpG[[j]],j] <- 1
}
fml <- "binomial"

#variables for gren
ind.vector <- rep(1,p)
for(j in 2:G){
  ind.vector[GroupsetCpG[[j]]] <- j
}

#Define folds and variables for the CV----
logname<-paste("logCVFarkas.txt",sep="")
write(paste(Sys.time(),logname),file=logname,append=T)

set.seed(3434) 
nfolds<-10
folds2<-produceFolds(n,nfolds,Y,balance=T)
#save(folds2,file="foldsFarkas.Rdata")
load("foldsFarkas.Rdata") #use the same folds for different methods

set.seed(3434)
B <- 10
Bfolds <- list()
for(b in 1:B){
  Bfolds[[b]] <- produceFolds(n,nfolds,Y,balance=T)
}
#save(Bfolds,file="BfoldsFarkas.Rdata")
load("BfoldsFarkas.Rdata")

alp <- c(0.01,0.3,0.8,1) #alpha variable used in elastic net
j<-2 #alpha alp[j]=0.3

set.seed(474)
subsampleVars <- sample(1:p,p/2,replace=F)
#save(subsampleVars,file="subsampleVars.Rdata")
load("subsampleVars.Rdata")


#p/4 still too large
set.seed(8686)
subsampleVars2 <- sample(1:p,size=floor(p/8),replace=T)
#save(subsampleVars2,file="subsampleVars2.Rdata")
load("subsampleVars2.Rdata")


#Fit squeezy on whole data for AIC results and normality check----
if(run_normalitycheck){
  fname <- paste(pathResults,"FarkasResAll",sep="")
  Res<-ecpc(Y,Xstd,groupsets=GroupsetsAll,hypershrinkage=hypershrinkage,
            postselection=F)
  res.squeezy <- squeezy(Y,Xstd,groupset=GroupsetsAll[[1]],
                         model="logistic",
                         alpha=alp[j], #elastic net parameter for which group-regularised EN has to be fit
                         fit.ecpc=Res,
                         method="MML",reCV=T)
  save(Res,res.squeezy,file=fname)
} 


#Do CV for ecpc+squeezy defaults----

if(run_squeezy){
  for(b in 1:B){
    folds2 <- Bfolds[[b]]
    fname <- paste(pathResults,"CVFarkassqueezyalpha",b,j,sep="")
    
    df<-data.frame() #data frame in which predicted values are stores for each of the samples in the left-out fold
    Res<-list() #list in which raw output of ecpc is stored (e.g. estimated regression coefficients)
    Res2<-list() #list in which raw output of ecpc is stored (e.g. estimated regression coefficients)
    res.squeezy <- list()
    res.squeezy2 <- list()
    res.squeezy3 <- list()
    res.squeezy4 <- list()
    res.squeezy5 <- list()
    #for(i in 1:nfolds){
    finalMatrix <- foreach(i=1:nfolds, .combine=rbind,
                           .packages = c("ecpc","squeezy")) %do% {
                       #fit ecpc multigroup for intialisation----
                       tic<-proc.time()[[3]]
                       Res[[i]]<-ecpc(Y[-folds2[[i]]],Xstd[-folds2[[i]],],
                                      groupsets=GroupsetsAll,hypershrinkage=hypershrinkage,
                                      Y2=Y[folds2[[i]]],X2=Xstd[folds2[[i]],],
                                      postselection=F)
                       Res[[i]]$timeGR <- proc.time()[[3]]-tic
                       
                       #fit squeezy multigroup, no reCV----
                       tic<-proc.time()[[3]]
                       res.squeezy[[i]] <- squeezy(Y[-folds2[[i]]],Xstd[-folds2[[i]],],groupset=GroupsetsAll[[1]],
                                                   Y2=Y[folds2[[i]]],X2=Xstd[folds2[[i]],],
                                                   model="logistic",
                                                   alpha=alp[j], #elastic net parameter for which group-regularised EN has to be fit
                                                   fit.ecpc=Res[[i]],
                                                   method="MML",reCV=F) 
                       res.squeezy[[i]]$timeGR <- proc.time()[[3]]-tic + Res[[i]]$timeGR
                       
                       #fit squeezy multigroup, with reCV----
                       tic <- proc.time()[[3]]
                       res.squeezy2[[i]] <- squeezy(Y[-folds2[[i]]],Xstd[-folds2[[i]],],groupset=GroupsetsAll[[1]],
                                                    Y2=Y[folds2[[i]]],X2=Xstd[folds2[[i]],],
                                                    model="logistic",
                                                    alpha=alp[j], #elastic net parameter for which group-regularised EN has to be fit
                                                    fit.ecpc=Res[[i]],
                                                    method="MML",reCV=T) 
                       res.squeezy2[[i]]$timeGR <- proc.time()[[3]]-tic + Res[[i]]$timeGR
                       
                       #fit ecpcEN with reCV----
                       tic <- proc.time()[[3]]
                       res.squeezy3[[i]] <- squeezy(Y[-folds2[[i]]],Xstd[-folds2[[i]],],groupset=GroupsetsAll[[1]],
                                                    Y2=Y[folds2[[i]]],X2=Xstd[folds2[[i]],],
                                                    model="logistic",
                                                    alpha=alp[j], #elastic net parameter for which group-regularised EN has to be fit
                                                    fit.ecpc=Res[[i]],
                                                    method="ecpcEN",reCV=T) 
                       res.squeezy3[[i]]$timeGR <- proc.time()[[3]]-tic + Res[[i]]$timeGR
                       
                       #fit ecpc one group for initialisation----
                       tic<-proc.time()[[3]]
                       temp <- list(list(1:p))
                       Res2[[i]]<-ecpc(Y[-folds2[[i]]],Xstd[-folds2[[i]],],
                                       groupsets=temp,hypershrinkage=hypershrinkage,
                                      Y2=Y[folds2[[i]]],X2=Xstd[folds2[[i]],],
                                      postselection=F)
                       Res2[[i]]$timeGR <- proc.time()[[3]]-tic
                       
                       #fit squeezy one group, no reCV----
                       tic<-proc.time()[[3]]
                       res.squeezy4[[i]] <- squeezy(Y[-folds2[[i]]],Xstd[-folds2[[i]],],groupset=list(1:p),
                                                    Y2=Y[folds2[[i]]],X2=Xstd[folds2[[i]],],
                                                    model="logistic",
                                                    alpha=alp[j], #elastic net parameter for which group-regularised EN has to be fit
                                                    fit.ecpc=Res2[[i]],
                                                    method="MML",reCV=F) 
                       res.squeezy4[[i]]$timeGR <- proc.time()[[3]]-tic + Res2[[i]]$timeGR
                       
                       #fit squeezy one group, with reCV----
                       tic <- proc.time()[[3]]
                       res.squeezy5[[i]] <- squeezy(Y[-folds2[[i]]],Xstd[-folds2[[i]],],groupset=list(1:p),
                                                    Y2=Y[folds2[[i]]],X2=Xstd[folds2[[i]],],
                                                    model="logistic",
                                                    alpha=alp[j], #elastic net parameter for which group-regularised EN has to be fit
                                                    fit.ecpc=Res2[[i]],
                                                    method="MML",reCV=T) 
                       res.squeezy5[[i]]$timeGR <- proc.time()[[3]]-tic + Res2[[i]]$timeGR
                       
                       
                       #Summarise results----
                       df2<-data.frame("Ypred"=c(res.squeezy[[i]]$YpredApprox,res.squeezy2[[i]]$YpredApprox,
                                                 res.squeezy3[[i]]$YpredApprox,res.squeezy4[[i]]$YpredApprox,
                                                 res.squeezy5[[i]]$YpredApprox))
                       df2$Method <- rep(c("squeezy.m","squeezy.m.reCV","ecpcEN","squeezy.s","squeezy.s.reCV"),each=length(folds2[[i]]))
                       df2$NumberSelectedVars <- rep(c(sum(res.squeezy[[i]]$betaApprox!=0),sum(res.squeezy2[[i]]$betaApprox!=0),
                                                       sum(res.squeezy3[[i]]$betaApprox!=0),sum(res.squeezy4[[i]]$betaApprox!=0),
                                                       sum(res.squeezy5[[i]]$betaApprox!=0)),each=length(folds2[[i]]))
                       df2$Fold <- i
                       df2$Sample <- rep(folds2[[i]],5)
                       df2$Time <-  rep(c(res.squeezy[[i]]$timeGR,res.squeezy2[[i]]$timeGR, 
                                          res.squeezy3[[i]]$timeGR,res.squeezy4[[i]]$timeGR,
                                          res.squeezy5[[i]]$timeGR),each=length(folds2[[i]]))
                       df2$Truth <- rep(Y[folds2[[i]]],5)
                       df2$Alpha <- alp[j]
                       
                       write(paste(Sys.time(),"fold",i,"of",nfolds,"done"),file=logname,append=T)
                       
                       list("Res"=Res,"df"=df2,"res.ecpc"=res.squeezy,"res.ecpc2"=res.squeezy2,
                            "res.ecpc3"=res.squeezy3,"res.ecpc4"=res.squeezy4,"res.ecpc5"=res.squeezy5)
                     }
    
    Res <- lapply(1:nfolds,function(i) finalMatrix[i,1][[1]][[i]])
    #Res2 <- lapply(1:nfolds,function(i) finalMatrix[i,2][[1]][[i]])
    df2 <- lapply(1:nfolds,function(i) finalMatrix[i,2][[1]])
    df <- df2[[1]]; for(i in 2:nfolds) df <- rbind(df,df2[[i]])
    res.squeezy <- lapply(1:nfolds,function(i) finalMatrix[i,3][[1]][[i]])
    res.squeezy2 <- lapply(1:nfolds,function(i) finalMatrix[i,4][[1]][[i]])
    res.squeezy3 <- lapply(1:nfolds,function(i) finalMatrix[i,5][[1]][[i]])
    res.squeezy4 <- lapply(1:nfolds,function(i) finalMatrix[i,6][[1]][[i]])
    res.squeezy5 <- lapply(1:nfolds,function(i) finalMatrix[i,7][[1]][[i]])
    #save(Res,df,res.squeezy,res.squeezy2,res.squeezy3,
    #     res.squeezy4,res.squeezy5,file=paste(pathResults,"CVFarkassqueezy",sep=""))
    #stopCluster(cl); rm(cl)
    
    #data frame with summary statistics
    Summdf <- df %>% dplyr::group_by(Method,Alpha) %>% summarise(Brier = mean((Ypred-Truth)^2),
                                                                 CVLL = sum(log(Ypred)*Truth+log(1-Ypred)*(1-Truth)),
                                                                 NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
    
    df$Method<-as.factor(df$Method)
    
    dfROC<-data.frame()
    for(i in levels(df$Method)){
      temp<-data.frame()
      cutoffs<-rev(seq(0,1,by=0.001))
      rocGR <- GRridge::roc(probs=df$Ypred[df$Method==i],true=df$Truth[df$Method==i],cutoffs=cutoffs)
      temp<-data.frame("FPR"=rocGR[1,],"TPR"=rocGR[2,],"Accuracy"=rocGR[3,])
      temp$Method <- i
      temp$Alpha <- alp[j]
      temp$AUC<-c(GRridge::auc(rocGR))
      temp$NumberSelectedVars<-mean(df$NumberSelectedVars[df$Method==i])
      dfROC<-rbind(dfROC,temp)
    }
    dfAUC <- dfROC %>% group_by(Method,Alpha) %>% summarise(AUC=mean(AUC),
                                                            NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
    
    
    if(all(Summdf$NumberSelectedVars==dfAUC$NumberSelectedVars)){ #add AUC to the summary stastistics data frame
      Summdf$AUC <- dfAUC$AUC
    }
    Summdf
    save(df,Summdf,dfROC,dfAUC,Res,Res2,res.squeezy,res.squeezy2,res.squeezy3,
         res.squeezy4,res.squeezy5,file=fname) #this comparison
  }
}

#Do CV for elastic net-----

if(runEN){
  for(b in 1:B){
    folds2 <- Bfolds[[b]]
    fname <- paste(pathResults,"CVFarkasEN",b,sep="")
    
    dfElNet <- data.frame()
    dfBetaElNet <- data.frame()
    ResElNet <- list()
    for(j in 1:length(alp)){
      for(i in 1:nfolds){
        tic <- proc.time()[[3]]
        ResElNet[[i]]<-cv.glmnet(Xstd[-folds2[[i]],],Y[-folds2[[i]]],alpha=alp[j],family="binomial",
                                 intercept=T) 
        ResElNet[[i]]$time <- proc.time()[[3]] - tic
        Ypred<-predict(ResElNet[[i]],newx=Xstd[folds2[[i]],], type="response",s=ResElNet[[i]]$lambda.min)
        betaEN <- coef(ResElNet[[i]],s=ResElNet[[i]]$lambda.min,exact=T)
        
        df2<-data.frame("Ypred"=c(Ypred))
        df2$Method <- "EN"
        df2$NumberSelectedVars <- sum(betaEN!=0)
        df2$Fold <- i
        df2$Sample <- folds2[[i]]
        df2$Time <-  ResElNet[[i]]$time
        df2$Truth <- Y[folds2[[i]]]
        df2$Alpha <- alp[j]
        dfElNet<-rbind(dfElNet,df2)
        
        write(paste(Sys.time(),"fold",i,"of",nfolds,"done, elastic net"),file=logname,append=T)
        
      }
    }
    
    #data frame with summary statistics
    SummdfElNet <- dfElNet %>% group_by(Alpha,Method) %>% summarise(Brier = mean((Ypred-Truth)^2),
                                                                    CVLL = sum(log(Ypred)*Truth+log(1-Ypred)*(1-Truth)),
                                                                    NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
    
    dfElNet$Method<-as.factor(dfElNet$Method)
    
    dfROCElNet<-data.frame()
    for(i in levels(dfElNet$Method)){
      for(j in unique(dfElNet$Alpha)){
        temp<-data.frame()
        cutoffs<-rev(seq(0,1,by=0.001))
        rocGR <- GRridge::roc(probs=dfElNet$Ypred[dfElNet$Method==i & dfElNet$Alpha==j],
                              true=dfElNet$Truth[dfElNet$Method==i & dfElNet$Alpha==j],
                              cutoffs=cutoffs)
        temp<-data.frame("FPR"=rocGR[1,],"TPR"=rocGR[2,],"Accuracy"=rocGR[3,])
        temp$Method <- i
        temp$Alpha <- j
        temp$AUC<-c(GRridge::auc(rocGR))
        temp$NumberSelectedVars<-mean(dfElNet$NumberSelectedVars[dfElNet$Method==i & dfElNet$Alpha==j])
        dfROCElNet<-rbind(dfROCElNet,temp)
      }
    }
    dfAUCElNet <- dfROCElNet %>% group_by(Method,Alpha) %>% 
      summarise(AUC=mean(AUC),NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
    
    #plot(dfAUCElNet$NumberSelectedVars,dfAUCElNet$AUC)
    
    if(all(SummdfElNet$NumberSelectedVars==dfAUCElNet$NumberSelectedVars)){ #add AUC to the summary stastistics data frame
      SummdfElNet$AUC <- dfAUCElNet$AUC
    }
    
    save(dfElNet,SummdfElNet,dfROCElNet,dfAUCElNet,ResElNet,dfBetaElNet,file=fname) 
  }
}

#Do CV for fwelnet----

if(run_fwelnet){ #set to 1 to perform analysis or run code inside block manually
  for(b in 1:B){
    folds2 <- Bfolds[[b]]
    fname <- paste(pathResults,"CVFarkasfwENalpha",b,j,sep="")
    
    df<-data.frame() #data frame in which predicted values are stores for each of the samples in the left-out fold
    res.fwEN <- list()
    #for(i in 1:nfolds){
    finalMatrix <- foreach(i=1:nfolds, .combine=rbind,
                           .packages = c("glmnet","fwelnet")) %dopar% {
               tic<-proc.time()[[3]]
               res.fwEN[[i]] <-  cv.fwelnet(x=Xstd[-folds2[[i]],],y=Y[-folds2[[i]]],z=Z,
                                            alpha=alp[j],family=fml,standardize=F)
               res.fwEN[[i]]$time <- proc.time()[[3]]-tic
               
               ind.minlam <- which(res.fwEN[[i]]$lambda==res.fwEN[[i]]$lambda.min)
               
               betafwEN <- res.fwEN[[i]]$glmfit$beta[,ind.minlam]
               a0fwEN <- res.fwEN[[i]]$glmfit$a0[ind.minlam]
               lp <- Xstd[folds2[[i]],]%*%betafwEN+a0fwEN
               Ypred <- 1/(1+exp(-lp))
               
               df2<-data.frame("Ypred"=Ypred)
               df2$Method <- "fwEN"
               df2$NumberSelectedVars <- sum(betafwEN!=0)
               df2$Fold <- i
               df2$Sample <- folds2[[i]]
               df2$Time <-  res.fwEN[[i]]$time
               df2$Truth <- Y[folds2[[i]]]
               df2$Alpha <- alp[j]
               
               write(paste(Sys.time(),"fold",i,"of",nfolds,"done"),file=logname,append=T)
               
               list("res.fwEN"=res.fwEN,"df"=df2)
             }
    
    res.fwEN <- lapply(1:nfolds,function(i) finalMatrix[i,1][[1]][[i]])
    df2 <- lapply(1:nfolds,function(i) finalMatrix[i,2][[1]])
    df <- df2[[1]]; for(i in 2:nfolds) df <- rbind(df,df2[[i]])
    
    #save(res.fwEN,df,file=paste(pathResults,"CVFarkasfwENalpha",j,sep=""))
    #stopCluster(cl); rm(cl)
    
    #data frame with summary statistics
    Summdf <- df %>% group_by(Method,Alpha) %>% summarise(Brier = mean((Ypred-Truth)^2),
                                                          CVLL = sum(log(Ypred)*Truth+log(1-Ypred)*(1-Truth)),
                                                          NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
    
    df$Method<-as.factor(df$Method)
    
    dfROC<-data.frame()
    for(i in levels(df$Method)){
      temp<-data.frame()
      cutoffs<-rev(seq(0,1,by=0.001))
      rocGR <- GRridge::roc(probs=df$Ypred[df$Method==i],true=df$Truth[df$Method==i],cutoffs=cutoffs)
      temp<-data.frame("FPR"=rocGR[1,],"TPR"=rocGR[2,],"Accuracy"=rocGR[3,])
      temp$Method <- i
      temp$Alpha <- alp[j]
      temp$AUC<-c(GRridge::auc(rocGR))
      temp$NumberSelectedVars<-mean(df$NumberSelectedVars[df$Method==i])
      dfROC<-rbind(dfROC,temp)
    }
    dfAUC <- dfROC %>% group_by(Method,Alpha) %>% summarise(AUC=mean(AUC),
                                                            NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
    
    #plot(dfAUC$NumberSelectedVars,dfAUC$AUC)
    
    if(all(Summdf$NumberSelectedVars==dfAUC$NumberSelectedVars)){ #add AUC to the summary stastistics data frame
      Summdf$AUC <- dfAUC$AUC
    }
    Summdf
    save(df,Summdf,dfROC,dfAUC,res.fwEN,file=fname) #this comparison 
  }
}

#Do CV for ipflasso (small grid)----

if(run_ipf){ 
  #make list with proposal lambdas in the form of ratios 1:1:..:1, 0.5:..:1, etc.
  #weights vary from 2^-x to 2^x
  ratios <- 2^(0:1)
  A <- matrix(ratios,c(length(ratios),1))
  for(i in 2:G){
    l<-sapply(ratios,function(x){
      cbind(A,rep(x,length(ratios)))
    },simplify=F)
    A<-do.call(rbind,l)
  }
  
  pflist <- lapply(1:dim(A)[1],function(x)c(A[x,]))
  
  
  GroupsetIPF <- GroupsetsAll[[1]] #all variables
  
  for(b in 1:B){
    folds2 <- Bfolds[[b]]
    fname <- paste(pathResults,"CVFarkasipfalpha",b,j,sep="")
    
    df<-data.frame() #data frame in which predicted values are stores for each of the samples in the left-out fold
    #res.ipf <- list()
    #use sequential loop instead of parallel loop, else memory issues
    for(i in 1:nfolds){
      # finalMatrix <- foreach(i=1:nfolds, .combine=rbind,
      #                        .packages = c("glmnet","ipflasso")) %dopar% {
      tic<-proc.time()[[3]]
      res.ipf <- cvr2.ipflasso(X=Xstd[-folds2[[i]],],Y=Y[-folds2[[i]]],
                               blocks=GroupsetIPF,alpha=alp[j],
                               standardize = F,family="binomial",
                               type.measure="deviance",nfolds=10,ncv=1,pflist=pflist)
      res.ipf$timeGR <- proc.time()[[3]]-tic
      a0.ipf <- res.ipf$coeff[1,res.ipf$ind.bestlambda]
      beta.ipf <- c(res.ipf$coeff[-1,res.ipf$ind.bestlambda])
      
      Ypred <- 1/(1+exp(-Xstd[folds2[[i]],]%*%beta.ipf - a0.ipf))
      
      df2<-data.frame("Ypred"=Ypred)
      df2$Method <- "ipflasso"
      df2$NumberSelectedVars <- sum(beta.ipf!=0)
      df2$Fold <- i
      df2$Sample <- folds2[[i]]
      df2$Time <-  res.ipf$timeGR
      df2$Truth <- Y[folds2[[i]]]
      df2$Alpha <- alp[j]
      
      write(paste(Sys.time(),"fold",i,"of",nfolds,"done"),file=logname,append=T)
      
      df <- rbind(df,df2)
      #list("df"=df2,"df2"=df2)
    }
    
    # df2 <- lapply(1:nfolds,function(i) finalMatrix[i,1][[1]])
    # df <- df2[[1]]; for(i in 2:nfolds) df <- rbind(df,df2[[i]])
    # #res.ipf <- lapply(1:nfolds,function(i) finalMatrix[i,2][[1]][[i]])
    # save(df,file=paste(pathResults,"CVFarkasipfalpha",j,sep=""))
    #stopCluster(cl); rm(cl)
    
    #data frame with summary statistics
    Summdf <- df %>% group_by(Method,Alpha) %>% summarise(Brier = mean((Ypred-Truth)^2),
                                                          CVLL = sum(log(Ypred)*Truth+log(1-Ypred)*(1-Truth)),
                                                          NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
    
    df$Method<-as.factor(df$Method)
    
    dfROC<-data.frame()
    for(i in levels(df$Method)){
      temp<-data.frame()
      cutoffs<-rev(seq(0,1,by=0.001))
      rocGR <- GRridge::roc(probs=df$Ypred[df$Method==i],true=df$Truth[df$Method==i],cutoffs=cutoffs)
      temp<-data.frame("FPR"=rocGR[1,],"TPR"=rocGR[2,],"Accuracy"=rocGR[3,])
      temp$Method <- i
      temp$Alpha <- alp[j]
      temp$AUC<-c(GRridge::auc(rocGR))
      temp$NumberSelectedVars<-mean(df$NumberSelectedVars[df$Method==i])
      dfROC<-rbind(dfROC,temp)
    }
    dfAUC <- dfROC %>% group_by(Method,Alpha) %>% summarise(AUC=mean(AUC),
                                                            NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
    
    #plot(dfAUC$NumberSelectedVars,dfAUC$AUC)
    
    if(all(Summdf$NumberSelectedVars==dfAUC$NumberSelectedVars)){ #add AUC to the summary stastistics data frame
      Summdf$AUC <- dfAUC$AUC
    }
    Summdf
    save(df,Summdf,dfROC,dfAUC,file=fname) #this comparison 
  }
}



#Do CV for ipflasso (large grid)----

if(run_ipflarge){
  #make list with proposal lambdas in the form of ratios 1:1:..:1, 0.5:..:1, etc.
  #weights vary from 2^-x to 2^x
  ratios <- 2^(0:2)
  A <- matrix(ratios,c(length(ratios),1))
  for(i in 2:G){
    l<-sapply(ratios,function(x){
      cbind(A,rep(x,length(ratios)))
    },simplify=F)
    A<-do.call(rbind,l)
  }
  
  pflist <- lapply(1:dim(A)[1],function(x)c(A[x,]))
  
  
  #ipflasso: fit on subsample, else memory issues
  GroupsetIPF <- createGroupset(as.factor(ind.vector[subsampleVars2])) #subsample variables
  
  for(b in 1:B){
    folds2 <- Bfolds[[b]]
    fname <- paste(pathResults,"CVFarkasipfLargealpha",b,j,sep="")
    
    df<-data.frame() #data frame in which predicted values are stores for each of the samples in the left-out fold
    #res.ipf <- list()
    #use sequential loop instead of parallel loop, else memory issues
    for(i in 1:nfolds){
      # finalMatrix <- foreach(i=1:nfolds, .combine=rbind,
      #                        .packages = c("glmnet","penalized","mvtnorm","gglasso","mgcv",
      #                                      "CVXR","GRridge","expm","Rsolnp","ipflasso")) %dopar% {
      tic<-proc.time()[[3]]
      res.ipf <- cvr2.ipflasso(X=Xstd[-folds2[[i]],subsampleVars2],Y=Y[-folds2[[i]]],
                               blocks=GroupsetIPF,alpha=alp[j],
                               standardize = F,family="binomial",
                               type.measure="deviance",nfolds=10,ncv=1,pflist=pflist)
      res.ipf$timeGR <- proc.time()[[3]]-tic
      a0.ipf <- res.ipf$coeff[1,res.ipf$ind.bestlambda]
      beta.ipf <- c(res.ipf$coeff[-1,res.ipf$ind.bestlambda])
      
      Ypred <- 1/(1+exp(-Xstd[folds2[[i]],subsampleVars2]%*%beta.ipf - a0.ipf))
      
      df2<-data.frame("Ypred"=Ypred)
      df2$Method <- "ipflasso"
      df2$NumberSelectedVars <- sum(beta.ipf!=0)
      df2$Fold <- i
      df2$Sample <- folds2[[i]]
      df2$Time <-  res.ipf$timeGR
      df2$Truth <- Y[folds2[[i]]]
      df2$Alpha <- alp[j]
      
      write(paste(Sys.time(),"fold",i,"of",nfolds,"done"),file=logname,append=T)
      
      df <- rbind(df,df2)
      #list("df"=df2,"df2"=df2)
    }
    
    # df2 <- lapply(1:nfolds,function(i) finalMatrix[i,1][[1]])
    # df <- df2[[1]]; for(i in 2:nfolds) df <- rbind(df,df2[[i]])
    # #res.ipf <- lapply(1:nfolds,function(i) finalMatrix[i,2][[1]][[i]])
    # save(df,file=paste(pathResults,"CVFarkasipfalpha",j,sep=""))
    #stopCluster(cl); rm(cl)
    
    #data frame with summary statistics
    Summdf <- df %>% group_by(Method,Alpha) %>% summarise(Brier = mean((Ypred-Truth)^2),
                                                          CVLL = sum(log(Ypred)*Truth+log(1-Ypred)*(1-Truth)),
                                                          NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
    
    df$Method<-as.factor(df$Method)
    
    dfROC<-data.frame()
    for(i in levels(df$Method)){
      temp<-data.frame()
      cutoffs<-rev(seq(0,1,by=0.001))
      rocGR <- GRridge::roc(probs=df$Ypred[df$Method==i],true=df$Truth[df$Method==i],cutoffs=cutoffs)
      temp<-data.frame("FPR"=rocGR[1,],"TPR"=rocGR[2,],"Accuracy"=rocGR[3,])
      temp$Method <- i
      temp$Alpha <- alp[j]
      temp$AUC<-c(GRridge::auc(rocGR))
      temp$NumberSelectedVars<-mean(df$NumberSelectedVars[df$Method==i])
      dfROC<-rbind(dfROC,temp)
    }
    dfAUC <- dfROC %>% group_by(Method,Alpha) %>% summarise(AUC=mean(AUC),
                                                            NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
    
    #plot(dfAUC$NumberSelectedVars,dfAUC$AUC)
    
    if(all(Summdf$NumberSelectedVars==dfAUC$NumberSelectedVars)){ #add AUC to the summary stastistics data frame
      Summdf$AUC <- dfAUC$AUC
    }
    Summdf
    save(df,Summdf,dfROC,dfAUC,file=fname) #this comparison 
  }
}



#Do CV for gren----

if(run_gren){
  for(b in 1:B){
    folds2 <- Bfolds[[b]]
    fname <- paste(pathResults,"CVFarkasgrenalpha",b,j,sep="")
    
    #standardize in gren cannot be set to false, remove couple of genes that are only 0 in subsamples (else error in gren)
    ind.remove <- c()
    for(i in 1:nfolds){
      remove <- which(is.nan(apply(scale(Xstd[-folds2[[i]],]),2,mean)))
      ind.remove <- unique(c(ind.remove,remove))
    }
    length(ind.remove)
    
    df<-data.frame() #data frame in which predicted values are stores for each of the samples in the left-out fold
    res.gren <- list()
    #for(i in 1:nfolds){
    #Note: loop executed sequentially as parallel runs into memory issue
    finalMatrix <- foreach(i=1:nfolds, .combine=rbind,
                           .packages = c("gren")) %do% {
                       tic<-proc.time()[[3]]
                       res.gren[[i]] <- gren(x=Xstd[-folds2[[i]],subsampleVars],y=Y[-folds2[[i]]],
                                             partitions=list(ind.vector[subsampleVars]),
                                             alpha=alp[j],trace=F)
                       res.gren[[i]]$timeGR <- proc.time()[[3]]-tic
                       
                       betagren <- as.vector(coef(res.gren[[i]], type="groupreg", s=res.gren[[i]]$lambda))[-1]
                       a0gren <- as.vector(coef(res.gren[[i]], type="groupreg", s=res.gren[[i]]$lambda))[1]
                       
                       Ypred <- 1/(1+exp(-Xstd[folds2[[i]],subsampleVars]%*%betagren - a0gren))
                       
                       df2<-data.frame("Ypred"=Ypred)
                       df2$Method <- "gren"
                       df2$NumberSelectedVars <- sum(betagren!=0)
                       df2$Fold <- i
                       df2$Sample <- folds2[[i]]
                       df2$Time <-  res.gren[[i]]$timeGR
                       df2$Truth <- Y[folds2[[i]]]
                       df2$Alpha <- alp[j]
                       
                       write(paste(Sys.time(),"fold",i,"of",nfolds,"done"),file=logname,append=T)
                       
                       list("df"=df2,"res.gren"=res.gren)
                     }
    
    df2 <- lapply(1:nfolds,function(i) finalMatrix[i,1][[1]])
    df <- df2[[1]]; for(i in 2:nfolds) df <- rbind(df,df2[[i]])
    res.gren <- lapply(1:nfolds,function(i) finalMatrix[i,2][[1]][[i]])
    #save(df,res.gren,file=paste(pathResults,"CVFarkasgren",j,sep=""))
    #stopCluster(cl); rm(cl)
    
    #data frame with summary statistics
    Summdf <- df %>% group_by(Method,Alpha) %>% summarise(Brier = mean((Ypred-Truth)^2),
                                                          CVLL = sum(log(Ypred)*Truth+log(1-Ypred)*(1-Truth)),
                                                          NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
    
    df$Method<-as.factor(df$Method)
    
    dfROC<-data.frame()
    for(i in levels(df$Method)){
      temp<-data.frame()
      cutoffs<-rev(seq(0,1,by=0.001))
      rocGR <- GRridge::roc(probs=df$Ypred[df$Method==i],true=df$Truth[df$Method==i],cutoffs=cutoffs)
      temp<-data.frame("FPR"=rocGR[1,],"TPR"=rocGR[2,],"Accuracy"=rocGR[3,])
      temp$Method <- i
      temp$Alpha <- alp[j]
      temp$AUC<-c(GRridge::auc(rocGR))
      temp$NumberSelectedVars<-mean(df$NumberSelectedVars[df$Method==i])
      dfROC<-rbind(dfROC,temp)
    }
    dfAUC <- dfROC %>% group_by(Method,Alpha) %>% summarise(AUC=mean(AUC),
                                                            NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
    
    
    if(all(Summdf$NumberSelectedVars==dfAUC$NumberSelectedVars)){ #add AUC to the summary stastistics data frame
      Summdf$AUC <- dfAUC$AUC
    }
    Summdf
    save(df,Summdf,dfROC,dfAUC,res.gren,file=fname) #this comparison 
  }
}

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
colpal <- "Dark2"
colsfill <- brewer.pal(3,"Dark2")[1:2]
colsAUC <- brewer.pal(10,"RdYlBu")
colsAUC <- seq_gradient_pal(low="grey50",high="white")((0:9)/10)
colsAUC <- seq_gradient_pal(low="black",high="grey50")((0:9)/10)


#Load results, given in files:----
Summdf2 <- data.frame()
df2 <- data.frame()
SummGrps <- data.frame()

for(b in 1:B){
  #squeezy defaults
  j<-2 #alpha
  fname <- paste(pathResults,"CVFarkassqueezyalpha",b,j,sep="")
  load(fname)
  Summdf$B <- b
  df$B <- b
  Summdf2 <- rbind(Summdf2,Summdf)
  df2 <- rbind(df2,df)
  
  #gather multigroup estimates
  for(i in 1:nfolds){
    temp <- data.frame("Taus"=res.squeezy2[[i]]$tauMR,
                       "TausApprox"=1/res.squeezy2[[i]]$lambdaApprox)
    temp$Tauglobal <- 1/res.squeezy2[[i]]$lambdaglobal
    temp$TausRelative <- temp$Taus/temp$Tauglobal
    temp$Fold <- i
    temp$B <- b
    temp$Group <- names(GroupsetCpG)
    
    SummGrps <- rbind(SummGrps,temp)
  }
  
  #elastic net
  fname <- paste(pathResults,"CVFarkasEN",b,sep="")
  load(fname)
  SummdfElNet$B <- b
  dfElNet$B <- b
  Summdf2 <- rbind(Summdf2,SummdfElNet[SummdfElNet$Alpha==0.3,])
  df2 <- rbind(df2,dfElNet[dfElNet$Alpha==0.3,])
  
  #fwelnet
  fname <- paste(pathResults,"CVFarkasfwENalpha",b,j,sep="")
  load(fname)
  Summdf$Method <- "fwen.groups"
  df$Method <- "fwen.groups"
  Summdf$B <- b
  df$B <- b
  Summdf2 <- rbind(Summdf2,Summdf)
  df2 <- rbind(df2,df)
  
  #ipflasso 
  fname <- paste(pathResults,"CVFarkasipfalpha",b,j,sep="")
  load(fname)
  Summdf$B <- b
  df$B <- b
  Summdf2 <- rbind(Summdf2,Summdf)
  df2 <- rbind(df2,df)
  
  #ipflasso 
  fname <- paste(pathResults,"CVFarkasipfLargealpha",b,j,sep="")
  load(fname)
  Summdf$B <- b
  df$B <- b
  Summdf$Method <- "ipflasso2"
  df$Method <- "ipflasso2"
  Summdf2 <- rbind(Summdf2,Summdf)
  df2 <- rbind(df2,df)
  
  #gren (only half of all variables used)
  j<-2 #alpha
  fname <- paste(pathResults,"CVFarkasgrenalpha",b,j,sep="")
  load(fname)
  Summdf$B <- b
  df$B <- b
  Summdf2 <- rbind(Summdf2,Summdf)
  df2 <- rbind(df2,df)
  
}
SummGrps$Group <- factor(SummGrps$Group,levels=names(GroupsetCpG)[c(1,3,4,2,6,5)],
                         labels=names(GroupsetCpG)[c(1,3,4,2,6,5)])

#gren: memory problems
#ipflasso2: memory problems


df2$Method <- factor(df2$Method,levels=unique(df2$Method)[c(6,7,8,9,10,3,4,1,5,2)],
                     labels=c("EN","fwEN","ipf","ipf2","gren",
                              "ecpcEN squeezy","squeezy (single)","squeezy (multi)",
                              "squeezy (single+reCV)","squeezy (multi+reCV)"))

df2$MethodNum <- as.factor(as.numeric(as.factor(df2$Method)))
df2$MethodBoth <- as.factor(paste(df2$MethodNum,". ",df2$Method,sep=""))

Summdf2$Method <- factor(Summdf2$Method,levels=unique(Summdf2$Method)[c(6,7,8,9,10,1,4,2,5,3)],
                         labels=c("EN","fwEN","ipf","ipf2","gren",
                                  "ecpcEN squeezy","squeezy (single)","squeezy (multi)",
                                  "squeezy (single+reCV)","squeezy (multi+reCV)"))
Summdf2$MethodNum <- as.factor(as.numeric(as.factor(Summdf2$Method)))
temp <- paste(Summdf2$MethodNum,". ",Summdf2$Method,sep="")
Summdf2$MethodBoth <- factor(temp,levels=unique(temp)[c(6,7,8,9,10,1,4,2,5,3)])



#Plot AUC for different methods----

#Summdf2$reCV <- as.factor(Summdf2$reCV)

figname<-paste(pathFigures,"FigFarkasAUC.pdf",sep="")
pdf(width = wdthpdf*1.2, height = hghtpdf,
    file = figname)
p1<-ggplot(Summdf2)+
  geom_line(aes(x=Alpha,y=AUC,col=MethodBoth),
            linetype=2,size=ls,alpha=0.6)+
  geom_point(aes(x=Alpha,y=AUC,col=MethodBoth),size=ps)+
  geom_text_repel(aes(x=Alpha,y=AUC,label=MethodNum),
                  min.segment.length = 0.1,size=ts/3)+
  #geom_line(aes(x=Alpha,y=AUC,col=Method,linetype=reCV))+
  scale_color_manual(values=colsAUC,name="Method")+
  coord_cartesian(ylim=lims)+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
p1
dev.off()
p1

#Table average computing times+sd, AUC+sd and number of selected variables+sd----
#Time
timeTable <- df2[df2$Alpha!=0,] %>% group_by(Method) %>% 
  summarise(meanTime=mean(Time),sdTime=sd(Time)) %>% ungroup()
# library("writexl")
# fname <- paste(pathFigures,"TimeTableRawFarkas.xlsx")
# write_xlsx(timeTable,fname)

#AUC
AUCtable <- Summdf2[Summdf2$Alpha!=0,] %>% group_by(Method) %>% 
  summarise(meanAUC=mean(AUC),sdAUC=sd(AUC)) %>% ungroup()
# fname <- paste(pathFigures,"AUCTableRawFarkas.xlsx")
# write_xlsx(AUCtable,fname)

#Number of variables
varsTable <- df2[df2$Alpha!=0,] %>% group_by(Method) %>% 
  summarise(meanVars=mean(NumberSelectedVars),sdVars=sd(NumberSelectedVars)) %>% ungroup()
# fname <- paste(pathFigures,"VarsTableRawFarkas.xlsx")
# write_xlsx(varsTable,fname)

#Cross-validated likelihood
CVLLTable <- Summdf2[Summdf2$Alpha!=0,] %>% group_by(Method) %>% 
  summarise(meanCVLL=mean(CVLL),sdAUC=sd(CVLL)) %>% ungroup()
# fname <- paste(pathFigures,"CVLLTableRawFarkas.xlsx")
# write_xlsx(CVLLTable,fname)


#Plot qq-plot normality check----
#data have been centered, remove first row for check
fname <- paste(pathResults,"FarkasResAll",sep="")
load(fname)
p1 <- normalityCheckQQ(Xstd[-1,], groupset=GroupsetsAll[[1]], fit.squeezy = res.squeezy)

figname<-paste(pathFigures,"FigFarkasQQ.pdf",sep="")
pdf(width = wdthpdf, height = hghtpdf,
    file = figname)
p1<-p1+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
p1
dev.off()
p1

#Plot group weights----
figname<-paste(pathFigures,"FigFarkasgroupweights.pdf",sep="")
pdf(width = wdthpdf, height = hghtpdf,
    file = figname)
p2<-ggplot(SummGrps)+aes(x=Group,y=TausApprox)+
  geom_jitter(size=ps,alpha=0.6,height=0)+
  #geom_boxplot()+
  labs(y="Inverse elastic net group penalties")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
p2
dev.off()
p2

