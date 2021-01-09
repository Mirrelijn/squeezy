##miRNA expression data (n=88, p=2114)
#NOTE: in Rstudio; press ALT+o for overview of code section headings

#Load libraries and set paths----
setwd("C:\\Users\\VNOB-0728\\Documents\\Server\\Analysis squeezy") 
pathData <- "" #path to data and co-data
library(ecpc)
pathResults <- "./miRNAResults/" #results are saved in Folder Results (should exist in working directory), or set to "" to save results in working directory
pathFigures <- "./miRNAFigures/"
# load libraries required for ecpc
library(MASS)
library(penalized)
library(glmnet)
library(mvtnorm)
library(gglasso)
library(mgcv)
library(CVXR)
library(GRridge)
library(randomForest)
library(expm)
library(Rsolnp)
library(foreach)
library(doParallel)
library(graper)

library(fwelnet)
library(ipflasso)
library(gren)

source("squeezy.R")
source("multiridge.R")

#load libraries needed for storing and plotting results
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggraph)
library(igraph)
library(RColorBrewer)
library(scales)

#optional: 
if(0){ #set to 1 to setup parallel backend to use many processors
  cores=detectCores()
  cl <- makeCluster(7) #not to overload your computer
  registerDoParallel(cl)
}

#Load data----
load(paste(pathData,"forMagnusN88.Rdata",sep=""))

means <- apply(as.matrix(mirnormcen_resp),1,mean) #abundance: average gene expression
Xcen <- as.matrix(mirnormcen_resp)-means
dim(Xcen) #centered data
sds <- apply(Xcen,1,sd) #standard deviations
Xstd<- t(Xcen/sds) #standardized gene expression

n <- dim(Xstd)[1] #number of samples
p <- dim(Xstd)[2] #number of covariates

Y <- resp
Y<-as.numeric(Y)-1 #response: numeric values of 0/1

#Define co-data groupings and types of hypershrinkage----
#Groupings based on FDR (missing values are miRNAs with FDR>=0.5)
codata <- read.delim(paste(pathData,"results_all.txt",sep=""))
#Needs some cleaning first as names codata miRs and data miRs are not exactly the same
#First remove double miRNAs in co-data, set BFDRs to geometric average
namesCodata <- paste(as.character(codata$miRNA)," ",sep="") #names of miRs in codata, add extra space to find less ambiguous
namesmiR <- colnames(Xstd) #names miRs in data
whichFDR <- sapply(1:length(namesCodata),function(x)grep(namesCodata[x],namesmiR) ) #data index of miRs in list
double <- which(sapply(whichFDR,function(x){length(x)>1})) #couple of miRs are twice in the data, but once in the codata
codata[double,c("miRNA","precursor","BFDR_MNminM","BFDR_PNminP")]
namesmiR[unlist(whichFDR[double])] #8 double miR in data 
table(namesCodata)[table(namesCodata)>1]
codata2 <- codata[-double,c("miRNA","BFDR_MNminM","BFDR_PNminP")]
for(i in 1:4){
  doublename <- codata$miRNA[double[i]]
  temp<-codata[codata$miRNA==doublename,c("miRNA","BFDR_MNminM","BFDR_PNminP")]
  temp[1,c(2,3)] <- sqrt(temp[1,c(2,3)]*temp[2,c(2,3)]) #take geometric average of the double miRNAs
  codata2 <- rbind(codata2,temp[1,])
}
#then match codata names with data names
namesCodata <- paste(as.character(codata2$miRNA)," ",sep="") #names of miRs in codata, add extra space to find less ambiguous
whichFDR <- sapply(1:length(namesCodata),function(x)grep(namesCodata[x],namesmiR) ) #data index of miRs in list
double <- which(sapply(whichFDR,function(x){length(x)>1})) #couple of miRs are twice in the data, but once in the codata
codata2[double,]
namesmiR[unlist(whichFDR[double])] #8 miR in data are double 
table(namesCodata)[table(namesCodata)>1] #check: no double miRNAs in codata
FDRselected <- unlist(whichFDR) #data index of miRs in vector
codataInd <- unlist(sapply(1:length(whichFDR),function(x) rep(x,length(whichFDR[[x]])))) #codata index of miRs

#Grouping FDR1: based on metastatic versus normal (BFDR_MNminM)
#ecpc requires a discretised version of the continuous co-data;
# we use an adaptive discretisation, by using hierarchical groups of varying grid size,
# and using hierarchical lasso shrinkage on group level to select hierarchical groups.
#First create a list with the groups of covariates varying in size;
# splitMedian splits continuous co-data recursively at the median to form two new groups, 
# split="lower" splits only the lower half group
FDR1 <- rep(NaN,p)
FDR1[FDRselected] <- codata2$BFDR_MNminM[codataInd]
GroupingFDR1 <- splitMedian(values=FDR1[FDRselected],index=FDRselected,minGroupSize=50,split="lower")
GroupingFDR1 <- c(GroupingFDR1,list(which(!((1:p)%in%FDRselected)))) #add group with miRs which have no FDR (in this case, that means FDR>0.5)
#Then define the hierarchy by forming groups on group level
HierarchyFDR1 <- obtainHierarchy(GroupingFDR1) 

#Grouping FDR2: based on primary tumor tissue versus normal (BFDR_PNminP)
FDR2 <- rep(NaN,p)
FDR2[FDRselected] <- codata2$BFDR_PNminP[codataInd]

values <- FDR2
values[is.nan(values)] <- 0.5
G<-8
GroupingFDR2 <- CreatePartition(values,ngroup=G,uniform=T,decreasing = F) #8 groups (~10 samples per group parameter)

GroupingsAll <- list(GroupingFDR2)
hypershrinkage <- "none"

#fwelnet matrix
Z <- matrix(rep(0,p*G),c(p,G))
for(j in 1:G){
  Z[GroupingFDR2[[j]],j] <- 1
}
fml <- "binomial"

Z2 <- matrix(values,nrow = p)

#gren
ind.vector <- rep(1,p)
for(j in 2:G){
  ind.vector[GroupingFDR2[[j]]] <- j
}

#Define folds and variables for the CV----
logname<-paste("logCVmiRNA.txt",sep="")
write(paste(Sys.time(),logname),file=logname,append=T)

set.seed(4739) 
nfolds<-10
folds2<-produceFolds(n,nfolds,Y,balance=T)
#save(folds2,file="foldsmiRNA")
load("foldsmiRNA") #use the same folds for different methods

alp <- c(0.01,0.3,0.8,1) #alpha variable used in elastic net

#Do CV for ecpc+squeezy----
j<-1 #alpha
#load MML for same global lambda and ecpc fit
method <- "MML"
fname <- paste(pathResults,"CVmiRNAsqueezy",method,"alpha",j,sep="")
load(fname)

j<-4
ecpcinit <- T
reCV <- T
method <- "ecpcEN"
selectAIC <- T
fname <- paste(pathResults,"CVmiRNAsqueezy",method,"alpha",j,"reCV",reCV,"ecpcinit",sep="")
if(selectAIC) fname <- paste(pathResults,"CVmiRNAsqueezyAIC",method,"alpha",j,"reCV",reCV,"ecpcinit",sep="")


if(0){ #set to 1 to perform analysis or run code inside block manually
  df<-data.frame() #data frame in which predicted values are stores for each of the samples in the left-out fold
  #Res<-list() #list in which raw output of ecpc is stored (e.g. estimated regression coefficients)
  res.squeezy <- list()
  #for(i in 1:nfolds){
  finalMatrix <- foreach(i=1:nfolds, .combine=rbind,
                         .packages = c("glmnet","penalized","mvtnorm","gglasso","mgcv",
                                       "CVXR","GRridge","expm","Rsolnp","ecpc")) %dopar% {
             tic<-proc.time()[[3]]
             # Res[[i]]<-ecpc(Y[-folds2[[i]]],Xstd[-folds2[[i]],],GroupingsAll,hypershrinkage=hypershrinkage,
             #                Y2=Y[folds2[[i]]],X2=Xstd[folds2[[i]],],#lambda=lambdas[i],
             #                postselection=F)
             # Res[[i]]$timeGR <- proc.time()[[3]]-tic
             res.squeezy[[i]] <- squeezy(Y[-folds2[[i]]],Xstd[-folds2[[i]],],grouping=GroupingsAll[[1]],
                                        Y2=Y[folds2[[i]]],X2=Xstd[folds2[[i]],],
                                        model="logistic",
                                        alpha=alp[j], #elastic net parameter for which group-regularised EN has to be fit
                                        fit.ecpc=Res[[i]],ecpcinit=ecpcinit,
                                        method=method,reCV=reCV,selectAIC=selectAIC) #sigmahat found in ecpc
             #res.squeezy[[i]]$timeGR <- proc.time()[[3]]-tic
             res.squeezy[[i]]$timeGR <- proc.time()[[3]]-tic + Res[[i]]$timeGR
             
             # lp <- Xstd[folds2[[i]],]%*%res.squeezy[[i]]$betaApprox 
             # res.squeezy[[i]]$YpredApprox <- 1/(1+exp(-lp))
             
             df2<-data.frame("Ypred"=c(res.squeezy[[i]]$YpredApprox,res.squeezy[[i]]$YpredMR,
                                       Res[[i]]$Ypred,Res[[i]]$Ypredridge))
             df2$Method <- rep(c(paste(c("squeezy","multiridge"),method),"ecpc","ordinary.ridge"),each=length(folds2[[i]]))
             df2$NumberSelectedVars <- rep(c(sum(res.squeezy[[i]]$betaApprox!=0),
                                             sum(res.squeezy[[i]]$betaMR!=0),
                                             sum(Res[[i]]$beta!=0),p),each=length(folds2[[i]]))
             df2$Fold <- i
             df2$Sample <- rep(folds2[[i]],4)
             df2$Time <-  rep(c(res.squeezy[[i]]$timeGR,res.squeezy[[i]]$timeGR, Res[[i]]$timeGR,NaN),each=length(folds2[[i]]))
             df2$Truth <- rep(Y[folds2[[i]]],4)
             df2$Alpha <- rep(c(alp[j],0,0,0),each=length(folds2[[i]]))
             
             write(paste(Sys.time(),"fold",i,"of",nfolds,"done"),file=logname,append=T)
             
             list("Res"=Res,"df"=df2,"res.ecpc"=res.squeezy)
           }
  
  Res <- lapply(1:nfolds,function(i) finalMatrix[i,1][[1]][[i]])
  df2 <- lapply(1:nfolds,function(i) finalMatrix[i,2][[1]])
  df <- df2[[1]]; for(i in 2:nfolds) df <- rbind(df,df2[[i]])
  res.squeezy <- lapply(1:nfolds,function(i) finalMatrix[i,3][[1]][[i]])
  save(Res,df,res.squeezy,file=paste(pathResults,"CVmiRNAsqueezy",method,sep=""))
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
  save(df,Summdf,dfROC,dfAUC,Res,res.squeezy,file=fname) #this comparison 
}

#Do CV for ecpc+squeezy defaults----
j<-1 #alpha
fname <- paste(pathResults,"CVmiRNAsqueezyalpha",j,sep="")
load(fname)

j<-1#alpha
fname <- paste(pathResults,"CVmiRNAsqueezyalpha",j,sep="")


if(0){ #set to 1 to perform analysis or run code inside block manually
  df<-data.frame() #data frame in which predicted values are stores for each of the samples in the left-out fold
  #Res<-list() #list in which raw output of ecpc is stored (e.g. estimated regression coefficients)
  #Res2<-list() #list in which raw output of ecpc is stored (e.g. estimated regression coefficients)
  res.squeezy <- list()
  res.squeezy2 <- list()
  res.squeezy3 <- list()
  res.squeezy4 <- list()
  res.squeezy5 <- list()
  #for(i in 1:nfolds){
  finalMatrix <- foreach(i=1:nfolds, .combine=rbind,
                         .packages = c("glmnet","penalized","mvtnorm","gglasso","mgcv",
                                       "CVXR","GRridge","expm","Rsolnp","ecpc")) %dopar% {
               #fit ecpc multigroup----
               # tic<-proc.time()[[3]]
               # Res[[i]]<-ecpc(Y[-folds2[[i]]],Xstd[-folds2[[i]],],GroupingsAll,hypershrinkage=hypershrinkage,
               #                Y2=Y[folds2[[i]]],X2=Xstd[folds2[[i]],],#lambda=lambdas[i],
               #                postselection=F)
               # Res[[i]]$timeGR <- proc.time()[[3]]-tic
               
               #fit squeezy multigroup, no reCV----
               tic<-proc.time()[[3]]
               res.squeezy[[i]] <- squeezy(Y[-folds2[[i]]],Xstd[-folds2[[i]],],grouping=GroupingsAll[[1]],
                                           Y2=Y[folds2[[i]]],X2=Xstd[folds2[[i]],],
                                           model="logistic",
                                           alpha=alp[j], #elastic net parameter for which group-regularised EN has to be fit
                                           fit.ecpc=Res[[i]],
                                           method="MML",reCV=F) #sigmahat found in ecpc
               #res.squeezy[[i]]$timeGR <- proc.time()[[3]]-tic
               res.squeezy[[i]]$timeGR <- proc.time()[[3]]-tic + Res[[i]]$timeGR
               
               #fit squeezy multigroup, with reCV----
               tic <- proc.time()[[3]]
               res.squeezy2[[i]] <- squeezy(Y[-folds2[[i]]],Xstd[-folds2[[i]],],grouping=GroupingsAll[[1]],
                                           Y2=Y[folds2[[i]]],X2=Xstd[folds2[[i]],],
                                           model="logistic",
                                           alpha=alp[j], #elastic net parameter for which group-regularised EN has to be fit
                                           fit.ecpc=Res[[i]],
                                           method="MML",reCV=T) #sigmahat found in ecpc
               #res.squeezy[[i]]$timeGR <- proc.time()[[3]]-tic
               res.squeezy2[[i]]$timeGR <- proc.time()[[3]]-tic + Res[[i]]$timeGR
               
               #fit ecpcEN with reCV----
               tic <- proc.time()[[3]]
               res.squeezy3[[i]] <- squeezy(Y[-folds2[[i]]],Xstd[-folds2[[i]],],grouping=GroupingsAll[[1]],
                                            Y2=Y[folds2[[i]]],X2=Xstd[folds2[[i]],],
                                            model="logistic",
                                            alpha=alp[j], #elastic net parameter for which group-regularised EN has to be fit
                                            fit.ecpc=Res[[i]],
                                            method="ecpcEN",reCV=T) #sigmahat found in ecpc
               #res.squeezy[[i]]$timeGR <- proc.time()[[3]]-tic
               res.squeezy3[[i]]$timeGR <- proc.time()[[3]]-tic + Res[[i]]$timeGR
               
               #fit ecpc one group----
               # tic<-proc.time()[[3]]
               # Res2[[i]]<-ecpc(Y[-folds2[[i]]],Xstd[-folds2[[i]],],list(list(1:p)),hypershrinkage=hypershrinkage,
               #                Y2=Y[folds2[[i]]],X2=Xstd[folds2[[i]],],#lambda=lambdas[i],
               #                postselection=F)
               # Res2[[i]]$timeGR <- proc.time()[[3]]-tic
               
               #fit squeezy one group, no reCV----
               tic<-proc.time()[[3]]
               res.squeezy4[[i]] <- squeezy(Y[-folds2[[i]]],Xstd[-folds2[[i]],],grouping=list(1:p),
                                           Y2=Y[folds2[[i]]],X2=Xstd[folds2[[i]],],
                                           model="logistic",
                                           alpha=alp[j], #elastic net parameter for which group-regularised EN has to be fit
                                           fit.ecpc=Res2[[i]],
                                           method="MML",reCV=F) #sigmahat found in ecpc
               #res.squeezy[[i]]$timeGR <- proc.time()[[3]]-tic
               res.squeezy4[[i]]$timeGR <- proc.time()[[3]]-tic + Res2[[i]]$timeGR
               
               #fit squeezy one group, with reCV----
               tic <- proc.time()[[3]]
               res.squeezy5[[i]] <- squeezy(Y[-folds2[[i]]],Xstd[-folds2[[i]],],grouping=list(1:p),
                                            Y2=Y[folds2[[i]]],X2=Xstd[folds2[[i]],],
                                            model="logistic",
                                            alpha=alp[j], #elastic net parameter for which group-regularised EN has to be fit
                                            fit.ecpc=Res2[[i]],
                                            method="MML",reCV=T) #sigmahat found in ecpc
               #res.squeezy[[i]]$timeGR <- proc.time()[[3]]-tic
               res.squeezy5[[i]]$timeGR <- proc.time()[[3]]-tic + Res2[[i]]$timeGR
               
               # lp <- Xstd[folds2[[i]],]%*%res.squeezy[[i]]$betaApprox 
               # res.squeezy[[i]]$YpredApprox <- 1/(1+exp(-lp))
               
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
               
               list("Res"=Res,"Res2"=Res2,"df"=df2,"res.ecpc"=res.squeezy,"res.ecpc2"=res.squeezy2,
                    "res.ecpc3"=res.squeezy3,"res.ecpc4"=res.squeezy4,"res.ecpc5"=res.squeezy5)
             }
  
  Res <- lapply(1:nfolds,function(i) finalMatrix[i,1][[1]][[i]])
  Res2 <- lapply(1:nfolds,function(i) finalMatrix[i,2][[1]][[i]])
  df2 <- lapply(1:nfolds,function(i) finalMatrix[i,3][[1]])
  df <- df2[[1]]; for(i in 2:nfolds) df <- rbind(df,df2[[i]])
  res.squeezy <- lapply(1:nfolds,function(i) finalMatrix[i,4][[1]][[i]])
  res.squeezy2 <- lapply(1:nfolds,function(i) finalMatrix[i,5][[1]][[i]])
  res.squeezy3 <- lapply(1:nfolds,function(i) finalMatrix[i,6][[1]][[i]])
  res.squeezy4 <- lapply(1:nfolds,function(i) finalMatrix[i,7][[1]][[i]])
  res.squeezy5 <- lapply(1:nfolds,function(i) finalMatrix[i,8][[1]][[i]])
  save(Res,Res2,df,res.squeezy,res.squeezy2,res.squeezy3,
       res.squeezy4,res.squeezy5,file=paste(pathResults,"CVmiRNAsqueezy",sep=""))
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
  save(df,Summdf,dfROC,dfAUC,Res,Res2,df,res.squeezy,res.squeezy2,res.squeezy3,
       res.squeezy4,res.squeezy5,file=fname) #this comparison 
}

#Do CV for ecpc+squeezy for different values of alpha----
method <- "MML"
reCV <- T
j<-1 #alpha
fname <- paste(pathResults,"CVmiRNAsqueezy",method,"alpha",j,sep="")
load(fname) #load Res, res.squeezy

fname <- paste(pathResults,"CVmiRNAsqueezy",method,"seqalphareCV",reCV,sep="")
alp <- seq(0,1,length.out=100)

if(0){ #set to 1 to perform analysis or run code inside block manually
  dfAlp<-data.frame() #data frame in which predicted values are stores for each of the samples in the left-out fold
  SummdfAlp<-data.frame()
  dfGrpsAlp <- data.frame()
  #for(i in 1:nfolds){
  for(j in 1:length(alp)){
    finalMatrix <- foreach(i=1:nfolds, .combine=rbind,
                           .packages = c("glmnet","penalized","mvtnorm","gglasso","mgcv",
                                         "CVXR","GRridge","expm","Rsolnp","ecpc")) %dopar% {
         tic<-proc.time()[[3]]
         temp <- squeezy(Y[-folds2[[i]]],Xstd[-folds2[[i]],],grouping=GroupingsAll[[1]],
                                     Y2=Y[folds2[[i]]],X2=Xstd[folds2[[i]],],
                                     model="logistic",
                                     alpha=alp[j], #elastic net parameter for which group-regularised EN has to be fit
                                     lambdas = res.squeezy[[i]]$lambdaMR,
                                     lambdaglobal=res.squeezy[[i]]$lambdaglobal,
                                     method=method,reCV=reCV) #sigmahat found in ecpc
         temp$timeGR <- proc.time()[[3]]-tic
         
         
         df2<-data.frame("Ypred"=c(temp$YpredApprox))
         df2$Method <- paste("squeezy",method)
         df2$NumberSelectedVars <- sum(temp$betaApprox!=0)
         df2$Fold <- i
         df2$Sample <- folds2[[i]]
         df2$Time <-  temp$timeGR
         df2$Truth <- Y[folds2[[i]]]
         df2$Alpha <- alp[j]
         
         dfGrps2 <- data.frame(lambdaApprox = c(temp$lambdaApprox,use.names=F)) 
         dfGrps2$Alpha <- alp[j]
         dfGrps2$Fold <- i
         dfGrps2$Group <- 1:G
         
         write(paste(Sys.time(),"fold",i,"of",nfolds,"done"),file=logname,append=T)
         
         list("dfGrps"=dfGrps2,"df"=df2)
       }
    dfGrps2 <- lapply(1:nfolds,function(i) finalMatrix[i,1][[1]])
    dfGrps <- dfGrps2[[1]]; for(i in 2:nfolds) dfGrps <- rbind(dfGrps,dfGrps2[[i]])
    df2 <- lapply(1:nfolds,function(i) finalMatrix[i,2][[1]])
    df <- df2[[1]]; for(i in 2:nfolds) df <- rbind(df,df2[[i]])
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
    
    SummdfAlp <- rbind(SummdfAlp,Summdf)
    dfAlp <- rbind(dfAlp,df)
    dfGrpsAlp <- rbind(dfGrpsAlp,dfGrps)
  }
  
  save(dfAlp,SummdfAlp,dfGrpsAlp,file=fname) #this comparison 
}

#Do CV for elastic net-----
fname <- paste(pathResults,"CVmiRNAEN",sep="")

if(0){
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

#Do CV for fwelnet----
j<-4 #alpha
fname <- paste(pathResults,"CVmiRNAfwENalpha",j,sep="")


if(0){ #set to 1 to perform analysis or run code inside block manually
  df<-data.frame() #data frame in which predicted values are stores for each of the samples in the left-out fold
  res.fwEN <- list()
  #for(i in 1:nfolds){
  finalMatrix <- foreach(i=1:nfolds, .combine=rbind,
                         .packages = c("glmnet","penalized","mvtnorm","gglasso","mgcv",
                                       "CVXR","GRridge","expm","Rsolnp","fwelnet")) %dopar% {
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
  
  save(res.fwEN,df,file=paste(pathResults,"CVmiRNAfwENalpha",j,sep=""))
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

#Do CV for fwelnet continuous FDR2----
j<-4 #alpha
fname <- paste(pathResults,"CVmiRNAfwEN2alpha",j,sep="")


if(0){ #set to 1 to perform analysis or run code inside block manually
  df<-data.frame() #data frame in which predicted values are stores for each of the samples in the left-out fold
  res.fwEN <- list()
  #for(i in 1:nfolds){
  finalMatrix <- foreach(i=1:nfolds, .combine=rbind,
                         .packages = c("glmnet","penalized","mvtnorm","gglasso","mgcv",
                                       "CVXR","GRridge","expm","Rsolnp","fwelnet")) %dopar% {
           tic<-proc.time()[[3]]
           res.fwEN[[i]] <-  cv.fwelnet(x=Xstd[-folds2[[i]],],y=Y[-folds2[[i]]],z=Z2,
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
  
  save(res.fwEN,df,file=paste(pathResults,"CVmiRNAfwENalpha",j,sep=""))
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

#Do CV for ipflasso----
j<-4 #alpha
fname <- paste(pathResults,"CVmiRNAipfalpha",j,sep="")

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

if(0){ #set to 1 to perform analysis or run code inside block manually
  df<-data.frame() #data frame in which predicted values are stores for each of the samples in the left-out fold
  #res.ipf <- list()
  #for(i in 1:3){
    finalMatrix <- foreach(i=1:nfolds, .combine=rbind,
                           .packages = c("glmnet","penalized","mvtnorm","gglasso","mgcv",
                                         "CVXR","GRridge","expm","Rsolnp","ipflasso")) %dopar% {
    tic<-proc.time()[[3]]
    res.ipf <- cvr2.ipflasso(X=Xstd[-folds2[[i]],],Y=Y[-folds2[[i]]],
                                  blocks=GroupingsAll[[1]],alpha=alp[j],
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
    
    #df <- rbind(df,df2)
    list("df"=df2,"df2"=df2)
  }
  
  df2 <- lapply(1:nfolds,function(i) finalMatrix[i,1][[1]])
  df <- df2[[1]]; for(i in 2:nfolds) df <- rbind(df,df2[[i]])
  #res.ipf <- lapply(1:nfolds,function(i) finalMatrix[i,2][[1]][[i]])
  save(df,file=paste(pathResults,"CVmiRNAipfalpha",j,sep=""))
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


#Do CV for gren----
j<-4 #alpha
fname <- paste(pathResults,"CVmiRNAgrenalpha",j,sep="")

#standardize in gren cannot be set to false, remove couple of genes that are only 0 in subsamples
ind.remove <- c()
for(i in 1:nfolds){
  remove <- which(is.nan(apply(scale(Xstd[-folds2[[i]],]),2,mean)))
  ind.remove <- unique(c(ind.remove,remove))
}
length(ind.remove)

if(0){ #set to 1 to perform analysis or run code inside block manually
  df<-data.frame() #data frame in which predicted values are stores for each of the samples in the left-out fold
  res.gren <- list()
  #for(i in 1:nfolds){
  finalMatrix <- foreach(i=1:nfolds, .combine=rbind,
                         .packages = c("glmnet","penalized","mvtnorm","gglasso","mgcv",
                                       "CVXR","GRridge","expm","Rsolnp","gren")) %dopar% {
           tic<-proc.time()[[3]]
           res.gren[[i]] <- gren(x=Xstd[-folds2[[i]],-ind.remove],y=Y[-folds2[[i]]],
                                 partitions=list(ind.vector[-ind.remove]),
                                 alpha=alp[j],trace=F)
           res.gren[[i]]$timeGR <- proc.time()[[3]]-tic
           
           betagren <- rep(0,p)
           betagren[-ind.remove] <- c(coef(res.gren[[i]], type="groupreg"))[-1]
           a0gren <- c(coef(res.gren[[i]], type="groupreg"))[1]
           Ypred <- 1/(1+exp(-Xstd[folds2[[i]],]%*%betagren - a0gren))
           
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
  save(df,res.gren,file=paste(pathResults,"CVmiRNAgren",j,sep=""))
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
display.brewer.all(3,colorblindFriendly=T)
colpal <- "Dark2"
colsfill <- brewer.pal(3,"Dark2")[1:2]
colsAUC <- seq_gradient_pal(low="grey50",high="white")((0:9)/10)
colsAUC <- seq_gradient_pal(low="black",high="grey50")((0:9)/10)

#Load results, given in files:----
Summdf2 <- data.frame()
df2 <- data.frame()
for(j in 1:length(alp)){
  fname <- paste(pathResults,"CVmiRNAsqueezyalpha",j,sep="")
  load(fname)
  Summdf2 <- rbind(Summdf2,Summdf)
  df2 <- rbind(df2,df)
  
  #fwelnet
  fname <- paste(pathResults,"CVmiRNAfwENalpha",j,sep="")
  load(fname)
  # Summdf$reCV <- NaN
  # df$reCV <- NaN
  Summdf$Method <- "fwen.groups"
  df$Method <- "fwen.groups"
  Summdf2 <- rbind(Summdf2,Summdf)
  df2 <- rbind(df2,df)
  
  #fwelnet, continuous
  fname <- paste(pathResults,"CVmiRNAfwEN2alpha",j,sep="")
  load(fname)
  # Summdf$reCV <- NaN
  # df$reCV <- NaN
  Summdf$Method <- "fwen.continuous"
  df$Method <- "fwen.continuous"
  Summdf2 <- rbind(Summdf2,Summdf)
  df2 <- rbind(df2,df)

  #ipflasso
  fname <- paste(pathResults,"CVmiRNAipfalpha",j,sep="")
  load(fname)
  # Summdf$reCV <- NaN
  # df$reCV <- NaN
  Summdf2 <- rbind(Summdf2,Summdf)
  df2 <- rbind(df2,df)
  
  #gren
  fname <- paste(pathResults,"CVmiRNAgrenalpha",j,sep="")
  load(fname)
  # Summdf$reCV <- NaN
  # df$reCV <- NaN
  Summdf2 <- rbind(Summdf2,Summdf)
  df2 <- rbind(df2,df)
}
#elastic net
fname <- paste(pathResults,"CVmiRNAEN",sep="")
load(fname)
# SummdfElNet$reCV <- NaN
# dfElNet$reCV <- NaN
Summdf2 <- rbind(Summdf2,SummdfElNet)
df2 <- rbind(df2,dfElNet)

df2$Method <- factor(df2$Method,levels=unique(df2$Method)[c(10,6,7,8,9,3,4,1,5,2)],
                     labels=c("EN","fwEN","fwEN (continuous)","ipf","gren",
                     "ecpcEN squeezy","squeezy (single)","squeezy (multi)",
                     "squeezy (single+reCV)","squeezy (multi+reCV)"))

df2$MethodNum <- as.factor(as.numeric(as.factor(df2$Method)))
df2$MethodBoth <- as.factor(paste(df2$MethodNum,". ",df2$Method,sep=""))

Summdf2$Method <- factor(Summdf2$Method,levels=unique(Summdf2$Method)[c(10,6,7,8,9,1,4,2,5,3)],
                     labels=c("EN","fwEN","fwEN (continuous)","ipf","gren",
                              "ecpcEN squeezy","squeezy (single)","squeezy (multi)",
                              "squeezy (single+reCV)","squeezy (multi+reCV)"))
Summdf2$MethodNum <- as.factor(as.numeric(as.factor(Summdf2$Method)))
temp <- paste(Summdf2$MethodNum,". ",Summdf2$Method,sep="")
Summdf2$MethodBoth <- factor(temp,levels=unique(temp)[c(10,6,7,8,9,1,4,2,5,3)])

#alpha vs AUC plots 
method <- "MML"
# fname <- paste(pathResults,"CVmiRNAsqueezy",method,"seqalpha",sep="")
reCV <- T
fname <- paste(pathResults,"CVmiRNAsqueezy",method,"seqalphareCV",reCV,sep="")
load(fname)
SummdfAlp$reCV <- reCV
SummdfAlp2 <- SummdfAlp
#dfGrpsAlp2 <- dfGrpsAlp

reCV <- F
fname <- paste(pathResults,"CVmiRNAsqueezy",method,"seqalphareCV",reCV,sep="")
load(fname)
SummdfAlp$reCV <- reCV
SummdfAlp2 <- rbind(SummdfAlp2,SummdfAlp)
#dfGrpsAlp <- rbind(dfGrpsAlp,dfGrps)

lims <- range(Summdf2$AUC,SummdfAlp2$AUC)

#Plot AUC for different methods----

#Summdf2$reCV <- as.factor(Summdf2$reCV)

figname<-paste(pathFigures,"FigmiRNAAUC.pdf",sep="")
pdf(width = wdthpdf*1.2, height = hghtpdf,
    file = figname)
p1<-ggplot(Summdf2)+
  geom_line(aes(x=Alpha,y=AUC,col=MethodBoth),
            linetype=2,size=ls,alpha=0.2)+
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

#Plot computing time for different methods----
summTime <- df2 %>% group_by(Method,Alpha) %>% summarise(meanTime=mean(Time)) %>% ungroup()
summTime2 <- df2 %>% group_by(Method,Alpha,Fold) %>% summarise(Time=mean(Time)) %>% ungroup()
# figname<-paste(pathFigures,"FigmiRNAtime.pdf",sep="")
# pdf(width = wdthpdf, height = hghtpdf,
#     file = figname)
# p1<-ggplot(summTime2[summTime2$Alpha!=0&!is.nan(summTime2$Time),])+
#   geom_jitter(aes(x=Alpha,y=Time,col=Method),width=0.01,height=0,alpha=0.6,size=ps)+
#   geom_line(data=summTime[summTime$Alpha!=0&!is.nan(summTime$meanTime),],
#             aes(x=Alpha,y=meanTime,col=Method),size=ls,alpha=0.4,linetype=2)+
#   scale_y_log10()+
#   theme_bw()+
#   theme(axis.text.x=element_text(size=ts),
#         axis.text.y=element_text(size=ts),
#         axis.title.x=element_text(size=ts+2),
#         axis.title.y=element_text(size=ts+2),
#         legend.text=element_text(size=ts),
#         legend.title=element_text(size=ts+2),
#         strip.text=element_text(size=ts))#,
# p1
# dev.off()
# p1

#Table average computing times + standard deviation----
timeTable <- df2[df2$Alpha!=0,] %>% group_by(MethodBoth) %>% 
  summarise(meanTime=mean(Time),sdTime=sd(Time)) %>% ungroup()
library("writexl")
fname <- paste(pathFigures,"TimeTableRaw.xlsx")
write_xlsx(timeTable,fname)

#Plot AUC for different alpha----
figname<-paste(pathFigures,"FigmiRNAAUCalp.pdf",sep="")
pdf(width = wdthpdf, height = hghtpdf,
    file = figname)
p2<-ggplot(SummdfAlp2)+
  geom_line(aes(x=Alpha,y=AUC,linetype=reCV,col=reCV),alpha=0.2,size=ls)+
  geom_point(aes(x=Alpha,y=AUC,shape=reCV,col=reCV),size=ps)+
  coord_cartesian(ylim=lims)+
  scale_linetype_manual(values=c(2,1))+
  scale_color_manual(values=c(rev(colsAUC)[1],colsAUC[1]))+
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

#Plot both auc reCV and all methods----
figname<-paste(pathFigures,"FigmiRNAAUCAll.pdf",sep="")
pdf(width = wdthpdf, height = hghtpdf,
    file = figname)
p4<-ggarrange(p2,p1,common.legend=F,nrow=1,widths=c(1.5,2),heights = c(1,1))+
  #theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
p4
dev.off()
p4
ggsave(figname,width= wdthpdf*2.5, height = hghtpdf)

#Plot AUC for different squeezy default methods----
lims <- range(Summdf2$AUC)
Summdf2$reCV <- as.factor(Summdf2$reCV)

# figname<-paste("FigmiRNAAUC.pdf",sep="")
# pdf(width = wdthpdf, height = hghtpdf,
#     file = figname
ggplot(Summdf2[Summdf2$Alpha!=0&grepl("squeezy",Summdf2$Method),])+
  geom_point(aes(x=Alpha,y=AUC,col=Method))+
  geom_line(aes(x=Alpha,y=AUC,col=Method,linetype=reCV))+
  #coord_cartesian(ylim=lims)+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
# dev.off()

#Plot number of selected variables(+intercept) for different methods----
SummVars <- df2 %>% group_by(Method,Alpha) %>% summarise(meanVars=mean(NumberSelectedVars+1)) %>% ungroup()

# figname<-paste("FigmiRNAvars.pdf",sep="")
# pdf(width = wdthpdf, height = hghtpdf,
#     file = figname)
ggplot(df2[df2$Alpha!=0,])+
  geom_point(aes(x=Alpha,y=NumberSelectedVars+1,col=Method),
             size=ps,alpha=0.6)+
  geom_line(data=SummVars[SummVars$Alpha!=0,],
            aes(x=Alpha,y=meanVars,col=Method),
            size=ls,alpha=0.4,linetype=2)+
  #coord_cartesian(ylim=c(0,750))+
  #xlim(0.2,1)+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
# dev.off()


#Plot lambdaapprox for different alpha----
# figname<-paste("FigmiRNAAUCalp.pdf",sep="")
# pdf(width = wdthpdf, height = hghtpdf,
#     file = figname
summGrps <- dfGrpsAlp %>% group_by(Alpha,Group) %>% summarise(meanLam=mean(lambdaApprox)) %>% ungroup()
dfGrpsAlp$Group <- as.factor(dfGrpsAlp$Group)
summGrps$Group <- as.factor(summGrps$Group)

ggplot(dfGrpsAlp)+
  #geom_point(aes(x=Alpha,y=lambdaApprox,col=Group))+
  geom_line(data=summGrps,aes(x=Alpha,y=meanLam,col=Group))+
  #coord_cartesian(ylim=lims)+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
# dev.off()

