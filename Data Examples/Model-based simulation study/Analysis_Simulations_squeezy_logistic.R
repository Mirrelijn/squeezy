###Simulation study - logistic regression
#NOTE: in Rstudio; press ALT+o for overview of code section headings
#This file may be used to replicate the analysis for Supplementary Section B.1.1

#make sure working directory is set to where this file is
setwd("")
pathResults <- "" #if desired, separate folder can be specified to save results..
pathFigures <- "" #..and figures

#optional: set up parallel back-end to use multiple processors
runParallel <- F

#Load libraries----
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
if(!requireNamespace("pROC")) install.packages("pROC")

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
# cores=detectCores()
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
alphaEN <- 1 #elastic net alpha variable 
n2<-300 #sample size of each test data set
Nblock <- 10     #correlation blocks in X
CorX <- 0.2 #correlation strength
q<-0 #set q*100% of the betas to 0
nSim <- 50 #number of simulated training&test data sets (reset inside for-loop)
n<-150 #number of training samples (reset inside for-loop)
p<-1200 #number of covariates (reset inside for-loop)

#data simulation settings to run (one/multiple of "1"-"3")
Settings <- as.character(1:3)

#set T/F to run the following methods
run_squeezy <- T #ecpc+squeezy, squeezy

dist <- "laplace"
model <- "logistic"

#Do simulations----
#set.seed(101) #used for initial values comparison
#setting <- "1"
for(setting in as.character(c(1,2,3))){
  #Set specific simulation setting variables----
  switch(setting,
         "1"={ #setting 1: no signal (single lambda setting)
           taugrp <- rep(0.1,5)  #single lambda setting
           sigmasq <- 1
           n<-150
           p<-1200 
           nSim <- 50
           G<-5 #number of groups for which co-data is generated
           
           set.seed(seeds[1])
         },
         "2"={ #setting 2: medium signal (multigroup setting with small differences)
           taugrp <- c(0.01,0.05,0.1,0.4,0.8)
           taugrp <- mean(taugrp) + (taugrp-mean(taugrp))/2
           sigmasq <- 1
           n<-150
           p<-1200
           nSim <- 50
           G<-5 #number of groups for which co-data is generated
           
           set.seed(seeds[2])
         },
         "3"={ #setting 3: high signal (multigroup setting with large differences)
           taugrp <- c(0.01,0.05,0.1,0.4,0.8)
           sigmasq <- 1
           n<-150
           p<-1200 
           nSim <- 50
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
  
  #Generate correlated data----
  betas <- list()
  X <- list()
  X2 <- list()
  Y <- list()
  Y2 <- list()
  sigma2_emp <- c() #empirical residual variance in simulated linear regression
  lambda_emp <- list() #empirical residual variance in simulated linear regression with lasso prior
  lambdaR_emp <- list() #empirical ridge lambda
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
      
      sigma2_emp[i] <- n/(n-1)*mean((Y[[i]] - as.numeric(X[[i]]%*%betas[[i]]))^2) #residual variance
    }else if(model=="logistic"){
      probs <- 1/(1+exp(-as.numeric(X[[i]]%*%betas[[i]])))
      probs2 <- 1/(1+exp(-as.numeric(X2[[i]]%*%betas[[i]])))
      Y[[i]] <- sapply(probs,function(x) rbinom(n=1,size=1,prob=x))
      Y2[[i]] <- sapply(probs2,function(x) rbinom(n=1,size=1,prob=x))
      
      fml <- "binomial"
      
      sigma2_emp[i] <- 1 #no variance parameter in logistic regression
      var_y <- mean(probs*(1-probs))
    }
    
    #empirical 'true' estimates of the simulations for lasso prior
    tau2_emp[[i]] <- sapply(indRank[[i]], function(x) 1/2/length(x)*sum(abs(betas[[i]][x]))) #0.5*average |\beta_j| in each group
    lambda_emp[[i]] <- sigma2_emp[i]/tau2_emp[[i]]
    lambdaR_emp[[i]] <- lambda_emp[[i]]^2/8
    dfGrps2 <- data.frame("Group"=1:G,"Lambda"=lambda_emp[[i]],"Sigma2"=rep(sigma2_emp[i],G),
                        "Method"=rep("True_emp",G),"Tau2"=tau2_emp[[i]])
    dfGrps2$Dataset <- i
    dfGrps2$n <- n
    dfGrps2$G <- G
    dfGrps2$p <- p
    dfGrps2$Setting <- setting
    dfGrps <- rbind(dfGrps,dfGrps2)
  }
  save(dfGrps,file=paste(pathResults,paste("simresSqueezy_GroupsEmpTrue",model,alphaEN,setting,sep="_"),".Rdata",sep=""))
  


  #Fit squeezy (set run_squeezy=T to run)----
  if(run_squeezy){
    logname<-paste("logSimulations.txt",sep="")
    write(paste(Sys.time(),logname),file=logname,append=T)
    
    simpara <- list(n=n,ntest=n2,p=p,G=G,CorX=CorX,Nblock=Nblock,
                    alphaEN=alphaEN, model=model, dist= dist,taugrp=taugrp,sigmasq=sigmasq)
    fname <- paste(pathResults,paste("simresSqueezy",model,alphaEN,setting,sep="_"),".Rdata",sep="")
    
    nSim2 <- nSim
    #nSim2 <- 5
    
    df <- data.frame()
    dfGrps <- data.frame()
    
    #for(i in 1:nSim2){
    finalMatrix <- foreach(i=1:nSim2, .combine=rbind,
                           .packages = c("ecpc","squeezy")) %dopar% {
             write(paste(Sys.time(),i,logname),file=logname,append=T)
             
             #Fit ecpc+squeezy defaults----
             
             #first fit ecpc to find group-ridge penalty estimates for initialisation
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
             
             #multigroup, empirical
             pmt <- proc.time()[[3]]
             res.ecpcsqueezy_emp <- squeezy(Y[[i]],X[[i]],groupset=indRank[[i]],Y2=Y2[[i]],X2=X2[[i]],
                                         model=model,intrcpt=T,
                                         alpha=alphaEN, #elastic net parameter for which group-regularised EN has to be fit
                                         method="MML", #one of c("ecpcEN","MML","MML.noDeriv","CV")
                                         lambdas=lambdaR_emp[[i]],reCV=F,lambdaglobal=mean(lambdaR_emp[[i]])) #fit of ecpc function
             ctsqueezy_emp <- proc.time()[[3]]-pmt;
             ctsqueezy_emp
             
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
                                         model=model,intrcpt=T,reCV=F,
                                         alpha=alphaEN, #elastic net parameter for which group-regularised EN has to be fit
                                         method="MML", #one of c("ecpcEN","MML","MML.noDeriv","CV")
                                         fit.ecpc=res.ecpc2) #fit of ecpc function
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
             
             #default ecpcEN
             pmt <- proc.time()[[3]]
             res.ecpcsqueezy5 <- squeezy(Y[[i]],X[[i]],groupset=indRank[[i]],Y2=Y2[[i]],X2=X2[[i]],
                                         model=model,intrcpt=T,
                                         alpha=alphaEN, #elastic net parameter for which group-regularised EN has to be fit
                                         method="ecpcEN", #one of c("ecpcEN","MML","MML.noDeriv","CV")
                                         fit.ecpc=res.ecpc) #fit of ecpc function
             ctsqueezy5 <- proc.time()[[3]]-pmt +ctecpc;
             ctsqueezy5
             
             #Brier score (MSE)
             MSEmulti <- res.ecpcsqueezy$MSEApprox #e.g. MSE on independent test set of EN-approximation
             MSEmultireCV <- res.ecpcsqueezy2$MSEApprox #e.g. MSE on independent test set of EN-approximation
             MSEsingle <- res.ecpcsqueezy3$MSEApprox #e.g. MSE on independent test set of EN-approximation
             MSEsinglereCV <- res.ecpcsqueezy4$MSEApprox #e.g. MSE on independent test set of EN-approximation
             MSEecpcEN <- res.ecpcsqueezy5$MSEApprox #e.g. MSE on independent test set of EN-approximation
             MSEecpc_emp <- res.ecpcsqueezy_emp$MSEApprox
             
             mses_ecpcsq <- c(MSEmulti,MSEmultireCV, MSEsingle, MSEsinglereCV,MSEecpcEN,MSEecpc_emp)
             MLs_ecpcsq <- c(res.ecpcsqueezy$MLfinal,res.ecpcsqueezy2$MLfinal, 
                             res.ecpcsqueezy3$MLfinal, res.ecpcsqueezy4$MLfinal,
                             res.ecpcsqueezy5$MLfinal,res.ecpcsqueezy_emp$MLfinal)
             
             #AUC
             rocpROC_multi <- pROC::roc(predictor = c(res.ecpcsqueezy$YpredApprox), 
                                       response = Y2[[i]], 
                                       smooth = F, auc = T, levels = c(0, 1), direction = "<")
             rocpROC_multireCV <- pROC::roc(predictor = c(res.ecpcsqueezy2$YpredApprox), 
                                         response = Y2[[i]], 
                                         smooth = F, auc = T, levels = c(0, 1), direction = "<")
             rocpROC_single <- pROC::roc(predictor = c(res.ecpcsqueezy3$YpredApprox), 
                                         response = Y2[[i]], 
                                         smooth = F, auc = T, levels = c(0, 1), direction = "<")
             rocpROC_singlereCV <- pROC::roc(predictor = c(res.ecpcsqueezy4$YpredApprox), 
                                         response = Y2[[i]], 
                                         smooth = F, auc = T, levels = c(0, 1), direction = "<")
             rocpROC_ecpcEN <- pROC::roc(predictor = c(res.ecpcsqueezy5$YpredApprox), 
                                         response = Y2[[i]], 
                                         smooth = F, auc = T, levels = c(0, 1), direction = "<")
             rocpROC_emp <- pROC::roc(predictor = c(res.ecpcsqueezy_emp$YpredApprox), 
                                         response = Y2[[i]], 
                                         smooth = F, auc = T, levels = c(0, 1), direction = "<")

             
             AUCs_ecpcsq <- c(rocpROC_multi$auc[1], rocpROC_multireCV$auc[1],
                              rocpROC_single$auc[1], rocpROC_singlereCV$auc[1],
                              rocpROC_ecpcEN$auc[1],rocpROC_emp$auc[1])
             
             #mean log likelihood
             LL_multi <- mean(Y2[[i]]*log(res.ecpcsqueezy$YpredApprox) + (1-Y2[[i]])*log(1-res.ecpcsqueezy$YpredApprox))
             LL_multireCV <- mean(Y2[[i]]*log(res.ecpcsqueezy2$YpredApprox) + (1-Y2[[i]])*log(1-res.ecpcsqueezy2$YpredApprox))
             LL_single <- mean(Y2[[i]]*log(res.ecpcsqueezy3$YpredApprox) + (1-Y2[[i]])*log(1-res.ecpcsqueezy3$YpredApprox))
             LL_singlereCV <- mean(Y2[[i]]*log(res.ecpcsqueezy4$YpredApprox) + (1-Y2[[i]])*log(1-res.ecpcsqueezy4$YpredApprox))
             LL_ecpcEN <- mean(Y2[[i]]*log(res.ecpcsqueezy5$YpredApprox) + (1-Y2[[i]])*log(1-res.ecpcsqueezy5$YpredApprox))
             LL_emp <- mean(Y2[[i]]*log(res.ecpcsqueezy_emp$YpredApprox) + (1-Y2[[i]])*log(1-res.ecpcsqueezy_emp$YpredApprox))
             
             LL_ecpcsq <- c(LL_multi, LL_multireCV, LL_single, LL_singlereCV, LL_ecpcEN, LL_emp)
             
             #summarise results----
             
             df2 <- data.frame("time"=c(ctsqueezy, ctsqueezy2,ctsqueezy3,ctsqueezy4,ctsqueezy5,ctsqueezy_emp),
                               "AUC"=AUCs_ecpcsq,
                               "ML"=MLs_ecpcsq,
                               "logLik"=LL_ecpcsq,
                               "Brier"=mses_ecpcsq,
                               "Method"=c("squeezy.m","squeezy.m.reCV","squeezy.s","squeezy.s.reCV","ecpcEN","true_empirical"))
             df2$Dataset <- i
             df2$n <- n
             df2$G <- G
             df2$p <- p
             
             dfGrps <- data.frame("Group"=rep(1:G,4),
                                  "Lambda"=c(res.ecpcsqueezy$lambdaApprox,res.ecpcsqueezy5$lambdaApprox,
                                             lambdaApproxTrue,res.ecpcsqueezy_emp$lambdaApprox),
                                  "Sigma2"=rep(c(res.ecpcsqueezy$sigmahat,res.ecpcsqueezy5$sigmahat,
                                                 1,1),each=G),
                                  "Method"=rep(c("squeezy.m","ecpcEN","Truth","True_emp"),each=G))
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
#display.brewer.all(3,colorblindFriendly=T)
colpal <- "Dark2"
colsfill <- brewer.pal(3,"Dark2")[1:2]
colsAUC <- brewer.pal(8,"Dark2")
colsMSE <- seq_gradient_pal(low="grey50",high="white")((0:8)/9)
colsTime <- seq_gradient_pal(low="black",high="grey50")((0:8)/9)

#load data----
n <- 150
p <- 1200
#load data for plots performance AUC & likelihood
alphaEN <- 1
df2 <- data.frame()
for(alphaEN in c(1)){
  for(setting in as.character(c(1,2,3))){
    fname <- paste(pathResults,paste("simresSqueezy",model,alphaEN,setting,sep="_"),".Rdata",sep="")
    load(fname)
    df$Setting <- setting
    df$Alpha <- alphaEN
    df2 <- rbind(df2,df)
  }
}

df2$Method <- factor(df2$Method,levels=unique(df2$Method)[c(5,3,1,4,2,6)],
                     labels=c("ecpcEN squeezy","squeezy (single)",
                              "squeezy (multi)","squeezy (single+reCV)","squeezy (multi+reCV)","Sample estimate"))

df2$Setting2 <- as.numeric(df2$Setting)
df2$Setting <- factor(df2$Setting,levels=c("1","2","3"),
                      labels=c("No groups","Weakly informative groups","Informative groups"))


#load data for plot group estimates
alphaEN <- 1
dfGrps2 <- data.frame()
n <- 150
p <- 1200
for(setting in as.character(c(1,2,3))){
  fname <- paste(pathResults,paste("simresSqueezy",model,alphaEN,setting,sep="_"),".Rdata",sep="")
  load(fname)
  dfGrps2 <- rbind(dfGrps2,dfGrps)
  
  #empirical truth
  fname <- paste(pathResults,paste("simresSqueezy_GroupsEmpTrue",model,alphaEN,setting,sep="_"),".Rdata",sep="")
  load(fname)
  dfGrps2 <- rbind(dfGrps2,dfGrps)
}
dfGrps2$Setting <- factor(dfGrps2$Setting,levels=c("1","2","3"),
                          labels=c("No groups","Weakly informative groups","Informative groups"))
dfGrps2$Method <- factor(dfGrps2$Method ,levels=unique(dfGrps2$Method)[c(2,1,4,3)],
                         labels=c("ecpcEN squeezy", "squeezy (multi)", "sample estimate", "Truth"))



#Plot: AUC prediction all with outliers----
alphaEN <- 1
figname<-paste(pathFigures,"FigSimsAUC",model,alphaEN,".pdf",sep="")
pdf(width = wdthpdf*1.5, height = hghtpdf*1.1,
    file = figname)
p1<-ggplot(df2[df2$Alpha==alphaEN,])+aes(x=Method,y=AUC)+
  geom_hline(yintercept=0.5, linetype=2)+
  geom_boxplot(aes(fill=Method))+
  facet_wrap(.~Setting,scales="free")+
  #scale_fill_manual(values=colsMSE,name="Method")+
  lims(y=c(0.4,0.85))+
  labs(y="AUC",x="Method")+
  theme_bw()+
  theme(axis.text.x=element_blank(),#element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts),
        legend.position ="bottom")#,
p1
dev.off()
p1

#Plot: log likelihood prediction all with outliers----
alphaEN <- 1
figname<-paste(pathFigures,"FigSimslogLike",model,alphaEN,".pdf",sep="")
pdf(width = wdthpdf*1.5, height = hghtpdf*1.1,
    file = figname)
p1<-ggplot(df2[df2$Alpha==alphaEN,])+aes(x=Method,y=logLik)+
  geom_hline(yintercept=log(0.5), linetype=2)+
  geom_boxplot(aes(fill=Method))+
  facet_wrap(.~Setting)+
  #scale_fill_manual(values=colsMSE,name="Method")+
  #lims(y=c(0.4,0.85))+
  labs(y="Mean log likelihood",x="Method")+
  theme_bw()+
  theme(axis.text.x=element_blank(),#element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts),
        legend.position ="bottom")#,
p1
dev.off()
p1

#Plot: group lambda estimates versus truth----

lb <- min(1/dfGrps2$Lambda[dfGrps2$Lambda<Inf])
dfGrps2$Lambdacutoff <- pmin(dfGrps2$Lambda,10^7)
figname<-paste(pathFigures,"FigSimsEstimatesLambda",model,alphaEN,".pdf",sep="")
pdf(width = wdthpdf*1.5, height = hghtpdf*1.1,
    file = figname)
p1<-ggplot(dfGrps2[dfGrps2$Method!="Truth",])+aes(x=factor(Group),y=Lambdacutoff)+
  geom_boxplot(aes(fill=Method))+
  facet_wrap(.~Setting,scales="free")+
  scale_y_log10()+
  #scale_fill_manual(values=colsMSE[1:3],name="Method")+
  labs(y=expression(lambda[g]), x="Group")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts),
        legend.position ="bottom")#,
p1
dev.off()
p1

figname<-paste(pathFigures,"FigSimsEstimatesTau",model,alphaEN,".pdf",sep="")
pdf(width = wdthpdf*1.5, height = hghtpdf*1.1,
    file = figname)
p1<-ggplot(dfGrps2[dfGrps2$Method!="Truth",])+aes(x=factor(Group),y=Tau2)+
  geom_boxplot(aes(fill=Method))+
  facet_wrap(.~Setting,scales="free")+
  #scale_y_log10()+
  #scale_fill_manual(values=colsMSE[1:3],name="Method")+
  labs(y=expression(tau[g]^2),x="Group")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts),
        legend.position ="bottom")#,
p1
dev.off()
p1