squeezy <- function(Y,X,grouping,alpha=1,model=NULL,
                      X2=NULL,Y2=NULL,unpen=NULL,intrcpt=T,
                      method=c("ecpcEN","MML","MML.noDeriv","CV"),
                      fold=10,compareMR = T,selectAIC=F,
                      fit.ecpc=NULL,
                      lambdas=NULL,lambdaglobal=NULL,
                      lambdasinit=NULL,sigmasq=NULL,
                      ecpcinit=T,SANN=F,minlam=10^-3,
                      standardise_Y=NULL,reCV=NULL,opt.sigma=NULL,
                      resultsAICboth=F){
  #Y: response
  #X: observed data
  #grouping: list with index of covariates of co-data groups
  #alpha: elastic net penalty parameter (as in glmnet)
  #X2: independent observed data
  #Y2: independent response data
  #unpen: vector with index of unpenalised variables
  #intrcpt: should intercept be included?
  #fold: number of folds used in global lambda cross-validation
  #sigmasq: (linear regression) noise level
  #method: "CV","MML.noDeriv", or "MML"
  #compareMR: T/F to return betas for multiridge estimate, and predictions for Y2 if X2 is given
  #selectAIC: T/F to compare AIC of multiridge model and ordinary ridge model. Return best one.
  
  #Set-up variables ---------------------------------------------------------------------------
  groupings <- list(grouping)
  n <- dim(X)[1] #number of samples
  p <- dim(X)[2] #number of covariates 
  if(!is.null(X2)) n2<-dim(X2)[1] #number of samples in independent data set x2 if given
  
  #settings default methods
  if(length(method)==4){
    method <- "MML"
  }
  switch(method,
         "ecpcEN"={
           if(is.null(fit.ecpc)) stop("provide ecpc fit results")
           if(is.null(lambdas)) lambdas <- fit.ecpc$sigmahat/(fit.ecpc$gamma*fit.ecpc$tauglobal)
           if(is.null(lambdaglobal)) lambdaglobal <- fit.ecpc$sigmahat/fit.ecpc$tauglobal
           if(is.null(sigmasq)) sigmasq <- fit.ecpc$sigmahat
           if(is.null(standardise_Y)) standardise_Y <- T
           if(is.null(reCV)) reCV <- T
           if(is.null(opt.sigma)) opt.sigma <- F
         },
         "MML"={
           if(!is.null(fit.ecpc)){
             if(is.null(lambdasinit)) lambdasinit <- fit.ecpc$sigmahat/(fit.ecpc$gamma*fit.ecpc$tauglobal)
             if(is.null(lambdaglobal)) lambdaglobal <- fit.ecpc$sigmahat/fit.ecpc$tauglobal
             if(is.null(sigmasq)) sigmasq <- fit.ecpc$sigmahat
           } 
           if(is.null(standardise_Y)) standardise_Y <- F
           if(is.null(opt.sigma)) opt.sigma <- T
           if(is.null(reCV)){
             reCV <- F; if(!opt.sigma) reCV <- T
           } 
           
         },
         "MML.noDeriv"={
           if(!is.null(fit.ecpc)){
             if(is.null(lambdasinit)) lambdasinit <- fit.ecpc$sigmahat/(fit.ecpc$gamma*fit.ecpc$tauglobal)
             if(is.null(lambdaglobal)) lambdaglobal <- fit.ecpc$sigmahat/fit.ecpc$tauglobal
             if(is.null(sigmasq)) sigmasq <- fit.ecpc$sigmahat
           } 
           if(is.null(standardise_Y)) standardise_Y <- F
           if(is.null(opt.sigma)) opt.sigma <- T
           if(is.null(reCV)){
             reCV <- F; if(!opt.sigma) reCV <- T
           } 
         },
         "CV"={
           if(!is.null(fit.ecpc)){
             if(is.null(lambdasinit)) lambdasinit <- fit.ecpc$sigmahat/(fit.ecpc$gamma*fit.ecpc$tauglobal)
             if(is.null(lambdaglobal)) lambdaglobal <- fit.ecpc$sigmahat/fit.ecpc$tauglobal
             if(is.null(sigmasq)) sigmasq <- fit.ecpc$sigmahat
           } 
           if(is.null(standardise_Y)) standardise_Y <- F
           if(is.null(opt.sigma)) opt.sigma <- F
           if(is.null(reCV)) reCV <- T
         },
  )
  
  if(is.null(model)){
    if(all(is.element(Y,c(0,1))) || is.factor(Y)){
      model <- "logistic" 
    } else if(all(is.numeric(Y)) & !(is.matrix(Y) && dim(Y)[2]==2)){
      model <- "linear"
    }else{
      model <- "cox"
    }
  }

  switch(model,
         'linear'={
           fml <- 'gaussian'
           sd_y <- sqrt(var(Y)*(n-1)/n)[1]
           if(standardise_Y){
             Y <- Y/sd_y
             sd_y_former <- sd_y
             sd_y <- 1
             if(!is.null(sigmasq)) sigmasq <- sigmasq/sd_y_former
           }
           if(method=="MML") minlam <- max(minlam,10^-4*var(Y))
         },
         'logistic'={
           fml <- 'binomial'
           opt.sigma <- F
           standardise_Y <- F
           sd_y <- 1 #do not standardise y in logistic setting
           sd_y_former <- sd_y
           
           #set response to numeric 0,1
           levelsY<-cbind(c(0,1),c(0,1))
           if(!all(is.element(Y,c(0,1)))){
             oldLevelsY<-levels(Y)
             levels(Y)<-c("0","1")
             Y<-as.numeric(Y)-1
             levelsY<-cbind(oldLevelsY,c(0,1))
             colnames(levelsY)<-c("Old level names","New level names")
             if(!is.null(Y2)){
               levels(Y2)<-c("0","1")
               Y2<-as.numeric(Y2)-1
             }
             if(!silent) print("Y is put in 0/1 format, see levelsY in output for new names")
           }
         },
         'cox'={
           fml <- 'cox'
           opt.sigma <- F
           standardise_Y <- F
           sd_y <- 1 #do not standardise y in cox regression setting
           sd_y_former <- sd_y
           intrcpt <- F #Cox does not use intercept
         }
  )
  # check whether unpenalised covariates are not in partition
  # and set penalty.factor of unpenalised covariates to 0 for glmnet
  penfctr <- rep(1,p) #factor=1 for penalised covariates
  if(length(unpen)>0){
    penfctr[unpen] <- 0 #factor=0 for unpenalised covariates
    if(any(unlist(groupings)%in%unpen)){
      warning("Unpenalised covariates removed from grouping")
      for(i in 1:length(groupings)){
        for(j in 1:length(groupings[[i]])){
          if(all(groupings[[i]][[j]]%in%unpen)){
            groupings[[i]][[j]] <- NULL #remove whole group if all covariates unpenalised
          }else{
            groupings[[i]][[j]] <- groupings[[i]][[j]][!(groupings[[i]][[j]]%in%unpen)]
          }
        }
      }
    }
  }
  
  G <- sapply(groupings,length) #1xm vector with G_i, number of groups in partition i
  m <- length(G) #number of partitions
  
  indGrpsGlobal <- list(1:G[1]) #global group index in case we have multiple partitions
  if(m>1){
    for(i in 2:m){
      indGrpsGlobal[[i]] <- (sum(G[1:(i-1)])+1):sum(G[1:i])
    }
  }
  Kg <- lapply(groupings,function(x)(sapply(x,length))) #m-list with G_i vector of group sizes in partition i
  #ind1<-ind
  
  #ind <- (matrix(1,G,1)%*%ind)==(1:G)#sparse matrix with ij element T if jth element in group i, otherwise F
  i<-unlist(sapply(1:sum(G),function(x){rep(x,unlist(Kg)[x])}))
  j<-unlist(unlist(groupings))
  ind <- sparseMatrix(i,j,x=1) #sparse matrix with ij element 1 if jth element in group i (global index), otherwise 0
  
  Ik <- lapply(1:m,function(i){
    x<-rep(0,sum(G))
    x[(sum(G[1:i-1])+1):sum(G[1:i])]<-1
    as.vector(x%*%ind)}) #list for each partition with px1 vector with number of groups beta_k is in
  #sparse matrix with ij element 1/Ij if beta_j in group i
  
  #make co-data matrix Z (Zt transpose of Z as in paper, with co-data matrices stacked for multiple groupings)
  Zt<-ind; 
  if(G[1]>1){
    Zt[1:G[1],]<-t(t(ind[1:G[1],])/apply(ind[1:G[1],],2,sum))
  }
  if(m>1){
    for(i in 2:m){
      if(G[i]>1){
        Zt[indGrpsGlobal[[i]],]<-t(t(ind[indGrpsGlobal[[i]],])/
                                     apply(ind[indGrpsGlobal[[i]],],2,sum))
      }
    }
  }
  
  if(dim(Zt)[2]<p) Zt <- cbind(Zt,matrix(rep(NaN,(p-dim(Zt)[2])*sum(G)),c(sum(G),p-dim(Zt)[2])))
  
  
  #Extend data to make artifical non-overlapping groups----
  Xxtnd <- do.call(cbind,lapply(groupings[[1]],function(group){t(t(X[,group])/sqrt(Ik[[1]][group]))}))
  #create new group indices for Xxtnd
  Kg2 <- c(1,Kg[[1]]) 
  G2 <- length(Kg2)-1
  groupxtnd <- lapply(2:length(Kg2),function(i){sum(Kg2[1:(i-1)]):(sum(Kg2[1:i])-1)}) #list of indices in each group
  groupxtnd2 <- unlist(sapply(1:G2,function(x){rep(x,Kg2[x+1])})) #vector with group number
 
  Xunpen <- NULL
  if(sum((1:p)%in%unpen)>0) Xunpen<-X[,(1:p)%in%unpen] #seperate matrix for unpenalised variables
  
  #datablockNo <- groupxtnd2
  #datablocks <- groupxtnd (or grouping if no overlapping groups)
  
  #datablocks: list with each element a data type containing indices of covariates with that data type 
  Xbl <- createXblocks(lapply(groupxtnd,function(ind) Xxtnd[,ind]))
  XXbl <- createXXblocks(lapply(groupxtnd,function(ind) Xxtnd[,ind]))
  
  #Find global lambda if not given for initial penalty and/or for AIC comparison----
  if(is.null(lambdaglobal)){
    #find initial global lambda
    cvperblock <- fastCV(list(Xxtnd),Y=Y,kfold=fold,fixedfolds = F,X1=Xunpen,intercept=intrcpt)
    lambda <- cvperblock$lambdas
    lambda[lambda==Inf] <- 10^6
  }else{
    lambda <- lambdaglobal
  }
    
  #Further optimise global lambda with same method as multi-group for AIC comparison
  if(selectAIC){
    if(method=="CV"){
      leftout <- CVfolds(Y=Y,kfold=fold,nrepeat=3,fixedfolds = F) #Create (repeated) CV-splits of the data
      lambda1groupfit <- optLambdasWrap(penaltiesinit=lambda, 
                                     XXblocks=list(apply(simplify2array(XXbl),c(1,2),sum)),
                                     Y=Y,folds=leftout,
                                     X1=Xunpen,intercept=intrcpt,
                                     score=ifelse(model == "linear", "mse", "loglik"),model=model,
                                     maxItropt2 = 500,reltol = 10^-3)
      lambda <- lambda1groupfit$optpen
    }else if(method=="MML.noDeriv"){
      lambda1groupfit <- optLambdas_mgcv(penaltiesinit=lambda,Y=Y,
                                         XXblocks=list(apply(simplify2array(XXbl),c(1,2),sum)),
                                         model=model,reltol=1e-3,
                                         maxItropt=500,tracescore=F,fixedseed =F,
                                         optmethod = "Nelder-Mead",
                                         sigmasq=sigmasq) #TD: intercept?
      lambda <- lambda1groupfit$optpen
    }else if(method=="MML"){
      sigmahat<-1 
      if(!is.null(sigmasq) & model=="linear") sigmahat <- sigmasq
      
      if(opt.sigma){
        lambda1groupfit <- optim(par=c(log(sigmahat),log(lambda)), fn=minML.LA.ridgeGLM,
                                 gr=dminML.LA.ridgeGLM, method="BFGS",minlam=minlam,
                                 XXblocks = list(apply(simplify2array(XXbl),c(1,2),sum)) ,
                                 Y=Y,sigmasq=sigmahat,model=model,
                                 intrcpt=intrcpt, Xunpen=Xunpen,opt.sigma=opt.sigma)
        sigmasq <- exp(lambda1groupfit$par[1])+minlam
        lambda <- exp(lambda1groupfit$par[-1])+minlam
      }else{
        lambda1groupfit <- optim(par=log(lambda), fn=minML.LA.ridgeGLM,
                                 gr=dminML.LA.ridgeGLM, method="BFGS",minlam=minlam,
                                 XXblocks = list(apply(simplify2array(XXbl),c(1,2),sum)) ,
                                 Y=Y,sigmasq=sigmahat,model=model,
                                 intrcpt=intrcpt, Xunpen=Xunpen,opt.sigma=opt.sigma)
        lambda <- exp(lambda1groupfit$par)+minlam
      }
      
    }
  }
  
  #find estimate for sigma for 1 group for AIC comparison
  sigmahat <- 1
  if(model=="linear"){
    if(!is.null(sigmasq)) sigmahat <- sigmasq
    else{
      XXT1grp <- SigmaFromBlocks(XXblocks = list(apply(simplify2array(XXbl),c(1,2),sum)),lambda)
      if(length(unpen)>0 | intrcpt){
        Xunpen2 <- Xunpen
        if(intrcpt) Xunpen2 <- cbind(Xunpen,rep(1,n))
        if(intrcpt && length(unpen)==0){
          betaunpenML1grp <- sum(Y)/n
        }else{
          temp <- solve(XXT1grp+diag(rep(1,n)),Xunpen2)
          betaunpenML1grp <- solve(t(Xunpen2)%*%temp , t(temp)%*%Y)
        }
        sigmahat <- c(t(Y-Xunpen2%*%betaunpenML1grp)%*%
                              solve(XXT1grp+diag(rep(1,n)),Y-Xunpen2%*%betaunpenML1grp)/n)
      }else{
        sigmahat <- c(t(Y)%*%solve(XXT1grp+diag(rep(1,n)),Y)/n)
      }
    }
  }
  if(selectAIC){
    sigmahat1group <- sigmahat
  }
  
  #Compute ridge penalties and corresponding group variance estimates----
  if(is.null(lambdas)){
    #If not given, set initial lambdas to global lambda
    if(is.null(lambdasinit)|!ecpcinit){
      lambdasinit <- rep(lambda,G)
    }
    if(any(lambdasinit==Inf)){
      lambdasinit[lambdasinit>2*lambda] <- 2*lambda #truncate at 2 times the global lambda
    }
    
    #Find joint lambdas:
    if(method=="CV"){
      leftout <- CVfolds(Y=Y,kfold=fold,nrepeat=3,fixedfolds = F) #Create (repeated) CV-splits of the data
      jointlambdas <- optLambdasWrap(penaltiesinit=lambdasinit, XXblocks=XXbl,Y=Y,folds=leftout,
                                     X1=Xunpen,intercept=intrcpt,
                                     score=ifelse(model == "linear", "mse", "loglik"),model=model,
                                     maxItropt2 = 500,reltol = 10^-3,traceCV=F)
      lambdas <- jointlambdas$optpen
    }else if(method=="MML.noDeriv"){
      #browser()
      if(ecpcinit){
        if(SANN){
          jointlambdas <- optLambdas_mgcvWrap(penaltiesinit=lambdasinit, XXblocks=XXbl,Y=Y,
                                              model=model,reltol=1e-4,
                                              maxItropt2=1000,tracescore=F,fixedseed =F,
                                              optmethod2 = "Nelder-Mead",
                                              sigmasq=sigmahat,opt.sigma=opt.sigma) #TD: intercept?
        }else{
          jointlambdas <- optLambdas_mgcv(penaltiesinit=lambdasinit, XXblocks=XXbl,Y=Y,
                                          model=model,reltol=1e-4,
                                          maxItropt=1000,tracescore=F,fixedseed =F,
                                          optmethod = "Nelder-Mead",
                                          sigmasq=sigmahat,opt.sigma=opt.sigma) #TD: intercept?
        }
        
      }else{
        if(SANN){
          jointlambdas <- optLambdas_mgcvWrap(penaltiesinit=rep(lambda,G), XXblocks=XXbl,Y=Y,
                                              model=model,reltol=1e-4,
                                              maxItropt2=1000,tracescore=F,fixedseed =F,
                                              optmethod2 = "Nelder-Mead",
                                              sigmasq=sigmahat,opt.sigma=opt.sigma) #TD: intercept?
        }else{
          jointlambdas <- optLambdas_mgcv(penaltiesinit=rep(lambda,G), XXblocks=XXbl,Y=Y,
                                          model=model,reltol=1e-4,
                                          maxItropt=1000,tracescore=F,fixedseed =F,
                                          optmethod = "Nelder-Mead",
                                          sigmasq=sigmahat,opt.sigma=opt.sigma) #TD: intercept?
        }
      }
      
      # print("-log(ML) for initial penalties")
      # if(ecpcinit){
      #   MLinit <- mgcv_lambda(penalties=lambdasinit, XXblocks=XXbl,Y=Y, model=model, sigmasq = sigmasq) #-log(ml) for ecpc
      # }else{
      #   MLinit <- mgcv_lambda(penalties=rep(lambda,G), XXblocks=XXbl,Y=Y, model=model, sigmasq = sigmasq) #-log(ml) for global
      # }

      if(opt.sigma){
        sigmasq <- jointlambdas$optpen[1]
        lambdas <- jointlambdas$optpen[-1]
      }else{
        lambdas <- jointlambdas$optpen
      }
      # print("-log(ML) for final penalties")
      # MLfinal <- mgcv_lambda(penalties=lambdas, XXblocks=XXbl,Y=Y, model=model, sigmasq = sigmasq) #-log(ml) for ecpc

      
    }else if(method=="MML"){
      if(opt.sigma){
        jointlambdas <- optim(par=c(log(sigmahat),log(lambdasinit)), 
                              fn=minML.LA.ridgeGLM, 
                              gr=dminML.LA.ridgeGLM, method="BFGS",minlam=minlam,
                              XXblocks = XXbl , Y=Y,opt.sigma=opt.sigma,model=model,
                              intrcpt=intrcpt, Xunpen=Xunpen)
        sigmasq <- exp(jointlambdas$par[1])+minlam
        lambdas <- exp(jointlambdas$par[-1])+minlam
        
        # MLinit <- minML.LA.ridgeGLM(loglambdas=c(log(sigmahat),log(lambdasinit)),opt.sigma = T,
        #                             XXblocks = XXbl , Y=Y,model=model,intrcpt=intrcpt,minlam=0)
        # MLfinal <- minML.LA.ridgeGLM(loglambdas=c(log(sigmasq),log(lambdas)),
        #                              opt.sigma = T,
        #                              XXblocks = XXbl , Y=Y,model=model,intrcpt=intrcpt,minlam=0)
        
      }else{
        jointlambdas <- optim(par=log(lambdasinit), fn=minML.LA.ridgeGLM, 
                              gr=dminML.LA.ridgeGLM, method="BFGS",minlam=minlam,
                              XXblocks = XXbl , Y=Y,sigmasq=sigmahat,model=model,
                              intrcpt=intrcpt, Xunpen=Xunpen)
        lambdas <- exp(jointlambdas$par)+minlam
        
        # MLinit <- minML.LA.ridgeGLM(loglambdas=c(log(lambdasinit)),opt.sigma = F,
        #                             sigmasq=sigmahat,
        #                             XXblocks = XXbl , Y=Y,model=model,intrcpt=intrcpt,minlam=0)
        # MLfinal <- minML.LA.ridgeGLM(loglambdas=c(log(lambdas)),
        #                              opt.sigma = F,sigmasq=sigmasq,
        #                              XXblocks = XXbl , Y=Y,model=model,intrcpt=intrcpt,minlam=0)
      }
        
      
    }
  }else{
    lambdasinit <- lambdas
  }
  
  sigmahat <- 1 #sigma not in model for logistic: set to 1
  if(model=="linear"){
    if(!is.null(sigmasq)) sigmahat <- sigmasq
    else{
      XXT <- SigmaFromBlocks(XXblocks = XXbl,lambdas)
      if(length(unpen)>0 | intrcpt){
        Xunpen2 <- Xunpen
        if(intrcpt) Xunpen2 <- cbind(Xunpen,rep(1,n))
        if(intrcpt && length(unpen)==0){
          betaunpenML <- sum(Y)/n
        }else{
          temp <- solve(XXT+diag(rep(1,n)),Xunpen2)
          betaunpenML <- solve(t(Xunpen2)%*%temp , t(temp)%*%Y)
        }
        sigmahat <- c(t(Y-Xunpen2%*%betaunpenML)%*%solve(XXT+diag(rep(1,n)),Y-Xunpen2%*%betaunpenML)/n)
      }else{
        sigmahat <- c(t(Y)%*%solve(XXT+diag(rep(1,n)),Y)/n)
      }
    }
  }
  
  MLinit <- minML.LA.ridgeGLM(loglambdas=log(lambdasinit),opt.sigma = F,sigmasq=sigmahat,
                              XXblocks = XXbl , Y=Y,model=model,intrcpt=intrcpt,minlam=0)
  MLfinal <- minML.LA.ridgeGLM(loglambdas=log(lambdas),
                               opt.sigma = F, sigmasq=sigmahat,
                               XXblocks = XXbl , Y=Y,model=model,intrcpt=intrcpt,minlam=0)
  # print("-log(ML) for initial penalties"); print(MLinit)
  # print("-log(ML) for final penalties"); print(MLfinal)
  
  #May compare 1 group versus multiple groups with AIC
  if(selectAIC){
    # res1group <- squeezy(Y=Y,X=X,grouping=list(1:p),alpha=0,model=model,
    #                        unpen=unpen,intrcpt=intrcpt,
    #                        fold=fold,sigmasq=sigmasq,method=method,
    #                        compareMR = F,selectAIC=F)
    # 
    # lambda1group <- res1group$lambdaMR
    # sigmahat1group <- res1group$sigmahat
    lambda1group <- lambda
    if(sum((1:p)%in%unpen)>0){
      AIC1group <- AIC.LA.ridgeGLM(log(lambda1group),XXblocks=list(apply(simplify2array(XXbl),c(1,2),sum)),
                                   Y=Y,sigmasq=sigmahat1group,Xunpen=Xunpen,intrcpt=intrcpt,model=model)
      AICmultigroup <- AIC.LA.ridgeGLM(log(lambdas),XXblocks=XXbl,
                                       Y=Y,sigmasq=sigmahat,Xunpen=Xunpen,intrcpt=intrcpt,model=model)
    }else{
      AIC1group <- AIC.LA.ridgeGLM(log(lambda1group),
                                   XXblocks=list(apply(simplify2array(XXbl),c(1,2),sum)),
                                   Y=Y,sigmasq=sigmahat1group,intrcpt=intrcpt,model=model)
      AICmultigroup <- AIC.LA.ridgeGLM(log(lambdas),XXblocks=XXbl,
                                       Y=Y,sigmasq=sigmahat,intrcpt=intrcpt,model=model)
    }
    #If model with only one group has lower AIC, select that model
    if(AIC1group <= AICmultigroup){
      lambdasNotOptimalAIC <- lambdas
      sigmahatNotOptimalAIC <- sigmahat
      lambdas <- rep(lambda1group,G)
      sigmahat <- sigmahat1group
      modelbestAIC <- "onegroup"
    }else{
      lambdasNotOptimalAIC <- lambda1group
      sigmahatNotOptimalAIC <- sigmahat1group
      modelbestAIC <- "multigroup"
    }
  }
  
  tauglobal<- sigmahat/lambda #set global group variance
  gamma <- lambda/lambdas #= (sigmahat/lambdas)/tauglobal
  lambdap<- lambda/(as.vector(gamma%*%Zt)) #=sigmahat/(tauglobal*as.vector(gamma%*%Zt)) #specific ridge penalty for beta_k
  lambdap[lambdap<0]<-Inf 
  
  
  #Compute regression model----
  if(compareMR||alpha==0){
    #Compute ridge betas with multiridge
    XXT <- SigmaFromBlocks(XXbl,penalties=lambdas) #create nxn Sigma matrix = sum_b [lambda_b)^{-1} X_b %*% t(X_b)]
    if(sum((1:p)%in%unpen)>0){
      fit <- IWLSridge(XXT,Y=Y, model=model,intercept=intrcpt,X1=X[,(1:p)%in%unpen]) #Fit. fit$etas contains the n linear predictors
    }else{
      fit <- IWLSridge(XXT,Y=Y, model=model,intercept=intrcpt) #Fit. fit$etas contains the n linear predictors
    }
    betafit <- betasout(fit, Xblocks=Xbl, penalties=lambdas) #Find betas.
    a0MR <- 0
    if(intrcpt) a0MR <- c(betafit[[1]][1]) #intercept
    betaMR <- rep(0,p) 
    betaMR[(1:p)%in%unpen] <- betafit[[1]][-1] #unpenalised variables
    for(i in 1:length(groupings[[1]])){
      betaMR[groupings[[1]][[i]]] <- betaMR[groupings[[1]][[i]]] + betafit[[1+i]]/sqrt(Ik[[1]][groupings[[1]][[i]]])
    }  
    rm(betafit)
  }
  
  #Transform group variance estimates to desired group prior estimates----
  
  pen <- which(!((1:p)%in%unpen))
  if(any(is.nan(sqrt(lambdap[pen])))){browser()}
  
  # #Compute ridge betas with glmnet
  # 
  # penfctr2 <- penfctr
  # penfctr2[pen] <- lambdap[pen]
  # not0 <- which(penfctr2!=Inf)
  # lambdaoverall <- sum(penfctr2[not0])/length(penfctr2[not0]) #global penalty parameter
  # penfctr2 <- penfctr2/lambdaoverall #scale penalty factor such that sums to p
  # #all(penfctr2*lambdaoverall==lambdap)
  # 
  # Xacc <- X
  # Xacc[,pen] <- as.matrix(X[,pen] %*% sparseMatrix(i=1:length(pen),j=1:length(pen),
  #                                                  x=c(1/sqrt(penfctr2[pen]))))
  # if(model=="cox"){
  #   glmGR <- glmnet(Xacc[,not0],Y,alpha=0,#lambda = lambdaoverall/n*sd_y,
  #                   family=fml,standardize = F,
  #                   penalty.factor=penfctr[not0], thresh=10^-10)
  #   temp <- coef(glmGR,s=lambdaoverall/n*sd_y,exact=T,x=Xacc[,not0],y=Y,
  #                penalty.factor=penfctr[not0],family=fml,
  #                thresh=10^-10)
  # }else{
  #   glmGR <- glmnet(Xacc[,not0],Y,alpha=0,#lambda = lambdaoverall/n*sd_y,
  #                   family=fml, intercept = intrcpt, standardize = F,
  #                   penalty.factor=penfctr[not0], thresh=10^-10)
  #   #betaGLM <- rep(0,p); betaGLM[not0] <- as.vector(glmGR$beta) 
  #   temp <- coef(glmGR,s=lambdaoverall/n*sd_y,exact=T,x=Xacc[,not0],y=Y,
  #                penalty.factor=penfctr[not0],family=fml,
  #                intercept=intrcpt, thresh=10^-10)
  # }
  # betaGLM <- rep(0,p); betaGLM[not0] <- temp[-1]
  # betaGLM[pen] <- c(1/sqrt(penfctr2[pen])) * betaGLM[pen]
  # a0GLM <- temp[1]
  # 
  # plot(betaMR,betaGLM); abline(0,1)
  
  # #with penfctr and glmnet: does not always converge
  # penfctr2 <- penfctr
  # penfctr2[pen] <- lambdap[pen]
  # not0 <- which(penfctr2!=Inf)
  # lambdaoverall <- sum(penfctr2[not0])/length(penfctr2[not0]) #global penalty parameter
  # penfctr2 <- penfctr2/lambdaoverall #scale penalty factor such that sums to p
  # if(model=="cox"){
  #   glmGR <- glmnet(X[,not0],Y,alpha=0,family=fml,standardize = F,
  #                   penalty.factor=penfctr2[not0], thresh=10^-10)
  #   temp <- coef(glmGR,s=lambdaoverall/n*sd_y,exact=T,x=X[,not0],y=Y,
  #                penalty.factor=penfctr2[not0],family=fml,
  #                thresh=10^-10)
  # }else{
  #   glmGR <- glmnet(X[,not0],Y,alpha=0,family=fml,#lambda=lambdaoverall/n*sd_y,
  #                   intercept = intrcpt, standardize = F,
  #                   penalty.factor=penfctr2[not0], thresh=10^-10)
  #   #betaGLM2 <- rep(0,p); betaGLM2[not0] <- as.vector(glmGR$beta)
  #   
  #   temp <- coef(glmGR,s=lambdaoverall/n*sd_y,exact=T,x=X[,not0],y=Y,
  #                penalty.factor=penfctr2[not0],family=fml,
  #                intercept=intrcpt, thresh=10^-12)
  # }
  # betaGLM2 <- rep(0,p); betaGLM2[not0] <- temp[-1]
  # plot(betaGLM2,betaGLM); abline(0,1); 
  # points(betaGLM2[groupings[[1]][[1]]],betaGLM[groupings[[1]][[1]]],col="red")
  
  #Transform penalties and compute elastic net betas with glmnet
  if(alpha <= 1){
    varFunc <- function(lam,alpha=0.5,tausq){
      t1 <- - alpha/(1-alpha)^(3/2)/sqrt(lam)*
        exp(dnorm(alpha*sqrt(lam)/sqrt(1-alpha),log=T) -
              pnorm(-alpha*sqrt(lam)/sqrt(1-alpha),log.p=T))
      varBeta <- 1/lam/(1-alpha) + alpha^2/(1-alpha)^2 + t1
      f <- varBeta - tausq
      return(f)
    }
    lamEN <- function(alpha,tausq){
      if(alpha==0){
        lamEN <- 1/tausq
      }else if(alpha==1){
        lamEN <- sqrt(2/tausq)
      }else if(tausq>10^6){
        lamEN <- 1/tausq
        if(alpha>0.2){
          lb2 <- 10^-6
          ub2 <- sqrt(2/10^6)
          temp <- try(uniroot(varFunc,c(0.9*lb2,1.1*ub2),alpha=alpha,tausq=10^6,tol=10^-6))
          if(class(temp)[1]!="try-error"){
            if(temp$root<lamEN) lamEN <- temp$root
          }
        }
      }else if(tausq<10^-6){
        lamEN <- sqrt(2/tausq)
        if(alpha<0.8){
          lb2 <- sqrt(2/10^-7)
          ub2 <- 10^7
          temp <- try(uniroot(varFunc,c(0.9*lb2,1.1*ub2),alpha=alpha,tausq=10^-7,tol=10^-6))
          if(class(temp)[1]!="try-error"){
            if(temp$root>lamEN) lamEN <- temp$root
          }
        }
      }else{
        lb <- min(1/tausq,sqrt(2/tausq))
        ub <- max(1/tausq,sqrt(2/tausq))
        #if(sign(varFunc(0.9*lb,alpha=alpha,tausq=tausq))==sign(varFunc(1.1*ub,alpha=alpha,tausq=tausq))) browser()
        temp <- try(uniroot(varFunc,c(0.9*lb,1.1*ub),alpha=alpha,tausq=tausq,tol=10^-6))
        if(class(temp)[1]=="try-error"){
          if(tausq<0.5) lamEN <- sqrt(2/tausq)
          if(tausq>0.5) lamEN <- 1/tausq
        }else{
          lamEN <- temp$root
        }
      }
      return(lamEN)
    }
    alphat <- 1/(1+sd_y*(1-alpha)/alpha)
    lambdasEN <- sapply(sigmahat/lambdas,function(tau) lamEN(tausq = tau,alpha=alpha))
    tauEN <- 1/lambdasEN
    lambdasENhat <- lambdasEN*sigmahat*((1-alpha)+ alpha/sd_y)
    
    uniqueTaus <- unique(as.vector(c(sigmahat/lambdas)%*%Zt))
    if(dim(X)[2]==dim(Xxtnd)[2]){
      lambdap<- (as.vector(lambdasENhat%*%Zt)) 
    }else{ #some groups are overlapping
      lambdasENunique <- sapply(uniqueTaus,function(tau) lamEN(tausq = tau,alpha=alpha))
      lambdasENuniquehat <- lambdasENunique*sigmahat*((1-alpha)+ alpha/sd_y)
      indUnique <- sapply(as.vector(c(sigmahat/lambdas)%*%Zt),function(x)which(uniqueTaus==x))
      lambdap <- lambdasENuniquehat[indUnique]
    }
    
    if(alpha==0){
      beta <- betaMR
      a0 <- a0MR
    }else{
      penfctr2 <- penfctr
      penfctr2[pen] <- lambdap[pen] 
      not0 <- which(penfctr2!=Inf)
      lambdaEN <- sum(penfctr2[not0])/length(penfctr2[not0])
      penfctr2 <- penfctr2/lambdaEN
      if(model=="cox"){
        glmGR <- glmnet(X[,not0],Y,alpha=alphat,family=fml,standardize = F,
                        penalty.factor=penfctr2[not0], thresh=10^-10)
      }else{
        glmGR <- glmnet(X[,not0],Y,alpha=alphat,family=fml,
                        intercept = intrcpt, standardize = F,
                        penalty.factor=penfctr2[not0], thresh=10^-10)
      }
      
      if(reCV){
        glmGR.cv <- cv.glmnet(X[,not0],Y,alpha=alphat,family=fml,
                              intercept = intrcpt, standardize = F,
                              penalty.factor=penfctr2[not0], thresh=10^-10)
        sopt <- glmGR.cv$lambda.min
      } else sopt <- lambdaEN/n*sd_y
      
      temp <- coef(glmGR,s=sopt,exact=T,x=X[,not0],y=Y,alpha=alphat,
                   penalty.factor=penfctr2[not0],family=fml,intercept=intrcpt) 
      beta <- rep(0,p); beta[not0] <- temp[-1]
      a0 <- temp[1]
      #plot(beta)
      #print(lambdaEN/n*sd_y)
    }
    
  # }else if(alpha==1){ #Compute lasso estimate directly
  #   tauEN <-  sqrt(sigmahat/2/lambdas) #Laplace prior parameter estimates
  #   lambdasEN <- sqrt(sigmahat)/sd_y * sqrt(2*lambdas) #=1/tauL*sigmahat/sd_y #note not same notation as above
  #   lambdap<- lambda/(as.vector(gamma%*%Zt)) #=sigmahat/(tauglobal*as.vector(gamma%*%Zt)) #specific ridge penalty for beta_k
  #   lambdap[lambdap<0]<-Inf
  #   
  #   penfctr2 <- penfctr
  #   penfctr2[pen] <- sqrt(sigmahat)/sd_y * sqrt(2*lambdap[pen]) 
  #   not0 <- which(penfctr2!=Inf)
  #   lambda1 <- sum(penfctr2[not0])/length(penfctr2[not0])
  #   penfctr2 <- penfctr2/lambda1
  #   if(model=="cox"){
  #     glmGR <- glmnet(X[,not0],Y,alpha=1,family=fml,standardize = F,
  #                     penalty.factor=penfctr2[not0], thresh=10^-10)
  #   }else{
  #     glmGR <- glmnet(X[,not0],Y,alpha=1,family=fml,
  #                     intercept = intrcpt, standardize = F,
  #                     penalty.factor=penfctr2[not0], thresh=10^-10)
  #   }
  #   
  #   if(reCV){
  #     glmGR.cv <- cv.glmnet(X[,not0],Y,alpha=1,family=fml,
  #                           intercept = intrcpt, standardize = F,
  #                           penalty.factor=penfctr2[not0], thresh=10^-10)
  #     sopt <- glmGR.cv$lambda.min
  #   } else sopt <- lambda1/n*sd_y
  #   
  #   temp <- coef(glmGR,s=sopt,exact=T,x=X[,not0],y=Y,
  #                penalty.factor=penfctr2[not0],family=fml,intercept=intrcpt) 
  #   beta <- rep(0,p); beta[not0] <- temp[-1]
  #   a0 <- temp[1]
  #   #plot(beta)
  #   #print(lambda1/n*sd_y)
  #   
  #   # browser()
  #   # #or with penalized:
  #   # penfctr2 <- penfctr
  #   # penfctr2[pen] <- sqrt(sigmahat)/sd_y * sqrt(2*lambdap[pen])
  #   # penLasso <- penalized(penalized=X[,not0],response=Y/sd_y,model=model,
  #   #                       lambda1 = penfctr2,
  #   #                       standardize = F)
  #   # betaPEN <- rep(0,p); betaPEN[not0] <- penLasso@penalized*sd_y
  #   # a0PEN <- penLasso@unpenalized
  #   # plot(betaPEN,beta); abline(0,1)
  }else{stop("alpha should be between 0 and 1")}
  
  if(standardise_Y){
    beta <- beta*sd_y_former
    a0 <- a0*sd_y_former
  }  
  
  #Compute predictions for independent data if given----
  if(!is.null(X2)){
    #Ypredridge <- predict(glmR,newx=X2)
    #browser()
    
    if(compareMR){
      if(model=="linear"){
        X2c <- cbind(X2,rep(1,n2))
        YpredMR <- X2c %*% c(betaMR,a0MR)
        MSEMR <- sum((YpredMR-Y2)^2)/n2
      } 
      if(model=='logistic'){
        X2c <- cbind(X2,rep(1,n2))
        YpredMR <- 1/(1+exp(-X2c %*% c(betaMR,a0MR)))
        MSEMR <- sum((YpredMR-Y2)^2)/n2
      }else if(model=="cox"){
        expXb<-exp(X %*% c(betaMR))
        h0 <- sapply(1:length(Y[,1]),function(i){Y[i,2]/sum(expXb[Y[,1]>=Y[i,1]])})#updated baseline hazard in censored times for left out samples
        H0 <- sapply(Y2[,1],function(Ti){sum(h0[Y[,1]<=Ti])})
        YpredMR <- H0*exp(X2 %*% betaMR)
        MSEMR<- sum((YpredMR-Y2[,2])^2)/n2
      }
    }
    
    #Compute test error for elastic net/lasso estimate
    if(model=="linear"){
      X2c <- cbind(X2,rep(1,n2))
      YpredApprox <- X2c %*% c(beta,a0)
      MSEApprox <- sum((YpredApprox-Y2)^2)/n2
    } 
    if(model=='logistic'){
      X2c <- cbind(X2,rep(1,n2))
      YpredApprox <- 1/(1+exp(-X2c %*% c(beta,a0)))
      MSEApprox <- sum((YpredApprox-Y2)^2)/n2
    }else if(model=="cox"){
      expXb<-exp(X %*% c(beta))
      h0 <- sapply(1:length(Y[,1]),function(i){Y[i,2]/sum(expXb[Y[,1]>=Y[i,1]])})#updated baseline hazard in censored times for left out samples
      H0 <- sapply(Y2[,1],function(Ti){sum(h0[Y[,1]<=Ti])})
      YpredApprox <- H0*exp(X2 %*% beta)
      MSEApprox<- sum((YpredApprox-Y2[,2])^2)/n2
    }
    
  }
  
  #Return output----
  output <- list(betaApprox = beta, #lasso
                 a0Approx = a0,
                 #tauApprox = tauEN, #prior parameter estimates
                 lambdaApprox = lambdasEN, #elastic net group penalty
                 tauMR = sigmahat/lambdas, #prior parameter estimates
                 lambdaMR = lambdas,
                 lambdaglobal = lambda,
                 sigmahat = sigmahat,
                 MLinit=MLinit,
                 MLfinal=MLfinal)
                 #betaMR2 = betaGLM, #multiridge penalty with glmnet
                 #a0MR2 = a0GLM)
  if(compareMR){
    output$betaMR <- betaMR #multiridge package
    output$a0MR <- a0MR
  }
  if(!is.null(X2)){
    output$YpredApprox<-YpredApprox #predictions for test set
    output$MSEApprox <- MSEApprox #MSE on test set
    if(compareMR){
      output$YpredMR <-YpredMR #predictions for test set
      output$MSEMR <- MSEMR #MSE on test set
    }
  }
  if(selectAIC){
    output$AICmodels <- list(
      "multigroup"=list(
        "lambdas"=lambdas,
        "sigmahat"=sigmahat,
        "AIC"=AICmultigroup
      ),
      "onegroup"=list(
        "lambdas"=lambda1group,
        "sigmahat"=sigmahat1group,
        "AIC"=AIC1group
      )
    )
    if(modelbestAIC!="multigroup"){
      output$AICmodels$multigroup$lambdas <- lambdasNotOptimalAIC
      output$AICmodels$multigroup$sigmahat <- sigmahatNotOptimalAIC
    }
    if(resultsAICboth){
        output$AICmodels$onegroup$fit <- squeezy(Y,X,grouping,alpha=alpha,model=model,
                X2=X2,Y2=Y2,unpen=unpen,intrcpt=intrcpt,
                method="MML",fold=fold,
                compareMR = compareMR,selectAIC=F,
                fit.ecpc=NULL,
                lambdas=rep(lambda1group,G),lambdaglobal=lambda1group,
                sigmasq=sigmahat1group,
                standardise_Y=standardise_Y,reCV=standardise_Y,resultsAICboth=F)
        output$AICmodels$multigroup$fit <- squeezy(Y,X,grouping,alpha=alpha,model=model,
                       X2=X2,Y2=Y2,unpen=unpen,intrcpt=intrcpt,
                       method="MML",fold=fold,
                       compareMR = compareMR,selectAIC=F,
                       fit.ecpc=NULL,
                       lambdas=output$AICmodels$multigroup$lambdas,
                       lambdaglobal=lambda1group,
                       sigmasq=output$AICmodels$multigroup$sigmahat,
                       standardise_Y=standardise_Y,reCV=standardise_Y,resultsAICboth=F)
    }
    
    output$modelbestAIC <- modelbestAIC
  }
  
  return(output)
}

#auxiliary function
#log-determinant; needed below
.logdet <- function(mat) return(determinant(mat,logarithm=T)$modulus[1])

.SigmaBlocks <- function(XXblocks,lambdas){ #computes X Laminv X^T from block cross-products
  nblocks <- length(XXblocks)
  if(nblocks != length(lambdas)){
    print("Error: Number of penalty parameters should equal number of blocks")
    return(NULL)
  } else {
    Sigma<-Reduce('+', lapply(1:nblocks,function(i) XXblocks[[i]] * 1/lambdas[i]))
    return(Sigma)
  }
}

minML.LA.ridgeGLM<- function(loglambdas,XXblocks,Y,sigmasq=1,
                             Xunpen=NULL,intrcpt=F,model,minlam=0,opt.sigma=F){
  #compute Laplace approximation (LA) of the minus log marginal likelihood of ridge penalised GLMs
  #(NOTE: for now only implemented for linear and logistic)
  #Input:
  #loglambdas: logarithm of group penalties
  #XXblocks: nxn matrices for each group
  #Y: response vector
  #sigmasq: noise level in linear regression (1 for logistic)
  #Xunpen: unpenalised variables
  #intercept: T/F to include intercept
  #model: "linear" or "logistic" 
  #opt.sigma: (linear case) T/F if log(sigma^2) is given as first argument of 
  #           loglambdas for optimisation purposes
  if(model!="linear"){
    opt.sigma <- F
    sigmasq <- 1
  } 
  
  n<-length(Y)
  lambdas <- exp(loglambdas)+minlam
  if(opt.sigma){
    sigmasq <- lambdas[1]
    lambdas <- lambdas[-1]
  }
  Lam0inv <- .SigmaBlocks(XXblocks,lambdas) 
  if(!is.null(Xunpen)){
    fit <- IWLSridge(XXT=Lam0inv,Y=Y, model=model,intercept=intrcpt,
                     X1=Xunpen,maxItr = 500,eps=10^-12) #Fit. fit$etas contains the n linear predictors
  }else{
    fit <- IWLSridge(XXT=Lam0inv,Y=Y, model=model,intercept=intrcpt,
                     maxItr = 500,eps=10^-12) #Fit. fit$etas contains the n linear predictors
  }
  eta <- fit$etas + fit$eta0
  
  if(model=="linear"){
    mu <- eta
    W <- diag(rep(1,n))
    Hpen <- Lam0inv - Lam0inv %*% solve(solve(W)+Lam0inv, Lam0inv)
    
    t1 <- dmvnorm(c(Y),mean=eta,sigma=sigmasq*diag(rep(1,n)),log=T) #log likelihood in linear predictor
  }else{
    mu <- 1/(1+exp(-eta))
    W <- diag(c(mu)*c(1-mu))
    
    #Hpen <- Lam0inv - Lam0inv %*% solve(diag(1/diag(W))+Lam0inv, Lam0inv)
    Hpen <- Lam0inv - Lam0inv %*% solve(diag(1,n)+W%*%Lam0inv, W%*%Lam0inv)
    
    t1 <- sum(Y*log(mu)+(1-Y)*log(1-mu)) #log likelihood in linear predictor
  }
  
  t2 <- -1/2/sigmasq*sum(c(Y-mu)*eta)
  t3 <- 1/2 * .logdet(diag(rep(1,n))-W%*%Hpen)
  #return(-t2)
  return(-(t1+t2+t3))
}

dminML.LA.ridgeGLM<- function(loglambdas,XXblocks,Y,sigmasq=1,
                             Xunpen=NULL,intrcpt=F,model,minlam=0,opt.sigma=F){
  #compute Laplace approximation (LA) of the first derivative of the minus log 
  #  marginal likelihood of ridge penalised GLMs
  #(NOTE: for now only implemented for linear and logistic)
  #Input:
  #loglambdas: logarithm of group penalties 
  #XXblocks: nxn matrices for each group
  #Y: response vector
  #sigmasq: noise level in linear regression (1 for logistic)
  #Xunpen: unpenalised variables
  #intercept: T/F to include intercept
  #model: "linear" or "logistic" 
  #minlam: minimum of lambda that is added to the value given (useful in optimisation settings)
  #opt.sigma: (linear case) T/F if derivative to log(sigma^2) should be given for optimisation purposes
  #           Note that sigma^2 should then be given as the first argument of loglambdas     
  
  if(model!="linear") opt.sigma <- F
  
  n<-length(Y)
  lambdas <- exp(loglambdas)+minlam
  if(opt.sigma){
    sigmasq <- lambdas[1]
    lambdas <- lambdas[-1]
  }
  rhos<-log(lambdas)
  Lam0inv <- .SigmaBlocks(XXblocks,lambdas) 
  if(!is.null(Xunpen)){
    fit <- IWLSridge(XXT=Lam0inv,Y=Y, model=model,intercept=intrcpt,
                     X1=Xunpen,maxItr = 500,eps=10^-12) #Fit. fit$etas contains the n linear predictors
  }else{
    fit <- IWLSridge(XXT=Lam0inv,Y=Y, model=model,intercept=intrcpt,
                     maxItr = 500,eps=10^-12) #Fit. fit$etas contains the n linear predictors
  }
  eta <- fit$etas + fit$eta0
  
  if(intrcpt) Xunpen <- cbind(Xunpen,rep(1,n))
  #eta0 <- fit$eta0 #TD: add unpenalised variables/intercept
  
  #GLM-type specific terms
  if(model=="linear"){
    mu <- eta
    W <- diag(rep(1,n))
    dwdeta <- rep(0,n)
  }else if(model=="logistic"){
    mu <- 1/(1+exp(-eta))
    W <- diag(c(mu)*c(1-mu))
    dwdeta <- c(mu*(1-mu)*(1-2*mu))
  }else{
    stop(print(paste("Only model type linear and logistic supported")))
  }
  
  Hpen <- Lam0inv - Lam0inv %*% solve(solve(W)+Lam0inv, Lam0inv)
  if(!is.null(Xunpen)){
    WP1W <- diag(rep(1,n)) - Xunpen%*%solve(t(Xunpen)%*%W%*%Xunpen,t(Xunpen)%*%W)
  }else{
    WP1W <- diag(rep(1,n))
  }
  #L2 <- WP1W%*%(diag(rep(1,n))-t(solve(solve(W)+WP1W%*%Lam0inv,WP1W%*%Lam0inv)))
  L2 <- WP1W%*%(diag(rep(1,n))-t(solve(diag(1,n)+W%*%WP1W%*%Lam0inv,W%*%WP1W%*%Lam0inv)))
  Hj <- lapply(1:length(XXblocks),function(i){
    L2%*%XXblocks[[i]]/lambdas[i]
  })
  
  if(is.null(Xunpen)) {Hres <- .Hpen(diag(W),Lam0inv);Hmat <- Hres$Hmat} else {
    Hres <- .Hunpen(diag(W),Lam0inv,Xunpen); Hmat <- Hres$Hmat
  }
  # detadrho <- lapply(1:length(rhos),function(j){
  #   -(diag(rep(1,n))-L2%*%Lam0inv%*%W)%*%t(Hj[[j]])%*%c(Y-mu+W%*%eta)})
  detadrho <- lapply(1:length(rhos),function(j){
    -exp(rhos[j])*(diag(rep(1,n))-Hmat%*%W)%*%t(Hj[[j]]/lambdas[j])%*%c(Y-mu+W%*%eta)})
  dWdrho <- lapply(1:length(rhos),function(j){dwdeta*c(detadrho[[j]])})
  
  t1 <- sapply(1:length(rhos),function(j){-1/sigmasq*c((Y-mu))%*%detadrho[[j]]})
  t2 <- sapply(1:length(rhos),function(j){1/2/sigmasq*c(-W%*%eta+Y-mu)%*%detadrho[[j]]})
  #t3 <- sapply(1:length(rhos),function(j){-0.5*sum(diag(solve(solve(W)+Lam0inv,XXblocks[[j]]/sigmasq/lambdas[j])))})
  t3 <- sapply(1:length(rhos),function(j){-0.5*sum(diag(solve(diag(1,n)+W%*%Lam0inv,W%*%XXblocks[[j]]/lambdas[j])))})
  t4 <- sapply(1:length(rhos),function(j){0.5*sum(diag(Hpen)*dWdrho[[j]])})
  
  #return(t2)
  if(!opt.sigma) return((t1+t2+t3+t4)) #derivative to rho
  
  #for linear model, may want to return derivative to log(sigma^2) as well
  # browser()
  # invMat <- solve(diag(rep(sigmasq,n))+Lam0inv*sigmasq)
  # ts.1 <- -0.5*sum(diag(invMat))
  # ts.2 <- -0.5*t(Y-mu)%*%invMat%*%invMat%*%(Y-mu)
  
  ts.1 <- n/2/sigmasq
  ts.2 <- -1/2/sigmasq^2*t(Y-eta)%*%Y
  
  return(c(c(ts.1+ts.2)*sigmasq,(t1+t2+t3+t4))) #derivative to rho
}

AIC.LA.ridgeGLM<- function(loglambdas,XXblocks,Y,sigmasq=1,
                             Xunpen=NULL,intrcpt=F,model,minlam=0){
  
  minLL <- minML.LA.ridgeGLM(loglambdas,XXblocks,Y,sigmasq=sigmasq,
                             Xunpen=Xunpen,intrcpt=intrcpt,model=model,minlam=minlam)
  k <- length(XXblocks) #number of parameters
  if(model=="linear") k <- k+1 #sigma estimated in linear model as well
  AIC <- 2*minLL + 2*k
  return(AIC)
}

.cv.squeezy <- function(Y,X,folds,type.measure="MSE",...){
  #Input:
  #X: nxp observed data matrix
  #Y: n-dimensional vector for response
  #folds: number of folds to evaluate response, or list with samples in each fold
  #type.measure: type of performance measure for evaluation (MSE or AUC)
  more.args <- list(...) #more arguments needed for squeezy
  
  n <- dim(X)[1] #number of samples
  p <- dim(X)[2] #number of covariates
  
  if(!is.element("model",names(more.args))){
    if(all(is.element(Y,c(0,1))) || is.factor(Y)){
      model <- "logistic" 
    } else if(all(is.numeric(Y)) & !(is.matrix(Y) && dim(Y)[2]==2)){
      model <- "linear"
    }else{
      model <- "cox"
    }
  }else{
    model <- more.args$model
  }
  
  if(is.numeric(folds)){ #number of folds given
    folds2<-.produceFolds(n,folds,Y,balance=balance,model=model) #produce folds balanced in response
  }else{
    folds2 <- folds
  }
  nfolds <- length(folds2)
  
  
  grpsno <- 1:length(more.args$grouping) #vector with group numbers of the grouping
  
  Res <- list() #list for raw data
  df <- data.frame() #data frame for prediction results
  dfGrps <- data.frame() #data frame for hyperparameter estimation results
  for(i in 1:nfolds){
    tic <- proc.time()[[3]]
    Res[[i]] <- do.call(squeezy,args=c(list(X=X[-folds2[[i]],],Y=Y[-folds2[[i]]],
                                            X2=X[folds2[[i]],],Y2=Y[folds2[[i]]]),
                                            more.args)
    )
    toc <- proc.time()[[3]]
    Res[[i]]$time <- toc - tic
    
    #update dfPred
    df2<-data.frame("Ypred"=c(Res[[i]]$YpredMR,Res[[i]]$YpredApprox))
    df2$Method <- rep(c("multiridge","squeezy_lasso"),each=length(folds2[[i]]))
    df2$NumberSelectedVars <- rep(c(sum(Res[[i]]$betaMR1!=0),
                                    sum(Res[[i]]$betaApprox!=0)),each=length(folds2[[i]]))
    df2$Fold <- i
    df2$Sample <- rep(folds2[[i]],2)
    df2$Time <-  Res[[i]]$time
    df2$Truth <- rep(Y[folds2[[i]]],2)
    df <- rbind(df,df2)
    
    #update dfGrps
    df3<-data.frame("Group"=rep(grpsno,2),
                    "Group parameter"=c(Res[[i]]$tauMR,Res[[i]]$tauApprox),
                    "Penalty parameter"=c(Res[[i]]$lambdaMR,Res[[i]]$lambdaApprox) #note: penalty parameter on glmnet-scale
                    )
    df3$Method <- rep(c("multiridge","squeezy_lasso"),each=length(grpsno))
    df3$Fold <- i
    dfGrps<-rbind(dfGrps,df3)
  }
  
  #data frame with performance measure
  if(is.factor(df$Truth)){
    warning("Response Y given as factor, transformed to numeric to compute AUC")
    if(!silent) print(levels(df$Truth)[1],"transformed to",0)
    if(!silent) print(levels(df$Truth)[2],"transformed to",1)
    df$Truth <- as.numeric(df$Truth)-1
  }
  if(type.measure=="MSE"){
    dfCVM <- df %>% group_by(Method,Fold) %>% summarise(CVM = mean((Ypred-Truth)^2),Type="MSE",
                                                        NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  }
  else if(type.measure=="AUC"){
    dfROC<-data.frame()
    for(i in levels(df$Method)){
      temp<-data.frame()
      cutoffs<-rev(seq(0,1,by=0.001))
      rocGR <- roc(probs=df$Ypred[df$Method==i],true=df$Truth[df$Method==i],cutoffs=cutoffs)
      temp<-data.frame("FPR"=rocGR[1,],"TPR"=rocGR[2,],"Accuracy"=rocGR[3,])
      temp$Method <- i
      temp$AUC<-c(auc(rocGR))
      temp$NumberSelectedVars<-mean(df$NumberSelectedVars[df$Method==i])
      dfROC<-rbind(dfROC,temp)
    }
    dfCVM <- dfROC %>% group_by(Method) %>% summarise(CVM=mean(AUC),Type="AUC",
                                                      NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  }else{
    warning(paste("The type of measure",type.measure,"is not yet supported."))
  }
  
  return(list("Res"=Res,"dfPred"=df,"dfGrps"=dfGrps,"dfCVM"=dfCVM))
}

.traintest.squeezy <- function(Y,X,Y2,X2,type.measure="MSE",
                                multi_grouping=F,grouping,
                                args.ecpc=NULL,ecpcinit=T,
                              ncores=1,
                                ...){
  #Input:
  #X: nxp observed data matrix
  #Y: n-dimensional vector for response
  #X2: independent test observed data
  #Y2: independent test response data
  #multi_grouping: is a separate grouping given for train/test data sets?
  #grouping: as in squeezy, or list of groupings for each train/test data set
  #type.measure: type of performance measure for evaluation (MSE or AUC)
  #args.ecpc: list with arguments used for ecpc function
  #ecpcinit: T/F should ecpc be fit for initialisation?
  
  args.squeezy <- list(...) #more arguments needed for squeezy
  
  if(is.list(X)) nSim <- length(X)
  else{
    nSim <- 1
    X <- list(X)
    Y <- list(Y)
    X2 <- list(X2)
    Y2 <- list(Y2)
  } 
  
  if(!is.element("model",c(names(args.ecpc),names(args.squeezy)))){
    if(all(is.element(Y[[1]],c(0,1))) || is.factor(Y[[1]])){
      model <- "logistic" 
    } else if(all(is.numeric(Y[[1]])) & !(is.matrix(Y[[1]]) && dim(Y[[1]])[2]==2)){
      model <- "linear"
    }else{
      model <- "cox"
    }
  }else if(is.element("model",names(args.squeezy))){
    model <- args.squeezy$model
  }else{
    model <- args.ecpc$model
  }
  
  Res <- list() #list for raw data
  df <- data.frame() #data frame for prediction results
  dfGrps <- data.frame() #data frame for hyperparameter estimation results
  if(ncores==1){
    for(i in 1:nSim){
      n <- dim(X[[i]])[1] #number of samples
      p <- dim(X[[i]])[2] #number of covariates
      
      if(!multi_grouping) groupingTemp <- grouping
      else groupingTemp <- grouping[[i]]
      grpsno <- 1:length(groupingTemp) #vector with group numbers of the grouping
      
      if(ecpcinit){
        tic <- proc.time()[[3]]
        res <- do.call(ecpc,args=c(list(X=X[[i]],Y=Y[[i]],
                                        X2=X2[[i]],Y2=Y2[[i]],
                                        groupings=list(groupingTemp),
                                        hypershrinkage="none",postselection = F,model=model),
                                   args.ecpc)
        )
        Res[[i]] <- do.call(squeezy,args=c(list(X=X[[i]],Y=Y[[i]],
                                                X2=X2[[i]],Y2=Y2[[i]],
                                                grouping=groupingTemp,fit.ecpc=res),
                                           args.squeezy)
        )
        toc <- proc.time()[[3]]
        Res[[i]]$time <- toc - tic
      }else{
        tic <- proc.time()[[3]]
        Res[[i]] <- do.call(squeezy,args=c(list(X=X[[i]],Y=Y[[i]],
                                                X2=X2[[i]],Y2=Y2[[i]],
                                                grouping=groupingTemp),
                                           args.squeezy)
        )
        toc <- proc.time()[[3]]
        Res[[i]]$time <- toc - tic
      }
      
      #update dfPred
      df2<-data.frame("Ypred"=c(Res[[i]]$YpredMR,Res[[i]]$YpredApprox))
      df2$Method <- rep(c("multiridge","squeezy_EN"),each=length(Y2[[i]]))
      df2$NumberSelectedVars <- rep(c(sum(Res[[i]]$betaMR!=0),
                                      sum(Res[[i]]$betaApprox!=0)),each=length(Y2[[i]]))
      df2$Dataset <- i
      df2$Sample <- rep(1:length(Y2[[i]]),2)
      df2$Time <-  Res[[i]]$time
      df2$Truth <- rep(Y2[[i]],2)
      df <- rbind(df,df2)
      
      #update dfGrps
      df3<-data.frame("Group"=rep(grpsno,2),
                      "Group parameter"=c(Res[[i]]$tauMR,Res[[i]]$tauApprox),
                      "Penalty parameter"=c(Res[[i]]$lambdaMR,Res[[i]]$lambdaApprox) #note: penalty parameter on glmnet-scale
      )
      df3$Method <- rep(c("multiridge","squeezy_EN"),each=length(grpsno))
      df3$Dataset <- i
      dfGrps<-rbind(dfGrps,df3)
    }
  }else{
    cl <- makeCluster(ncores) #set up parallel cluster
    registerDoParallel(cl)
    finalMatrix <- foreach(i=1:nfolds, .combine=rbind, 
                           .packages = c("glmnet","penalized","mvtnorm","gglasso",
                                         "Matrix","Rsolnp","ecpc")) %dopar% {
           n <- dim(X[[i]])[1] #number of samples
           p <- dim(X[[i]])[2] #number of covariates
           
           if(!multi_grouping) groupingTemp <- grouping
           else groupingTemp <- grouping[[i]]
           grpsno <- 1:length(groupingTemp) #vector with group numbers of the grouping
           
           if(ecpcinit){
             tic <- proc.time()[[3]]
             res <- do.call(ecpc,args=c(list(X=X[[i]],Y=Y[[i]],
                                             X2=X2[[i]],Y2=Y2[[i]],
                                             groupings=list(groupingTemp),
                                             hypershrinkage="none",postselection = F,model=model),
                                        args.ecpc)
             )
             Res[[i]] <- do.call(squeezy,args=c(list(X=X[[i]],Y=Y[[i]],
                                                     X2=X2[[i]],Y2=Y2[[i]],
                                                     grouping=groupingTemp,fit.ecpc=res),
                                                args.squeezy)
             )
             toc <- proc.time()[[3]]
             Res[[i]]$time <- toc - tic
           }else{
             tic <- proc.time()[[3]]
             Res[[i]] <- do.call(squeezy,args=c(list(X=X[[i]],Y=Y[[i]],
                                                     X2=X2[[i]],Y2=Y2[[i]],
                                                     grouping=groupingTemp),
                                                args.squeezy)
             )
             toc <- proc.time()[[3]]
             Res[[i]]$time <- toc - tic
           }
           
           #update dfPred
           df2<-data.frame("Ypred"=c(Res[[i]]$YpredMR,Res[[i]]$YpredApprox))
           df2$Method <- rep(c("multiridge","squeezy_EN"),each=length(Y2[[i]]))
           df2$NumberSelectedVars <- rep(c(sum(Res[[i]]$betaMR!=0),
                                           sum(Res[[i]]$betaApprox!=0)),each=length(Y2[[i]]))
           df2$Dataset <- i
           df2$Sample <- rep(1:length(Y2[[i]]),2)
           df2$Time <-  Res[[i]]$time
           df2$Truth <- rep(Y2[[i]],2)
           #df <- rbind(df,df2)
           
           #update dfGrps
           df3<-data.frame("Group"=rep(grpsno,2),
                           "Group parameter"=c(Res[[i]]$tauMR,Res[[i]]$tauApprox),
                           "Penalty parameter"=c(Res[[i]]$lambdaMR,Res[[i]]$lambdaApprox) #note: penalty parameter on glmnet-scale
           )
           df3$Method <- rep(c("multiridge","squeezy_EN"),each=length(grpsno))
           df3$Dataset <- i
           #dfGrps<-rbind(dfGrps,df3)
           
           list("Res"=Res,"df"=df2,"dfGrps"=df3)
    }
    Res <- lapply(1:nfolds,function(i) finalMatrix[i,1][[1]][[i]])
    df2 <- lapply(1:nfolds,function(i) finalMatrix[i,2][[1]])
    dfGrps2 <- lapply(1:nfolds,function(i) finalMatrix[i,3][[1]])
    df <- df2[[1]]; for(i in 2:nfolds) df <- rbind(df,df2[[i]])
    dfGrps <- dfGrps2[[1]]; for(i in 2:nfolds) dfGrps <- rbind(dfGrps,dfGrps2[[i]])
    stopCluster(cl); rm(cl)
  }
  
  #data frame with performance measure
  if(is.factor(df$Truth)){
    warning("Response Y given as factor, transformed to numeric to compute AUC")
    if(!silent) print(levels(df$Truth)[1],"transformed to",0)
    if(!silent) print(levels(df$Truth)[2],"transformed to",1)
    df$Truth <- as.numeric(df$Truth)-1
  }
  if(type.measure=="MSE"){
    dfCVM <- df %>% group_by(Method,Dataset) %>% summarise(CVM = mean((Ypred-Truth)^2),Type="MSE",
                                                        NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  }
  else if(type.measure=="AUC"){
    dfROC<-data.frame()
    for(i in levels(df$Method)){
      temp<-data.frame()
      cutoffs<-rev(seq(0,1,by=0.001))
      rocGR <- roc(probs=df$Ypred[df$Method==i],true=df$Truth[df$Method==i],cutoffs=cutoffs)
      temp<-data.frame("FPR"=rocGR[1,],"TPR"=rocGR[2,],"Accuracy"=rocGR[3,])
      temp$Method <- i
      temp$AUC<-c(auc(rocGR))
      temp$NumberSelectedVars<-mean(df$NumberSelectedVars[df$Method==i])
      dfROC<-rbind(dfROC,temp)
    }
    dfCVM <- dfROC %>% group_by(Method) %>% summarise(CVM=mean(AUC),Type="AUC",
                                                      NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  }else{
    warning(paste("The type of measure",type.measure,"is not yet supported."))
  }
  
  return(list("Res"=Res,"dfPred"=df,"dfGrps"=dfGrps,"dfCVM"=dfCVM))
}

.outercv.glmnet <- function(Y,X,folds,type.measure="MSE",...){
  #Input:
  #X: nxp observed data matrix
  #Y: n-dimensional vector for response
  #folds: number of folds to evaluate response, or list with samples in each fold
  #type.measure: type of performance measure for evaluation (MSE or AUC)
  more.args <- list(...) #more arguments needed for squeezy
  
  n <- dim(X)[1] #number of samples
  p <- dim(X)[2] #number of covariates
  
  if(!is.element("model",names(more.args))){
    if(all(is.element(Y,c(0,1))) || is.factor(Y)){
      model <- "logistic" 
    } else if(all(is.numeric(Y)) & !(is.matrix(Y) && dim(Y)[2]==2)){
      model <- "linear"
    }else{
      model <- "cox"
    }
  }else{
    model <- more.args$model
  }
  
  if(is.numeric(folds)){ #number of folds given
    folds2<-.produceFolds(n,folds,Y,balance=balance,model=model) #produce folds balanced in response
  }else{
    folds2 <- folds
  }
  nfolds <- length(folds2)
  
  Res <- list() #list for raw data
  df <- data.frame() #data frame for prediction results
  dfGrps <- data.frame() #data frame for hyperparameter estimation results
  for(i in 1:nfolds){
    #first compute optimal lambda
    tic <- proc.time()[[3]]
    Res[[i]] <- do.call(cv.glmnet,args=c(list(x=X[-folds2[[i]],],y=Y[-folds2[[i]]]),
                                         more.args)
    )
    toc <- proc.time()[[3]]
    Res[[i]]$time <- toc - tic
    
    temp <- coef.glmnet(Res[[i]],s=Res[[i]]$lambda.min,exact=T,
                        x=X[-folds2[[i]],],y=Y[-folds2[[i]]])
    Res[[i]]$betaL <- temp[-1]
    #a0 <- temp[1]
    fit <- do.call(glmnet,args=c(list(x=X[-folds2[[i]],],y=Y[-folds2[[i]]]),
                                    more.args))
    Ypred <- predict.glmnet(fit,newx=X[folds2[[i]],],
                            s=Res[[i]]$lambda.min,exact=T,
                            x=X[-folds2[[i]],],y=Y[-folds2[[i]]])
    #update dfPred
    df2<-data.frame("Ypred"=c(Ypred))
    df2$Method <- "lasso"
    df2$NumberSelectedVars <- sum(Res[[i]]$betaL!=0)
    df2$Fold <- i
    df2$Sample <- folds2[[i]]
    df2$Time <-  Res[[i]]$time
    df2$Truth <- Y[folds2[[i]]]
    df <- rbind(df,df2)
    
    #update dfGrps
    df3<-data.frame("Penalty parameter"=Res[[i]]$lambda.min,
                    "Group"=1) #note: penalty parameter on glmnet-scale
    df3$Method <- "lasso"
    df3$Fold <- i
    dfGrps<-rbind(dfGrps,df3)
  }
  
  #data frame with performance measure
  if(is.factor(df$Truth)){
    warning("Response Y given as factor, transformed to numeric to compute AUC")
    if(!silent) print(levels(df$Truth)[1],"transformed to",0)
    if(!silent) print(levels(df$Truth)[2],"transformed to",1)
    df$Truth <- as.numeric(df$Truth)-1
  }
  if(type.measure=="MSE"){
    dfCVM <- df %>% group_by(Method,Fold) %>% summarise(CVM = mean((Ypred-Truth)^2),Type="MSE",
                                                        NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  }
  else if(type.measure=="AUC"){
    dfROC<-data.frame()
    for(i in levels(df$Method)){
      temp<-data.frame()
      cutoffs<-rev(seq(0,1,by=0.001))
      rocGR <- roc(probs=df$Ypred[df$Method==i],true=df$Truth[df$Method==i],cutoffs=cutoffs)
      temp<-data.frame("FPR"=rocGR[1,],"TPR"=rocGR[2,],"Accuracy"=rocGR[3,])
      temp$Method <- i
      temp$AUC<-c(auc(rocGR))
      temp$NumberSelectedVars<-mean(df$NumberSelectedVars[df$Method==i])
      dfROC<-rbind(dfROC,temp)
    }
    dfCVM <- dfROC %>% group_by(Method) %>% summarise(CVM=mean(AUC),Type="AUC",
                                                      NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  }else{
    warning(paste("The type of measure",type.measure,"is not yet supported."))
  }
  
  return(list("Res"=Res,"dfPred"=df,"dfGrps"=dfGrps,"dfCVM"=dfCVM))
}

.traintest.glmnet <- function(Y,X,Y2,X2,type.measure="MSE",
                               ...){
  #Input:
  #X: nxp observed data matrix
  #Y: n-dimensional vector for response
  #X2: independent test observed data
  #Y2: independent test response data
  #type.measure: type of performance measure for evaluation (MSE or AUC)
  more.args <- list(...) #more arguments needed for glmnet
  
  if(is.list(X)) nSim <- length(X)
  else{
    nSim <- 1
    X <- list(X)
    Y <- list(Y)
    X2 <- list(X2)
    Y2 <- list(Y2)
  } 
  
  if(!is.element("model",names(more.args))){
    if(all(is.element(Y[[1]],c(0,1))) || is.factor(Y[[1]])){
      model <- "logistic" 
    } else if(all(is.numeric(Y[[1]])) & !(is.matrix(Y[[1]]) && dim(Y[[1]])[2]==2)){
      model <- "linear"
    }else{
      model <- "cox"
    }
  }else{
    model <- more.args$model
  }
  
  Res <- list() #list for raw data
  df <- data.frame() #data frame for prediction results
  dfGrps <- data.frame() #data frame for hyperparameter estimation results
  for(i in 1:nSim){
    n <- dim(X[[i]])[1] #number of samples
    p <- dim(X[[i]])[2] #number of covariates
    
    tic <- proc.time()[[3]]
    Res[[i]] <- do.call(cv.glmnet,args=c(list(x=X[[i]],y=Y[[i]]),
                                         more.args)
    )
    toc <- proc.time()[[3]]
    Res[[i]]$time <- toc - tic
    
    temp <- coef.glmnet(Res[[i]],s=Res[[i]]$lambda.min,exact=T,
                        x=X[[i]],y=Y[[i]])
    Res[[i]]$betaL <- temp[-1]
    #a0 <- temp[1]
    fit <- do.call(glmnet,args=c(list(x=X[[i]],y=Y[[i]]),
                                 more.args))
    Ypred <- predict.glmnet(fit,newx=X2[[i]],
                            s=Res[[i]]$lambda.min,exact=T,
                            x=X[[i]],y=Y[[i]])
    
    #update dfPred
    df2<-data.frame("Ypred"=c(Ypred))
    df2$Method <- "lasso"
    df2$NumberSelectedVars <- sum(Res[[i]]$betaL!=0)
    df2$Dataset <- i
    df2$Sample <- 1:length(Y2[[i]])
    df2$Time <-  Res[[i]]$time
    df2$Truth <- Y2[[i]]
    df <- rbind(df,df2)
    
    #update dfGrps
    df3<-data.frame("Penalty parameter"=Res[[i]]$lambda.min,
                    "Group"=1) #note: penalty parameter on glmnet-scale
    df3$Method <- "lasso"
    df3$Dataset <- i
    dfGrps<-rbind(dfGrps,df3)
  }
  
  #data frame with performance measure
  if(is.factor(df$Truth)){
    warning("Response Y given as factor, transformed to numeric to compute AUC")
    if(!silent) print(levels(df$Truth)[1],"transformed to",0)
    if(!silent) print(levels(df$Truth)[2],"transformed to",1)
    df$Truth <- as.numeric(df$Truth)-1
  }
  if(type.measure=="MSE"){
    dfCVM <- df %>% group_by(Method,Dataset) %>% summarise(CVM = mean((Ypred-Truth)^2),Type="MSE",
                                                           NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  }
  else if(type.measure=="AUC"){
    dfROC<-data.frame()
    for(i in levels(df$Method)){
      temp<-data.frame()
      cutoffs<-rev(seq(0,1,by=0.001))
      rocGR <- roc(probs=df$Ypred[df$Method==i],true=df$Truth[df$Method==i],cutoffs=cutoffs)
      temp<-data.frame("FPR"=rocGR[1,],"TPR"=rocGR[2,],"Accuracy"=rocGR[3,])
      temp$Method <- i
      temp$AUC<-c(auc(rocGR))
      temp$NumberSelectedVars<-mean(df$NumberSelectedVars[df$Method==i])
      dfROC<-rbind(dfROC,temp)
    }
    dfCVM <- dfROC %>% group_by(Method) %>% summarise(CVM=mean(AUC),Type="AUC",
                                                      NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  }else{
    warning(paste("The type of measure",type.measure,"is not yet supported."))
  }
  
  return(list("Res"=Res,"dfPred"=df,"dfGrps"=dfGrps,"dfCVM"=dfCVM))
}

.Hunpen <- function(WV,XXT,X1){
  #Function from multiridge-package written by Mark van de Wiel
  #WV:weigths as vector (n)
  #XXT: sample cross-product (nxn) penalized variables
  #X1 (p1xn) design matrix unpenalized variables
  n <- length(WV)
  WsqrtV <- sqrt(WV)
  WinvsqrtV <- 1/WsqrtV
  X1W <- WsqrtV * X1
  X1aux <- solve(t(X1W) %*% X1W) %*% t(X1W)
  X1Wproj <- X1W %*% X1aux
  P1W <- diag(n) - X1Wproj
  GammaW <- t(t(WsqrtV * XXT) * WsqrtV) #faster
  P1GammaW <- P1W %*% GammaW
  # cons <- mean(abs(P1GammaW))
  # solve(diag(n)/cons+P1GammaW/cons)
  invforH2<-try(solve(diag(n) + P1GammaW))
  if(class(invforH2)[1]=="try-error"){
    svdXXT <- svd(diag(n) + P1GammaW)
    svdd <- svdXXT$d
    #reci <- 1/svdd[1:n]
    reci <- c(1/svdd[1:n-1],0)
    invforH2 <- svdXXT$v %*% (reci * t(svdXXT$u))
  }
  #invforH2 <- solve(diag(n) + P1GammaW)
  #H2 <- Winvsqrt %*% GammaW %*% (diag(n) -invforH2 %*% P1GammaW) %*% P1W %*% Winvsqrt
  MW <- WsqrtV * t(t((diag(n) - invforH2 %*% P1GammaW) %*% P1W) * WinvsqrtV)
  H20 <- GammaW %*% (WinvsqrtV * MW)
  H2 <- t(t(WinvsqrtV * H20)) #faster
  #Hboth <- Winvsqrt %*% X1Wproj %*% (diag(n) - Wsqrt %*% H2 %*% Wsqrt) %*% Winvsqrt + H2
  KW <- t(t(X1aux %*% (diag(n) - t(t(WsqrtV * H2) * WsqrtV))) * WinvsqrtV)
  Hmat <-(WinvsqrtV * X1W) %*% KW  + H2
  return(list(Hmat=Hmat,MW=MW,KW=KW))
}

#Produce balanced folds----
.produceFolds <- function(nsam,outerfold,response,model="logistic",balance=TRUE,fixedfolds=F){
  if(fixedfolds) set.seed(3648310) #else set.seed(NULL)
  if(model=="linear") balance=F
  if(!balance){
    rand<-sample(1:nsam)
    grs1 <- floor(nsam/outerfold)
    grs2 <- grs1+1
    ngr1 <- outerfold*grs2 - nsam
    folds <- lapply(1:outerfold,function(xg) {
      if(xg <= ngr1) els <- rand[(1+(xg-1)*grs1):(xg*grs1)] else els <- rand[(ngr1*grs1 + 1+(xg-ngr1-1)*grs2):(ngr1*grs1 + (xg-ngr1)*grs2)]
      return(els)
    }
    )} else {  #balanced folds
      if(model=="logistic") if(class(response)=="factor") nev <- which((as.numeric(response)-1)==1) else nev <- which(response==1)  
      if(model=="survival") nev <- which(response[,1]==1)    
      nsamev <- length(nev) 
      randev<-sample(nev)
      grs1 <- floor(nsamev/outerfold)
      grs2 <- grs1+1
      ngr1 <- outerfold*grs2 - nsamev
      foldsev <- lapply(1:outerfold,function(xg) {
        if(xg <= ngr1) els <- randev[(1+(xg-1)*grs1):(xg*grs1)] else els <- randev[(ngr1*grs1 + 1+(xg-ngr1-1)*grs2):(ngr1*grs1 + (xg-ngr1)*grs2)]
        return(els)
      }
      )
      nonev <- setdiff(1:nsam,nev)
      nsamnonev <- length(nonev) 
      randnonev<-sample(nonev)
      grs1 <- floor(nsamnonev/outerfold)
      grs2 <- grs1+1
      ngr1 <- outerfold*grs2 - nsamnonev
      foldsnonev <- lapply(1:outerfold,function(xg) {
        if(xg <= ngr1) els <- randnonev[(1+(xg-1)*grs1):(xg*grs1)] else els <- randnonev[(ngr1*grs1 + 1+(xg-ngr1-1)*grs2):(ngr1*grs1 + (xg-ngr1)*grs2)]
        return(els)
      }
      )
      folds <- lapply(1:outerfold,function(i) c(foldsev[[i]],foldsnonev[[i]]))
    }
  return(folds)
}
