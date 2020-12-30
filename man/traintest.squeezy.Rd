\name{traintest.squeezy}
\alias{traintest.squeezy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
traintest.squeezy(Y, X, Y2, X2, type.measure = "MSE", multi_grouping = F, grouping, args.ecpc = NULL, ecpcinit = T, ncores = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Y2}{
%%     ~~Describe \code{Y2} here~~
}
  \item{X2}{
%%     ~~Describe \code{X2} here~~
}
  \item{type.measure}{
%%     ~~Describe \code{type.measure} here~~
}
  \item{multi_grouping}{
%%     ~~Describe \code{multi_grouping} here~~
}
  \item{grouping}{
%%     ~~Describe \code{grouping} here~~
}
  \item{args.ecpc}{
%%     ~~Describe \code{args.ecpc} here~~
}
  \item{ecpcinit}{
%%     ~~Describe \code{ecpcinit} here~~
}
  \item{ncores}{
%%     ~~Describe \code{ncores} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Y, X, Y2, X2, type.measure = "MSE", multi_grouping = F, 
    grouping, args.ecpc = NULL, ecpcinit = T, ncores = 1, ...) 
{
    args.squeezy <- list(...)
    if (is.list(X)) 
        nSim <- length(X)
    else {
        nSim <- 1
        X <- list(X)
        Y <- list(Y)
        X2 <- list(X2)
        Y2 <- list(Y2)
    }
    if (!is.element("model", c(names(args.ecpc), names(args.squeezy)))) {
        if (all(is.element(Y[[1]], c(0, 1))) || is.factor(Y[[1]])) {
            model <- "logistic"
        }
        else if (all(is.numeric(Y[[1]])) & !(is.matrix(Y[[1]]) && 
            dim(Y[[1]])[2] == 2)) {
            model <- "linear"
        }
        else {
            model <- "cox"
        }
    }
    else if (is.element("model", names(args.squeezy))) {
        model <- args.squeezy$model
    }
    else {
        model <- args.ecpc$model
    }
    Res <- list()
    df <- data.frame()
    dfGrps <- data.frame()
    if (ncores == 1) {
        for (i in 1:nSim) {
            n <- dim(X[[i]])[1]
            p <- dim(X[[i]])[2]
            if (!multi_grouping) 
                groupingTemp <- grouping
            else groupingTemp <- grouping[[i]]
            grpsno <- 1:length(groupingTemp)
            if (ecpcinit) {
                tic <- proc.time()[[3]]
                res <- do.call(ecpc, args = c(list(X = X[[i]], 
                  Y = Y[[i]], X2 = X2[[i]], Y2 = Y2[[i]], groupings = list(groupingTemp), 
                  hypershrinkage = "none", postselection = F, 
                  model = model), args.ecpc))
                Res[[i]] <- do.call(squeezy, args = c(list(X = X[[i]], 
                  Y = Y[[i]], X2 = X2[[i]], Y2 = Y2[[i]], grouping = groupingTemp, 
                  fit.ecpc = res), args.squeezy))
                toc <- proc.time()[[3]]
                Res[[i]]$time <- toc - tic
            }
            else {
                tic <- proc.time()[[3]]
                Res[[i]] <- do.call(squeezy, args = c(list(X = X[[i]], 
                  Y = Y[[i]], X2 = X2[[i]], Y2 = Y2[[i]], grouping = groupingTemp), 
                  args.squeezy))
                toc <- proc.time()[[3]]
                Res[[i]]$time <- toc - tic
            }
            df2 <- data.frame(Ypred = c(Res[[i]]$YpredMR, Res[[i]]$YpredApprox))
            df2$Method <- rep(c("multiridge", "squeezy_EN"), 
                each = length(Y2[[i]]))
            df2$NumberSelectedVars <- rep(c(sum(Res[[i]]$betaMR != 
                0), sum(Res[[i]]$betaApprox != 0)), each = length(Y2[[i]]))
            df2$Dataset <- i
            df2$Sample <- rep(1:length(Y2[[i]]), 2)
            df2$Time <- Res[[i]]$time
            df2$Truth <- rep(Y2[[i]], 2)
            df <- rbind(df, df2)
            df3 <- data.frame(Group = rep(grpsno, 2), `Group parameter` = c(Res[[i]]$tauMR, 
                Res[[i]]$tauApprox), `Penalty parameter` = c(Res[[i]]$lambdaMR, 
                Res[[i]]$lambdaApprox))
            df3$Method <- rep(c("multiridge", "squeezy_EN"), 
                each = length(grpsno))
            df3$Dataset <- i
            dfGrps <- rbind(dfGrps, df3)
        }
    }
    else {
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
        finalMatrix <- foreach(i = 1:nfolds, .combine = rbind, 
            .packages = c("glmnet", "penalized", "mvtnorm", "gglasso", 
                "Matrix", "Rsolnp", "ecpc")) \%dopar\% {
            n <- dim(X[[i]])[1]
            p <- dim(X[[i]])[2]
            if (!multi_grouping) 
                groupingTemp <- grouping
            else groupingTemp <- grouping[[i]]
            grpsno <- 1:length(groupingTemp)
            if (ecpcinit) {
                tic <- proc.time()[[3]]
                res <- do.call(ecpc, args = c(list(X = X[[i]], 
                  Y = Y[[i]], X2 = X2[[i]], Y2 = Y2[[i]], groupings = list(groupingTemp), 
                  hypershrinkage = "none", postselection = F, 
                  model = model), args.ecpc))
                Res[[i]] <- do.call(squeezy, args = c(list(X = X[[i]], 
                  Y = Y[[i]], X2 = X2[[i]], Y2 = Y2[[i]], grouping = groupingTemp, 
                  fit.ecpc = res), args.squeezy))
                toc <- proc.time()[[3]]
                Res[[i]]$time <- toc - tic
            }
            else {
                tic <- proc.time()[[3]]
                Res[[i]] <- do.call(squeezy, args = c(list(X = X[[i]], 
                  Y = Y[[i]], X2 = X2[[i]], Y2 = Y2[[i]], grouping = groupingTemp), 
                  args.squeezy))
                toc <- proc.time()[[3]]
                Res[[i]]$time <- toc - tic
            }
            df2 <- data.frame(Ypred = c(Res[[i]]$YpredMR, Res[[i]]$YpredApprox))
            df2$Method <- rep(c("multiridge", "squeezy_EN"), 
                each = length(Y2[[i]]))
            df2$NumberSelectedVars <- rep(c(sum(Res[[i]]$betaMR != 
                0), sum(Res[[i]]$betaApprox != 0)), each = length(Y2[[i]]))
            df2$Dataset <- i
            df2$Sample <- rep(1:length(Y2[[i]]), 2)
            df2$Time <- Res[[i]]$time
            df2$Truth <- rep(Y2[[i]], 2)
            df3 <- data.frame(Group = rep(grpsno, 2), `Group parameter` = c(Res[[i]]$tauMR, 
                Res[[i]]$tauApprox), `Penalty parameter` = c(Res[[i]]$lambdaMR, 
                Res[[i]]$lambdaApprox))
            df3$Method <- rep(c("multiridge", "squeezy_EN"), 
                each = length(grpsno))
            df3$Dataset <- i
            list(Res = Res, df = df2, dfGrps = df3)
        }
        Res <- lapply(1:nfolds, function(i) finalMatrix[i, 1][[1]][[i]])
        df2 <- lapply(1:nfolds, function(i) finalMatrix[i, 2][[1]])
        dfGrps2 <- lapply(1:nfolds, function(i) finalMatrix[i, 
            3][[1]])
        df <- df2[[1]]
        for (i in 2:nfolds) df <- rbind(df, df2[[i]])
        dfGrps <- dfGrps2[[1]]
        for (i in 2:nfolds) dfGrps <- rbind(dfGrps, dfGrps2[[i]])
        stopCluster(cl)
        rm(cl)
    }
    if (is.factor(df$Truth)) {
        warning("Response Y given as factor, transformed to numeric to compute AUC")
        if (!silent) 
            print(levels(df$Truth)[1], "transformed to", 0)
        if (!silent) 
            print(levels(df$Truth)[2], "transformed to", 1)
        df$Truth <- as.numeric(df$Truth) - 1
    }
    if (type.measure == "MSE") {
        dfCVM <- df \%>\% group_by(Method, Dataset) \%>\% summarise(CVM = mean((Ypred - 
            Truth)^2), Type = "MSE", NumberSelectedVars = mean(NumberSelectedVars)) \%>\% 
            ungroup()
    }
    else if (type.measure == "AUC") {
        dfROC <- data.frame()
        for (i in levels(df$Method)) {
            temp <- data.frame()
            cutoffs <- rev(seq(0, 1, by = 0.001))
            rocGR <- roc(probs = df$Ypred[df$Method == i], true = df$Truth[df$Method == 
                i], cutoffs = cutoffs)
            temp <- data.frame(FPR = rocGR[1, ], TPR = rocGR[2, 
                ], Accuracy = rocGR[3, ])
            temp$Method <- i
            temp$AUC <- c(auc(rocGR))
            temp$NumberSelectedVars <- mean(df$NumberSelectedVars[df$Method == 
                i])
            dfROC <- rbind(dfROC, temp)
        }
        dfCVM <- dfROC \%>\% group_by(Method) \%>\% summarise(CVM = mean(AUC), 
            Type = "AUC", NumberSelectedVars = mean(NumberSelectedVars)) \%>\% 
            ungroup()
    }
    else {
        warning(paste("The type of measure", type.measure, "is not yet supported."))
    }
    return(list(Res = Res, dfPred = df, dfGrps = dfGrps, dfCVM = dfCVM))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
