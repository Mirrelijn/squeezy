\name{traintest.glmnet}
\alias{traintest.glmnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
traintest.glmnet(Y, X, Y2, X2, type.measure = "MSE", ...)
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
function (Y, X, Y2, X2, type.measure = "MSE", ...) 
{
    more.args <- list(...)
    if (is.list(X)) 
        nSim <- length(X)
    else {
        nSim <- 1
        X <- list(X)
        Y <- list(Y)
        X2 <- list(X2)
        Y2 <- list(Y2)
    }
    if (!is.element("model", names(more.args))) {
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
    else {
        model <- more.args$model
    }
    Res <- list()
    df <- data.frame()
    dfGrps <- data.frame()
    for (i in 1:nSim) {
        n <- dim(X[[i]])[1]
        p <- dim(X[[i]])[2]
        tic <- proc.time()[[3]]
        Res[[i]] <- do.call(cv.glmnet, args = c(list(x = X[[i]], 
            y = Y[[i]]), more.args))
        toc <- proc.time()[[3]]
        Res[[i]]$time <- toc - tic
        temp <- coef.glmnet(Res[[i]], s = Res[[i]]$lambda.min, 
            exact = T, x = X[[i]], y = Y[[i]])
        Res[[i]]$betaL <- temp[-1]
        fit <- do.call(glmnet, args = c(list(x = X[[i]], y = Y[[i]]), 
            more.args))
        Ypred <- predict.glmnet(fit, newx = X2[[i]], s = Res[[i]]$lambda.min, 
            exact = T, x = X[[i]], y = Y[[i]])
        df2 <- data.frame(Ypred = c(Ypred))
        df2$Method <- "lasso"
        df2$NumberSelectedVars <- sum(Res[[i]]$betaL != 0)
        df2$Dataset <- i
        df2$Sample <- 1:length(Y2[[i]])
        df2$Time <- Res[[i]]$time
        df2$Truth <- Y2[[i]]
        df <- rbind(df, df2)
        df3 <- data.frame(`Penalty parameter` = Res[[i]]$lambda.min, 
            Group = 1)
        df3$Method <- "lasso"
        df3$Dataset <- i
        dfGrps <- rbind(dfGrps, df3)
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
