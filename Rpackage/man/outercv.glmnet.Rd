\name{outercv.glmnet}
\alias{outercv.glmnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
outercv.glmnet(Y, X, folds, type.measure = "MSE", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{folds}{
%%     ~~Describe \code{folds} here~~
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
function (Y, X, folds, type.measure = "MSE", ...) 
{
    more.args <- list(...)
    n <- dim(X)[1]
    p <- dim(X)[2]
    if (!is.element("model", names(more.args))) {
        if (all(is.element(Y, c(0, 1))) || is.factor(Y)) {
            model <- "logistic"
        }
        else if (all(is.numeric(Y)) & !(is.matrix(Y) && dim(Y)[2] == 
            2)) {
            model <- "linear"
        }
        else {
            model <- "cox"
        }
    }
    else {
        model <- more.args$model
    }
    if (is.numeric(folds)) {
        folds2 <- produceFolds(n, folds, Y, balance = balance, 
            model = model)
    }
    else {
        folds2 <- folds
    }
    nfolds <- length(folds2)
    Res <- list()
    df <- data.frame()
    dfGrps <- data.frame()
    for (i in 1:nfolds) {
        tic <- proc.time()[[3]]
        Res[[i]] <- do.call(cv.glmnet, args = c(list(x = X[-folds2[[i]], 
            ], y = Y[-folds2[[i]]]), more.args))
        toc <- proc.time()[[3]]
        Res[[i]]$time <- toc - tic
        temp <- coef.glmnet(Res[[i]], s = Res[[i]]$lambda.min, 
            exact = T, x = X[-folds2[[i]], ], y = Y[-folds2[[i]]])
        Res[[i]]$betaL <- temp[-1]
        fit <- do.call(glmnet, args = c(list(x = X[-folds2[[i]], 
            ], y = Y[-folds2[[i]]]), more.args))
        Ypred <- predict.glmnet(fit, newx = X[folds2[[i]], ], 
            s = Res[[i]]$lambda.min, exact = T, x = X[-folds2[[i]], 
                ], y = Y[-folds2[[i]]])
        df2 <- data.frame(Ypred = c(Ypred))
        df2$Method <- "lasso"
        df2$NumberSelectedVars <- sum(Res[[i]]$betaL != 0)
        df2$Fold <- i
        df2$Sample <- folds2[[i]]
        df2$Time <- Res[[i]]$time
        df2$Truth <- Y[folds2[[i]]]
        df <- rbind(df, df2)
        df3 <- data.frame(`Penalty parameter` = Res[[i]]$lambda.min, 
            Group = 1)
        df3$Method <- "lasso"
        df3$Fold <- i
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
        dfCVM <- df \%>\% group_by(Method, Fold) \%>\% summarise(CVM = mean((Ypred - 
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
