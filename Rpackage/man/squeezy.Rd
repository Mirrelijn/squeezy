\name{squeezy}
\alias{squeezy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
squeezy(Y, X, grouping, alpha = 1, model = NULL, X2 = NULL, Y2 = NULL, unpen = NULL, intrcpt = T, method = c("fixed.ecpc", "MML", "MML.noDeriv", "CV"), fold = 10, compareMR = T, selectAIC = F, fit.ecpc = NULL, lambdas = NULL, lambdaglobal = NULL, lambdasinit = NULL, sigmasq = NULL, ecpcinit = T, SANN = F, minlam = 10^-3, standardise_Y = NULL, reCV = NULL, opt.sigma = NULL, resultsAICboth = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{grouping}{
%%     ~~Describe \code{grouping} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{X2}{
%%     ~~Describe \code{X2} here~~
}
  \item{Y2}{
%%     ~~Describe \code{Y2} here~~
}
  \item{unpen}{
%%     ~~Describe \code{unpen} here~~
}
  \item{intrcpt}{
%%     ~~Describe \code{intrcpt} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{fold}{
%%     ~~Describe \code{fold} here~~
}
  \item{compareMR}{
%%     ~~Describe \code{compareMR} here~~
}
  \item{selectAIC}{
%%     ~~Describe \code{selectAIC} here~~
}
  \item{fit.ecpc}{
%%     ~~Describe \code{fit.ecpc} here~~
}
  \item{lambdas}{
%%     ~~Describe \code{lambdas} here~~
}
  \item{lambdaglobal}{
%%     ~~Describe \code{lambdaglobal} here~~
}
  \item{lambdasinit}{
%%     ~~Describe \code{lambdasinit} here~~
}
  \item{sigmasq}{
%%     ~~Describe \code{sigmasq} here~~
}
  \item{ecpcinit}{
%%     ~~Describe \code{ecpcinit} here~~
}
  \item{SANN}{
%%     ~~Describe \code{SANN} here~~
}
  \item{minlam}{
%%     ~~Describe \code{minlam} here~~
}
  \item{standardise_Y}{
%%     ~~Describe \code{standardise_Y} here~~
}
  \item{reCV}{
%%     ~~Describe \code{reCV} here~~
}
  \item{opt.sigma}{
%%     ~~Describe \code{opt.sigma} here~~
}
  \item{resultsAICboth}{
%%     ~~Describe \code{resultsAICboth} here~~
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
function (Y, X, grouping, alpha = 1, model = NULL, X2 = NULL, 
    Y2 = NULL, unpen = NULL, intrcpt = T, method = c("fixed.ecpc", 
        "MML", "MML.noDeriv", "CV"), fold = 10, compareMR = T, 
    selectAIC = F, fit.ecpc = NULL, lambdas = NULL, lambdaglobal = NULL, 
    lambdasinit = NULL, sigmasq = NULL, ecpcinit = T, SANN = F, 
    minlam = 10^-3, standardise_Y = NULL, reCV = NULL, opt.sigma = NULL, 
    resultsAICboth = F) 
{
    groupings <- list(grouping)
    n <- dim(X)[1]
    p <- dim(X)[2]
    if (!is.null(X2)) 
        n2 <- dim(X2)[1]
    if (length(method) == 4) {
        method <- "MML"
    }
    switch(method, fixed.ecpc = {
        if (is.null(fit.ecpc)) stop("provide ecpc fit results")
        if (is.null(lambdas)) lambdas <- fit.ecpc$sigmahat/(fit.ecpc$gamma * 
            fit.ecpc$tauglobal)
        if (is.null(lambdaglobal)) lambdaglobal <- fit.ecpc$sigmahat/fit.ecpc$tauglobal
        if (is.null(sigmasq)) sigmasq <- fit.ecpc$sigmahat
        if (is.null(standardise_Y)) standardise_Y <- T
        if (is.null(reCV)) reCV <- T
        if (is.null(opt.sigma)) opt.sigma <- F
    }, MML = {
        if (!is.null(fit.ecpc)) {
            if (is.null(lambdasinit)) lambdasinit <- fit.ecpc$sigmahat/(fit.ecpc$gamma * 
                fit.ecpc$tauglobal)
            if (is.null(lambdaglobal)) lambdaglobal <- fit.ecpc$sigmahat/fit.ecpc$tauglobal
            if (is.null(sigmasq)) sigmasq <- fit.ecpc$sigmahat
        }
        if (is.null(standardise_Y)) standardise_Y <- F
        if (is.null(opt.sigma)) opt.sigma <- T
        if (is.null(reCV)) {
            reCV <- F
            if (!opt.sigma) reCV <- T
        }
    }, MML.noDeriv = {
        if (!is.null(fit.ecpc)) {
            if (is.null(lambdasinit)) lambdasinit <- fit.ecpc$sigmahat/(fit.ecpc$gamma * 
                fit.ecpc$tauglobal)
            if (is.null(lambdaglobal)) lambdaglobal <- fit.ecpc$sigmahat/fit.ecpc$tauglobal
            if (is.null(sigmasq)) sigmasq <- fit.ecpc$sigmahat
        }
        if (is.null(standardise_Y)) standardise_Y <- F
        if (is.null(opt.sigma)) opt.sigma <- T
        if (is.null(reCV)) {
            reCV <- F
            if (!opt.sigma) reCV <- T
        }
    }, CV = {
        if (!is.null(fit.ecpc)) {
            if (is.null(lambdasinit)) lambdasinit <- fit.ecpc$sigmahat/(fit.ecpc$gamma * 
                fit.ecpc$tauglobal)
            if (is.null(lambdaglobal)) lambdaglobal <- fit.ecpc$sigmahat/fit.ecpc$tauglobal
            if (is.null(sigmasq)) sigmasq <- fit.ecpc$sigmahat
        }
        if (is.null(standardise_Y)) standardise_Y <- F
        if (is.null(opt.sigma)) opt.sigma <- F
        if (is.null(reCV)) reCV <- T
    }, )
    if (is.null(model)) {
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
    switch(model, linear = {
        fml <- "gaussian"
        sd_y <- sqrt(var(Y) * (n - 1)/n)[1]
        if (standardise_Y) {
            Y <- Y/sd_y
            sd_y_former <- sd_y
            sd_y <- 1
        }
        if (method == "MML") minlam <- max(minlam, 10^-4 * var(Y))
    }, logistic = {
        fml <- "binomial"
        opt.sigma <- F
        standardise_Y <- F
        sd_y <- 1
        sd_y_former <- sd_y
        levelsY <- cbind(c(0, 1), c(0, 1))
        if (!all(is.element(Y, c(0, 1)))) {
            oldLevelsY <- levels(Y)
            levels(Y) <- c("0", "1")
            Y <- as.numeric(Y) - 1
            levelsY <- cbind(oldLevelsY, c(0, 1))
            colnames(levelsY) <- c("Old level names", "New level names")
            if (!is.null(Y2)) {
                levels(Y2) <- c("0", "1")
                Y2 <- as.numeric(Y2) - 1
            }
            if (!silent) print("Y is put in 0/1 format, see levelsY in output for new names")
        }
    }, cox = {
        fml <- "cox"
        opt.sigma <- F
        standardise_Y <- F
        sd_y <- 1
        sd_y_former <- sd_y
        intrcpt <- F
    })
    penfctr <- rep(1, p)
    if (length(unpen) > 0) {
        penfctr[unpen] <- 0
        if (any(unlist(groupings) \%in\% unpen)) {
            warning("Unpenalised covariates removed from grouping")
            for (i in 1:length(groupings)) {
                for (j in 1:length(groupings[[i]])) {
                  if (all(groupings[[i]][[j]] \%in\% unpen)) {
                    groupings[[i]][[j]] <- NULL
                  }
                  else {
                    groupings[[i]][[j]] <- groupings[[i]][[j]][!(groupings[[i]][[j]] \%in\% 
                      unpen)]
                  }
                }
            }
        }
    }
    G <- sapply(groupings, length)
    m <- length(G)
    indGrpsGlobal <- list(1:G[1])
    if (m > 1) {
        for (i in 2:m) {
            indGrpsGlobal[[i]] <- (sum(G[1:(i - 1)]) + 1):sum(G[1:i])
        }
    }
    Kg <- lapply(groupings, function(x) (sapply(x, length)))
    i <- unlist(sapply(1:sum(G), function(x) {
        rep(x, unlist(Kg)[x])
    }))
    j <- unlist(unlist(groupings))
    ind <- sparseMatrix(i, j, x = 1)
    Ik <- lapply(1:m, function(i) {
        x <- rep(0, sum(G))
        x[(sum(G[1:i - 1]) + 1):sum(G[1:i])] <- 1
        as.vector(x \%*\% ind)
    })
    Zt <- ind
    if (G[1] > 1) {
        Zt[1:G[1], ] <- t(t(ind[1:G[1], ])/apply(ind[1:G[1], 
            ], 2, sum))
    }
    if (m > 1) {
        for (i in 2:m) {
            if (G[i] > 1) {
                Zt[indGrpsGlobal[[i]], ] <- t(t(ind[indGrpsGlobal[[i]], 
                  ])/apply(ind[indGrpsGlobal[[i]], ], 2, sum))
            }
        }
    }
    if (dim(Zt)[2] < p) 
        Zt <- cbind(Zt, matrix(rep(NaN, (p - dim(Zt)[2]) * sum(G)), 
            c(sum(G), p - dim(Zt)[2])))
    Xxtnd <- do.call(cbind, lapply(groupings[[1]], function(group) {
        t(t(X[, group])/sqrt(Ik[[1]][group]))
    }))
    Kg2 <- c(1, Kg[[1]])
    G2 <- length(Kg2) - 1
    groupxtnd <- lapply(2:length(Kg2), function(i) {
        sum(Kg2[1:(i - 1)]):(sum(Kg2[1:i]) - 1)
    })
    groupxtnd2 <- unlist(sapply(1:G2, function(x) {
        rep(x, Kg2[x + 1])
    }))
    Xunpen <- NULL
    if (sum((1:p) \%in\% unpen) > 0) 
        Xunpen <- X[, (1:p) \%in\% unpen]
    Xbl <- createXblocks(lapply(groupxtnd, function(ind) Xxtnd[, 
        ind]))
    XXbl <- createXXblocks(lapply(groupxtnd, function(ind) Xxtnd[, 
        ind]))
    if (is.null(lambdaglobal)) {
        cvperblock <- fastCV(list(Xxtnd), Y = Y, kfold = fold, 
            fixedfolds = F, X1 = Xunpen, intercept = intrcpt)
        lambda <- cvperblock$lambdas
        lambda[lambda == Inf] <- 10^6
    }
    else {
        lambda <- lambdaglobal
    }
    if (selectAIC) {
        if (method == "CV") {
            leftout <- CVfolds(Y = Y, kfold = fold, nrepeat = 3, 
                fixedfolds = F)
            lambda1groupfit <- optLambdasWrap(penaltiesinit = lambda, 
                XXblocks = list(apply(simplify2array(XXbl), c(1, 
                  2), sum)), Y = Y, folds = leftout, X1 = Xunpen, 
                intercept = intrcpt, score = ifelse(model == 
                  "linear", "mse", "loglik"), model = model, 
                maxItropt2 = 500, reltol = 10^-3)
            lambda <- lambda1groupfit$optpen
        }
        else if (method == "MML.noDeriv") {
            lambda1groupfit <- optLambdas_mgcv(penaltiesinit = lambda, 
                Y = Y, XXblocks = list(apply(simplify2array(XXbl), 
                  c(1, 2), sum)), model = model, reltol = 0.001, 
                maxItropt = 500, tracescore = F, fixedseed = F, 
                optmethod = "Nelder-Mead", sigmasq = sigmasq)
            lambda <- lambda1groupfit$optpen
        }
        else if (method == "MML") {
            sigmahat <- 1
            if (!is.null(sigmasq) & model == "linear") 
                sigmahat <- sigmasq
            if (opt.sigma) {
                lambda1groupfit <- optim(par = c(log(sigmahat), 
                  log(lambda)), fn = minML.LA.ridgeGLM, gr = dminML.LA.ridgeGLM, 
                  method = "BFGS", minlam = minlam, XXblocks = list(apply(simplify2array(XXbl), 
                    c(1, 2), sum)), Y = Y, sigmasq = sigmahat, 
                  model = model, intrcpt = intrcpt, Xunpen = Xunpen, 
                  opt.sigma = opt.sigma)
                sigmasq <- exp(lambda1groupfit$par[1]) + minlam
                lambda <- exp(lambda1groupfit$par[-1]) + minlam
            }
            else {
                lambda1groupfit <- optim(par = log(lambda), fn = minML.LA.ridgeGLM, 
                  gr = dminML.LA.ridgeGLM, method = "BFGS", minlam = minlam, 
                  XXblocks = list(apply(simplify2array(XXbl), 
                    c(1, 2), sum)), Y = Y, sigmasq = sigmahat, 
                  model = model, intrcpt = intrcpt, Xunpen = Xunpen, 
                  opt.sigma = opt.sigma)
                lambda <- exp(lambda1groupfit$par) + minlam
            }
        }
    }
    sigmahat <- 1
    if (model == "linear") {
        if (!is.null(sigmasq)) 
            sigmahat <- sigmasq
        else {
            XXT1grp <- SigmaFromBlocks(XXblocks = list(apply(simplify2array(XXbl), 
                c(1, 2), sum)), lambda)
            if (length(unpen) > 0 | intrcpt) {
                Xunpen2 <- Xunpen
                if (intrcpt) 
                  Xunpen2 <- cbind(Xunpen, rep(1, n))
                if (intrcpt && length(unpen) == 0) {
                  betaunpenML1grp <- sum(Y)/n
                }
                else {
                  temp <- solve(XXT1grp + diag(rep(1, n)), Xunpen2)
                  betaunpenML1grp <- solve(t(Xunpen2) \%*\% temp, 
                    t(temp) \%*\% Y)
                }
                sigmahat <- c(t(Y - Xunpen2 \%*\% betaunpenML1grp) \%*\% 
                  solve(XXT1grp + diag(rep(1, n)), Y - Xunpen2 \%*\% 
                    betaunpenML1grp)/n)
            }
            else {
                sigmahat <- c(t(Y) \%*\% solve(XXT1grp + diag(rep(1, 
                  n)), Y)/n)
            }
        }
    }
    if (selectAIC) {
        sigmahat1group <- sigmahat
    }
    if (is.null(lambdas)) {
        if (is.null(lambdasinit) | !ecpcinit) {
            lambdasinit <- rep(lambda, G)
        }
        if (any(lambdasinit == Inf)) {
            lambdasinit[lambdasinit > 2 * lambda] <- 2 * lambda
        }
        if (method == "CV") {
            leftout <- CVfolds(Y = Y, kfold = fold, nrepeat = 3, 
                fixedfolds = F)
            jointlambdas <- optLambdasWrap(penaltiesinit = lambdasinit, 
                XXblocks = XXbl, Y = Y, folds = leftout, X1 = Xunpen, 
                intercept = intrcpt, score = ifelse(model == 
                  "linear", "mse", "loglik"), model = model, 
                maxItropt2 = 500, reltol = 10^-3, traceCV = F)
            lambdas <- jointlambdas$optpen
        }
        else if (method == "MML.noDeriv") {
            if (ecpcinit) {
                if (SANN) {
                  jointlambdas <- optLambdas_mgcvWrap(penaltiesinit = lambdasinit, 
                    XXblocks = XXbl, Y = Y, model = model, reltol = 1e-04, 
                    maxItropt2 = 1000, tracescore = F, fixedseed = F, 
                    optmethod2 = "Nelder-Mead", sigmasq = sigmahat, 
                    opt.sigma = opt.sigma)
                }
                else {
                  jointlambdas <- optLambdas_mgcv(penaltiesinit = lambdasinit, 
                    XXblocks = XXbl, Y = Y, model = model, reltol = 1e-04, 
                    maxItropt = 1000, tracescore = F, fixedseed = F, 
                    optmethod = "Nelder-Mead", sigmasq = sigmahat, 
                    opt.sigma = opt.sigma)
                }
            }
            else {
                if (SANN) {
                  jointlambdas <- optLambdas_mgcvWrap(penaltiesinit = rep(lambda, 
                    G), XXblocks = XXbl, Y = Y, model = model, 
                    reltol = 1e-04, maxItropt2 = 1000, tracescore = F, 
                    fixedseed = F, optmethod2 = "Nelder-Mead", 
                    sigmasq = sigmahat, opt.sigma = opt.sigma)
                }
                else {
                  jointlambdas <- optLambdas_mgcv(penaltiesinit = rep(lambda, 
                    G), XXblocks = XXbl, Y = Y, model = model, 
                    reltol = 1e-04, maxItropt = 1000, tracescore = F, 
                    fixedseed = F, optmethod = "Nelder-Mead", 
                    sigmasq = sigmahat, opt.sigma = opt.sigma)
                }
            }
            if (opt.sigma) {
                sigmasq <- jointlambdas$optpen[1]
                lambdas <- jointlambdas$optpen[-1]
            }
            else {
                lambdas <- jointlambdas$optpen
            }
        }
        else if (method == "MML") {
            if (opt.sigma) {
                jointlambdas <- optim(par = c(log(sigmahat), 
                  log(lambdasinit)), fn = minML.LA.ridgeGLM, 
                  gr = dminML.LA.ridgeGLM, method = "BFGS", minlam = minlam, 
                  XXblocks = XXbl, Y = Y, opt.sigma = opt.sigma, 
                  model = model, intrcpt = intrcpt, Xunpen = Xunpen)
                sigmasq <- exp(jointlambdas$par[1]) + minlam
                lambdas <- exp(jointlambdas$par[-1]) + minlam
            }
            else {
                jointlambdas <- optim(par = log(lambdasinit), 
                  fn = minML.LA.ridgeGLM, gr = dminML.LA.ridgeGLM, 
                  method = "BFGS", minlam = minlam, XXblocks = XXbl, 
                  Y = Y, sigmasq = sigmahat, model = model, intrcpt = intrcpt, 
                  Xunpen = Xunpen)
                lambdas <- exp(jointlambdas$par) + minlam
            }
        }
    }
    else {
        lambdasinit <- lambdas
    }
    sigmahat <- 1
    if (model == "linear") {
        if (!is.null(sigmasq)) 
            sigmahat <- sigmasq
        else {
            XXT <- SigmaFromBlocks(XXblocks = XXbl, lambdas)
            if (length(unpen) > 0 | intrcpt) {
                Xunpen2 <- Xunpen
                if (intrcpt) 
                  Xunpen2 <- cbind(Xunpen, rep(1, n))
                if (intrcpt && length(unpen) == 0) {
                  betaunpenML <- sum(Y)/n
                }
                else {
                  temp <- solve(XXT + diag(rep(1, n)), Xunpen2)
                  betaunpenML <- solve(t(Xunpen2) \%*\% temp, t(temp) \%*\% 
                    Y)
                }
                sigmahat <- c(t(Y - Xunpen2 \%*\% betaunpenML) \%*\% 
                  solve(XXT + diag(rep(1, n)), Y - Xunpen2 \%*\% 
                    betaunpenML)/n)
            }
            else {
                sigmahat <- c(t(Y) \%*\% solve(XXT + diag(rep(1, 
                  n)), Y)/n)
            }
        }
    }
    MLinit <- minML.LA.ridgeGLM(loglambdas = log(lambdasinit), 
        opt.sigma = F, sigmasq = sigmahat, XXblocks = XXbl, Y = Y, 
        model = model, intrcpt = intrcpt, minlam = 0)
    MLfinal <- minML.LA.ridgeGLM(loglambdas = log(lambdas), opt.sigma = F, 
        sigmasq = sigmahat, XXblocks = XXbl, Y = Y, model = model, 
        intrcpt = intrcpt, minlam = 0)
    if (selectAIC) {
        lambda1group <- lambda
        if (sum((1:p) \%in\% unpen) > 0) {
            AIC1group <- AIC.LA.ridgeGLM(log(lambda1group), XXblocks = list(apply(simplify2array(XXbl), 
                c(1, 2), sum)), Y = Y, sigmasq = sigmahat1group, 
                Xunpen = Xunpen, intrcpt = intrcpt, model = model)
            AICmultigroup <- AIC.LA.ridgeGLM(log(lambdas), XXblocks = XXbl, 
                Y = Y, sigmasq = sigmahat, Xunpen = Xunpen, intrcpt = intrcpt, 
                model = model)
        }
        else {
            AIC1group <- AIC.LA.ridgeGLM(log(lambda1group), XXblocks = list(apply(simplify2array(XXbl), 
                c(1, 2), sum)), Y = Y, sigmasq = sigmahat1group, 
                intrcpt = intrcpt, model = model)
            AICmultigroup <- AIC.LA.ridgeGLM(log(lambdas), XXblocks = XXbl, 
                Y = Y, sigmasq = sigmahat, intrcpt = intrcpt, 
                model = model)
        }
        if (AIC1group <= AICmultigroup) {
            lambdasNotOptimalAIC <- lambdas
            sigmahatNotOptimalAIC <- sigmahat
            lambdas <- rep(lambda1group, G)
            sigmahat <- sigmahat1group
            modelbestAIC <- "onegroup"
        }
        else {
            lambdasNotOptimalAIC <- lambda1group
            sigmahatNotOptimalAIC <- sigmahat1group
            modelbestAIC <- "multigroup"
        }
    }
    tauglobal <- sigmahat/lambda
    gamma <- lambda/lambdas
    lambdap <- lambda/(as.vector(gamma \%*\% Zt))
    lambdap[lambdap < 0] <- Inf
    if (compareMR) {
        XXT <- SigmaFromBlocks(XXbl, penalties = lambdas)
        if (sum((1:p) \%in\% unpen) > 0) {
            fit <- IWLSridge(XXT, Y = Y, model = model, intercept = intrcpt, 
                X1 = X[, (1:p) \%in\% unpen])
        }
        else {
            fit <- IWLSridge(XXT, Y = Y, model = model, intercept = intrcpt)
        }
        betafit <- betasout(fit, Xblocks = Xbl, penalties = lambdas)
        a0MR <- 0
        if (intrcpt) 
            a0MR <- c(betafit[[1]][1])
        betaMR <- rep(0, p)
        betaMR[(1:p) \%in\% unpen] <- betafit[[1]][-1]
        for (i in 1:length(groupings[[1]])) {
            betaMR[groupings[[1]][[i]]] <- betaMR[groupings[[1]][[i]]] + 
                betafit[[1 + i]]/sqrt(Ik[[1]][groupings[[1]][[i]]])
        }
        rm(betafit)
    }
    pen <- which(!((1:p) \%in\% unpen))
    if (any(is.nan(sqrt(lambdap[pen])))) {
        browser()
    }
    if (alpha < 1) {
        varFunc <- function(lam, alpha = 0.5, tausq) {
            t1 <- -alpha/(1 - alpha)^(3/2)/sqrt(lam) * exp(dnorm(alpha * 
                sqrt(lam)/sqrt(1 - alpha), log = T) - pnorm(-alpha * 
                sqrt(lam)/sqrt(1 - alpha), log.p = T))
            varBeta <- 1/lam/(1 - alpha) + alpha^2/(1 - alpha)^2 + 
                t1
            f <- varBeta - tausq
            return(f)
        }
        lamEN <- function(alpha, tausq) {
            if (alpha == 0) {
                lamEN <- 1/tausq
            }
            else if (alpha == 1) {
                lamEN <- sqrt(2/tausq)
            }
            else if (tausq < 0.5 & alpha > 0.95) {
                lamEN <- sqrt(2/tausq)
            }
            else if (tausq < 10^-7) {
                a <- tausq^2 * (1 - alpha)^3 - 2 * tausq * alpha^2 * 
                  (1 - alpha)
                b <- alpha^2 - 2 * tausq * (1 - alpha)^2
                c <- 1 - alpha
                D <- b^2 - 4 * a * c
                lamEN <- -b - sqrt(D)/2/a
                if (sign(lamEN) != 1) 
                  lamEN <- -b + sqrt(D)/2/a
            }
            else {
                lb <- min(1/tausq, sqrt(2/tausq))
                ub <- max(1/tausq, sqrt(2/tausq))
                temp <- try(uniroot(varFunc, c(0.9 * lb, 1.1 * 
                  ub), alpha = alpha, tausq = tausq, tol = 10^-6))
                if (class(temp)[1] == "try-error") {
                  lamEN <- sqrt(2/tausq)
                }
                else {
                  lamEN <- temp$root
                }
            }
            return(lamEN)
        }
        alphat <- 1/(1 + sd_y * (1 - alpha)/alpha)
        lambdasEN <- sapply(sigmahat/lambdas, function(tau) lamEN(tausq = tau, 
            alpha = alpha))
        tauEN <- 1/lambdasEN
        lambdasENhat <- lambdasEN * sigmahat * ((1 - alpha) + 
            alpha/sd_y)
        uniqueTaus <- unique(as.vector(c(sigmahat/lambdas) \%*\% 
            Zt))
        if (dim(X)[2] == dim(Xxtnd)[2]) {
            lambdap <- (as.vector(lambdasENhat \%*\% Zt))
        }
        else {
            lambdasENunique <- sapply(uniqueTaus, function(tau) lamEN(tausq = tau, 
                alpha = alpha))
            lambdasENuniquehat <- lambdasENunique * sigmahat * 
                ((1 - alpha) + alpha/sd_y)
            indUnique <- sapply(as.vector(c(sigmahat/lambdas) \%*\% 
                Zt), function(x) which(uniqueTaus == x))
            lambdap <- lambdasENuniquehat[indUnique]
        }
        penfctr2 <- penfctr
        penfctr2[pen] <- lambdap[pen]
        not0 <- which(penfctr2 != Inf)
        lambdaEN <- sum(penfctr2[not0])/length(penfctr2[not0])
        penfctr2 <- penfctr2/lambdaEN
        if (model == "cox") {
            glmGR <- glmnet(X[, not0], Y, alpha = alphat, family = fml, 
                standardize = F, penalty.factor = penfctr2[not0], 
                thresh = 10^-10)
        }
        else {
            glmGR <- glmnet(X[, not0], Y, alpha = alphat, family = fml, 
                intercept = intrcpt, standardize = F, penalty.factor = penfctr2[not0], 
                thresh = 10^-10)
        }
        if (reCV) {
            glmGR.cv <- cv.glmnet(X[, not0], Y, alpha = alphat, 
                family = fml, intercept = intrcpt, standardize = F, 
                penalty.factor = penfctr2[not0], thresh = 10^-10)
            sopt <- glmGR.cv$lambda.min
        }
        else sopt <- lambdaEN/n * sd_y
        temp <- coef(glmGR, s = sopt, exact = T, x = X[, not0], 
            y = Y, alpha = alphat, penalty.factor = penfctr2[not0], 
            family = fml, intercept = intrcpt)
        beta <- rep(0, p)
        beta[not0] <- temp[-1]
        a0 <- temp[1]
    }
    else if (alpha == 1) {
        tauEN <- sqrt(sigmahat/2/lambdas)
        lambdasEN <- sqrt(sigmahat)/sd_y * sqrt(2 * lambdas)
        lambdap <- lambda/(as.vector(gamma \%*\% Zt))
        lambdap[lambdap < 0] <- Inf
        penfctr2 <- penfctr
        penfctr2[pen] <- sqrt(sigmahat)/sd_y * sqrt(2 * lambdap[pen])
        not0 <- which(penfctr2 != Inf)
        lambda1 <- sum(penfctr2[not0])/length(penfctr2[not0])
        penfctr2 <- penfctr2/lambda1
        if (model == "cox") {
            glmGR <- glmnet(X[, not0], Y, alpha = 1, family = fml, 
                standardize = F, penalty.factor = penfctr2[not0], 
                thresh = 10^-10)
        }
        else {
            glmGR <- glmnet(X[, not0], Y, alpha = 1, family = fml, 
                intercept = intrcpt, standardize = F, penalty.factor = penfctr2[not0], 
                thresh = 10^-10)
        }
        if (reCV) {
            glmGR.cv <- cv.glmnet(X[, not0], Y, alpha = 1, family = fml, 
                intercept = intrcpt, standardize = F, penalty.factor = penfctr2[not0], 
                thresh = 10^-10)
            sopt <- glmGR.cv$lambda.min
        }
        else sopt <- lambda1/n * sd_y
        temp <- coef(glmGR, s = sopt, exact = T, x = X[, not0], 
            y = Y, penalty.factor = penfctr2[not0], family = fml, 
            intercept = intrcpt)
        beta <- rep(0, p)
        beta[not0] <- temp[-1]
        a0 <- temp[1]
    }
    else {
        stop("alpha should be between 0 and 1")
    }
    if (standardise_Y) {
        beta <- beta * sd_y_former
        a0 <- a0 * sd_y_former
    }
    if (!is.null(X2)) {
        if (compareMR) {
            if (model == "linear") {
                X2c <- cbind(X2, rep(1, n2))
                YpredMR <- X2c \%*\% c(betaMR, a0MR)
                MSEMR <- sum((YpredMR - Y2)^2)/n2
            }
            if (model == "logistic") {
                X2c <- cbind(X2, rep(1, n2))
                YpredMR <- 1/(1 + exp(-X2c \%*\% c(betaMR, a0MR)))
                MSEMR <- sum((YpredMR - Y2)^2)/n2
            }
            else if (model == "cox") {
                expXb <- exp(X \%*\% c(betaMR))
                h0 <- sapply(1:length(Y[, 1]), function(i) {
                  Y[i, 2]/sum(expXb[Y[, 1] >= Y[i, 1]])
                })
                H0 <- sapply(Y2[, 1], function(Ti) {
                  sum(h0[Y[, 1] <= Ti])
                })
                YpredMR <- H0 * exp(X2 \%*\% betaMR)
                MSEMR <- sum((YpredMR - Y2[, 2])^2)/n2
            }
        }
        if (model == "linear") {
            X2c <- cbind(X2, rep(1, n2))
            YpredApprox <- X2c \%*\% c(beta, a0)
            MSEApprox <- sum((YpredApprox - Y2)^2)/n2
        }
        if (model == "logistic") {
            X2c <- cbind(X2, rep(1, n2))
            YpredApprox <- 1/(1 + exp(-X2c \%*\% c(beta, a0)))
            MSEApprox <- sum((YpredApprox - Y2)^2)/n2
        }
        else if (model == "cox") {
            expXb <- exp(X \%*\% c(beta))
            h0 <- sapply(1:length(Y[, 1]), function(i) {
                Y[i, 2]/sum(expXb[Y[, 1] >= Y[i, 1]])
            })
            H0 <- sapply(Y2[, 1], function(Ti) {
                sum(h0[Y[, 1] <= Ti])
            })
            YpredApprox <- H0 * exp(X2 \%*\% beta)
            MSEApprox <- sum((YpredApprox - Y2[, 2])^2)/n2
        }
    }
    output <- list(betaApprox = beta, a0Approx = a0, tauApprox = tauEN, 
        lambdaApprox = lambdasEN, tauMR = sigmahat/lambdas, lambdaMR = lambdas, 
        lambdaglobal = lambda, sigmahat = sigmahat, MLinit = MLinit, 
        MLfinal = MLfinal)
    if (compareMR) {
        output$betaMR <- betaMR
        output$a0MR <- a0MR
    }
    if (!is.null(X2)) {
        output$YpredApprox <- YpredApprox
        output$MSEApprox <- MSEApprox
        if (compareMR) {
            output$YpredMR <- YpredMR
            output$MSEMR <- MSEMR
        }
    }
    if (selectAIC) {
        output$AICmodels <- list(multigroup = list(lambdas = lambdas, 
            sigmahat = sigmahat, AIC = AICmultigroup), onegroup = list(lambdas = lambda1group, 
            sigmahat = sigmahat1group, AIC = AIC1group))
        if (modelbestAIC != "multigroup") {
            output$AICmodels$multigroup$lambdas <- lambdasNotOptimalAIC
            output$AICmodels$multigroup$sigmahat <- sigmahatNotOptimalAIC
        }
        if (resultsAICboth) {
            output$AICmodels$onegroup$fit <- squeezy(Y, X, grouping, 
                alpha = alpha, model = model, X2 = X2, Y2 = Y2, 
                unpen = unpen, intrcpt = intrcpt, method = "MML", 
                fold = fold, compareMR = compareMR, selectAIC = F, 
                fit.ecpc = NULL, lambdas = rep(lambda1group, 
                  G), lambdaglobal = lambda1group, sigmasq = sigmahat1group, 
                standardise_Y = standardise_Y, reCV = standardise_Y, 
                resultsAICboth = F)
            output$AICmodels$multigroup$fit <- squeezy(Y, X, 
                grouping, alpha = alpha, model = model, X2 = X2, 
                Y2 = Y2, unpen = unpen, intrcpt = intrcpt, method = "MML", 
                fold = fold, compareMR = compareMR, selectAIC = F, 
                fit.ecpc = NULL, lambdas = output$AICmodels$multigroup$lambdas, 
                lambdaglobal = lambda1group, sigmasq = output$AICmodels$multigroup$sigmahat, 
                standardise_Y = standardise_Y, reCV = standardise_Y, 
                resultsAICboth = F)
        }
        output$modelbestAIC <- modelbestAIC
    }
    return(output)
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
