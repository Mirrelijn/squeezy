\name{dminML.LA.ridgeGLM}
\alias{dminML.LA.ridgeGLM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Partial derivatives of -log(ML) of ridge penalised GLMs
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Returns the partial derivatives (w.r.t. 'loglambdas') of the minus log Laplace approximation (LA) of the marginal likelihood of ridge penalised generalised linear models. Note: currently only implemented for linear and logistic regression.
}
\usage{
dminML.LA.ridgeGLM(loglambdas, XXblocks, Y, sigmasq = 1, 
                  Xunpen = NULL, intrcpt = F, model, minlam = 0, 
                  opt.sigma = F, offset = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loglambdas}{
%%     ~~Describe \code{loglambdas} here~~
Logarithm of the ridge penalties as returned by ecpc or squeezy; Gx1 vector.
}
  \item{XXblocks}{
%%     ~~Describe \code{XXblocks} here~~
List of sample covariance matrices X_g \%*\% t(X_g) for groups g = 1,..., G.
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
Response data; n-dimensional vector (n: number of samples) for linear and logistic outcomes.
}
  \item{sigmasq}{
%%     ~~Describe \code{sigmasq} here~~
(linear model only) Noise level (Y~N(X*beta,sd=sqrt(sigmasq))).
}
  \item{Xunpen}{
%%     ~~Describe \code{Xunpen} here~~
Unpenalised variables; nxp_1-dimensional matrix for p_1 unpenalised variables.
}
  \item{intrcpt}{
%%     ~~Describe \code{intrcpt} here~~
Should an intercept be included? Set to TRUE by default.
}
  \item{model}{
%%     ~~Describe \code{model} here~~
Type of model for the response; linear or logistic. 

}
  \item{minlam}{
%%     ~~Describe \code{minlam} here~~
Minimum value of lambda that is added to exp(loglambdas); set to 0 as default.
}
  \item{opt.sigma}{
%%     ~~Describe \code{opt.sigma} here~~
(linear model only) T/F if log(sigmasq) is given as first argument of 
  loglambdas for optimisation purposes
}
  \item{offset}{
%%     ~~Describe \code{offset} here~~t
(optional) n-dimensional vector with offset for linear predictors. 
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
Partial derivatives of the Laplace approximation of the minus log marginal likelihood to the model parameters 'loglambdas';

For opt.sigma=FALSE: Gx1-dimensional vector for the G log(group ridge penalties). 

For opt.sigma=TRUE (linear model only): (G+1)x1-dimensional vector for the partial derivative to log(sigmasq) (first element) and for the G log(group ridge penalties). 

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
function (loglambdas, XXblocks, Y, sigmasq = 1, Xunpen = NULL, 
    intrcpt = F, model, minlam = 0, opt.sigma = F, offset = NULL) 
{
    if (model != "linear") 
        opt.sigma <- F
    n <- length(Y)
    lambdas <- exp(loglambdas) + minlam
    if (opt.sigma) {
        sigmasq <- lambdas[1]
        lambdas <- lambdas[-1]
    }
    rhos <- log(lambdas)
    Lam0inv <- SigmaBlocks(XXblocks, lambdas)
    if (!is.null(Xunpen)) {
        fit <- IWLSridge(XXT = Lam0inv, Y = Y, model = model, 
            intercept = intrcpt, X1 = Xunpen, maxItr = 500, eps = 10^-12, 
            offset = offset)
    }
    else {
        fit <- IWLSridge(XXT = Lam0inv, Y = Y, model = model, 
            intercept = intrcpt, maxItr = 500, eps = 10^-12, 
            offset = offset)
    }
    eta <- fit$etas + fit$eta0
    if (intrcpt) 
        Xunpen <- cbind(Xunpen, rep(1, n))
    if (model == "linear") {
        mu <- eta
        W <- diag(rep(1, n))
        dwdeta <- rep(0, n)
    }
    else if (model == "logistic") {
        mu <- 1/(1 + exp(-eta))
        W <- diag(c(mu) * c(1 - mu))
        dwdeta <- c(mu * (1 - mu) * (1 - 2 * mu))
    }
    else {
        stop(print(paste("Only model type linear and logistic supported")))
    }
    Hpen <- Lam0inv - Lam0inv \%*\% solve(solve(W) + Lam0inv, Lam0inv)
    if (!is.null(Xunpen)) {
        WP1W <- diag(rep(1, n)) - Xunpen \%*\% solve(t(Xunpen) \%*\% 
            W \%*\% Xunpen, t(Xunpen) \%*\% W)
    }
    else {
        WP1W <- diag(rep(1, n))
    }
    L2 <- WP1W \%*\% (diag(rep(1, n)) - t(solve(diag(1, n) + W \%*\% 
        WP1W \%*\% Lam0inv, W \%*\% WP1W \%*\% Lam0inv)))
    Hj <- lapply(1:length(XXblocks), function(i) {
        L2 \%*\% XXblocks[[i]]/lambdas[i]
    })
    if (is.null(Xunpen)) {
        Hres <- .Hpen(diag(W), Lam0inv)
        Hmat <- Hres$Hmat
    }
    else {
        Hres <- .Hunpen(diag(W), Lam0inv, Xunpen)
        Hmat <- Hres$Hmat
    }
    detadrho <- lapply(1:length(rhos), function(j) {
        -exp(rhos[j]) * (diag(rep(1, n)) - Hmat \%*\% W) \%*\% t(Hj[[j]]/lambdas[j]) \%*\% 
            c(Y - mu + W \%*\% eta)
    })
    dWdrho <- lapply(1:length(rhos), function(j) {
        dwdeta * c(detadrho[[j]])
    })
    t1 <- sapply(1:length(rhos), function(j) {
        -1/sigmasq * c((Y - mu)) \%*\% detadrho[[j]]
    })
    t2 <- sapply(1:length(rhos), function(j) {
        1/2/sigmasq * c(-W \%*\% eta + Y - mu) \%*\% detadrho[[j]]
    })
    t3 <- sapply(1:length(rhos), function(j) {
        -0.5 * sum(diag(solve(diag(1, n) + W \%*\% Lam0inv, W \%*\% 
            XXblocks[[j]]/lambdas[j])))
    })
    t4 <- sapply(1:length(rhos), function(j) {
        0.5 * sum(diag(Hpen) * dWdrho[[j]])
    })
    if (!opt.sigma) 
        return((t1 + t2 + t3 + t4))
    ts.1 <- n/2/sigmasq
    ts.2 <- -1/2/sigmasq^2 * t(Y - eta) \%*\% Y
    return(c(c(ts.1 + ts.2) * sigmasq, (t1 + t2 + t3 + t4)))
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
