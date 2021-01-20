\name{minML.LA.ridgeGLM}
\alias{minML.LA.ridgeGLM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
-log(ML) of ridge penalised GLMs
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Returns the Laplace approximation (LA) of the minus log marginal likelihood of ridge penalised generalised linear models. Note: currently only implemented for linear and logistic regression.
}
\usage{
minML.LA.ridgeGLM(loglambdas, XXblocks, Y, sigmasq = 1, 
                  Xunpen = NULL, intrcpt = F, model, minlam = 0, 
                  opt.sigma = F)
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Laplace approximation of the minus log marginal likelihood for the ridge penalised GLM with model parameters 'loglambdas' and 'sigmasq' (for linear regression).
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
function (loglambdas, XXblocks, Y, sigmasq = 1, Xunpen = NULL, 
    intrcpt = F, model, minlam = 0, opt.sigma = F, offset = NULL) 
{
    if (model != "linear") {
        opt.sigma <- F
        sigmasq <- 1
    }
    n <- length(Y)
    lambdas <- exp(loglambdas) + minlam
    if (opt.sigma) {
        sigmasq <- lambdas[1]
        lambdas <- lambdas[-1]
    }
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
    if (model == "linear") {
        mu <- eta
        W <- diag(rep(1, n))
        Hpen <- Lam0inv - Lam0inv \%*\% solve(solve(W) + Lam0inv, 
            Lam0inv)
        t1 <- dmvnorm(c(Y), mean = eta, sigma = sigmasq * diag(rep(1, 
            n)), log = T)
    }
    else {
        mu <- 1/(1 + exp(-eta))
        W <- diag(c(mu) * c(1 - mu))
        Hpen <- Lam0inv - Lam0inv \%*\% solve(diag(1, n) + W \%*\% 
            Lam0inv, W \%*\% Lam0inv)
        t1 <- sum(Y * log(mu) + (1 - Y) * log(1 - mu))
    }
    t2 <- -1/2/sigmasq * sum(c(Y - mu) * eta)
    t3 <- 1/2 * logdet(diag(rep(1, n)) - W \%*\% Hpen)
    return(-(t1 + t2 + t3))
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
