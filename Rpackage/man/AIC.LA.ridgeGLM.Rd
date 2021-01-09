\name{AIC.LA.ridgeGLM}
\alias{AIC.LA.ridgeGLM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
ML AIC of a multi-group, ridge penalised GLMs
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Compute the AIC for the marginal likelihood (ML) of multi-group, ridge penalised generalised linear models. Note: currently only implemented for linear and logistic regression.
}
\usage{
AIC.LA.ridgeGLM(loglambdas, XXblocks, Y, sigmasq = 1, 
                Xunpen = NULL, intrcpt = F, model, minlam = 0)
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
  \item{AIC }{AIC of the model}
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
