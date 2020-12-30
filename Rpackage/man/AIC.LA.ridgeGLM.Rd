\name{AIC.LA.ridgeGLM}
\alias{AIC.LA.ridgeGLM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
AIC.LA.ridgeGLM(loglambdas, XXblocks, Y, sigmasq = 1, Xunpen = NULL, intrcpt = F, model, minlam = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loglambdas}{
%%     ~~Describe \code{loglambdas} here~~
}
  \item{XXblocks}{
%%     ~~Describe \code{XXblocks} here~~
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{sigmasq}{
%%     ~~Describe \code{sigmasq} here~~
}
  \item{Xunpen}{
%%     ~~Describe \code{Xunpen} here~~
}
  \item{intrcpt}{
%%     ~~Describe \code{intrcpt} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{minlam}{
%%     ~~Describe \code{minlam} here~~
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
function (loglambdas, XXblocks, Y, sigmasq = 1, Xunpen = NULL, 
    intrcpt = F, model, minlam = 0) 
{
    minLL <- minML.LA.ridgeGLM(loglambdas, XXblocks, Y, sigmasq = sigmasq, 
        Xunpen = Xunpen, intrcpt = intrcpt, model = model, minlam = minlam)
    k <- length(XXblocks)
    if (model == "linear") 
        k <- k + 1
    AIC <- 2 * minLL + 2 * k
    return(AIC)
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
