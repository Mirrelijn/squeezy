\name{squeezy}
\alias{squeezy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Group-adaptive elastic net
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function estimates group-specific elastic net penalties and fits a linear or logistic regression model.
}
\usage{
squeezy(Y, X, grouping, alpha = 1, model = NULL, X2 = NULL, 
        Y2 = NULL, unpen = NULL, intrcpt = T, 
        method = c("ecpcEN", "MML", "MML.noDeriv", "CV"), 
        fold = 10, compareMR = T, selectAIC = F, fit.ecpc = NULL, 
        lambdas = NULL, lambdaglobal = NULL, lambdasinit = NULL, 
        sigmasq = NULL, ecpcinit = T, SANN = F, minlam = 10^-3, 
        standardise_Y = NULL, reCV = NULL, opt.sigma = NULL, 
        resultsAICboth = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
Response data; n-dimensional vector (n: number of samples) for linear and logistic outcomes.
}
  \item{X}{
%%     ~~Describe \code{X} here~~
Observed data; (nxp)-dimensional matrix (p: number of covariates) with each row the observed high-dimensional feature vector of a sample.
}
  \item{grouping}{
%%     ~~Describe \code{grouping} here~~
Co-data grouping; list with G groups.  Each group is a vector containing the indices of the covariates in that group.
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
Elastic net penalty mixing parameter.
}
  \item{model}{
%%     ~~Describe \code{model} here~~
Type of model for the response; linear or logistic.
}
  \item{X2}{
%%     ~~Describe \code{X2} here~~
(optional) Independent observed data for which response is predicted.
}
  \item{Y2}{
%%     ~~Describe \code{Y2} here~~
(optional) Independent response data to compare with predicted response.
}
  \item{unpen}{
%%     ~~Describe \code{unpen} here~~
Unpenalised covariates; vector with indices of covariates that should not be penalised.
}
  \item{intrcpt}{
%%     ~~Describe \code{intrcpt} here~~
Should an intercept be included? Included by default for linear and logistic, excluded for Cox for which the baseline hazard is estimated.
}
  \item{method}{
%%     ~~Describe \code{method} here~~
Which method should be used to estimate the group-specific penalties? Default MML.
}
  \item{fold}{
%%     ~~Describe \code{fold} here~~
Number of folds used in inner cross-validation to estimate (initial) global ridge penalty lambda (if not given).
}
  \item{compareMR}{
%%     ~~Describe \code{compareMR} here~~
T/F to fit the multi-ridge model and return results for comparison.
}
  \item{selectAIC}{
%%     ~~Describe \code{selectAIC} here~~y
T/F to select the single-group model or multi-group model.
}
  \item{fit.ecpc}{
%%     ~~Describe \code{fit.ecpc} here~~
(optional) Model fit obtained by the function ecpc (from the ecpc R-package)
}
  \item{lambdas}{
%%     ~~Describe \code{lambdas} here~~
(optional) Group-specific ridge penalty parameters. If given, these are transformed to elastic net penalties.
}
  \item{lambdaglobal}{
%%     ~~Describe \code{lambdaglobal} here~~
(optional) Global ridge penalty parameter used for initialising the optimisation.
}
  \item{lambdasinit}{
%%     ~~Describe \code{lambdasinit} here~~
(optional) Group-specific ridge penalty parameters used for initialising the optimisation.
}
  \item{sigmasq}{
%%     ~~Describe \code{sigmasq} here~~
(linear model only) If given, noise level is fixed (Y~N(X*beta,sd=sqrt(sigmasq))).
}
  \item{ecpcinit}{
%%     ~~Describe \code{ecpcinit} here~~
T/F for using group-specific ridge penalties as given in `fit.ecpc' for initialising the optimisation.
}
  \item{SANN}{
%%     ~~Describe \code{SANN} here~~
('method'=MML.noDeriv only) T/F to use simulated annealing in optimisation of the ridge penalties.
}
  \item{minlam}{
%%     ~~Describe \code{minlam} here~~
Minimal value of group-specific ridge penalty used in the optimisation.
}
  \item{standardise_Y}{
%%     ~~Describe \code{standardise_Y} here~~
T/F should Y be standardised?
}
  \item{reCV}{
%%     ~~Describe \code{reCV} here~~
T/F should the elastic net penalties be recalibrated by cross-validation of a global rescaling penalty?
}
  \item{opt.sigma}{
%%     ~~Describe \code{opt.sigma} here~~
(linear model only) T/F to optimise sigmasq jointly with the ridge penalties.
}
  \item{resultsAICboth}{
%%     ~~Describe \code{resultsAICboth} here~~
(selectAIC=T only) T/F to return results of both the single-group and multi-group model.
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
}
