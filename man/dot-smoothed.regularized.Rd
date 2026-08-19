\name{.smoothed.regularized}
\alias{.smoothed.regularized}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Selects surrogate markers using an unweighted regularized approach
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Helper function used by select.surrogate with selection = "rr"; this fits an unweighted (ordinary) cross-validated lasso to the treated group to select a subset of candidate markers, refits an unpenalized (debiased) model on the selected markers, and estimates the residual treatment effect and proportion of treatment effect explained for the selected subset.
}
\usage{
.smoothed.regularized(data, surrogate_names, nlambda = 20,
lambda_rule = "lambda.1se", var = FALSE, boot.num = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
data frame of all observations (treated and control), including columns y, treat, and one column per candidate marker.
}
  \item{surrogate_names}{
%%     ~~Describe \code{surrogate_names} here~~
character vector of the candidate marker column names in data.
}
  \item{nlambda}{
%%     ~~Describe \code{nlambda} here~~
the number of penalty values considered along the lasso path, default is 20.
}
  \item{lambda_rule}{
%%     ~~Describe \code{lambda_rule} here~~
the cross-validation rule used to select the penalty; choices are "lambda.1se" or "lambda.min", default is "lambda.1se".
}
  \item{var}{
%%     ~~Describe \code{var} here~~
TRUE or FALSE; indicates whether a bootstrap variance estimate is requested, default is FALSE.
}
  \item{boot.num}{
%%     ~~Describe \code{boot.num} here~~
number of bootstrap resamples used when var = TRUE, default is 200.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list is returned:
 \item{selection }{a named numeric vector of the penalized lasso coefficients of the selected markers.}
 \item{delta }{the estimated overall treatment effect.}
 \item{delta.s }{the estimated residual treatment effect for the selected markers.}
 \item{pte }{the estimated proportion of treatment effect explained by the selected markers.}
 \item{var }{a list, as returned by .bootstrap_var, of bootstrap standard errors and confidence intervals for delta, delta.s, and pte; NULL if var = FALSE.}
}
\references{
Parast, L., Freestone, G., Tarr, J. (2026). Treatment-Aware Regularization via Multi-Objective Optimization for Surrogate Selection.
}
\author{
%%  ~~who you are~~
Layla Parast
}
\keyword{regression}
\keyword{smooth}
\keyword{robust}
\keyword{internal}
% __ONLY ONE__ keyword per line
