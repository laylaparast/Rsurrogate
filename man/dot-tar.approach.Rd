\name{.tar.approach}
\alias{.tar.approach}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Selects surrogate markers using treatment-aware regularization
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Helper function used by select.surrogate with selection = "TAR"; this fits a weighted lasso to the treated group, in which each candidate marker's penalty is weighted by the inverse of its estimated standardized treatment effect, and selects the penalty via cross-validation within the treated group to be the sparsest point along the weighted lasso path whose estimated proportion of treatment effect explained falls no more than epsilon below that of the full candidate set. The selected markers are then refit using an unpenalized (debiased) model, and the residual treatment effect and proportion of treatment effect explained are estimated for the selected subset.
}
\usage{
.tar.approach(data, folds = 5, epsilon = 0.05, nlambda = 20, surrogate_names,
c.constant = 1, small.n = 1e-6, var = FALSE, boot.num = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
data frame of all observations (treated and control), including columns y, treat, and one column per candidate marker.
}
  \item{folds}{
%%     ~~Describe \code{folds} here~~
the number of cross-validation folds, within the treated group, used to select the penalty, default is 5.
}
  \item{epsilon}{
%%     ~~Describe \code{epsilon} here~~
the tolerance used to select the penalty: the selected subset is the sparsest point along the weighted lasso path whose estimated proportion of treatment effect explained falls no more than epsilon below that of the full candidate set, default is 0.05.
}
  \item{nlambda}{
%%     ~~Describe \code{nlambda} here~~
the number of penalty values considered along the weighted lasso path, default is 20.
}
  \item{surrogate_names}{
%%     ~~Describe \code{surrogate_names} here~~
character vector of the candidate marker column names in data.
}
  \item{c.constant}{
%%     ~~Describe \code{c.constant} here~~
a constant controlling the strength of the treatment-aware penalty weighting, default is 1.
}
  \item{small.n}{
%%     ~~Describe \code{small.n} here~~
a small positive constant added to the estimated treatment effect magnitude of each candidate marker before inverting it to form the penalty weights, to avoid dividing by a value at or near zero, default is 1e-6.
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
 \item{selection }{a named numeric vector of the penalized (weighted) lasso coefficients of the selected markers.}
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
