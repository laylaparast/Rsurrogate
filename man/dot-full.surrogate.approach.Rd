\name{.full.surrogate.approach}
\alias{.full.surrogate.approach}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Estimates the proportion of treatment effect explained using all candidate markers
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Helper function used by select.surrogate with selection = "none"; this is a thin wrapper around R.s.estimate, run on all candidate markers, that adds a bootstrap variance option in place of R.s.estimate's perturbation-resampling.
}
\usage{
.full.surrogate.approach(data, surrogate_names, var = FALSE, boot.num = 200)
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
 \item{delta }{the estimated overall treatment effect.}
 \item{delta.s }{the estimated residual treatment effect using all candidate markers.}
 \item{pte }{the estimated proportion of treatment effect explained using all candidate markers.}
 \item{var }{a list, as returned by .bootstrap_var, of bootstrap standard errors and confidence intervals for delta, delta.s, and pte; NULL if var = FALSE.}
}
\references{
Parast, L., Freestone, G., Tarr, J. (2026). Treatment-Aware Regularization via Multi-Objective Optimization for Surrogate Selection.

Parast, L., McDermott, M., Tian, L. (2016). Robust estimation of the proportion of treatment effect explained by surrogate marker information. Statistics in Medicine, 35(10):1637-1653.
}
\author{
%%  ~~who you are~~
Layla Parast
}
\keyword{regression}
\keyword{smooth}
\keyword{robust}
\keyword{nonparametric}
\keyword{internal}
% __ONLY ONE__ keyword per line
