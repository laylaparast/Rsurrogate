\name{.make.delta.s}
\alias{.make.delta.s}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates the residual treatment effect for a selected set of markers
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Helper function used by select.surrogate with selection = "TAR" or selection = "rr"; given a (debiased) coefficient vector for a selected subset of markers, this forms the linear index implied by those coefficients and estimates the residual treatment effect using kernel smoothing, as in the multiple marker, robust approach described in the delta.s.estimate documentation. Unlike delta.s.estimate, values outside the observed support are always extrapolated to the closest non-missing value.
}
\usage{
.make.delta.s(coefs_vec, treat, control)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coefs_vec}{
%%     ~~Describe \code{coefs_vec} here~~
a named numeric vector of (debiased) coefficients for the selected markers; if of length zero (no markers selected), the residual treatment effect reduces to the overall treatment effect.
}
  \item{treat}{
%%     ~~Describe \code{treat} here~~
data frame of treated observations, including a column y and one column for each marker named in coefs_vec.
}
  \item{control}{
%%     ~~Describe \code{control} here~~
data frame of control observations, including a column y and one column for each marker named in coefs_vec.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
a numeric value; the estimated residual treatment effect for the selected markers.
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
