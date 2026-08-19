\name{.debias_p}
\alias{.debias_p}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Counts the number of selected markers
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Helper function used by select.surrogate with selection = "TAR"; this counts the number of nonzero (non-intercept) coefficients in a debiased coefficient vector, used to track the size of the selected set along the weighted lasso path.
}
\usage{
.debias_p(debias_coef, penalty)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{debias_coef}{
%%     ~~Describe \code{debias_coef} here~~
a one-column matrix of debiased coefficients, as returned by .debias_ols, with row names "Intercept" followed by the names of the selected markers.
}
  \item{penalty}{
%%     ~~Describe \code{penalty} here~~
not used; retained for compatibility when this function is applied across a list of debiased coefficient matrices, one per penalty value.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
an integer; the number of selected (non-intercept, nonzero) markers.
}
\references{
Parast, L., Freestone, G., Tarr, J. (2026). Treatment-Aware Regularization via Multi-Objective Optimization for Surrogate Selection.
}
\author{
%%  ~~who you are~~
Layla Parast
}
\keyword{regression}
\keyword{internal}
% __ONLY ONE__ keyword per line
