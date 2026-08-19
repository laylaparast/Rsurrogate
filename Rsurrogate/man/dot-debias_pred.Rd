\name{.debias_pred}
\alias{.debias_pred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Predicted outcome from a debiased refit
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Helper function used by select.surrogate with selection = "TAR"; this computes fitted values on a validation set from the debiased coefficients returned by .debias_ols, and is used to evaluate the cross-validated mean squared error along the weighted lasso path.
}
\usage{
.debias_pred(debias_coef, data_val)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{debias_coef}{
%%     ~~Describe \code{debias_coef} here~~
a one-column matrix of debiased coefficients, as returned by .debias_ols, with row names "Intercept" followed by the names of the selected markers.
}
  \item{data_val}{
%%     ~~Describe \code{data_val} here~~
a matrix or data frame of held-out (validation fold) observations, including a column named Intercept (all ones) and one column for each marker named in debias_coef.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
a numeric vector of predicted outcome values for the observations in data_val.
}
\references{
Parast, L., Freestone, G., Tarr, J. (2026). Treatment-Aware Regularization via Multi-Objective Optimization for Surrogate Selection.
}
\author{
%%  ~~who you are~~
Layla Parast
}
\keyword{regression}
\keyword{robust}
\keyword{internal}
% __ONLY ONE__ keyword per line
