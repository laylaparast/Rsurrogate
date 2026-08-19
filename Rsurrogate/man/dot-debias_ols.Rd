\name{.debias_ols}
\alias{.debias_ols}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Debiased ordinary least squares refit on lasso-selected markers
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Helper function used by select.surrogate with selection = "TAR"; given a vector of penalized lasso coefficients, this refits an unpenalized (ordinary least squares) model on the treated group using only the markers with nonzero coefficients, to remove the shrinkage bias of the lasso.
}
\usage{
.debias_ols(original_coefs, data_train)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{original_coefs}{
%%     ~~Describe \code{original_coefs} here~~
a named numeric vector of penalized lasso coefficients (excluding the intercept), for one value of the penalty.
}
  \item{data_train}{
%%     ~~Describe \code{data_train} here~~
data frame of treated observations, including a column y (the outcome) and one column for each marker named in original_coefs.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
a one-column matrix of the debiased (unpenalized) coefficients, with row names "Intercept" followed by the names of the selected markers. If the number of selected markers is at least as large as the number of treated observations, so that the debiased refit is not identified, a matrix of the same shape with all entries set to NA is returned instead.
}
\references{
Parast, L., Freestone, G., Tarr, J. (2026). Treatment-Aware Regularization via Multi-Objective Optimization for Surrogate Selection.
}
\author{
%%  ~~who you are~~
Layla Parast
}
\keyword{regression}
\keyword{models}
\keyword{robust}
\keyword{internal}
% __ONLY ONE__ keyword per line
