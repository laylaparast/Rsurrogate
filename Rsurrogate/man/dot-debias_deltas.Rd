\name{.debias_deltas}
\alias{.debias_deltas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Residual treatment effect from a debiased refit, for one cross-validation fold
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Helper function used by select.surrogate with selection = "TAR"; this computes the residual treatment effect, evaluated on a held-out cross-validation fold, implied by a debiased coefficient vector, used to check whether a given penalty value keeps the estimated proportion of treatment effect explained within the tolerance epsilon of the full candidate set.
}
\usage{
.debias_deltas(debias_coef, data_val, control_data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{debias_coef}{
%%     ~~Describe \code{debias_coef} here~~
a one-column matrix of debiased coefficients, as returned by .debias_ols, with row names "Intercept" followed by the names of the selected markers.
}
  \item{data_val}{
%%     ~~Describe \code{data_val} here~~
data frame of held-out (validation fold) treated observations, including a column y and one column for each marker named in debias_coef.
}
  \item{control_data}{
%%     ~~Describe \code{control_data} here~~
data frame of control observations, including a column y and one column for each marker named in debias_coef.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
a numeric value; the residual treatment effect for the selected markers, evaluated on data_val and control_data, or NA if debias_coef contains missing values (an unidentified debiased refit).
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
