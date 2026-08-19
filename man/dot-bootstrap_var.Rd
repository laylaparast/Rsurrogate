\name{.bootstrap_var}
\alias{.bootstrap_var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Bootstrap variance estimation for select.surrogate
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Helper function used by select.surrogate when var = TRUE; this resamples both treatment groups with replacement, reruns the entire estimation procedure on each resample, and summarizes the resulting bootstrap distribution of delta, delta.s, and the proportion of treatment effect explained.
}
\usage{
.bootstrap_var(data, fit_once, boot.num)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
data frame of all observations (treated and control), including columns y, treat, and one column per candidate marker.
}
  \item{fit_once}{
%%     ~~Describe \code{fit_once} here~~
a function of a single data frame argument (of the same form as data) that runs the point estimation procedure once and returns a list with elements delta, delta.s, and pte.
}
  \item{boot.num}{
%%     ~~Describe \code{boot.num} here~~
number of bootstrap resamples.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list is returned:
 \item{delta_se }{bootstrap standard error of delta.}
 \item{delta_ci }{a vector of size 2; the 2.5th and 97.5th percentiles of the bootstrap replicates of delta.}
 \item{delta.s_se }{bootstrap standard error of delta.s.}
 \item{delta.s_ci }{a vector of size 2; the 2.5th and 97.5th percentiles of the bootstrap replicates of delta.s.}
 \item{R.s_se }{bootstrap standard error of the proportion of treatment effect explained.}
 \item{R.s_ci }{a vector of size 2; the 2.5th and 97.5th percentiles of the bootstrap replicates of the proportion of treatment effect explained.}
 \item{boot.num }{the number of bootstrap resamples used.}
}
\references{
Parast, L., Freestone, G., Tarr, J. (2026). Treatment-Aware Regularization via Multi-Objective Optimization for Surrogate Selection.
}
\author{
%%  ~~who you are~~
Layla Parast
}
\keyword{robust}
\keyword{internal}
% __ONLY ONE__ keyword per line
