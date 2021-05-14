\name{fieller.ci}
\alias{fieller.ci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Constructs Fieller's confidence interval.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Constructs Fieller's confidence interval.
}
\usage{
fieller.ci(perturb.delta.s, perturb.delta, delta.s, delta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{perturb.delta.s}{
%%     ~~Describe \code{perturb.delta.s} here~~
numeric vector; the perturbed values for \eqn{\hat{\Delta}_S}, the residual treatment effect estimate, used in variance estimation and confidence interval construction.
}
  \item{perturb.delta}{
%%     ~~Describe \code{perturb.delta} here~~
numeric vector; the perturbed values for \eqn{\hat{\Delta}}, the treatment effect estimate, used in variance estimation and confidence interval construction.
}
  \item{delta.s}{
%%     ~~Describe \code{delta.s} here~~
the residual treatment effect, \eqn{\Delta_S}, estimate, \eqn{\hat{\Delta}_S}.
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
the treatment effect, \eqn{\Delta}, estimate, \eqn{\hat{\Delta}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
See documention for R.s.estimate for more detail.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a vector of length 2, lower bound of the 95\% confidence interval and upper bound of the 95\% confidence interval.
}

\references{
%% ~put references to the literature/web site here ~
Fieller, Edgar C. (1954). Some problems in interval estimation. Journal of the Royal Statistical Society. Series B (Methodological), 175-185.

Fieller, E. C. (1940). The biological standardization of insulin. Supplement to the Journal of the Royal Statistical Society, 1-64.

Parast, L., McDermott, M., Tian, L. (2016). Robust estimation of the proportion of treatment effect explained by surrogate marker information. Statistics in Medicine, 35(10):1637-1653.
}

\author{
%%  ~~who you are~~
Layla Parast
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
