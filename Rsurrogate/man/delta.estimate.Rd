\name{delta.estimate}
\alias{delta.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates treatment effect
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function calculates the treatment effect estimate, the difference in the average outcome in the treatment group minus the control group. This function is intended to be used for a fully observed continuous outcome. The user can also request a variance estimate, estimated using perturbating-resampling, and a 95\% confidence interval. If a confidence interval is requested two versions are provided: a normal approximation based interval and a quantile based interval, both use perturbation-resampling. 
}
\usage{
delta.estimate(yone,yzero, var = FALSE, conf.int = FALSE, weight = NULL, 
weight.perturb = NULL)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yone}{
%%     ~~Describe \code{yone} here~~
numeric vector; primary outcome for treated observations.
}
  \item{yzero}{
%%     ~~Describe \code{yzero} here~~
numeric vector; primary outcome for control observations.
}
  \item{var}{
%%     ~~Describe \code{var} here~~
TRUE or FALSE; indicates whether a variance estimate for delta is requested, default is FALSE.
}
  \item{conf.int}{
%%     ~~Describe \code{var} here~~
TRUE or FALSE; indicates whether a 95\% confidence interval for delta is requested, default is FALSE.
}
  \item{weight}{
%%     ~~Describe \code{weight} here~~
a n1+n0 by x matrix of weights where n1 = length of yone and n0 = length of yzero, default is null; generally not supplied by use but only used by other functions.
}
  \item{weight.perturb}{
%%     ~~Describe \code{weight} here~~
a n1+n0 by x matrix of weights where n1 = length of yone and n0 = length of yzero, default is null; generally used for confidence interval construction and may be supplied by user.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Let \eqn{Y^{(1)}} and \eqn{Y^{(0)}} denote the primary outcome under the treatment and primary outcome under the control,respectively. The treatment effect, \eqn{\Delta}, is the expected difference in \eqn{Y^{(1)}} compared to \eqn{Y^{(0)}}, \eqn{ \Delta=E(Y^{(1)}-Y^{(0)}).} We estimate \eqn{\Delta} as \deqn{\hat{\Delta}
= n_1^{-1} \sum_{i=1}^{n_1} Y_{1i}  - n_0^{-1} \sum_{i=1}^{n_0}  Y_{0i}} where \eqn{Y_{1i}} is the observed primary outcome for person \eqn{i} in the treated group, \eqn{Y_{0i}} is the observed primary outcome for person \eqn{i} in the control group, and \eqn{n_1} and \eqn{n_0} are the number of individuals in the treatment and control group, respectively. Randomized treatment assignment is assumed throughout this package.

Variance estimation and confidence interval construction are performed using perturbation-resampling. Specifically, let \eqn{\left \{ V^{(b)} = (V_{11}^{(b)}, ...V_{1n_1}^{(b)}, V_{01}^{(b)}, ...V_{0n_0}^{(b)})^T, b=1,....,D \right \}} be \eqn{n \times D} independent copies of a positive random variables \eqn{V} from a known distribution with unit mean and unit variance. Let
\deqn{\hat{\Delta}^{(b)}  = \frac{  \sum_{i=1}^{n_1} V_{1i}^{(b)} Y_{1i}}{ \sum_{i=1}^{n_1} V_{1i}^{(b)}}  - \frac{  \sum_{i=1}^{n_0} V_{0i}^{(b)} Y_{0i}}{ \sum_{i=1}^{n_0} V_{0i}^{(b)}}.}The variance of \eqn{\hat{\Delta}} is obtained as the empirical variance of \eqn{\{\hat{\Delta}^{(b)}, b = 1,...,D\}.} In this package, we use weights generated from an Exponential(1) distribution and use \eqn{D=500}.  We construct two versions of the \eqn{95\%} confidence interval for \eqn{\hat{\Delta}}: one based on a normal approximation confidence interval using the estimated variance and another taking the 2.5th and 97.5th empirical percentiles of \eqn{\hat{\Delta}^{(b)}}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list is returned:
 \item{delta }{the estimate, \eqn{\hat{\Delta}}, described above. }
 \item{var }{the variance estimate of \eqn{\hat{\Delta}}; if var = TRUE or conf.int = TRUE.}
 \item{conf.int.normal}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}} based on a normal approximation; if conf.int = TRUE.}
 \item{conf.int.quantile}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}} based on sample quantiles of the perturbed values, described above; if conf.int = TRUE.} 
}
\author{
%%  ~~who you are~~
Layla Parast
}
\examples{
data(d_example)
names(d_example)
delta.estimate(yone=d_example$y1, yzero=d_example$y0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}
