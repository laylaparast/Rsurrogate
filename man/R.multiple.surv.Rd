\name{R.multiple.surv}
\alias{R.multiple.surv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates the proportion of treatment effect explained by multiple surrogate markers measured at a specified time and primary outcome information up to that specified time
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function calculates the proportion of treatment effect on the primary outcome explained by multiple surrogate markers measured at \eqn{t_0} and primary outcome information up to \eqn{t_0}.  The user can also request a variance estimate, estimated using perturbating-resampling, and a 95\% confidence interval. If a confidence interval is requested three versions are provided: a normal approximation based interval, a quantile based interval and Fieller's confidence interval, all using perturbation-resampling. The user can also request an estimate of the incremental value of the surrogate marker information. 
}
\usage{
R.multiple.surv(xone, xzero, deltaone, deltazero, sone, szero, type =1, t, 
weight.perturb = NULL, landmark, extrapolate = FALSE, transform = FALSE, 
conf.int = FALSE, var = FALSE, incremental.value = FALSE, approx = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xone}{
%%     ~~Describe \code{xone} here~~
numeric vector, the observed event times in the treatment group, X = min(T,C) where T is the time of the primary outcome and C is the censoring time.
}
  \item{xzero}{
%%     ~~Describe \code{xzero} here~~
numeric vector, the observed event times in the control group, X = min(T,C) where T is the time of the primary outcome and C is the censoring time.
}
  \item{deltaone}{
%%     ~~Describe \code{deltaone} here~~
numeric vector, the event indicators for the treatment group, D = I(T<C) where T is the time of the primary outcome and C is the censoring time.

}
  \item{deltazero}{
%%     ~~Describe \code{deltazero} here~~
numeric vector, the event indicators for the control group, D = I(T<C) where T is the time of the primary outcome and C is the censoring time.
}
  \item{sone}{
%%     ~~Describe \code{sone} here~~
matrix of numeric values; surrogate marker measurements at \eqn{t_0} for treated observations. If \eqn{X_{1i}<t_0}, then the surrogate marker measurements should be NA. 

}
  \item{szero}{
%%     ~~Describe \code{szero} here~~
matrix of numeric values; surrogate marker measurements at \eqn{t_0} for control observations. If \eqn{X_{0i}<t_0}, then the surrogate marker measurements should be NA.
}
  \item{type}{
%%     ~~Describe \code{t} here~~
type of estimate; options are 1 = two-stage robust estimator, 2 = weighted two-stage robust estimator, 3 = double-robust estimator, 4 = two-stage model-based estimator, 5 = weighted estimator, 6 = double-robust model-bsed estimator; default is 1.
}  
\item{t}{
%%     ~~Describe \code{t} here~~
the time of interest.
}
  \item{weight.perturb}{
%%     ~~Describe \code{weight.perturb} here~~
weights used for perturbation resampling.
}
  \item{landmark}{
%%     ~~Describe \code{landmark} here~~
the landmark time \eqn{t_0} or time of surrogate marker measurement.
}
  \item{extrapolate}{
%%     ~~Describe \code{extrapolate} here~~
TRUE or FALSE; indicates whether the user wants to use extrapolation.
}
  \item{transform}{
%%     ~~Describe \code{transform} here~~
TRUE or FALSE; indicates whether the user wants to use a transformation for the surrogate marker pseudo-score. 
}
  \item{conf.int}{
%%     ~~Describe \code{var} here~~
TRUE or FALSE; indicates whether a 95\% confidence interval for delta is requested, default is FALSE.
}
  \item{var}{
%%     ~~Describe \code{var} here~~
TRUE or FALSE; indicates whether a variance estimate is requested, default is FALSE.
}
  \item{incremental.value}{
%%     ~~Describe \code{var} here~~
TRUE or FALSE; indicates whether the user would like to see the incremental value of the surrogate marker information, default is FALSE.
}
  \item{approx}{
TRUE or FALSE indicating whether an approximation should be used when calculating the probability of censoring; most relevant in settings where the survival time of interest for the primary outcome is greater than the last observed event but before the last censored case, default is TRUE.
}
}
\details{
Details can be found in Parast, L., Cai, T., Tian, L. (2020+). Evaluating Multiple Surrogate Markers with Censored Data. Under Review. 

Please email parast@rand.org if you would like a copy of this article. 

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list is returned:
 \item{delta }{the estimate, \eqn{\hat{\Delta}(t)}, described in delta.estimate documentation.}
 \item{delta.s }{the residual treatment effect estimate, \eqn{\hat{\Delta}_S(t,t_0)}.}
 \item{R.s }{the estimated proportion of treatment effect explained by the set of markers, \eqn{\hat{R}_S(t,t_0)}.}
 \item{delta.var }{the variance estimate of \eqn{\hat{\Delta}(t)}; if var = TRUE or conf.int = TRUE.}
 \item{delta.s.var }{the variance estimate of \eqn{\hat{\Delta}_S(t,t_0)}; if var = TRUE or conf.int = TRUE.}
 \item{R.s.var }{the variance estimate of \eqn{\hat{R}_S(t,t_0)}; if var = TRUE or conf.int = TRUE.}
 \item{conf.int.normal.delta}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}(t)} based on a normal approximation; if conf.int = TRUE.}
 \item{conf.int.quantile.delta}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}(t)} based on sample quantiles of the perturbed values; if conf.int = TRUE.} 
 \item{conf.int.normal.delta.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}_S(t,t_0)} based on a normal approximation; if conf.int = TRUE.}
 \item{conf.int.quantile.delta.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}_S(t,t_0)} based on sample quantiles of the perturbed values; if conf.int = TRUE.} 
 \item{conf.int.normal.R.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_S(t,t_0)} based on a normal approximation; if conf.int = TRUE.}
 \item{conf.int.quantile.R.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_S(t,t_0)} based on sample quantiles of the perturbed values; if conf.int = TRUE.} 
\item{conf.int.fieller.R.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_S(t,t_0)} based on Fieller's approach; if conf.int = TRUE.} 
  \item{delta.t }{the estimate, \eqn{\hat{\Delta}_T(t,t_0)}; if incremental.vaue = TRUE.}
 \item{R.t }{the estimated proportion of treatment effect explained by survival only, \eqn{\hat{R}_T(t,t_0)}; if incremental.vaue = TRUE.}
  \item{incremental.value}{the estimate of the incremental value of the surrogate markers, \eqn{\hat{IV}_S(t,t_0)}; if incremental.vaue = TRUE.}
 \item{delta.t.var }{the variance estimate of \eqn{\hat{\Delta}_T(t,t_0)}; if var = TRUE or conf.int = TRUE and  incremental.vaue = TRUE.}
 \item{R.t.var }{the variance estimate of \eqn{\hat{R}_T(t,t_0)}; if var = TRUE or conf.int = TRUE and incremental.vaue = TRUE.} 
  \item{incremental.value.var }{the variance estimate of \eqn{\hat{IV}_S(t,t_0)}; if var = TRUE or conf.int = TRUE and incremental.vaue = TRUE.}
 \item{conf.int.normal.delta.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}_T(t,t_0)} based on a normal approximation; if conf.int = TRUE and incremental.vaue = TRUE.}
 \item{conf.int.quantile.delta.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}_T(t,t_0)} based on sample quantiles of the perturbed values; if conf.int = TRUE and incremental.vaue = TRUE.} 
 \item{conf.int.normal.R.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_T(t,t_0)} based on a normal approximation; if conf.int = TRUE and incremental.vaue = TRUE.}
 \item{conf.int.quantile.R.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_T(t,t_0)} based on sample quantiles of the perturbed values; if conf.int = TRUE and incremental.vaue = TRUE.} 
\item{conf.int.fieller.R.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_T(t,t_0)} based on Fieller's approach; if conf.int = TRUE and incremental.vaue = TRUE.} 
\item{conf.int.normal.iv}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{IV}_S(t,t_0)} based on a normal approximation; if conf.int = TRUE and incremental.vaue = TRUE.}
 \item{conf.int.quantile.iv}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{IV}_S(t,t_0)} based on sample quantiles of the perturbed values; if conf.int = TRUE and incremental.vaue = TRUE.} 
}
\references{
%% ~put references to the literature/web site here ~
Parast, L., Cai, T., & Tian, L. (2017). Evaluating surrogate marker information using censored data. Statistics in Medicine, 36(11), 1767-1782.

Details can be found in Parast, L., Cai, T., Tian, L. (2020+). Evaluating Multiple Surrogate Markers with Censored Data. Under Review.
}
\author{
%%  ~~who you are~~
Layla Parast
}
\note{If the treatment effect is not significant, the user will receive the following message: "Warning: it looks like the treatment effect is not significant; may be difficult to interpret the residual treatment effect in this setting". If the treatment effect is negative, the user will receive the following message: "Warning: it looks like you need to switch the treatment groups" as this package assumes throughout that larger values of the event time are better.}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
data(d_example_multiple)
names(d_example_multiple)
\dontrun{
R.multiple.surv(xone = d_example_multiple$x1, xzero = d_example_multiple$x0, deltaone = 
d_example_multiple$delta1, deltazero = d_example_multiple$delta0, sone = 
as.matrix(d_example_multiple$s1), szero = as.matrix(d_example_multiple$s0), 
type =1, t = 1, landmark=0.5)
R.multiple.surv(xone = d_example_multiple$x1, xzero = d_example_multiple$x0, deltaone = 
d_example_multiple$delta1, deltazero = d_example_multiple$delta0, sone = 
as.matrix(d_example_multiple$s1), szero = as.matrix(d_example_multiple$s0), 
type =1, t = 1, landmark=0.5, conf.int = T)	
R.multiple.surv(xone = d_example_multiple$x1, xzero = d_example_multiple$x0, deltaone = 
d_example_multiple$delta1, deltazero = d_example_multiple$delta0, sone = 
as.matrix(d_example_multiple$s1), szero = as.matrix(d_example_multiple$s0), 
type =3, t = 1, landmark=0.5)
}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
\keyword{survival}
\keyword{robust}
\keyword{smooth}