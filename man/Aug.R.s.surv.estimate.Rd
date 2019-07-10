\name{Aug.R.s.surv.estimate}
\alias{Aug.R.s.surv.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the augmented estimator of the proportion of treatment effect explained by the surrogate marker information measured at a specified time and primary outcome information up to that specified time

}
\description{
This function calculates the augmented version of the proportion of treatment effect on the primary outcome explained by the surrogate marker information measured at \eqn{t_0} and primary outcome information up to \eqn{t_0}.  Variance estimates and 95 \% confidence intervals for the augmented estimates are provided automatically;  three versions of the confidence interval are provided: a normal approximation based interval, a quantile based interval and Fieller's confidence interval, all using perturbation-resampling. The user can also request an estimate of the incremental value of surrogate marker information. 
}
\usage{
Aug.R.s.surv.estimate(xone, xzero, deltaone, deltazero, sone, szero, t, 
weight.perturb = NULL, landmark, extrapolate = FALSE, transform = FALSE, 
basis.delta.one, basis.delta.zero, basis.delta.s.one = NULL, 
basis.delta.s.zero = NULL, incremental.value = FALSE)
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
numeric vector; surrogate marker measurement at \eqn{t_0} for treated observations, assumed to be continuous. If \eqn{X_{1i}<t_0}, then the surrogate marker measurement should be NA. 

}
  \item{szero}{
%%     ~~Describe \code{szero} here~~
numeric vector; surrogate marker measurement at \eqn{t_0} for control observations, assumed to be continuous. If \eqn{X_{1i}<t_0}, then the surrogate marker measurement should be NA.
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
TRUE or FALSE; indicates whether the user wants to use a transformation for the surrogate marker. 
}
  \item{basis.delta.one}{
either a vector of length \eqn{n_1} or a matrix with \eqn{n_1} rows; this is the basis transformation used for augmentation of \eqn{\hat{\Delta}(t)} for treated observations only, all values must be numeric
}
  \item{basis.delta.zero}{
either a vector of length \eqn{n_0} or a matrix with \eqn{n_0} rows; this is the basis transformation used for augmentation of \eqn{\hat{\Delta}(t)} for control observations only, all values must be numeric
}
  \item{basis.delta.s.one}{
either a vector of length \eqn{n_1} or a matrix with \eqn{n_1} rows; this is the basis transformation used for augmentation of \eqn{\hat{\Delta}_S(t,t_0)} for treated observations only, all values must be numeric; default is to assume this is the same as basis.delta.one i.e. that the same basis transformation is used for both quantities
}
  \item{basis.delta.s.zero}{
either a vector of length \eqn{n_0} or a matrix with \eqn{n_0} rows; this is the basis transformation used for augmentation of \eqn{\hat{\Delta}_S(t,t_0)} for control observations only, all values must be numeric; default is to assume this is the same as basis.delta.zero i.e. that the same basis transformation is used for both quantities
}
  \item{incremental.value}{
%%     ~~Describe \code{var} here~~
TRUE or FALSE; indicates whether the user would like to see the incremental value of the surrogate marker information, default is FALSE.
}
}
\details{
Please see R.s.surv.estimate documention for details about the estimates before augmentation is performed. Recent work has shown that augmentation can lead to improvements in efficiency by taking advantage of the association between baseline information, denoted here as \eqn{Z}, and the primary outcome. This function calculates the augmented estimates of the quantities of interest. For example, the augmented version of \eqn{ \hat{\Delta}(t)} is defined as:
\deqn{ \hat{\Delta}(t)^{AUG} =   \hat{\Delta}(t) + \gamma \{n_1^{-1}\sum_{i=1}^{n_1}h(Z_{1i})-n_0^{-1}\sum_{i=1}^{n_0}h(Z_{0i}) \} } 
where \eqn{Z_{gi}, i=1, 2, \cdots, n_g} are i.i.d. random vectors of baseline covariates from treatment group \eqn{g} and \eqn{h(\cdot)} is a basis transformation given a priori.  Due to treatment randomization, \eqn{ \{n_1^{-1}\sum_{i=1}^{n_1}h(Z_{1i})-n_0^{-1}\sum_{i=1}^{n_0}h(Z_{0i}) \}} converges to zero in probability as the sample size goes to infinity and thus the augmented estimator converges to the same limit as the original counterparts.  The quantity \eqn{\gamma} is selected such that the variance of \eqn{ \hat{\Delta}(t)^{AUG}} is minimized.  That is,  \eqn{\gamma = (\Xi_{12}) ( \Xi_{22} ) ^{-1}} where
\deqn{ \Xi_{12} = \mbox{cov}  \{ \hat{\Delta}(t), n_1^{-1}\sum_{i=1}^{n_1}h(Z_{1i})-n_0^{-1}\sum_{i=1}^{n_0}h(Z_{0i}) \},}
\deqn{  \Xi_{22} = \mbox{var}  \{n_1^{-1}\sum_{i=1}^{n_1}h(Z_{1i})-n_0^{-1}\sum_{i=1}^{n_0}h(Z_{0i})\}}
and thus we can obtain \eqn{\hat{\Delta}(t)^{AUG}} by replacing \eqn{\gamma} with a consistent estimator, \eqn{\hat{\gamma}} obtained using perturbation-resampling. A similar approach is used to obtain  \eqn{\hat{\Delta}_S(t)^{AUG}} and thus construct
\deqn{\hat{R}_S(t,t_0)^{AUG}=1-\frac{\hat{\Delta}_S(t,t_0)^{AUG}}{\hat{\Delta}(t)^{AUG}}.}

}
\value{
A list is returned:
 \item{aug.delta }{the estimate, \eqn{\hat{\Delta}(t)^{AUG}}.}
 \item{aug.delta.s }{the estimate, \eqn{\hat{\Delta}_S(t,t_0)^{AUG}}.}
 \item{aug.R.s }{the estimate, \eqn{\hat{R}_S(t,t_0)^{AUG}}.}
 \item{aug.delta.var }{the variance estimate of \eqn{\hat{\Delta}(t)^{AUG}}.}
 \item{aug.delta.s.var }{the variance estimate of \eqn{\hat{\Delta}_S(t,t_0)^{AUG}}.}
 \item{aug.R.s.var }{the variance estimate of \eqn{\hat{R}_S(t,t_0)^{AUG}}.}
 \item{conf.int.normal.aug.delta}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}(t)^{AUG}} based on a normal approximation.}
 \item{conf.int.quantile.aug.delta}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}(t)^{AUG}} based on sample quantiles of the perturbed values.} 
 \item{conf.int.normal.aug.delta.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}_S(t,t_0)^{AUG}} based on a normal approximation.}
 \item{conf.int.quantile.aug.delta.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}_S(t,t_0)^{AUG}} based on sample quantiles of the perturbed values.} 
 \item{conf.int.normal.R.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_S(t,t_0)^{AUG}} based on a normal approximation.}
 \item{conf.int.quantile.aug.R.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_S(t,t_0)^{AUG}} based on sample quantiles of the perturbed values..} 
\item{conf.int.fieller.aug.R.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_S(t,t_0)^{AUG}} based on Fieller's approach.} 
  \item{aug.delta.t }{the estimate, \eqn{\hat{\Delta}_T(t,t_0)^{AUG}}; if incremental.vaue = TRUE.}
 \item{aug.R.t }{the estimate, \eqn{\hat{R}_T(t,t_0)^{AUG}}; if incremental.vaue = TRUE.}
  \item{aug.incremental.value}{the estimate, \eqn{\hat{IV}_S(t,t_0)^{AUG}}; if incremental.vaue = TRUE.}
 \item{aug.delta.t.var }{the variance estimate of \eqn{\hat{\Delta}_T(t,t_0)^{AUG}}; if  incremental.vaue = TRUE.}
 \item{aug.R.t.var }{the variance estimate of \eqn{\hat{R}_T(t,t_0)^{AUG}}; if  incremental.vaue = TRUE.} 
  \item{aug.incremental.value.var }{the variance estimate of \eqn{\hat{IV}_S(t,t_0)^{AUG}}; if incremental.vaue = TRUE.}
 \item{aug.conf.int.normal.delta.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}_T(t,t_0)^{AUG}} based on a normal approximation; if incremental.vaue = TRUE.}
 \item{aug.conf.int.quantile.delta.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}_T(t,t_0)^{AUG}} based on sample quantiles of the perturbed values; if  incremental.vaue = TRUE.} 
 \item{aug.conf.int.normal.R.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_T(t,t_0)^{AUG}} based on a normal approximation; if  incremental.vaue = TRUE.}
 \item{aug.conf.int.quantile.R.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_T(t,t_0)^{AUG}} based on sample quantiles of the perturbed values; if  incremental.vaue = TRUE.} 
\item{aug.conf.int.fieller.R.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_T(t,t_0)^{AUG}} based on Fieller's approach, described above; if incremental.vaue = TRUE.} 
\item{aug.conf.int.normal.iv}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{IV}_S(t,t_0)^{AUG}} based on a normal approximation; if incremental.vaue = TRUE.}
 \item{aug.conf.int.quantile.iv}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{IV}_S(t,t_0)^{AUG}} based on sample quantiles of the perturbed values; if  incremental.vaue = TRUE.} 
}
\references{
Tian L, Cai T, Zhao L,Wei L. On the covariate-adjusted estimation for an overall treatment difference with data from a randomized comparative clinical trial. Biostatistics 2012; 13(2): 256-273.

Garcia TP, Ma Y, Yin G. Efficiency improvement in a class of survival models through model-free covariate incorporation. Lifetime Data Analysis 2011; 17(4): 552-565. 

Zhang M, Tsiatis AA, Davidian M. Improving efficiency of inferences in randomized clinical trials using auxiliary covariates. Biometrics 2008; 64(3): 707-715.

Parast L, Cai T and Tian L. Evaluating Surrogate Marker Information using Censored Data. Under Review.

}
\author{
Layla Parast
}
\note{If the treatment effect is not significant, the user will receive the following message: "Warning: it looks like the treatment effect is not significant; may be difficult to interpret the residual treatment effect in this setting". If the treatment effect is negative, the user will receive the following message: "Warning: it looks like you need to switch the treatment groups" as this package assumes throughout that larger values of the event time are better. If the observed support of the surrogate marker for the control group is outside the observed support of the surrogate marker for the treatment group, the user will receive the following message: "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation".
}
\examples{
#computationally intensive
#Aug.R.s.surv.estimate(xone = d_example_surv$x1, xzero = d_example_surv$x0,  
#deltaone = d_example_surv$delta1, deltazero = d_example_surv$delta0, 
#sone = d_example_surv$s1, szero = d_example_surv$s0, t=3, landmark = 1, 
#basis.delta.one = d_example_surv$z1 , basis.delta.zero = d_example_surv$z0)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
\keyword{survival}
\keyword{robust}
\keyword{smooth}
\keyword{augment}