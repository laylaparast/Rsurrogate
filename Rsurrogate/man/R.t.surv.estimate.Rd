\name{R.t.surv.estimate}
\alias{R.t.surv.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates the proportion of treatment effect explained by the primary outcome information up to a specified time

}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function calculates the proportion of treatment effect on the primary outcome explained by the treatment effect on the primary outcome up to \eqn{t_0}.  The user can also request a variance estimate, estimated using perturbating-resampling, and a 95\% confidence interval. If a confidence interval is requested three versions are provided: a normal approximation based interval, a quantile based interval and Fieller's confidence interval, all using perturbation-resampling.
}
\usage{
R.t.surv.estimate(xone, xzero, deltaone, deltazero, t, weight.perturb = NULL, 
landmark, var = FALSE, conf.int = FALSE, approx = T)
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
  \item{var}{
%%     ~~Describe \code{var} here~~
TRUE or FALSE; indicates whether a variance estimate for delta is requested, default is FALSE.
}
  \item{conf.int}{
%%     ~~Describe \code{var} here~~
TRUE or FALSE; indicates whether a 95\% confidence interval for delta is requested, default is FALSE.
}
  \item{approx}{
TRUE or FALSE indicating whether an approximation should be used when calculating the probability of censoring; most relevant in settings where the survival time of interest for the primary outcome is greater than the last observed event but before the last censored case, default is TRUE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Let \eqn{G} be the binary treatment indicator with \eqn{G=1} for treatment and \eqn{G=0} for control and we assume throughout that subjects are randomly assigned to a treatment group at baseline. Let \eqn{T} denote the time of the primary outcome of interest, death for example. We use potential outcomes notation such that \eqn{T^{(g)}} denotes the time of the primary outcome under treatment \eqn{G = g}. The proportion of treatment effect explained by T observed up to \eqn{t_0} only is \eqn{R_T(t,t_0) = 1-\Delta_T(t,t_0)/\Delta(t)} where
\deqn{\Delta_T(t, t_0) = P(T^{(0)}>t_0)P(T^{(1)}>t\mid T^{(1)}>t_0)-P(T^{(0)}>t).}  To estimate \eqn{R_T(t,t_0)}, we use the estimator \eqn{\hat{R}_T(t,t_0) = 1-\hat{\Delta}_T(t,t_0)/\hat{\Delta}(t)} where \eqn{\hat{\Delta}_T(t,t_0) =  \hat{\phi}_0(t_0)\hat{\phi}_1(t)/\hat{\phi}_1(t_0) - \hat{\phi}_0(t)} and \eqn{\hat{\phi}_g(u) = n_g^{-1} \sum_{i=1}^{n_g} \frac{I(X_{gi}>u)}{\hat{W}^C_g(u)}} for \eqn{g=1,0} where  \eqn{\widehat{W}^C_g(\cdot)} is the Kaplan-Meier estimator of survival for censoring for \eqn{g=1,0.} 

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list is returned:
 \item{delta }{the estimate, \eqn{\hat{\Delta}(t)}, described in delta.estimate documentation.}
 \item{delta.t }{the estimate, \eqn{\hat{\Delta}_T(t,t_0)}, described above.}
 \item{R.t }{the estimate, \eqn{\hat{R}_T(t,t_0)}, described above.}
 \item{delta.var }{the variance estimate of \eqn{\hat{\Delta}(t)}; if var = TRUE or conf.int = TRUE.}
 \item{delta.t.var }{the variance estimate of \eqn{\hat{\Delta}_T(t,t_0)}; if var = TRUE or conf.int = TRUE.}
 \item{R.t.var }{the variance estimate of \eqn{\hat{R}_T(t,t_0)}; if var = TRUE or conf.int = TRUE.}
 \item{conf.int.normal.delta}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}(t)} based on a normal approximation; if conf.int = TRUE.}
 \item{conf.int.quantile.delta}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}(t)} based on sample quantiles of the perturbed values, described above; if conf.int = TRUE.} 
 \item{conf.int.normal.delta.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}_T(t,t_0)} based on a normal approximation; if conf.int = TRUE.}
 \item{conf.int.quantile.delta.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}_T(t,t_0)} based on sample quantiles of the perturbed values, described above; if conf.int = TRUE.} 
 \item{conf.int.normal.R.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_T(t,t_0)} based on a normal approximation; if conf.int = TRUE.}
 \item{conf.int.quantile.R.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_T(t,t_0)} based on sample quantiles of the perturbed values, described above; if conf.int = TRUE.} 
\item{conf.int.fieller.R.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_T(t,t_0)} based on Fieller's approach, described above; if conf.int = TRUE.} 

}
\references{
%% ~put references to the literature/web site here ~
Parast, L., Cai, T., & Tian, L. (2017). Evaluating surrogate marker information using censored data. Statistics in Medicine, 36(11), 1767-1782.
}

\note{If the treatment effect is not significant, the user will receive the following message: "Warning: it looks like the treatment effect is not significant; may be difficult to interpret the residual treatment effect in this setting".}


\author{
%%  ~~who you are~~
Layla Parast
}
\examples{
data(d_example_surv)
names(d_example_surv)
% Functions have elapsed time longer than 10s
%R.t.surv.estimate(xone = d_example_surv$x1, xzero = d_example_surv$x0,  
%deltaone = d_example_surv$delta1, deltazero = d_example_surv$delta0, 
%t = 3, landmark = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
\keyword{survival}
\keyword{robust}