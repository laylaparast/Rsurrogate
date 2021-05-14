\name{delta.surv.estimate}
\alias{delta.surv.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates treatment effect in a survival setting
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function calculates the treatment effect in the survival setting i.e. the difference in survival at time t between the treatment group and the control group. The user can also request a variance estimate, estimated using perturbating-resampling, and a 95\% confidence interval. If a confidence interval is requested two versions are provided: a normal approximation based interval and a quantile based interval, both use perturbation-resampling. 

}
\usage{
delta.surv.estimate(xone, xzero, deltaone, deltazero, t, var = FALSE, conf.int 
= FALSE, weight = NULL, weight.perturb = NULL, approx = T)
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
  \item{var}{
%%     ~~Describe \code{var} here~~
TRUE or FALSE; indicates whether a variance estimate for delta is requested, default is FALSE.
}
  \item{conf.int}{
%%     ~~Describe \code{conf.int} here~~
TRUE or FALSE; indicates whether a 95\% confidence interval for delta is requested, default is FALSE.

}
  \item{weight}{
%%     ~~Describe \code{weight} here~~
a \eqn{n_1+n_0} by \eqn{x} matrix of weights where \eqn{n_1 =} sample size in treatment group and \eqn{n_0 =} sample size in the control group, default is null; generally not supplied by use but only used by other functions.

}
  \item{weight.perturb}{
%%     ~~Describe \code{weight} here~~
a \eqn{n_1+n_0} by \eqn{x} matrix of weights where \eqn{n_1 =}sample size in treatment group and \eqn{n_0 =} sample size in the control group, default is null; generally used for confidence interval construction and may be supplied by user.

}
  \item{approx}{
TRUE or FALSE indicating whether an approximation should be used when calculating the probability of censoring; most relevant in settings where the survival time of interest for the primary outcome is greater than the last observed event but before the last censored case, default is TRUE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Let \eqn{G} be the binary treatment indicator with \eqn{G=1} for treatment and \eqn{G=0} for control and we assume throughout that subjects are randomly assigned to a treatment group at baseline. Let \eqn{T} denote the time of the primary outcome of interest, death for example. We use potential outcomes notation such that \eqn{T^{(g)}} denotes the time of the primary outcome under treatment \eqn{G = g}. We define the treatment effect, \eqn{\Delta(t)}, as the difference in survival rates by time \eqn{t} under treatment versus control,
  \deqn{ \Delta(t)=E\{ I(T^{(1)}>t)\} - E\{I(T^{(0)}>t)\} = P(T^{(1)}>t) - P(T^{(0)}>t)}
where \eqn{t>t_0}

Due to censoring, our data consist of \eqn{n_1} observations \eqn{\{(X_{1i}, \delta_{1i}), i=1,...,n_1\}} from the treatment group \eqn{G=1} and \eqn{n_0} observations \eqn{\{(X_{0i}, \delta_{0i}), i=1,...,n_0\}} from the control group \eqn{G=0} where \eqn{X_{gi} = \min(T_{gi}, C_{ gi})},  \eqn{\delta_{gi} = I(T_{gi} < C_{gi})}, and \eqn{C_{gi}} denotes the censoring time for \eqn{g= 1,0}, for individual \eqn{i}. Throughout, we estimate the treatment effect \eqn{\Delta(t)} as
\deqn{ \hat{\Delta}(t)  = n_1^{-1} \sum_{i=1}^{n_1} \frac{I(X_{1i}>t)}{\hat{W}^C_1(t)} - n_0^{-1} \sum_{i=1}^{n_0} \frac{I(X_{0i}>t)}{\hat{W}^C_0(t)}}
where  \eqn{\hat{W}^C_g(\cdot)} is the Kaplan-Meier estimator of survival for censoring for \eqn{g=1,0.} 


Variance estimation and confidence interval construction are performed using perturbation-resampling. Specifically, let \eqn{\left \{ V^{(b)} = (V_{11}^{(b)}, ...V_{1n_1}^{(b)}, V_{01}^{(b)}, ...V_{0n_0}^{(b)})^T, b=1,....,D \right \}} be \eqn{n \times D} independent copies of a positive random variables \eqn{V} from a known distribution with unit mean and unit variance. Let
\deqn{\hat{\Delta}^{(b)} (t) = \frac{  \sum_{i=1}^{n_1} V_{1i}^{(b)} I(X_{1i}>t)}{ \sum_{i=1}^{n_1} V_{1i}^{(b)} \hat{W}_1^{C(b)}(t)}  -\frac{  \sum_{i=1}^{n_0} V_{0i}^{(b)} I(X_{0i}>t)}{ \sum_{i=1}^{n_0} V_{0i}^{(b)} \hat{W}_0^{C(b)}(t)}.} In this package, we use weights generated from an Exponential(1) distribution and use \eqn{D=500}. The variance of \eqn{\hat{\Delta}(t)} is obtained as the empirical variance of \eqn{\{\hat{\Delta}(t)^{(b)}, b = 1,...,D\}.} We construct two versions of the \eqn{95\%}confidence interval for \eqn{\hat{\Delta}(t)}: one based on a normal approximation confidence interval using the estimated variance and another taking the 2.5th and 97.5th empirical percentiles of \eqn{\hat{\Delta}(t)^{(b)}}.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list is returned:
 \item{delta }{the estimate, \eqn{\hat{\Delta}(t)}, described above. }
 \item{var }{the variance estimate of \eqn{\hat{\Delta}(t)}; if var = TRUE or conf.int = TRUE.}
 \item{conf.int.normal}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}(t)} based on a normal approximation; if conf.int = TRUE.}
 \item{conf.int.quantile}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}(t)} based on sample quantiles of the perturbed values, described above; if conf.int = TRUE.} 
}
\author{
%%  ~~who you are~~
Layla Parast
}
\examples{
data(d_example_surv)
names(d_example_surv)
delta.surv.estimate(xone = d_example_surv$x1, xzero = d_example_surv$x0,  
deltaone = d_example_surv$delta1, deltazero = d_example_surv$delta0, t = 3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}
\keyword{nonparametric}
\keyword{survival}
