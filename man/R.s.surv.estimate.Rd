\name{R.s.surv.estimate}
\alias{R.s.surv.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates the proportion of treatment effect explained by the surrogate marker information measured at a specified time and primary outcome information up to that specified time
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function calculates the proportion of treatment effect on the primary outcome explained by the surrogate marker information measured at \eqn{t_0} and primary outcome information up to \eqn{t_0}.  The user can also request a variance estimate, estimated using perturbating-resampling, and a 95\% confidence interval. If a confidence interval is requested three versions are provided: a normal approximation based interval, a quantile based interval and Fieller's confidence interval, all using perturbation-resampling. The user can also request an estimate of the incremental value of surrogate marker information. 
}
\usage{
R.s.surv.estimate(xone, xzero, deltaone, deltazero, sone, szero, t, 
weight.perturb = NULL, landmark, extrapolate = FALSE, transform = FALSE, 
conf.int = FALSE, var = FALSE, incremental.value = FALSE)
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
  \item{conf.int}{
%%     ~~Describe \code{var} here~~
TRUE or FALSE; indicates whether a 95\% confidence interval for delta is requested, default is FALSE.
}
  \item{var}{
%%     ~~Describe \code{var} here~~
TRUE or FALSE; indicates whether a variance estimate for delta is requested, default is FALSE.
}
  \item{incremental.value}{
%%     ~~Describe \code{var} here~~
TRUE or FALSE; indicates whether the user would like to see the incremental value of the surrogate marker information, default is FALSE.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
Let \eqn{G} be the binary treatment indicator with \eqn{G=1} for treatment and \eqn{G=0} for control and we assume throughout that subjects are randomly assigned to a treatment group at baseline. Let \eqn{T^{(1)}} and \eqn{T^{(0)}} denote the time of the primary outcome of interest, death for example, under the treatment and under the control, respectively. Let \eqn{S^{(1)}} and \eqn{S^{(0)}} denote the surrogate marker measured at time \eqn{t_0} under the treatment and the control, respectively.

 The residual treatment effect is defined as
\deqn{ \Delta_S(t,t_0) = P(T^{(0)} > t_0) \left \{\int \psi_1(t | s, t_0) dF_0(s | t_0)-P(T^{(0)}> t| T^{(0)}> t_0) \right \}} where \eqn{F_0(\cdot |  t_0)} is the cumulative distribution function of \eqn{S^{(0)}} conditional on \eqn{T^{(0)}> t_0} and \eqn{\psi_1(t |  s,t_0) = P(T^{(1)}> t |  S^{(1)}=s, T ^{(1)} > t_0)}. The proportion of treatment effect explained by the surrogate marker information measured at \eqn{t_0} and primary outcome information up to \eqn{t_0}, which we denote by \eqn{R_S(t,t_0)}, can be expressed using a contrast between \eqn{\Delta_S(t,t_0)} and \eqn{\Delta(t)}:
\deqn{R_S(t,t_0)=\{\Delta(t)-\Delta_S(t,t_0)\}/\Delta(t)=1-\Delta_S(t,t_0)/\Delta(t).}   The definition and estimation of \eqn{\Delta(t)} is described in the delta.surv.estimate documentation.

Due to censoring, our data consist of \eqn{n_1} observations \eqn{\{(X_{1i}, \delta_{1i}, S_{1i}), i=1,...,n_1\}} from the treatment group \eqn{G=1} and \eqn{n_0} observations \eqn{\{(X_{0i}, \delta_{0i}, S_{0i}), i=1,...,n_0\}} from the control group \eqn{G=0} where \eqn{X_{gi} = \min(T_{gi}, C_{ gi})},  \eqn{\delta_{gi} = I(T_{gi} < C_{gi})}, \eqn{C_{gi}} denotes the censoring time, and \eqn{S_{gi}} denotes the surrogate marker information measured at time \eqn{t_0}, for \eqn{g= 1,0}, for individual \eqn{i}. Note that if \eqn{X_{gi} < t_0}, then \eqn{S_{gi}} should be NA (not available).

To estimate \eqn{\Delta_S(t,t_0)}, we use a nonparametric kernel Nelson-Aalen estimator to estimate \eqn{\psi_1(t | s,t_0)} as \eqn{\hat{\psi}_1(t| s,t_0)  = \exp\{-\hat{\Lambda}_1(t| s,t_0) \}}, where
\deqn{\hat{\Lambda}_1(t| s,t_0) = \int_{t_0}^t \frac{\sum_{i=1}^{n_1} I(X_{1i}>t_0) K_h\{\gamma(S_{1i}) - \gamma(s)\}dN_{1i}(z)}{\sum_{i=1}^{n_1}  I(X_{1i}>t_0) K_h\{\gamma(S_{1i}) - \gamma(s)\} Y_{1i}(z)},}
is a consistent estimate of \eqn{\Lambda_1(t| s,t_0 ) =  -\log [\psi_1(t| s,t_0)], Y_{1i}(t) = I(X_{1i} \geq t),  N_{1i}(t) = I(X_{1i} \leq t) \delta_i,  K(\cdot)} is a smooth symmetric density function, \eqn{K_h(x) = K(x/h)/h, \gamma(\cdot)} is a given monotone transformation function, and \eqn{h} is a specified bandwidth. To obtain an appropriate \eqn{h} we first use bw.nrd to obtain \eqn{h_{opt}}; and then we let \eqn{h = h_{opt}n_1^{-c_0}}  with \eqn{c_0 = 0.11}.

Since \eqn{F_0(s | t_0) = P(S_{0i} \le s | X_{0i} > t_0)}, we empirically estimate \eqn{F_0(s | t_0)} using all subjects with \eqn{X_{0i} > t_0} as \deqn{
\hat{F}_0(s | t_0) = \frac{\sum_{i=1}^{n_0} I(S_{0i} \le s, X_{0i} > t_0)}{\sum_{i=1}^{n_0} I(X_{0i} > t_0)}.}
Subsequently, we construct an estimator for \eqn{\Delta_{S}(t,t_0)} as
\deqn{\hat{\Delta}_S(t,t_0) = n_0^{-1} \sum_{i=1}^{n_0} \left[\hat{\psi}_1(t| S_{0i},t_0) \frac{I(X_{0i} > t_0)}{\hat{W}^C_0(t_0)}  -   \frac{I(X_{0i} > t)}{\hat{W}^C_0(t)}\right]}
 where \eqn{\hat{W}^C_g(\cdot)} is the Kaplan-Meier estimator of survival for censoring for \eqn{g=1,0.} Finally, we estimate \eqn{R_S(t,t_0)} as \eqn{\hat{R}_S(t,t_0) =1- \hat{\Delta}_S(t,t_0)/\hat{\Delta}(t).} 

	Variance estimation and confidence interval construction are performed using perturbation-resampling. Specifically, let \eqn{\left \{ V^{(b)} = (V_{11}^{(b)}, ...V_{1n_1}^{(b)}, V_{01}^{(b)}, ...V_{0n_0}^{(b)})^T, b=1,....,D \right \}} be \eqn{n \times D} independent copies of a positive random variables \eqn{V} from a known distribution with unit mean and unit variance. Let
\deqn{\hat{\Delta}^{(b)}(t) = \frac{  \sum_{i=1}^{n_1} V_{1i}^{(b)} I(X_{1i}>t)}{ \sum_{i=1}^{n_1} V_{1i}^{(b)} \hat{W}_1^{C(b)}(t)}  -\frac{  \sum_{i=1}^{n_0} V_{0i}^{(b)} I(X_{0i}>t)}{ \sum_{i=1}^{n_0} V_{0i}^{(b)} \hat{W}_0^{C(b)}(t)}.} In this package, we use weights generated from an Exponential(1) distribution and use \eqn{D=500}. The variance of \eqn{\hat{\Delta}(t)} is obtained as the empirical variance of \eqn{\{\hat{\Delta}(t)^{(b)}, b = 1,...,D\}}. Variance estimates for \eqn{\hat{\Delta}_S(t,t_0)} and \eqn{\hat{R}_S(t,t_0)} are calculated similarly. We construct two versions of the \eqn{95\%} confidence interval for each estimate: one based on a normal approximation confidence interval using the estimated variance and another taking the 2.5th and 97.5th empirical percentile of the perturbed quantities.  In addition, we use Fieller's method to obtain a third confidence interval for \eqn{R_S(t,t_0)} as
\deqn{\left\{1-r:  \frac{(\hat{\Delta}_S(t,t_0)-r\hat{\Delta}(t))^2}{\hat{\sigma}_{11}-2r\hat\sigma_{12}+r^2\hat\sigma_{22}} \le c_{\alpha}\right\},}
where \eqn{\hat{\Sigma}=(\hat\sigma_{ij})_{1\le i,j\le 2}} and \eqn{c_\alpha} is the \eqn{(1-\alpha)}th percentile of
\deqn{\left\{\frac{\{\hat{\Delta}^{(b)}_S(t)-(1-\hat R_S(t,t_0))\hat{\Delta}(t)^{(b)}\}^2}{\hat{\sigma}_{11}-2(1-\hat R_S(t,t_0))\hat\sigma_{12}+(1-\hat R_S(t,t_0))^2\hat\sigma_{22}}, b=1, \cdots, C\right\}}
where \eqn{\alpha=0.05}.


Since the definition of \eqn{R_S(t,t_0)} considers the surrogate information as a combination of both \eqn{S} information and \eqn{T} information up to \eqn{t_0}, a logical inquiry would be how to assess the incremental value of the \eqn{S} information in terms of the proportion of treatment effect explained, when added to \eqn{T} information up to \eqn{t_0}. The proportion of treatment effect explained by \eqn{T} information up to \eqn{t_0} only is denoted as \eqn{R_T(t,t_0)} and is described in the documentation for R.t.surv.estimate. The incremental value of \eqn{S} information is defined as:
\deqn{
IV_S(t,t_0) = R_S(t,t_0) - R_T(t,t_0) = \frac{\Delta_T(t,t_0) - \Delta_S(t,t_0)}{\Delta (t)}.}

For estimation of \eqn{R_T(t_0)}, see documentation for R.t.surv.estimate. The quantity \eqn{IV_S(t,t_0)}  is then estimated by \eqn{\hat{IV}_S(t,t_0) = \hat{R}_S(t,t_0) - \hat{R}_T(t,t_0)}. Perturbation-resampling is used for variance estimation and confidence interval construction for this quantity, similar to the other quantities in this package.

Note that if the observed supports for S are not the same, then \eqn{\hat{\Lambda}_1(t| s,t_0)} for \eqn{S_{0i} = s} outside the support of \eqn{S_{1i}} may return NA (depending on the bandwidth). If extrapolation = TRUE, then the   \eqn{\hat{\Lambda}_1(t| s,t_0)} values for these surrogate values are set to the closest non-NA value. If transform = TRUE, then \eqn{S_{1i}} and \eqn{S_{0i}} are transformed such that the new transformed values, \eqn{S^{tr}_{1i}} and \eqn{S^{tr}_{0i}} are defined as: \eqn{S^{tr}_{gi} = F([S_{gi} - \mu]/\sigma)} for \eqn{g=0,1} where \eqn{F(\cdot)} is the cumulative distribution function for a standard normal random variable, and \eqn{\mu} and \eqn{\sigma} are the sample mean and standard deviation, respectively, of \eqn{\{(S_{1i}, S_{0i})^T, i \quad s.t. X_{gi} > t_0\}}.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list is returned:
 \item{delta }{the estimate, \eqn{\hat{\Delta}(t)}, described in delta.estimate documentation.}
 \item{delta.s }{the estimate, \eqn{\hat{\Delta}_S(t,t_0)}, described above.}
 \item{R.s }{the estimate, \eqn{\hat{R}_S(t,t_0)}, described above.}
 \item{delta.var }{the variance estimate of \eqn{\hat{\Delta}(t)}; if var = TRUE or conf.int = TRUE.}
 \item{delta.s.var }{the variance estimate of \eqn{\hat{\Delta}_S(t,t_0)}; if var = TRUE or conf.int = TRUE.}
 \item{R.s.var }{the variance estimate of \eqn{\hat{R}_S(t,t_0)}; if var = TRUE or conf.int = TRUE.}
 \item{conf.int.normal.delta}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}(t)} based on a normal approximation; if conf.int = TRUE.}
 \item{conf.int.quantile.delta}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}(t)} based on sample quantiles of the perturbed values, described above; if conf.int = TRUE.} 
 \item{conf.int.normal.delta.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}_S(t,t_0)} based on a normal approximation; if conf.int = TRUE.}
 \item{conf.int.quantile.delta.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}_S(t,t_0)} based on sample quantiles of the perturbed values, described above; if conf.int = TRUE.} 
 \item{conf.int.normal.R.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_S(t,t_0)} based on a normal approximation; if conf.int = TRUE.}
 \item{conf.int.quantile.R.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_S(t,t_0)} based on sample quantiles of the perturbed values, described above; if conf.int = TRUE.} 
\item{conf.int.fieller.R.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_S(t,t_0)} based on Fieller's approach, described above; if conf.int = TRUE.} 
  \item{delta.t }{the estimate, \eqn{\hat{\Delta}_T(t,t_0)}, described above; if incremental.vaue = TRUE.}
 \item{R.t }{the estimate, \eqn{\hat{R}_T(t,t_0)}, described above; if incremental.vaue = TRUE.}
  \item{incremental.value}{the estimate, \eqn{\hat{IV}_S(t,t_0)}, described above; if incremental.vaue = TRUE.}
 \item{delta.t.var }{the variance estimate of \eqn{\hat{\Delta}_T(t,t_0)}; if var = TRUE or conf.int = TRUE and  incremental.vaue = TRUE.}
 \item{R.t.var }{the variance estimate of \eqn{\hat{R}_T(t,t_0)}; if var = TRUE or conf.int = TRUE and incremental.vaue = TRUE.} 
  \item{incremental.value.var }{the variance estimate of \eqn{\hat{IV}_S(t,t_0)}; if var = TRUE or conf.int = TRUE and incremental.vaue = TRUE.}
 \item{conf.int.normal.delta.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}_T(t,t_0)} based on a normal approximation; if conf.int = TRUE and incremental.vaue = TRUE.}
 \item{conf.int.quantile.delta.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}_T(t,t_0)} based on sample quantiles of the perturbed values, described above; if conf.int = TRUE and incremental.vaue = TRUE.} 
 \item{conf.int.normal.R.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_T(t,t_0)} based on a normal approximation; if conf.int = TRUE and incremental.vaue = TRUE.}
 \item{conf.int.quantile.R.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_T(t,t_0)} based on sample quantiles of the perturbed values, described above; if conf.int = TRUE and incremental.vaue = TRUE.} 
\item{conf.int.fieller.R.t}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_T(t,t_0)} based on Fieller's approach, described above; if conf.int = TRUE and incremental.vaue = TRUE.} 
\item{conf.int.normal.iv}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{IV}_S(t,t_0)} based on a normal approximation; if conf.int = TRUE and incremental.vaue = TRUE.}
 \item{conf.int.quantile.iv}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{IV}_S(t,t_0)} based on sample quantiles of the perturbed values, described above; if conf.int = TRUE and incremental.vaue = TRUE.} 
}
\references{
%% ~put references to the literature/web site here ~
Parast L, Cai T and Tian L. Evaluating Surrogate Marker Information using Censored Data. Under Review.
}
\author{
%%  ~~who you are~~
Layla Parast
}
\note{If the treatment effect is not significant, the user will receive the following message: "Warning: it looks like the treatment effect is not significant; may be difficult to interpret the residual treatment effect in this setting". If the treatment effect is negative, the user will receive the following message: "Warning: it looks like you need to switch the treatment groups" as this package assumes throughout that larger values of the event time are better. If the observed support of the surrogate marker for the control group is outside the observed support of the surrogate marker for the treatment group, the user will receive the following message: "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation".}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
data(d_example_surv)
names(d_example_surv)
% Functions have elapsed time longer than 10s
%R.s.surv.estimate(xone = d_example_surv$x1, xzero = d_example_surv$x0,  
%deltaone = d_example_surv$delta1, deltazero = d_example_surv$delta0, 
%sone = d_example_surv$s1, szero = d_example_surv$s0, t=3, landmark = 1)
%R.s.surv.estimate(xone = d_example_surv$x1, xzero = d_example_surv$x0,  
%deltaone = d_example_surv$delta1, deltazero = d_example_surv$delta0, 
%sone = d_example_surv$s1, szero = d_example_surv$s0, t=3, landmark = 1, 
%incremental.value = TRUE)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
\keyword{survival}
\keyword{robust}
\keyword{smooth}