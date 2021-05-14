\name{R.s.estimate}
\alias{R.s.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates the proportion of treatment effect explained
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function calculates the proportion of treatment effect on the primary outcome explained by the treatment effect on the surrogate marker(s). This function is intended to be used for a fully observed continuous outcome. The user can also request a variance estimate and a 95\% confidence interval, both estimated using perturbating-resampling. If a confidence interval is requested three versions are provided: a normal approximation based interval, a quantile based interval, and Fieller's confidence interval.

}
\usage{
R.s.estimate(sone, szero, yone, yzero, var = FALSE, conf.int = FALSE, 
weight.perturb = NULL, number = "single", type = "robust",extrapolate = FALSE, 
transform = FALSE,warn.te = FALSE, warn.support = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sone}{
%%     ~~Describe \code{yone} here~~
numeric vector or matrix; surrogate marker for treated observations, assumed to be continuous. If there are multiple surrogates then this should be a matrix with \eqn{n_1} (number of treated observations) rows and n.s (number of surrogate markers) columns.
}
  \item{szero}{
%%     ~~Describe \code{yzero} here~~
numeric vector; surrogate marker for control observations, assumed to be continuous.If there are multiple surrogates then this should be a matrix with \eqn{n_0} (number of control observations) rows and n.s (number of surrogate markers) columns. 
}
  \item{yone}{
%%     ~~Describe \code{yone} here~~
numeric vector; primary outcome for treated observations, assumed to be continuous.
}
  \item{yzero}{
%%     ~~Describe \code{yzero} here~~
numeric vector; primary outcome for control observations, assumed to be continuous.
}
  \item{var}{
%%     ~~Describe \code{var} here~~
TRUE or FALSE; indicates whether a variance estimate is requested, default is FALSE.
}
  \item{conf.int}{
%%     ~~Describe \code{var} here~~
TRUE or FALSE; indicates whether a 95\% confidence interval is requested, default is FALSE
}
  \item{weight.perturb}{
%%     ~~Describe \code{weight} here~~
a \eqn{n_1+n_0} by \eqn{x} matrix of weights where \eqn{n_1 = } length of yone and \eqn{n_0 = }  length of yzero; used for perturbation-resampling, default is null.
}
  \item{number}{
%%     ~~Describe \code{number} here~~
specifies the number of surrogate markers; choices are "multiple" or "single", default is "single" 
}
  \item{type}{
%%     ~~Describe \code{type} here~~
specifies the type of estimation; choices are "robust" or "model" or "freedman", default is "robust"
}
  \item{extrapolate}{
%%     ~~Describe \code{extrapolate} here~~
TRUE or FALSE; indicates whether the user wants to use extrapolation.
}
  \item{transform}{
%%     ~~Describe \code{transform} here~~
TRUE or FALSE; indicates whether the user wants to use a transformation for the surrogate marker. 
}
 \item{warn.te}{
%%     ~~Describe \code{type} here~~
value to control warnings; user does not need to specify.
}
  \item{warn.support}{
%%     ~~Describe \code{type} here~~
value to control warnings; user does not need to specify.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
Let \eqn{Y^{(1)}} and \eqn{Y^{(0)}} denote the primary outcome under the treatment and primary outcome under the control,respectively. Let \eqn{S^{(1)}} and \eqn{S^{(0)}} denote the surrogate marker under the treatment and the surrogate marker under the control,respectively. The residual treatment effect is defined as
\deqn{ \Delta_S=\int_{-\infty}^{\infty}  E(Y^{(1)}|S^{(1)}=s) dF_0(s)  - \int_{-\infty}^{\infty} E(Y^{(0)}|S^{(0)}=s) dF_0(s),}
where \eqn{\Delta_S(s)= E(Y^{(1)}|S^{(1)}=s)-E(Y^{(0)}|S^{(0)}=s)}  and \eqn{F_0(\cdot)} is the marginal cumulative distribution function of \eqn{S^{(0)}}, the surrogate marker measure under the control. The proportion of treatment effect explained by the surrogate marker, which we denote by \eqn{R_S}, can be expressed using a contrast between \eqn{\Delta_S} and \eqn{\Delta}:
\deqn{R_S=\{\Delta-\Delta_S\}/\Delta=1-\Delta_S/\Delta.}   The definition and estimation of \eqn{\Delta} is described in the delta.estimate documentation.


A flexible model-based approach to estimate \eqn{\Delta_S} in the single marker setting is to specify:
\deqn{E(S^{(0)})=\alpha_0 \quad\mbox{and}\quad E(S^{(1)})-E(S^{(0)}) = \alpha_1,}
\deqn{ E(Y^{(0)} | S^{(0)}) = \beta_0 + \beta_1 S^{(0)}  \quad \mbox{and} \quad  E(Y^{(1)} | S^{(1)}) = (\beta_0 +\beta_2)+ (\beta_1+\beta_3) S^{(1)}. }
It can be shown that when these models hold, \eqn{\Delta_S = \beta_2 + \beta_3 \alpha_0}. Thus, reasonable estimates for \eqn{\Delta_S} and \eqn{R_S} using this approach would be \eqn{\hat{\Delta}_S =  \hat{\beta}_2 + \hat{\beta}_3 \hat{\alpha}_0} and \eqn{\hat{R}_S = 1-\hat{\Delta}_S / \hat{\Delta}.} 

For robust estimation of \eqn{\Delta_S} in the single marker setting, we estimate \eqn{\mu_1(s) = E(Y^{(1)}|S^{(1)}=s)} nonparametrically using kernel smoothing: 
\deqn{\hat{\mu}_1(s) = \frac{\sum_{i=1}^{n_1} K_h\left (S_{1i}-s \right ) Y_{1i} }{\sum_{i=1}^{n_1} K_h\left (S_{1i}-s \right )}}
where \eqn{S_{1i}} is the observed \eqn{S^{(1)}} for person \eqn{i}, \eqn{Y_{1i}} is the observed \eqn{Y^{(1)}} for person \eqn{i}, \eqn{K(\cdot)} is a smooth symmetric density function with finite support, \eqn{K_h(\cdot)=K(\cdot/h)/h} and \eqn{h} is a specified bandwidth. As in most nonparametric functional estimation procedures, the choice of the smoothing parameter \eqn{h}  is critical. To eliminate the impact of the bias of the conditional mean function on the resulting estimator, we require the standard undersmoothing assumption of \eqn{h=O(n_1^{-\delta})} with \eqn{\delta \in (1/4,1/3).}  To obtain an appropriate \eqn{h} we first use bw.nrd to obtain \eqn{h_{opt}}; and then we let \eqn{h = h_{opt}n_1^{-c_0}}  with \eqn{c_0 = 0.25}. We then estimate \eqn{\Delta_S} as
\deqn{ \hat{\Delta}_S= \sum_{i=1}^{n_0} \frac{\hat{\mu}_1(S_{0i})- Y_{0i}}{n_0}} where \eqn{S_{0i}} is the observed \eqn{S^{(0)}} for person \eqn{i} and \eqn{Y_{0i}} is the observed \eqn{Y^{(0)}} for person \eqn{i}. Lastly, we estimate \eqn{R_S} as \eqn{\hat{R}_S = 1-\hat{\Delta}_S/\hat{\Delta}}.

This function also allows for estimation of \eqn{R_S} using Freedman's approach.  Let \eqn{Y} denote the primary outcome, \eqn{S} denote the surrogate marker, and \eqn{G} denote the treatment group (0 for control, 1 for treatment). Freedman's approach to calculating the proportion of treatment effect explained by the surrogate marker is to fit the following two regression models:
\deqn{E(Y|G) = \gamma_0 + \gamma_1 I(G=1) \quad \mbox{and} \quad E(Y|G, S) = \gamma_{0S} + \gamma_{1S}I(G=1) + \gamma_{2S} S }
and estimating the proportion of treatment effect explained, denoted by \eqn{R_S}, as \eqn{ 1-\hat{\gamma}_{1S}/\hat{\gamma}_1}.
 

This function also estimates \eqn{R_S} in a multiple marking setting. A flexible model-based approach to estimate \eqn{\Delta_S} in the multiple marker setting is to specify models for \eqn{E(Y|G, S)} and \eqn{E(S_j | G)} for each \eqn{S_j} in \eqn{S = \{S_1,...S_p\}} (where p is the number of surrogate markers). Without loss of generality, consider the case where there are three surrogate markers, \eqn{S = \{S_1, S_2, S_3\}} and one specifies the following linear models:
\deqn{E(Y^{(0)} | S^{(0)}) = \beta_0 + \beta_1 S_1^{(0)} + \beta_2 S_2^{(0)} +  \beta_3 S_3^{(0)}}
\deqn{E(Y^{(1)} | S^{(1)}) = (\beta_0+\beta_4) + (\beta_1+\beta_5) S_1^{(1)} + (\beta_2+\beta_6) S_2^{(1)} +  (\beta_3+\beta_7) S_3^{(1)}}
\deqn{ E(S_j^{(0)}) = \alpha_j, ~~~~j=1,2,3.} 
It can be shown that when these models hold 
\deqn{\Delta_{S} = \beta_4 + \beta_5\alpha_1 + \beta_6 \alpha_2 + \beta_7 \alpha_3.} Thus, reasonable estimates for \eqn{\Delta_{S}} and \eqn{R_{S}} here would be easily obtained by replacing the unknown regression coefficients in the models above by their consistent estimators. 


For robust estimation of S \eqn{\Delta_S} in the multiple marker setting, we use a two-stage procedure combining the model-based approach and the nonparametric estimation procedure from the single marker setting. Specifically, we use a working semiparametric  model: \deqn{E(Y^{(1)}|S^{(1)}=S)=\beta_0 + \beta_1 S_1^{(1)} + \beta_2 S_2^{(1)} +  \beta_3 S_3^{(1)}} 
 and define \eqn{Q^{(1)} = \hat{\beta}_0 + \hat{\beta}_1 S_1^{(1)} + \hat{\beta}_2 S_2^{(1)} +  \hat{\beta}_3 S_3^{(1)}} and \eqn{Q^{(0)} = \hat{\beta}_0 + \hat{\beta}_1 S_1^{(0)} + \hat{\beta}_2 S_2^{(0)} +  \hat{\beta}_3 S_3^{(0)}} to reduce the dimension of \eqn{S} in the first stage  and in the second stage, we apply the robust approach used in the single marker setting to  estimate its surrogacy. 

To use Freedman's approach in the presence of multiple markers, the markers are simply additively entered into the second regression model.

Variance estimation and confidence interval construction are performed using perturbation-resampling. Specifically, let \eqn{\left \{ V^{(b)} = (V_{11}^{(b)}, ...V_{1n_1}^{(b)}, V_{01}^{(b)}, ...V_{0n_0}^{(b)})^T, b=1,....,D \right \}} be \eqn{n \times D} independent copies of a positive random variables \eqn{V} from a known distribution with unit mean and unit variance. Let
\deqn{\hat{\Delta}^{(b)}  = \frac{  \sum_{i=1}^{n_1} V_{1i}^{(b)} Y_{1i}}{ \sum_{i=1}^{n_1} V_{1i}^{(b)}}  - \frac{  \sum_{i=1}^{n_0} V_{0i}^{(b)} Y_{0i}}{ \sum_{i=1}^{n_0} V_{0i}^{(b)}}.} The variance of \eqn{\hat{\Delta}} is obtained as the empirical variance of \eqn{\{\hat{\Delta}^{(b)}, b = 1,...,D\}.} In this package, we use weights generated from an Exponential(1) distribution and use \eqn{D=500}. Variance estimates for \eqn{\hat{\Delta}_S} and \eqn{\hat{R}_S} are calculated similarly. We construct two versions of the \eqn{95\%} confidence interval for each estimate: one based on a normal approximation confidence interval using the estimated variance and another taking the 2.5th and 97.5th empirical percentile of the perturbed quantities.  In addition, we use Fieller's method to obtain a third confidence interval for \eqn{R_S} as
\deqn{\left\{1-r:  \frac{(\hat{\Delta}_S-r\hat{\Delta})^2}{\hat{\sigma}_{11}-2r\hat\sigma_{12}+r^2\hat\sigma_{22}} \le c_{\alpha}\right\},}
where \eqn{\hat{\Sigma}=(\hat\sigma_{ij})_{1\le i,j\le 2}} and \eqn{c_\alpha} is the \eqn{(1-\alpha)}th percentile of
\deqn{\left\{\frac{\{\hat{\Delta}^{(b)}_S-(1-\hat R_S)\hat{\Delta}^{(b)}\}^2}{\hat{\sigma}_{11}-2(1-\hat R_S)\hat\sigma_{12}+(1-\hat R_S)^2\hat\sigma_{22}}, b=1, \cdots, C\right\}}
where \eqn{\alpha=0.05}.

Note that if the observed supports for S are not the same, then \eqn{\hat{\mu}_1(s)} for \eqn{S_{0i} = s} outside the support of \eqn{S_{1i}} may return NA (depending on the bandwidth). If extrapolation = TRUE, then the   \eqn{\hat{\mu}_1(s)} values for these surrogate values are set to the closest non-NA value. If transform = TRUE, then \eqn{S_{1i}} and \eqn{S_{0i}} are transformed such that the new transformed values, \eqn{S^{tr}_{1i}} and \eqn{S^{tr}_{0i}} are defined as: \eqn{S^{tr}_{gi} = F([S_{gi} - \mu]/\sigma)} for \eqn{g=0,1} where \eqn{F(\cdot)} is the cumulative distribution function for a standard normal random variable, and \eqn{\mu} and \eqn{\sigma} are the sample mean and standard deviation, respectively, of \eqn{(S_{1i}, S_{0i})^T}.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list is returned:
 \item{R.s }{the estimate, \eqn{\hat{R}_S}, described above.}
 \item{R.s.var }{the variance estimate of \eqn{\hat{R}_S}; if var = TRUE or conf.int = TRUE.}
 \item{conf.int.normal.R.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_S} based on a normal approximation; if conf.int = TRUE.}
 \item{conf.int.quantile.R.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_S} based on sample quantiles of the perturbed values, described above; if conf.int = TRUE.} 
\item{conf.int.fieller.R.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_S} based on Fieller's approach, described above; if conf.int = TRUE.} 
For all options other then "freedman", the following are also returned:
 \item{delta }{the estimate, \eqn{\hat{\Delta}}, described in delta.estimate documentation.}
 \item{delta.s }{the estimate, \eqn{\hat{\Delta}_S}, described above.}
 \item{delta.var }{the variance estimate of \eqn{\hat{\Delta}}; if var = TRUE or conf.int = TRUE.}
 \item{delta.s.var }{the variance estimate of \eqn{\hat{\Delta}_S}; if var = TRUE or conf.int = TRUE.}
 \item{conf.int.normal.delta}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}} based on a normal approximation; if conf.int = TRUE.}
 \item{conf.int.quantile.delta}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}} based on sample quantiles of the perturbed values, described above; if conf.int = TRUE.} 
 \item{conf.int.normal.delta.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}_S} based on a normal approximation; if conf.int = TRUE.}
 \item{conf.int.quantile.delta.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}_S} based on sample quantiles of the perturbed values, described above; if conf.int = TRUE.} 
}
\references{
%% ~put references to the literature/web site here ~
Freedman, L. S., Graubard, B. I., & Schatzkin, A. (1992). Statistical validation of intermediate endpoints for chronic diseases. Statistics in medicine, 11(2), 167-178.

Parast, L., McDermott, M., Tian, L. (2016). Robust estimation of the proportion of treatment effect explained by surrogate marker information. Statistics in Medicine, 35(10):1637-1653.

Wang, Y., & Taylor, J. M. (2002). A measure of the proportion of treatment effect explained by a surrogate marker. Biometrics, 58(4), 803-812.

Fieller, Edgar C. (1954). Some problems in interval estimation. Journal of the Royal Statistical Society. Series B (Methodological), 175-185.

Fieller, E. C. (1940). The biological standardization of insulin. Supplement to the Journal of the Royal Statistical Society, 1-64.
}
\note{ If the treatment effect is not significant, the user will receive the following message: "Warning: it looks like the treatment effect is not significant; may be difficult to interpret the proportion of treatment effect explained in this setting". In the single marker case with the robust estimation approach, if the observed support of the surrogate marker for the control group is outside the observed support of the surrogate marker for the treatment group, the user will receive the following message: "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"}	
\author{
%%  ~~who you are~~
Layla Parast
}
\examples{
data(d_example)
names(d_example)
R.s.estimate(yone=d_example$y1, yzero=d_example$y0, sone=d_example$s1.a, szero=d_example$s0.a, 
number = "single", type = "robust")
R.s.estimate(yone=d_example$y1, yzero=d_example$y0, sone=cbind(d_example$s1.a,d_example$s1.b, 
d_example$s1.c), szero=cbind(d_example$s0.a, d_example$s0.b, d_example$s0.c), 
number = "multiple", type = "model")
}

\keyword{regression}
\keyword{models}
\keyword{smooth}
\keyword{robust}
\keyword{nonparametric}
