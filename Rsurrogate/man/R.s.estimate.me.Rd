\name{R.s.estimate.me}
\alias{R.s.estimate.me}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates the proportion of treatment effect explained correcting for measurement error in the surrogate marker

}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function calculates the proportion of treatment effect on the primary outcome explained by the treatment effect on a surrogate marker, correcting for measurement error in the surrogate marker. This function is intended to be used for a fully observed continuous outcome. The user must specify what type of estimation they would like (parametric or nonparametric estimation of the proportion explained, denoted by R) and what estimator they would like (see below for details). 

}
\usage{
R.s.estimate.me(sone, szero, yone, yzero, parametric = FALSE, estimator = "n", 
me.variance, extrapolate = TRUE, transform = FALSE, naive = FALSE, Ronly = TRUE)
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
  \item{parametric}{
%%     ~~Describe \code{parametric} here~~
TRUE or FALSE; indicates whether the user wants the parametric approach to be used (TRUE) or nonparametric (FALSE).
}
  \item{estimator}{
%%     ~~Describe \code{estimator} here~~
options are "d","q","n" for parametric and "q","n" for nonparametric; "d" stands for the disattenuated estimator, "q" stands for the SIMEX estimator with quadratic extrapolation, "n" stands for the SIMEX estimator with a nonlinear extrapolation. Note that the nonlinear extrapolation may have convergence issues with a small sample size; if this occurs, please consider using quadratic extrapolation instead.
}
  \item{me.variance}{
%%     ~~Describe \code{me.variance} here~~
the variance of the measurement error; must be provided.
}
  \item{extrapolate}{
%%     ~~Describe \code{extrapolate} here~~
TRUE or FALSE; indicates whether the user wants to use extrapolation.
}
  \item{transform}{
%%     ~~Describe \code{transform} here~~
TRUE or FALSE; indicates whether the user wants to use a transformation for the surrogate marker. 
}
\item{naive}{
%%     ~~Describe \code{naive} here~~
TRUE or FALSE; indicates whether the user wants the naive estimate (not correcting for measurement error) to also be calculated
}
  \item{Ronly}{
%%     ~~Describe \code{Ronly} here~~
TRUE or FALSE; indicates whether the user wants only R (and corresponding variance and confidence intervals) to be returned.
}
}
\details{
While there are many methods available to quantify the value of a surrogate marker, most assume that the marker is measured without error. This function calculates the proportion of treatment effect on the primary outcome explained by the treatment effect on a surrogate marker, correcting for measurement error in the surrogate marker. The user can choose either the parametric framework or nonparametric framework for estmation. Within the parametric framework there are three options for measurement error correction: the disattenuated estimator, the SIMEX estimator with quadratic extrapolation, and the SIMEX estimator with nonlinear extrapolation. Within the nonparametric framework there are two options for measurement error correction: the SIMEX estimator with quadratic extrapolation and the SIMEX estimator with nonlinear extrapolation. We describe each below. 

Let \eqn{G} be the binary treatment indicator with \eqn{G=1} indicating treatment and \eqn{G=0} indicating control (or placebo).  We assume throughout that subjects are randomly assigned to treatment or control at baseline. Let \eqn{Y} and \eqn{S} denote the continuous primary outcome and continuous surrogate marker, respectively, where \eqn{S} is measured post-baseline and is assumed to be a biomarker, clinical measurement, psychological test score, or other  physiological measurement. In the absence of measurement error, the observed data consists of \eqn{\{Y_i, S_i, G_i\}} for \eqn{i \in \{1,...,n\}}. With measurement error, instead of observing \eqn{S} we observe \eqn{W = S + U}, where \eqn{E(U|S) = 0} and the variance of \eqn{U} is \eqn{\sigma_u^2}. Such measurement error may be attributable to, for example, laboratory error. Thus, our observed data will consist of \eqn{\{Y_i, W_i, G_i\}} for \eqn{i \in \{1,...,n\}}. Throughout, we assume that  \eqn{\sigma_u^2} is known. Here, we are interested in estimating the proprtion of the treatment effect on the primary outcome that is explained by the treatment effect on the surrogate marker, denoted as \eqn{R_S}. 

To estimate \eqn{R_S} parametrically, we assume the following models \eqn{E(Y|G) = \beta_0 + \beta_1 G} and \eqn{E(Y|G,S) = \beta_0^* + \beta_1^*G + \beta_2^* S}. It can be shown that if these models hold,  \eqn{R_S=1-\beta_1^*/\beta_1}. When \eqn{W = S+U} is available instead of \eqn{S}, this measurement error does not affect estimation of \eqn{\beta_1}, but it does affect estimation of \eqn{\beta_1^*}, and \eqn{\beta_2^*}. Since estimation of \eqn{R_S} relies on estimation of \eqn{\beta_1} and \eqn{\beta_1^*}, we focus on the effect of measurement error on \eqn{\beta_1^*} estimation. The attenuation bias for \eqn{\hat \beta_1^*} and \eqn{\hat R} can be written out in closed form when the proportion of treatment effect is parametrically estimated as described above, when these specified models hold, and when the surrogate marker \eqn{S} is measured with error. There exist two methods to eliminate this bias when estimating \eqn{R_S}. Taking advantage of the fact that we can express the attenuation bias in closed form, the first is a straightforward disattenuated estimator: \eqn{\hat \beta _{1A} = \hat{\beta}_1^* -  \frac{  \hat{\beta}_2^* \{\Omega^2_{W} \Omega_{GW}-\Omega_{GW}(\Omega^2_{W} - \sigma_u^2)\}}{\Omega^2_{G}(\Omega^2_{W} - \sigma_u^2)-\Omega_{GW}\Omega_{GW}}} and \eqn{\hat{R}_{A} = 1- \left [ \hat{\beta}_1^* -  \frac{  \hat{\beta}_2^* \{\Omega^2_{W} \Omega_{GW}-\Omega_{GW}(\Omega^2_{W} - \sigma_u^2)\}}{\Omega^2_{G}(\Omega^2_{W} - \sigma_u^2)-\Omega_{GW}\Omega_{GW}}  \right] / \hat{\beta}_1} where \eqn{\Omega^2} denotes the sample variance or covariance.
		
The second method to eliminate this bias uses Simulation Extrapolation (SIMEX) estimation, which is a simulation-based method that involves first generating additional measurement error and observing how it affects the bias of the parameter estimate of interest, and then extrapolating this information to a setting with no measurement error. To incorporate SIMEX estimation within our surrogate marker framework, we define \eqn{W_{b,i}(\lambda) = W_i + \lambda^{1/2} \sigma_u \epsilon_{i,b}}
for \eqn{b=1,...,B} where \eqn{B=50}, \eqn{\epsilon_{i,b} \sim N(0,1)}, \eqn{\sigma_u} is assumed known, and \eqn{\lambda \in (0,0.25,0.5,0.75,1.0,} \eqn{1.25,1.5,1.75,2.0)} and for each iteration \eqn{b} and \eqn{\lambda} value, obtaining  \eqn{\hat \beta_{1b}^*(\lambda)} by fitting the regression model: \eqn{E(Y \mid W_b(\lambda),S) = \beta_{0b}^* + \beta_{1b}^* W_{b}(\lambda) + \beta_{2b}^* S.} We then calculate the average estimate for each quantity over the iterations \eqn{b=1,...,B} for each \eqn{\lambda} value, denoted as \eqn{\hat \beta^*_{1,S,\sigma^2_u(1+\lambda)} = \sum_{b=1}^B \hat \beta_{1b}^*(\lambda)}. The second step, extrapolation, takes these average estimates for each \eqn{\lambda} value and extrapolates using a function \eqn{G(\Gamma, \lambda)} to obtain the estimated quantity if \eqn{\lambda=-1}. For the extrapolation step, we use both a quadratic extrapolation and nonlinear extrapolation i.e., we solve for \eqn{\Gamma = (\alpha_0, \alpha_1, \alpha_2)^T} in  \eqn{\hat \beta^*_{1,S,\sigma^2_u(1+\lambda)} = \alpha_0 + \alpha_1 \lambda + \alpha_2 \lambda^2} and \eqn{\hat \beta^*_{1,S,\sigma^2_u(1+\lambda)}= \alpha_0 + \alpha_1 /( \alpha_2 + \lambda)}, respectively. Using the estimates of \eqn{\alpha_0, \alpha_1, \alpha_2}, we calculate the predicted \eqn{\hat \beta^*_{1,S,\sigma^2_u(1+\lambda)}} when \eqn{\lambda = -1}. In essence, the simulations add successively larger measurement errors of size \eqn{(1+\lambda)\sigma^2_u} and then extrapolate to the case when \eqn{\lambda = -1} such that the measurement error is 0. We denote the resulting estimator of \eqn{\beta_1^*} as \eqn{\hat{\beta}^*_{1,SIMEX} = G(\hat \Gamma, -1)} and define \eqn{\hat{R}_{SIMEX} = 1- \hat{\beta}^*_{1,SIMEX}/ \hat \beta_1.} 
	
While the parametric approach to estimate the proportion of treatment effect explained by \eqn{S} is most commonly used in clinical practice, previous work has demonstrated biased results when the assumed models are not correctly specified. An alternative approach involves estimating the treatment effect, \eqn{\Delta}, and residual treatment effect, \eqn{\Delta_S}, as \eqn{R_S} is defined as \eqn{1-\Delta/\Delta_S}. The quantity \eqn{\Delta} can be estimated simply by \eqn{\hat{\Delta} = n_1^{-1}\sum_{i=1}^{n} Y_i I(G_i = 1) - n_0^{-1}\sum_{i=1}^{n} Y_i I(G_i = 0)}, where \eqn{n_1} and \eqn{n_0} denote the number of individuals in the treatment and control groups, respectively. The quantity \eqn{\Delta_S} can be estimated nonparametrically using kernel smoothing as \eqn{\hat{\Delta}_S = n_0^{-1} \sum_{i: G_i = 0}\hat{\mu}_1(S_i) - n_0^{-1}\sum_{i=1}^{n} Y_i I(G_i = 0)}
where \eqn{\hat{\mu}_1(s) = \{ \sum_{j: G_j = 1} K_h(S_j - s)Y_j \}/ \{\sum_{j:G_j = 1} K_h(S_j - s)\}}, \eqn{K(\cdot)} is a smooth symmetric density function with finite support, \eqn{K_h(\cdot)=K(\cdot/h)/h} and \eqn{h} is a specified bandwidth such that \eqn{h=O(n_1^{-\nu})} with \eqn{\nu \in (1/4,1/2).} 

When \eqn{W = S + U} is available instead of \eqn{S}, estimation of \eqn{\Delta} is not affected whereas estimation of \eqn{\Delta_S} is affected and thus, the nonparametric estimation procedure described above results in a biased estimate of \eqn{R_S}. Unlike the parametric approach, the attenuation bias cannot be expressed in closed form. Within this nonparametric framework, SIMEX estimation can be used to correct for measurement error. We implement the estimation procedure as described above where we first generate additional measurement error to obtain \eqn{W_{b,i}(\lambda)} and for each iteration \eqn{b} and \eqn{\lambda} values obtain \eqn{\hat{\Delta}_{S,b}(\lambda) = n_0^{-1} \sum_{i: G_i = 0} \left \{ \frac{\sum_{j: G_j = 1} K_h(W_{b,j}(\lambda) - W_{b,i}(\lambda))Y_j}{\sum_{j:G_j = 1} K_h(W_{b,j}(\lambda)- W_{b,i}(\lambda))} \right \} - n_0^{-1}\sum_{i=1}^{n} Y_i I(G_i = 0).} We then calculate the average estimate for each quantity over the iterations \eqn{b=1,...,B} for each \eqn{\lambda} value, denoted as \eqn{\hat{\Delta}_{S,\sigma_u^2(1+\lambda)} = \sum_{b=1}^B \hat{\Delta}_{S,b}(\lambda)} and extrapolate using a function \eqn{G(\Gamma, \lambda)}; we specifically use the quadratic and nonlinear functions as in the parametric setting. We denote the resulting estimator of \eqn{\Delta_S} as \eqn{\hat{\Delta}_{S,SIMEX} = G(\hat \Gamma, -1)} and define \eqn{\hat{R}_{S,SIMEX} = 1- \hat{\Delta}_{S,SIMEX} / \hat \Delta.} 

In this function, parametric estimation is equivalent to Freedman's approach in the R.s.estimate documentation; nonparametric estimation is equivalent to the robust approach in the R.s.estimate documentation. Variance estimates for all estimators are calculated in this function based on derived closed form variance expressions. For all approaches, confidence intervals for \eqn{\Delta_S} can be constructed using a normal approximation; confidence intervals for \eqn{R_S} can be constructed using either a normal approximation or using Fieller's method, all of which are provided in this function. Details regarding the asymptotic properties of these estimators and closed form variance calculation can be found in: Parast, L., Garcia, TP, Prentice, RL, Carroll, RJ (2021). Robust Methods to Correct for Measurement Error when Evaluating a Surrogate Marker. Biometrics, In press. 
}
\value{
A list is returned:

\item{R.naive}{the naive estimate of the proportion of treatment effect explained by the surrogate marker; only if naive = TRUE}
\item{R.naive.var}{the estimated variance of the naive estimate of the proportion of treatment effect explained by the surrogate marker; only if naive = TRUE}
\item{R.naive.CI.normal}{the 95\% confidence interval using the normal approximation for the naive estimate of the proportion of treatment effect explained by the surrogate marker; only if naive = TRUE}
\item{R.naive.CI.fieller}{the 95\% confidence interval using Fieller's approach for the naive estimate of the proportion of treatment effect explained by the surrogate marker; only if naive = TRUE}
\item{B1star.naive}{the naive estimate of the adjusted regression coefficient for treatment; only if naive = TRUE and Ronly = FALSE and parametric = TRUE}
\item{B1star.naive.var}{the estimated variance of the naive estimate of the adjusted regression coefficient for treatment; only if naive = TRUE and Ronly = FALSE and parametric = TRUE}
\item{B1star.naive.CI.normal}{the 95\% confidence interval using the normal approximation for the naive estimate of the adjusted regression coefficient for treatment; only if naive = TRUE and Ronly = FALSE and parametric = TRUE}
\item{deltas.naive}{the naive estimate of the residual treatment effect; only if naive = TRUE and Ronly = FALSE and parametric = FALSE}
\item{deltas.naive.var}{the estimated variance of the naive estimate of the residual treatment effect; only if naive = TRUE and Ronly = FALSE and parametric = FALSE}
\item{deltas.naive.CI.normal}{the 95\% confidence interval using the normal approximation for the naive estimate of the residual treatment effect; only if naive = TRUE and Ronly = FALSE and parametric = FALSE}
\item{R.corrected.dis}{the corrected disattenuated estimate of the proportion of treatment effect explained by the surrogate marker; only if parametric = TRUE and estimator ="d"}
\item{R.corrected.var.dis}{the estimated variance of the corrected disattenuated estimate of the proportion of treatment effect explained by the surrogate marker; only if naive = TRUE}
\item{R.corrected.CI.normal.dis}{the 95\% confidence interval using the normal approximation for the corrected disattenuated estimate of the proportion of treatment effect explained by the surrogate marker; only if parametric = TRUE and estimator ="d"}
\item{R.corrected.CI.fieller.dis}{the 95\% confidence interval using Fieller's approach for the corrected disattenuated estimate of the proportion of treatment effect explained by the surrogate marker; only if parametric = TRUE and estimator ="d"}
\item{B1star.corrected.dis}{the corrected disattenuated estimate of the adjusted regression coefficient for treatment; only if parametric = TRUE and estimator = "d" and Ronly = FALSE}
\item{B1star.corrected.var.dis}{the estimated variance of the corrected disattenuated estimate of the adjusted regression coefficient for treatment; only if parametric = TRUE and estimator = "d" and Ronly = FALSE}
\item{B1star.corrected.CI.normal.dis}{the 95\% confidence interval using the normal approximation for the corrected disattenuated estimate of the adjusted regression coefficient for treatment; only if parametric = TRUE and estimator = "d" and Ronly = FALSE}
\item{R.corrected.q}{the corrected SIMEX (quadratic) estimate of the proportion of treatment effect explained by the surrogate marker; only if estimator = "q"}
\item{R.corrected.var.q}{the estimated variance of the corrected SIMEX (quadratic) estimate of the proportion of treatment effect explained by the surrogate marker; only if estimator = "q"}
\item{R.corrected.CI.normal.q}{the 95\% confidence interval using the normal approximation for the corrected SIMEX (quadratic) estimate of the proportion of treatment effect explained by the surrogate marker; only if estimator = "q"}
\item{R.corrected.CI.fieller.q}{the 95\% confidence interval using Fieller's approach for the corrected SIMEX (quadratic) estimate of the proportion of treatment effect explained by the surrogate marker; only if estimator = "q"}
\item{B1star.corrected.q}{the corrected SIMEX (quadratic) estimate of the adjusted regression coefficient for treatment; only if estimator = "q" and Ronly = FALSE and parametric = TRUE}
\item{B1star.corrected.var.q}{the estimated variance of the corrected SIMEX (quadratic) estimate of the adjusted regression coefficient for treatment; only if estimator = "q" and Ronly = FALSE and parametric = TRUE}
\item{B1star.corrected.CI.normal.q}{the 95\% confidence interval using the normal approximation for the corrected SIMEX (quadratic) estimate of the adjusted regression coefficient for treatment; only if estimator = "q" and Ronly = FALSE and parametric = TRUE}
\item{deltas.corrected.q}{the corrected SIMEX (quadratic) estimate of the residual treatment effect; only if estimator = "q" and Ronly = FALSE and parametric = FALSE}
\item{deltas.corrected.var.q}{the estimated variance of the corrected SIMEX (quadratic) estimate of the residual treatment effect; only if estimator = "q" and Ronly = FALSE and parametric = FALSE}
\item{deltas.corrected.CI.normal.q}{the 95\% confidence interval using the normal approximation for the corrected SIMEX (quadratic) estimate of the residual treatment effect; only if estimator = "q" and Ronly = FALSE and parametric = FALSE}
\item{R.corrected.nl}{the corrected SIMEX (nonlinear) estimate of the proportion of treatment effect explained by the surrogate marker; only if estimator = "q"}
\item{R.corrected.var.nl}{the estimated variance of the corrected SIMEX (nonlinear) estimate of the proportion of treatment effect explained by the surrogate marker; only if estimator = "q"}
\item{R.corrected.CI.normal.nl}{the 95\% confidence interval using the normal approximation for the corrected SIMEX (nonlinear) estimate of the proportion of treatment effect explained by the surrogate marker; only if estimator = "q"}
\item{R.corrected.CI.fieller.nl}{the 95\% confidence interval using Fieller's approach for the corrected SIMEX (nonlinear) estimate of the proportion of treatment effect explained by the surrogate marker; only if estimator = "q"}
\item{B1star.corrected.nl}{the corrected SIMEX (nonlinear) estimate of the adjusted regression coefficient for treatment; only if estimator = "q" and Ronly = FALSE and parametric = TRUE}
\item{B1star.corrected.var.nl}{the estimated variance of the corrected SIMEX (nonlinear) estimate of the adjusted regression coefficient for treatment; only if estimator = "q" and Ronly = FALSE and parametric = TRUE}
\item{B1star.corrected.CI.normal.nl}{the 95\% confidence interval using the normal approximation for the corrected SIMEX (nonlinear) estimate of the adjusted regression coefficient for treatment; only if estimator = "q" and Ronly = FALSE and parametric = TRUE}
\item{deltas.corrected.nl}{the corrected SIMEX (nonlinear) estimate of the residual treatment effect; only if estimator = "q" and Ronly = FALSE and parametric = FALSE}
\item{deltas.corrected.var.nl}{the estimated variance of the corrected SIMEX (nonlinear) estimate of the residual treatment effect; only if estimator = "q" and Ronly = FALSE and parametric = FALSE}
\item{deltas.corrected.CI.normal.nl}{the 95\% confidence interval using the normal approximation for the corrected SIMEX (nonlinear) estimate of the residual treatment effect; only if estimator = "q" and Ronly = FALSE and parametric = FALSE} 
}
\references{
%% ~put references to the literature/web site here ~
Parast, L., Garcia, TP, Prentice, RL, Carroll, RJ (2021). Robust Methods to Correct for Measurement Error when Evaluating a Surrogate Marker. Biometrics, In press. 
}
\author{
%%  ~~who you are~~
Layla Parast
}

\examples{
data(d_example_me)
names(d_example_me)
R.s.estimate.me(yone=d_example_me$y1, yzero=d_example_me$y0, sone=d_example_me$s1, 
szero=d_example_me$s0, parametric = TRUE, estimator = "d", me.variance = 0.5, 
naive= TRUE, Ronly = FALSE)
R.s.estimate.me(yone=d_example_me$y1, yzero=d_example_me$y0, sone=d_example_me$s1, 
szero=d_example_me$s0, parametric = TRUE, estimator = "q", me.variance = 0.5, 
naive= FALSE, Ronly = TRUE)

#estimating measurement error variance with replicates
replicates = rbind(cbind(d_example_me$s1_rep1, d_example_me$s1_rep2, 
d_example_me$s1_rep3), cbind(d_example_me$s0_rep1, d_example_me$s0_rep2, 
d_example_me$s0_rep3))
mean.i = apply(replicates,1,mean, na.rm = TRUE)
num.i = apply(replicates,1,function(x) sum(!is.na(x)))
var.u = sum((replicates-mean.i)^2, na.rm = TRUE)/sum(num.i)
var.u
R.s.estimate.me(yone=d_example_me$y1, yzero=d_example_me$y0, sone=d_example_me$s1, 
szero=d_example_me$s0, parametric = TRUE, estimator = "d", me.variance = var.u, 
naive= TRUE, Ronly = FALSE)

\donttest{R.s.estimate.me(yone=d_example_me$y1, yzero=d_example_me$y0, 
sone=d_example_me$s1, szero=d_example_me$s0, parametric = FALSE, estimator = "q", 
me.variance = 0.5, naive= FALSE, Ronly = TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
\keyword{models}
\keyword{smooth}
\keyword{robust}
\keyword{nonparametric}

