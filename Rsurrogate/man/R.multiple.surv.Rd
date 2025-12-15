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
type of estimate; options are 1 = two-stage robust estimator, 2 = weighted two-stage robust estimator, 3 = double-robust estimator, 4 = two-stage model-based estimator, 5 = weighted estimator, 6 = double-robust model-based estimator; default is 1.
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
Let \eqn{G \in \{A, B\}} be the binary treatment indicator and we assume that subjects are randomly assigned to either treatment group \eqn{A} or \eqn{B} at baseline. Let \eqn{T} denote the time to the occurrence of the primary outcome, death for example, and \eqn{S = (S_1,S_2,...,S_k)^{T}} denote the vector of \eqn{k} surrogate markers measured at a given time \eqn{t_0}. Let \eqn{T^{(g)}} and \eqn{S^{(g)}} denote the counterfactual  event time and surrogate marker measurements under treatment \eqn{G = g} for \eqn{g \in \{A, B\}}.  In practice, we only observe \eqn{(T, S)=(T^{(A)}, S^{(A)})} or \eqn{(T^{(B)}, S^{(B)})} depending on whether \eqn{G=A} or \eqn{B.} The treatment effect, \eqn{\Delta(t)}, is the treatment difference in survival rates at time \eqn{t > t_0},
\eqn{\Delta(t)=E\{ I(T^{(A)}>t)\} - E\{I(T^{(B)}>t)\} = P(T^{(A)}>t) - P(T^{(B)}>t)}
where \eqn{I(\cdot)} is the indicator function. For individuals who are censored or experience the primary outcome before \eqn{t_0}, we assume that their \eqn{S} information is not available.

The surrogate marker information at time \eqn{t_0} is defined as a combination of the observed information on \eqn{I(T >t_0)} and the observed \eqn{S} at \eqn{t_0}, denoted by \eqn{Q_{t_0}  = \{Q_{t_0}^{(g)}, g = A, B\}}, where 
\eqn{Q_{t_0} ^{(g)} = \{ T ^{(g)} \wedge  t_0, S ^{(g)} I(T ^{(g)} > t_0)\}}. With information on \eqn{Q_{t_0}}, the residual treatment effect is defined as:
\eqn{\Delta_{S}(t,t_0) = E\{ I(T ^{(A)} > t) - I(T ^{(B)}>t)  \mid Q_{t_0}^{(A)} = Q_{t_0}^{(B)} \}
 =P(T^{(B)} > t_0) \int_{S}  \psi_A(t \mid S, t_0) dF_{B}(S \mid t_0) -P(T^{(B)}> t)}
where \eqn{S = (s_1, ..., s_k)^{T}}, 
\eqn{F_{g}(S \mid t_0)  = P(S^{(g)}\le S \mid T^{(g)}>t_0),}  
\eqn{\psi_g(t\mid S,t_0) = P(T^{(g)}> t\mid S^{(g)}=S, T^{(g)}> t_0).} The proportion of the treatment effect on the primary outcome that is explained by the treatment effect on \eqn{Q_{t_0}}  is \eqn{R_{S}(t,t_0)=1-\Delta_{S}(t,t_0)/\Delta}. This function provides 6 different estimators for \eqn{R_{S}(t,t_0)} using censored data. 

Due to censoring, the observed data consist of \eqn{n} observations \eqn{\{(G_i, X_{i}, \delta_{i}, S_{i}I(X_{i} > t_0)), i=1,...,n\}} from the two treatment groups, where \eqn{X_{i} = \min(T_{i}, C_{i})}, \eqn{\delta_{i} = I(T_{i} < C_{i})},  and \eqn{C_{i}} denotes the censoring time for the \eqn{i}th subject. We assume independent censoring i.e., \eqn{(T_i, S_i) \perp C_i \mid G_i}. For ease of notation, we also let \eqn{\{(X_{gi}, \delta_{gi}, S_{gi}I(X_{gi} > t_0)), i=1,...,n_g\}} denote \eqn{n_g = \sum_{i=1}^n I(G_i = g)} observations from treatment group \eqn{g \in \{A,B\},} where \eqn{X_{gi}=\min(T_{gi}, C_{gi})} and \eqn{\delta_{gi}=I(T_{gi}<C_{gi}).} Without loss of generality, we assume that \eqn{\bar{\pi}_g = n_g/n \to \pi_g \in (0,1)} as \eqn{n\to \infty}. Throughout, we estimate the treatment effect \eqn{\Delta(t) =P(T^{(A)}>t) - P(T^{(B)}>t)} as
\eqn{ \hat{\Delta}(t)  = n_A^{-1} \sum_{i=1}^{n_A} \frac{I(X_{Ai}>t)}{\hat{W}^C_A(t)} - n_B^{-1} \sum_{i=1}^{n_B} \frac{I(X_{Bi}>t)}{\hat{W}^C_B(t)}}
where  \eqn{\hat{W}^C_g(t)} is the Kaplan-Meier estimator of \eqn{W^C_g(t) = P(C_{g} > t)}. 

We first describe the two-stage robust estimator which involves a two-stage procedure combining the use of a working model and a nonparametric estimation procedure for \eqn{\Delta_{S}(t,t_0)}. The idea is simply to summarize \eqn{S} into a univariate score \eqn{U} and then construct a nonparametric estimator for \eqn{R_S(t,t_0)} treating \eqn{U} as \eqn{S}. To construct \eqn{U}, we approximate the conditional distribution of \eqn{T^{(A)} \mid S^{(A)}, T^{(A)} > t_0} by using a working semiparametric  model such as the landmark proportional hazards model \eqn{q_A(S) = P(T^{(A)} >t \mid T^{(A)} >t_0, S^{(A)}) = \exp\{-\Lambda_0(t|t_0) \exp(\beta ^{T} S^{(A)}) \},  t>t_0, }
where \eqn{ \Lambda_0(t|t_0)} is the unspecified baseline cumulative hazard function for \eqn{T^{(A)}} conditional on \eqn{\{T^{(A)} > t_0\}} and \eqn{\beta} is an unknown vector of coefficients. Let \eqn{\hat{\beta}} be the maximizer of the corresponding log partial likelihood function and \eqn{\hat{\Lambda}_0(t|t_0)} be the Breslow-type estimate of baseline hazard. If one were to assume that this working model is correctly specified, then a consistent estimate of \eqn{\Delta_{S}(t,t_0)} would simply be: \eqn{\hat{\Delta}_{S}^M=n_B^{-1} \sum_{i=1}^{n_B} [ \exp \{ -\hat{\Lambda}_0(t|t_0)\exp(\hat{\beta} ^{T} S_{Bi}) \} \frac{I(X_{Bi} > t_0)}{\hat{W}^C_B(t_0)}  -   \frac{I(X_{Bi} > t)}{\hat{W}^C_B(t)}].} We refer to this estimate as the two-stage model-based estimator (option 4 for type).
Instead of relying on correct specification of this model, we use the resulting  score \eqn{U = \beta_0^{T}S} as a univariate ``pseudo-marker" to summarize the \eqn{k} surrogates. In the second stage, to estimate \eqn{\Delta_{S}(t,t_0)}, we apply a nonparametric approach with \eqn{S} represented by the univariate marker \eqn{U}.  Specifically, we use a  kernel Nelson-Aalen estimator to nonparametrically estimate  \eqn{\phi_A(t|u, t_0)=P(T^{(A)}> t\mid U^{(A)}=u, T^{(A)}> t_0) = \exp\{-\Lambda_A(t|u, t_0 )\}} as \eqn{\hat \phi_A(t|u, t_0)  = \exp\{-\hat{\Lambda}_A(t|u, t_0) \}}, where
\eqn{\hat{\Lambda}_A(t|u, t_0) = \int_{t_0}^t \frac{\sum_{i=1}^{n_A} I(X_{Ai}>t_0) K_h\{\gamma(\hat{U}_{Ai}) - \gamma(u)\}dN_{Ai}(z)}{\sum_{i=1}^{n_A}  I(X_{Ai}>t_0) K_h\{\gamma(\hat{U}_{Ai}) - \gamma(u)\} Y_{Ai}(z)},}
\eqn{\hat{U}_{Ai} = \hat{\beta} ^{T} S_{Ai}}, \eqn{\hat{U}_{Bi} = \hat{\beta} ^{T} S_{Bi}}, \eqn{Y_{Ai} = I(X_{Ai} \geq t)}, \eqn{N_{Ai}(t) = I(X_{Ai} \leq t) \delta_{Ai},  K(\cdot)} is a smooth symmetric density function, \eqn{K_h(x) = K(x/h)/h,} and \eqn{\gamma(\cdot)} is a given monotone transformation function. We then estimate \eqn{\Delta_{S}(t,t_0)} as
\eqn{\hat{\Delta}_{S}(t,t_0) = n_B^{-1} \sum_{i=1}^{n_B} [\hat{\phi}_A(t|\hat{U}_{Bi},t_0) \frac{I(X_{Bi} > t_0)}{\hat{W}^C_B(t_0)}  -   \frac{I(X_{Bi} > t)}{\hat{W}^C_B(t)}]} and \eqn{\hat{R}_{S}(t,t_0) =1- \hat{\Delta}_{S}(t,t_0)/\hat{\Delta}(t).} We refer to this estimate as the two-stage robust estimator (option 1 for type).

The next estimator borrows ideas from the extensive causal inference literature focusing on double robust estimators  two-stage weighted estimator with a propensity score weight explicitly balancing the two treatment groups with respect to the distribution of \eqn{S}. The weighting enables us to ``adjust" the distribution of \eqn{S^{(A)}} before constructing the conditional survival estimate \eqn{\hat \phi_A(t|u, t_0).} This approach results in a double-robust estimator of \eqn{\Delta_{S}(t, t_0)}, which is consistent when either \eqn{U^{(A)}} captures all the information about the relationship between \eqn{I(T^{(A)}\ge t)} and \eqn{S^{(A)}} or the propensity score model for \eqn{\pi(S,t_0)=P(G_i=B|S_i=S, T_{i} > t_0)} is correctly specified. While \eqn{\pi(S,t_0)} depends on \eqn{t_0}, for simplicity, we drop \eqn{t_0} from our notation and simply use \eqn{\pi(S)}.

Regression models can be imposed to obtain estimates for \eqn{\pi(S)}. For example, a simple logistic regression model can be imposed for 
\eqn{\tilde{\pi}(S)=P(G_i=B|S_i=S, X_{i}>t_0)} with
%%
\eqn{logit\{\tilde{\pi}(S)\} = \alpha_0 + \alpha_1 ^{T} S,}
where \eqn{\alpha_0} and \eqn{\alpha_1} are estimated only among those with   \eqn{X_{gi} > t_0} to account for censoring. The propensity score of interest, \eqn{\pi(S),} can be derived from \eqn{\tilde{\pi}(S)} directly since 
\eqn{logit\{\pi(S)\}=logit\{\tilde{\pi}(S)\} + \log\{W_A^C(t_0)/W_B^C(t_0)\},}
which follows from the assumption that \eqn{(T_{gi}, S_{gi}) \perp C_{gi}.} We then modify  the above expression by weighting observations with the estimated \eqn{L(S_{Ai})=\pi(S_{Ai})/\{1-\pi(S_{Ai})\}}  and obtain
\deqn{\hat{\Lambda}^w_A(t|u, t_0)  = \int_{t_0}^t \frac{\sum_{i=1}^{n_A} \hat{L}(S_{Ai}) I(X_{Ai}>t_0) K_h\{\gamma(\hat{U}_{Ai}) - \gamma(u)\}dN_{Ai}(z)}{\sum_{i=1}^{n_A}  \hat{L}(S_{Ai}) I(X_{Ai}>t_0) K_h\{\gamma(\hat{U}_{Ai}) - \gamma(u)\} Y_{Ai}(z)},},
 where \eqn{\hat{L}(S_{gi}) = \exp(\hat{\alpha}_0+\hat{\alpha}_1^{T} S_{gi})\hat{W}^C_B(t_0)/\hat{W}^C_A(t_0).}

Subsequently, we define
\eqn{\hat{\Delta}^w_{S}(t,t_0) = n_B^{-1} \sum_{i=1}^{n_B} [\hat{\phi}^w_A(t|\hat{U}_{Bi},t_0) \frac{I(X_{Bi} > t_0)}{\hat{W}^C_B(t_0)}  -   \frac{I(X_{Bi} > t)}{\hat{W}^C_B(t)}]} and \eqn{\hat{R}^w_{S}(t,t_0) =1- \hat{\Delta}^w_{S}(t,t_0)/\hat{\Delta}(t)}
where \eqn{\hat \phi^w_A(t|t_0,u)  = \exp\{-\hat{\Lambda}^w_A(t|t_0,u) \}.} We refer to this estimate as the weighted two-stage robust estimator (option 2 for type).

While the two-stage weighted estimator reflects one way to enhance the robustness of an initial estimator, the idea of combining a propensity-score type model and a regression-type model has certainly been extensively studied in the causal inference literature and a more familiar double-robust estimator can be constructed as: \eqn{\hat{\Delta}_{S}^{DR}(t,t_0)  = n^{-1} [\sum_{i=1}^{n_A}\frac{I(X_{Ai}>t)}{\hat{W}_{A}^C(t)\bar{\pi}_B} \hat{L}(S_{Ai}) - \sum_{i=1}^{n_B} \frac{I(X_{Bi}>t)}{\hat{W}_{B}^C(t)\bar{\pi}_B}  ]  - n^{-1} [ \sum_{i=1}^{n_A} \frac{ \hat{\phi}_A(t\mid\hat{U}_{Ai},t_0) I(X_{Ai} > t_0) }{\hat{W}_{A}^C(t_0)\bar{\pi}_B}  \hat{L}(S_{Ai}) - \sum_{i=1}^{n_B}\frac{ \hat{\phi}_A(t\mid \hat{U}_{Bi},t_0) I(X_{Bi} > t_0) }{\hat{W}_{B}^C(t_0)\bar{\pi}_B} ] }
and \eqn{\hat{R}_{S}^{DR}(t,t_0) =1- \hat{\Delta}_{S}^{DR}(t,t_0)/\hat{\Delta}(t)}, where  \eqn{\hat{\phi}_A(t|\hat{U}_{gi},t_0)} is the (unweighted) estimate of \eqn{\phi_A(t\mid u,t_0)} used in \eqn{\hat{\Delta}_{S}(t,t_0)}. We refer to this estimate as the double robust estimator (option 3 for type).

The weighted estimator (option type 5) is defined as: 
\eqn{ \hat{\Delta}^{PS}_{S}(t,t_0) =  n^{-1} \sum_{i=1}^n  \{\frac{I(X_{i}>t)}{\hat{W}_{G_i}^C(t)\bar{\pi}_B}  [ I(G_i = A)\hat{L}(S_{Ai})   - I(G_i = B )  ]  \} }
and \eqn{\hat{R}^{PS}_{S}(t,t_0) =1- \hat{\Delta}^{PS}_{S}(t,t_0)/\hat{\Delta}(t).} This estimator  completely relies on the correct specification of \eqn{\pi(S)}. The double-robust model-based estimator (option 6 for type) is defined as \eqn{\hat{\Delta}_{S}^{DR2}(t,t_0)}  and \eqn{\hat{R}_{S}^{DR2}(t,t_0) =1- \hat{\Delta}_{S}^{DR2}(t,t_0)/\hat{\Delta}(t)} which are constructed parallel to the construction of  \eqn{\hat{\Delta}_{S}^{DR}(t,t_0)} i.e., a combination of  \eqn{\hat{\Delta}_{S}^M(t,t_0)} and \eqn{\hat{R}^{PS}_{S}(t,t_0)}.
 
Variance estimates are obtained using perturbation resampling. If a confidence interval is requested three versions are provided: a normal approximation based interval, a quantile based interval and Fieller's confidence interval, all using perturbation-resampling. An estimate of the incremental value of the surrogate marker information can also be requested; this essentially compared the proportion explained by the surrogate information vs. the proportion explained by \eqn{T} alone up to \eqn{t_0}. Details can be found in Parast, L., Cai, T., & Tian, L. (2021). Evaluating multiple surrogate markers with censored data. Biometrics, 77(4), 1315-1327.
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

Parast, L., Cai, T., & Tian, L. (2021). Evaluating multiple surrogate markers with censored data. Biometrics, 77(4), 1315-1327.
}
\author{
%%  ~~who you are~~
Layla Parast
}
\note{If the treatment effect is not significant, the user will receive the following message: "Warning: it looks like the treatment effect is not significant; may be difficult to interpret the residual treatment effect in this setting".}

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