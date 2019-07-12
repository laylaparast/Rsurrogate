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
options are "d","q","n" for parametric and "q","n" for nonparametric; "d" stands for the disattenuated estimator, "q" stands for the SIMEX estimator with quadratic extrapolation, "n" stands for the SIMEX estimator with a nonlinear extrapolation.
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
%%  ~~ If necessary, more details than the description above ~~
Details can be found in Parast, L., Garcia, TP, Prentice, RL, Carroll, RJ (2019+). Robust Methods to Correct for Measurement Error when Evaluating a Surrogate Marker. Under Review. 

Please email parast@rand.org if you would like a copy of this article. 
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
Parast, L., Garcia, TP, Prentice, RL, Carroll, RJ (2019+). Robust Methods to Correct for Measurement Error when Evaluating a Surrogate Marker. Under Review.
}
\author{
%%  ~~who you are~~
Layla Parast
}

\examples{
data(d_example)
names(d_example)
R.s.estimate.me(yone=d_example$y1, yzero=d_example$y0, sone=d_example$s1.a, szero=d_example$s0.a, 
parametric = TRUE, estimator = "d", me.variance = 0.5, naive= TRUE, Ronly = FALSE)
R.s.estimate.me(yone=d_example$y1, yzero=d_example$y0, sone=d_example$s1.a, szero=d_example$s0.a, 
parametric = TRUE, estimator = "q", me.variance = 0.5, naive= FALSE, Ronly = TRUE)
\donttest{R.s.estimate.me(yone=d_example$y1, yzero=d_example$y0, sone=d_example$s1.a, szero=d_example$s0.a, 
parametric = FALSE, estimator = "q", me.variance = 0.5, naive= FALSE, Ronly = TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
\keyword{models}
\keyword{smooth}
\keyword{robust}
\keyword{nonparametric}

