\name{delta.multiple.surv}
\alias{delta.multiple.surv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates robust residual treatment effect accounting for multiple surrogate markers at a specified time and primary outcome information up to that specified time
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function calculates the robust estimate of the residual treatment effect accounting for multiple surrogate markers measured at \eqn{t_0} and primary outcome information up to \eqn{t_0} i.e. the hypothetical treatment effect if both the surrogate marker distributions at \eqn{t_0} and survival up to \eqn{t_0} in the treatment group look like the surrogate marker distributions and survival up to \eqn{t_0} in the control group. Ideally this function is only used as a helper function and is not directly called. 

}
\usage{
delta.multiple.surv(xone, xzero, deltaone, deltazero, sone, szero, type =1, t, 
weight.perturb = NULL, landmark, extrapolate = FALSE, transform = FALSE,
approx = T)
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
TRUE or FALSE; indicates whether the user wants to use a transformation for the surrogate marker psuedo-score. 
}
  \item{approx}{
TRUE or FALSE indicating whether an approximation should be used when calculating the probability of censoring; most relevant in settings where the survival time of interest for the primary outcome is greater than the last observed event but before the last censored case, default is TRUE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
	Details are included in the documentation for R.multiple.surv.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\eqn{\hat{\Delta}_S(t,t_0)}, the residual treatment effect estimate accounting for multiple surrogarte markers measured at \eqn{t_0} and primary outcome information up to \eqn{t_0}.
}
\references{
%% ~put references to the literature/web site here ~
Details can be found in Parast, L., Cai, T., Tian, L. (2021). Evaluating Multiple Surrogate Markers with Censored Data. Biometrics, In press.
}

\note{If the treatment effect is not significant, the user will receive the following message: "Warning: it looks like the treatment effect is not significant; may be difficult to interpret the residual treatment effect in this setting".}




\author{
%%  ~~who you are~~
Layla Parast
}
\examples{
data(d_example_multiple)
names(d_example_multiple)
\dontrun{
delta.multiple.surv(xone = d_example_multiple$x1, xzero = d_example_multiple$x0, deltaone =
 d_example_multiple$delta1, deltazero = d_example_multiple$delta0, sone = 
 as.matrix(d_example_multiple$s1), szero = as.matrix(d_example_multiple$s0), 
 type =1, t = 1, landmark=0.5)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
\keyword{survival}
\keyword{robust}
\keyword{smooth}