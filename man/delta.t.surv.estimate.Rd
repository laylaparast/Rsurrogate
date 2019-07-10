\name{delta.t.surv.estimate}
\alias{delta.t.surv.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates robust residual treatment effect accounting only for primary outcome information up to a specified time

}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function calculates the robust estimate of the residual treatment effect accounting only for primary outcome information up to \eqn{t_0} i.e. the hypothetical treatment effect if survival up to \eqn{t_0} in the treatment group looks like survival up to \eqn{t_0} in the control group. Ideally this function is only used as a helper function and is not directly called. 

}
\usage{
delta.t.surv.estimate(xone, xzero, deltaone, deltazero, t, weight.perturb = NULL,
landmark)
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
}
\details{
	Details are included in the documentation for R.t.surv.estimate.
}
\references{
%% ~put references to the literature/web site here ~
Parast L, Cai T and Tian L. Evaluating Surrogate Marker Information using Censored Data. Under Review.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\eqn{\hat{\Delta}_T(t,t_0)}, the robust residual treatment effect estimate accounting only for survival up to \eqn{t_0}.
}
\author{
%%  ~~who you are~~
Layla Parast
}

\note{If the treatment effect is not significant, the user will receive the following message: "Warning: it looks like the treatment effect is not significant; may be difficult to interpret the residual treatment effect in this setting". If the treatment effect is negative, the user will receive the following message: "Warning: it looks like you need to switch the treatment groups" as this package assumes throughout that larger values of the event time are better.}



\examples{
data(d_example_surv)
names(d_example_surv)
% Functions have elapsed time longer than 10s
%delta.t.surv.estimate(xone = d_example_surv$x1, xzero = d_example_surv$x0,  
%deltaone = d_example_surv$delta1, deltazero = d_example_surv$delta0, t = 3, 
%landmark = 1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
\keyword{survival}
\keyword{robust}
