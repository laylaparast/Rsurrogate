\name{d_example_surv}
\alias{d_example_surv}
\docType{data}
\title{
Hypothetical survival data
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Hypothetical survival data to be used in examples.
}
\usage{data(d_example_surv)}
\format{
 A list with 8 elements representing 500 observations from a control group and 500 observations from a treatment group:
  \describe{
    \item{\code{s1}}{Surrogate marker measurement for treated observations; this marker is measured at time = 0.5. For observations that experience the primary outcome or are censored before 0.5, this value is NA.}
    \item{\code{x1}}{The observed event or censoring time for treated observations; X  = min(T, C) where T is the time of the primary outcome and C is the censoring time.}
    \item{\code{delta1}}{The indicator identifying whether the treated observation was observed to have the event or was censored; D =1*(T<C) where T is the time of the primary outcome and C is the censoring time.}
   \item{\code{s0}}{Surrogate marker measurement for control observations; this marker is measured at time = 0.5. For observations that experience the primary outcome or are censored before 0.5, this value is NA.}
    \item{\code{x0}}{The observed event or censoring time for control observations; X  = min(T, C) where T is the time of the primary outcome and C is the censoring time.}
    \item{\code{delta0}}{The indicator identifying whether the control observation was observed to have the event or was censored; D =1*(T<C) where T is the time of the primary outcome and C is the censoring time.}
    \item{\code{z1}}{A baseline covariate value for treated observations.}
    \item{\code{z0}}{A baseline covariate value for control observations.}
}}
\examples{
data(d_example_surv)
names(d_example_surv)
}
\keyword{datasets}
