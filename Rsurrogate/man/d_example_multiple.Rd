\name{d_example_multiple}
\alias{d_example_multiple}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Hypothetical survival data with multiple surrogate markers
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Hypothetical survival data with multiple surrogate markers to be used in examples.
}
\usage{data(d_example_multiple)}
\format{
 A list with 6 elements representing 1000 observations from a control group and 1000 observations from a treatment group:
  \describe{
    \item{\code{s1}}{Surrogate marker measurements for treated observations; these markers are measured at time = 0.5. For observations that experience the primary outcome or are censored before 0.5, the surrogate values are NA.}
    \item{\code{x1}}{The observed event or censoring time for treated observations; X  = min(T, C) where T is the time of the primary outcome and C is the censoring time.}
    \item{\code{delta1}}{The indicator identifying whether the treated observation was observed to have the event or was censored; D =1*(T<C) where T is the time of the primary outcome and C is the censoring time.}
   \item{\code{s0}}{Surrogate marker measurements for control observations; these markers are measured at time = 0.5. For observations that experience the primary outcome or are censored before 0.5, the surrogate values are NA.}
    \item{\code{x0}}{The observed event or censoring time for control observations; X  = min(T, C) where T is the time of the primary outcome and C is the censoring time.}
    \item{\code{delta0}}{The indicator identifying whether the control observation was observed to have the event or was censored; D =1*(T<C) where T is the time of the primary outcome and C is the censoring time.}
}}
\examples{
data(d_example_multiple)
names(d_example_multiple)
}
\keyword{datasets}
