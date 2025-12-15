\name{d_example_me}
\alias{d_example_me}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Hypothetical data with replicate measurements
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Hypothetical data to be used in measurement error example.
}
\usage{data(d_example_me)}
\format{
  A list with 10 elements representing 500 observations from a control group and 500 observations from a treatment group:
  \describe{
    \item{\code{y1}}{Primary outcome for treated observations.}
    \item{\code{s1}}{Surrogate marker for treated observations.}
    \item{\code{s1_rep1}}{Replicate measurement of the surrogate marker for treated observations.}
    \item{\code{s1_rep2}}{Replicate measurement of the surrogate marker for treated observations.}
    \item{\code{s1_rep3}}{Replicate measurement of the surrogate marker for treated observations.}
    \item{\code{y0}}{Primary outcome for control observations.}
    \item{\code{s0}}{Surrogate marker for control observations.}
    \item{\code{s0_rep1}}{Replicate measurement of the surrogate marker for control observations.}
    \item{\code{s0_rep2}}{Replicate measurement of the surrogate marker for control observations.}
    \item{\code{s0_rep3}}{Replicate measurement of the surrogate marker for control observations.}
  }
}
\examples{
data(d_example_me)
names(d_example_me)
}
\keyword{datasets}
