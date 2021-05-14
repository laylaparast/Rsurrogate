\name{d_example}
\alias{d_example}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Hypothetical data
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Hypothetical data to be used in examples.
}
\usage{data(d_example)}
\format{
 A list with 8 elements representing 500 observations from a control group and 500 observations from a treatment group:
  \describe{
    \item{\code{s1.a}}{First surrogate marker measurement for treated observations.}
    \item{\code{s1.b}}{Second surrogate marker measurement for treated observations.}
    \item{\code{s1.c}}{Third surrogate marker measurement for treated observations.}
    \item{\code{y1}}{Primary outcome for treated observations.}
    \item{\code{s0.a}}{First surrogate marker measurement for control observations.}
    \item{\code{s0.b}}{Second surrogate marker measurement for control observations.}
    \item{\code{s0.c}}{Third surrogate marker measurement for control observations.}
    \item{\code{y0}}{Primary outcome for control observations.}
}
}
\examples{
data(d_example)
names(d_example)
}
\keyword{datasets}
