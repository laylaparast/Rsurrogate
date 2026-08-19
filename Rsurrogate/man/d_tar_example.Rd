\name{d_tar_example}
\alias{d_tar_example}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Example data for treatment-aware regularization
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
One simulated dataset, generated from Setting 1 of the simulation study of Parast, Freestone, and Tarr (2026), to be used in the examples for select.surrogate. Nine correlated candidate surrogate markers are generated for each observation; the first three (s1, s2, s3) both shift between treatment groups and drive the outcome, while the remaining six are not associated with the treatment effect on the outcome.
}
\usage{data(d_tar_example)}
\format{
 A list with 4 elements representing 500 observations from a treatment group and 500 observations from a control group, each with 9 candidate surrogate markers:
  \describe{
    \item{\code{sone}}{a 500 by 9 numeric matrix of candidate surrogate markers for the treated observations, with columns named s1, ..., s9.}
    \item{\code{szero}}{a 500 by 9 numeric matrix of candidate surrogate markers for the control observations, with columns named s1, ..., s9.}
    \item{\code{yone}}{a numeric vector of length 500; the primary outcome for the treated observations.}
    \item{\code{yzero}}{a numeric vector of length 500; the primary outcome for the control observations.}
}
}
\references{
Parast, L., Freestone, G., Tarr, J. (2026). Treatment-Aware Regularization via Multi-Objective Optimization for Surrogate Selection.
}
\examples{
data(d_tar_example)
names(d_tar_example)
}
\keyword{datasets}
