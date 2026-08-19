\name{.te.calculate}
\alias{.te.calculate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates standardized treatment effect magnitude for one candidate marker
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Helper function used by select.surrogate with selection = "TAR"; this calculates the absolute standardized difference in means between the treated and control groups for one candidate surrogate marker, used to construct the treatment-aware penalty weights.
}
\usage{
.te.calculate(s, treat, control)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{
%%     ~~Describe \code{s} here~~
character; the name of the candidate marker column in treat and control.
}
  \item{treat}{
%%     ~~Describe \code{treat} here~~
data frame of treated observations, containing a column named s.
}
  \item{control}{
%%     ~~Describe \code{control} here~~
data frame of control observations, containing a column named s.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
the absolute value of the difference in means between treat and control for marker s, divided by the pooled standard deviation.
}
\references{
Parast, L., Freestone, G., Tarr, J. (2026). Treatment-Aware Regularization via Multi-Objective Optimization for Surrogate Selection.
}
\author{
%%  ~~who you are~~
Layla Parast
}
\keyword{regression}
\keyword{robust}
\keyword{internal}
% __ONLY ONE__ keyword per line
