\name{calculate.var.np}
\alias{calculate.var.np}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
calculates closed form variance estimate for R
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
calculates closed form variance estimate for R; used in R function that corrects for measurement error
}
\usage{
calculate.var.np(s1, s0, y1, y0, extrapolate = TRUE)
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
  \item{extrapolate}{
%%     ~~Describe \code{extrapolate} here~~
TRUE or FALSE; indicates whether the user wants to use extrapolation.

}
}

\value{
  \item{total}{matrix needed for variance calculation}
  \item{psionly}{matrix needed for variance calculation}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Parast, L., McDermott, M., Tian, L. (2016). Robust estimation of the proportion of treatment effect explained by surrogate marker information. Statistics in Medicine, 35(10):1637-1653.

Parast, L., Garcia, TP, Prentice, RL, Carroll, RJ (2021). Robust Methods to Correct for Measurement Error when Evaluating a Surrogate Marker. Biometrics, In press.

}
\author{
%%  ~~who you are~~
Layla Parast
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}% use one of  RShowDoc("KEYWORDS")
\keyword{htest}% __ONLY ONE__ keyword per line
