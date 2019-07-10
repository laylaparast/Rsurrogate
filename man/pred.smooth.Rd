\name{pred.smooth}
\alias{pred.smooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates expected outcome for control group surrogate values
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Helper function; calculates the expected outcome for control group surrogate values using treatment group surrogate and outcome information.
}
\usage{
pred.smooth(zz, zi.one, bw = NULL, y1, weight = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zz}{
%%     ~~Describe \code{zz} here~~
surrogate values in the control group
}
  \item{zi.one}{
%%     ~~Describe \code{zi.one} here~~
surrogate values in the treatment group; used to estimate conditional mean function
}
  \item{bw}{
%%     ~~Describe \code{bw} here~~
bandwidth
}
  \item{y1}{
%%     ~~Describe \code{y1} here~~
outcome in the treatment group
}
  \item{weight}{
%%     ~~Describe \code{weight} here~~
weight
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
Details are included in the documentation for R.s.estimate

}


\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
expected outcome for each surrogate value in the control group
}
\references{
%% ~put references to the literature/web site here ~
Parast, L., McDermott, M., Tian, L. (2016). Robust estimation of the proportion of treatment effect explained by surrogate marker information. Statistics in medicine, 35(10):1637-1653.
}
\author{
%%  ~~who you are~~
Layla Parast
}
\keyword{smooth}
\keyword{robust}
\keyword{nonparametric}
\keyword{internal}
% __ONLY ONE__ keyword per line
