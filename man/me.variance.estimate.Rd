\name{me.variance.estimate}
\alias{me.variance.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates measurement error variance given replicate data. 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Estimates measurement error variance given replicate data using a simple components of variance analysis. 
}
\usage{
me.variance.estimate(replicates)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{replicates}{
%%     ~~Describe \code{replicates} here~~
matrix of data where each row indicates a subject and each column is a replicated measurement; columns can have NAs when subjects have different numbers of measurements. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Estimates measurement error variance given replicate data using a simple components of variance analysis. 
}
\value{
estimate of measurement error variance
}
\references{
%% ~put references to the literature/web site here ~
Carroll, R. J., Ruppert, D., Crainiceanu, C. M., and Stefanski, L. A. (2006). Measurementerror in nonlinear models: a modern perspective. Chapman and Hall/CRC.

Parast, L., Garcia, T. P., Prentice, R. L., & Carroll, R. J. (2022). Robust methods to correct for measurement error when evaluating a surrogate marker. Biometrics, 78(1), 9-23.
}
\author{
%%  ~~who you are~~
Layla Parast
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}% use one of  RShowDoc("KEYWORDS")
\keyword{htest}% __ONLY ONE__ keyword per line
