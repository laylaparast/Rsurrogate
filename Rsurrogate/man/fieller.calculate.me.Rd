\name{fieller.calculate.me}
\alias{fieller.calculate.me}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Constructs Fieller's confidence interval.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Constructs Fieller's confidence interval - closed form, not using perturbation.
}
\usage{
fieller.calculate.me(a, B, var.mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
value of numerator
}
  \item{B}{
value of denominator
}
  \item{var.mat}{
variance-covariance matrix for a and B
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Produces a confidence interval for 1-a/B when closed form variance and covariance estimates are available. See documention for R.s.estimate.me for more detail.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a vector of length 2, lower bound of the 95\% confidence interval and upper bound of the 95\% confidence interval.
}
\references{
%% ~put references to the literature/web site here ~
Parast, L., Garcia, T. P., Prentice, R. L., & Carroll, R. J. (2022). Robust methods to correct for measurement error when evaluating a surrogate marker. Biometrics, 78(1), 9-23.

Fieller, Edgar C. (1954). Some problems in interval estimation. Journal of the Royal Statistical Society. Series B (Methodological), 175-185.

Fieller, E. C. (1940). The biological standardization of insulin. Supplement to the Journal of the Royal Statistical Society, 1-64.

}
\author{
Layla Parast
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}% 
