\name{VTM}
\alias{VTM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Repeats a row.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Helper function; this function creates a matrix that repeats vc, dm times where each row is equal to the vc vector. 
}
\usage{
VTM(vc, dm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vc}{
%%     ~~Describe \code{vc} here~~
the vector to repeat.
}
  \item{dm}{
%%     ~~Describe \code{dm} here~~
number of rows.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
a matrix that repeats vc, dm times where each row is equal to the vc vector
}
\keyword{array}
\keyword{internal}