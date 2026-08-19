\name{select.surrogate}
\alias{select.surrogate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Selects an optimal subset of surrogate markers and estimates the proportion of treatment effect explained
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function estimates the proportion of treatment effect on the primary outcome explained by the treatment effect on a set of surrogate markers, optionally first selecting a subset of the candidate markers using treatment-aware regularization (TAR) or an unweighted regularized approach. The user can also request a bootstrap variance estimate and 95\% confidence intervals.
}
\usage{
select.surrogate(sone, szero, yone, yzero, selection = "none", var = FALSE,
boot.num = 200, nlambda = 20, folds = 5, epsilon = 0.05, c.constant = 1,
small.n = 1e-6, lambda_rule = "lambda.1se")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sone}{
%%     ~~Describe \code{sone} here~~
numeric matrix or data frame; candidate surrogate markers for treated observations, with \eqn{n_1} (number of treated observations) rows and \eqn{k} (number of candidate surrogate markers) columns. If the columns are not named, they are labeled s1, ..., sk.
}
  \item{szero}{
%%     ~~Describe \code{szero} here~~
numeric matrix or data frame; candidate surrogate markers for control observations, with \eqn{n_0} (number of control observations) rows and \eqn{k} columns, in the same order as sone.
}
  \item{yone}{
%%     ~~Describe \code{yone} here~~
numeric vector; primary outcome for treated observations.
}
  \item{yzero}{
%%     ~~Describe \code{yzero} here~~
numeric vector; primary outcome for control observations.
}
  \item{selection}{
%%     ~~Describe \code{selection} here~~
specifies how the subset of surrogate markers is chosen; choices are "none" (no selection; the proportion of treatment effect explained is estimated using all candidate markers, as in R.s.estimate with number = "multiple"), "TAR" (treatment-aware regularization), or "rr" (an unweighted, regularized approach), default is "none".
}
  \item{var}{
%%     ~~Describe \code{var} here~~
TRUE or FALSE; indicates whether a bootstrap variance estimate and 95\% confidence intervals are requested, default is FALSE.
}
  \item{boot.num}{
%%     ~~Describe \code{boot.num} here~~
number of bootstrap resamples used when var = TRUE; on each resample, both treatment groups are resampled with replacement and the entire selection procedure, including re-selection of the penalty, is rerun, default is 200.
}
  \item{nlambda}{
%%     ~~Describe \code{nlambda} here~~
the number of penalty values considered along the lasso path used by selection = "TAR" or selection = "rr", default is 20.
}
  \item{folds}{
%%     ~~Describe \code{folds} here~~
the number of cross-validation folds, within the treated group, used by selection = "TAR" to select the penalty, default is 5. Not used if selection is "none" or "rr".
}
  \item{epsilon}{
%%     ~~Describe \code{epsilon} here~~
the tolerance used by selection = "TAR": the selected subset is the sparsest point along the weighted lasso path whose estimated proportion of treatment effect explained falls no more than epsilon below that of the full candidate set, default is 0.05. Not used if selection is "none" or "rr".
}
  \item{c.constant}{
%%     ~~Describe \code{c.constant} here~~
a constant controlling the strength of the treatment-aware penalty weighting in selection = "TAR", default is 1. Not used if selection is "none" or "rr".
}
  \item{small.n}{
%%     ~~Describe \code{small.n} here~~
a small positive constant added to the estimated treatment effect magnitude of each candidate marker before inverting it to form the treatment-aware penalty weights in selection = "TAR", to avoid dividing by a value at or near zero, default is 1e-6. Not used if selection is "none" or "rr".
}
  \item{lambda_rule}{
%%     ~~Describe \code{lambda_rule} here~~
the cross-validation rule used to select the penalty in selection = "rr"; choices are "lambda.1se" or "lambda.min", default is "lambda.1se". Not used if selection is "none" or "TAR".
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Let \eqn{Y^{(1)}} and \eqn{Y^{(0)}} denote the primary outcome under treatment and control, respectively, and let \eqn{\bold{S}^{(1)}} and \eqn{\bold{S}^{(0)}} denote the vector of \eqn{k} candidate surrogate markers under treatment and control, respectively. As in R.s.estimate, the proportion of treatment effect explained by a set of markers \eqn{\bold{S}} is \eqn{R_{\bold{S}}=1-\Delta_{\bold{S}}/\Delta}, where \eqn{\Delta} is the overall treatment effect and \eqn{\Delta_{\bold{S}}} is the residual treatment effect after accounting for \eqn{\bold{S}}, estimated using a working model to reduce the dimension of \eqn{\bold{S}} followed by kernel smoothing, exactly as described in the R.s.estimate documentation for number = "multiple", type = "robust".

When selection = "none", all \eqn{k} candidate markers are used, and the estimate of \eqn{R_{\bold{S}}} agrees with that of R.s.estimate run with number = "multiple".

When selection = "TAR" or selection = "rr", a subset \eqn{\widetilde{\bold{S}}} of the \eqn{k} candidate markers is first selected via a lasso fit to the treated group, and \eqn{R_{\widetilde{\bold{S}}}} is then estimated for the selected subset only, using a debiased (unpenalized) refit on the selected markers followed by the same kernel-smoothed approach. For selection = "rr", the lasso is unweighted and the penalty is chosen by ordinary cross-validation using lambda_rule. For selection = "TAR", each candidate marker's lasso penalty is instead weighted by the inverse of its estimated standardized treatment effect, so that markers with little or no treatment effect are penalized more heavily, and the penalty along the resulting weighted lasso path is selected via cross-validation, within the treated group, to be the sparsest point whose estimated proportion of treatment effect explained falls no more than epsilon below that of the full candidate set.

If a bootstrap variance estimate and confidence intervals are requested (var = TRUE), both treatment groups are resampled with replacement and the entire estimation procedure specified by selection is rerun on each resample, including re-selection of the penalty for selection = "TAR" or selection = "rr". This differs from the perturbation-resampling used elsewhere in this package. A 95\% confidence interval is provided for \eqn{\hat{\Delta}}, \eqn{\hat{\Delta}_{\bold{S}}}, and \eqn{\widehat{R}_{\bold{S}}} (or \eqn{\widehat{R}_{\widetilde{\bold{S}}}}), based on the 2.5th and 97.5th percentiles of the bootstrap replicates.

If the observed supports for the candidate markers (or the selected subset) are not the same across treatment groups, values outside the observed support are extrapolated to the closest non-missing value; unlike R.s.estimate, this is done automatically and cannot be turned off.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list is returned:
 \item{delta }{the estimate, \eqn{\hat{\Delta}}, described above.}
 \item{delta.s }{the estimate, \eqn{\hat{\Delta}_{\bold{S}}} (or \eqn{\hat{\Delta}_{\widetilde{\bold{S}}}} if selection is "TAR" or "rr"), described above.}
 \item{R.s }{the estimate, \eqn{\hat{R}_{\bold{S}}} (or \eqn{\hat{R}_{\widetilde{\bold{S}}}} if selection is "TAR" or "rr"), described above.}
 \item{selection }{a named numeric vector of the (penalized) lasso coefficients of the selected markers; only returned if selection is "TAR" or "rr".}
 \item{delta.var }{the bootstrap variance estimate of \eqn{\hat{\Delta}}; if var = TRUE.}
 \item{delta.s.var }{the bootstrap variance estimate of \eqn{\hat{\Delta}_{\bold{S}}}; if var = TRUE.}
 \item{R.s.var }{the bootstrap variance estimate of \eqn{\hat{R}_{\bold{S}}}; if var = TRUE.}
 \item{conf.int.delta}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}}, based on the 2.5th and 97.5th percentiles of the bootstrap replicates; if var = TRUE.}
 \item{conf.int.delta.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{\Delta}_{\bold{S}}}, based on the 2.5th and 97.5th percentiles of the bootstrap replicates; if var = TRUE.}
 \item{conf.int.R.s}{a vector of size 2; the 95\% confidence interval for \eqn{\hat{R}_{\bold{S}}}, based on the 2.5th and 97.5th percentiles of the bootstrap replicates; if var = TRUE.}
}
\references{
%% ~put references to the literature/web site here ~
Parast, L., Freestone, G., Tarr, J. (2026). Treatment-Aware Regularization via Multi-Objective Optimization for Surrogate Selection.

Parast, L., McDermott, M., Tian, L. (2016). Robust estimation of the proportion of treatment effect explained by surrogate marker information. Statistics in Medicine, 35(10):1637-1653.
}
\author{
%%  ~~who you are~~
Layla Parast
}
\examples{
data(d_tar_example)
names(d_tar_example)
select.surrogate(sone = d_tar_example$sone, szero = d_tar_example$szero,
yone = d_tar_example$yone, yzero = d_tar_example$yzero, selection = "none")
select.surrogate(sone = d_tar_example$sone, szero = d_tar_example$szero,
yone = d_tar_example$yone, yzero = d_tar_example$yzero, selection = "TAR")
\dontrun{
select.surrogate(sone = d_tar_example$sone, szero = d_tar_example$szero,
yone = d_tar_example$yone, yzero = d_tar_example$yzero, selection = "TAR",
var = TRUE, boot.num = 200)
}
}

\keyword{regression}
\keyword{models}
\keyword{smooth}
\keyword{robust}
\keyword{nonparametric}
