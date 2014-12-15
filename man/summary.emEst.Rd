% File src/library/base/man/summary.emEst.Rd
\name{summary.emEst}
\alias{summary.emEst}
\title{Summarize output from EM algorithm for structural equation models}
\description{
    Summarize data from object obtained from \code{\link{em}}.
}
\usage{
\method{summary}{emEst}(object, ...)
}
\arguments{
    \item{object}{estimated structural equation model of class `emEst'
    obtained from \code{\link{em}}.}
    \item{...}{additional arguments.}
}
\value{
    Returns a list that consists of the following components:
    \item{estimates}{table of estimated parameters with standard errors and
    t and p values.}
    \item{iterations}{iterations needed by EM algorithm till convergence.}
    \item{finallogLik}{final loglikelihood obtained by EM algorithm.}
    \item{loglikelihoods}{table of loglikelohoods for each iteration of EM
    algorithm with difference and relative change.}
}
\examples{
model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1, xi="x1-x3,x4-x6", eta="y1-y3")

# starting parameters
pars.start <- c(.6, .5, .4, .5, .4, .6, .5, .2, .6, .7, .3, .2, .5, .7, .3, .4,
                .6, .2, .3, .4, .3, .3, .4, .5, .2, .3, .3, .4, .5, .6, .8)

data(example_data)

\dontrun{
# estimate model
res <- em(model, pars.start, example_data)
summary(res)
}
}

