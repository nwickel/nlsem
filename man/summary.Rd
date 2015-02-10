% File src/library/base/man/summary.Rd
\name{summary}
\alias{summary.emEst}
\title{Summarize output from EM algorithm for structural equation models}
\description{
    Summarize data from object obtained from \code{\link{em}}.
}
\usage{
\method{summary}{emEst}(object, ...)
}
\arguments{
    \item{object}{estimated structural equation model of class \code{emEst}
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

