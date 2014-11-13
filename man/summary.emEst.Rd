% File src/library/base/man/summary.emEst.Rd
\name{summary.emEst}
\alias{summary.emEst}
\title{Summarize output from EM algorithm for structural equation model}
\description{
    Summarize data from object created with \code{\link{em}}.
}
\usage{
\method{summary}{emEst}(object, ...)
}
\arguments{
    \item{object}{estimated structural equation model of class 'emEst'
    obtained with \code{\link{em}}.}
    \item{...}{additional arguments.}
}
\examples{
\dontrun{
model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1, xi="x1-x3,x4-x6", eta="y1-y3")

# starting parameters
pars.start <- c(.6, .5, .4, .5, .4, .6, .5, .2, .6, .7, .3, .2, .5, .7, .3, .4,
                .6, .2, .3, .4, .3, .3, .4, .5, .2, .3, .3, .4, .5, .6, .8)

data(example_data)

# estimate model
res <- em(model, pars.start, example_data)
summary(res)
}
}

