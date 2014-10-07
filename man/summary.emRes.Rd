% File src/library/base/man/summary.emRes.Rd
\name{summary.emRes}
\alias{summary.emRes}
\title{Summarize data from LMS model.}
\description{
    Summarize data from object created from \code{\link{em}}.
}
\usage{
\method{summary}{emRes}(object, ...)
}
\arguments{
    \item{object}{Estimated structural equation model of class "emRes"
    obtained with \code{\link{em}}.}
    \item{...}{Additional arguments.}
}
\examples{
\dontrun{
specs <- specify_model(num.x=6, num.y=3, num.xi=2, num.eta=1, xi="x1-x3,x4-x6", eta="y1-y3")
model <- fill_matrices(specs1)

# starting parameters
pars.start <- c(.6, .5, .4, .5, .4, .6, .5, .2, .6, .7, .3, .2, .5, .7, .3, .4,
                .6, .2, .3, .4, .3, .3, .4, .5, .2, .3, .3, .4, .5, .6, .8)

data(example_data)

# estimate model
res <- em(model, pars.start, example_data)
summary(res)
}
}

