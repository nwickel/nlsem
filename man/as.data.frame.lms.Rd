% File src/library/base/man/as.data.frame.lms.Rd
\name{as.data.frame.lms}
\alias{as.data.frame.lms}
\alias{as.data.frame.semm}
\title{Coerce to a Data Frame}
\description{
     Functions to check if an object is a data frame, or coerce it if
     possible.
}
\usage{
\method{as.data.frame}{lms}(x, ...)
\method{as.data.frame}{semm}(x, ...)
}
\arguments{
    \item{x}{structural equation model of class 'lms' or 'semm'.}
    \item{...}{additional arguments.}
}
\examples{
# specify model
model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                     xi="x1-x3,x4-x6", eta="y1-y3", interaction="xi1:xi2")

# coerce to data frame
as.data.frame(model)
}

