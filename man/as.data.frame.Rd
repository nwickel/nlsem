% File src/library/base/man/as.data.frame.Rd
\name{as.data.frame}
\alias{as.data.frame.lms}
\alias{as.data.frame.semm}
\alias{as.data.frame.nsemm}
\title{Coerce to a Data Frame}
\description{
     Function to coerce an object created with \code{\link{specify_sem}} to
     a data frame.
}
\usage{
\method{as.data.frame}{lms}(x, ...)
\method{as.data.frame}{semm}(x, ...)
\method{as.data.frame}{nsemm}(x, ...)
}
\arguments{
    \item{x}{structural equation model of class \code{lms}, \code{semm}, or
    \code{nsemm}.}
    \item{...}{additional arguments.}
}
\value{
    Returns a data frame with first column \code{label} and one column for
    each latent class labeled \code{class1}, \code{class2} and so on.}
\seealso{\code{\link{specify_sem}}, \code{\link{create_sem}}}
\examples{
# specify model
model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                     xi="x1-x3,x4-x6", eta="y1-y3", interaction="xi1:xi2")

# coerce to data frame
as.data.frame(model)
}

