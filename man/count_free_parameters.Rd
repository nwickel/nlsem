% File src/library/base/man/count_free_parameters.Rd
\name{count_free_parameters}
\alias{count_free_parameters}
\title{Count free parameters of structural equation model}
\description{
    Counts free parameters of a structural equation model of class 'lms',
    'stemm', or 'nsemm'.
}
\usage{
    count_free_parameters(model)
}
\arguments{
    \item{model}{A model created with \code{\link{specify_sem}}.}
}
\value{
    Returns the number of free parameters in the \code{model} (numeric).
}
\examples{
model <- specify_sem(num.x=4, num.y=2, num.xi=2, num.eta=1, xi="x1-x2,x3-x4", eta="y1-y2")
count_free_parameters(model)
}

