% File src/library/base/man/free_parameters.Rd
\name{free_parameters}
\alias{free_parameters}
\title{Count free parameters of structural equation model}
\description{
    Counts free parameters of a structural equation model of class "lms".
}
\usage{
    free_parameters(model)
}
\arguments{
    \item{model}{A model created with \code{\link{specify_sem}} and \code{\link{fill_matrices}}.}
}
\examples{
specs <- specify_sem(num.x=4, num.y=2, num.xi=2, num.eta=1, xi="x1-x2,x3-x4", eta="y1-y2")
model <- fill_matrices(specs)
free_parameters(model)
}

