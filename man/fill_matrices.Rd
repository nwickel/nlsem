% File src/library/base/man/fill_matrices.Rd
\name{fill_matrices}
\alias{fill_matrices}
\title{Created model matrices for structural equation model.}
\description{
    Create model matrices with constraints from data frame that specifies
    structural euation model.
}
\usage{
    fill_matrices(dat)
}
\arguments{
    \item{dat}{Data frame with variables \code{label} and \code{ustart}
    usually created with \code{\link{specify_sem}}.}
}
\examples{
specs <- specify_sem(num.x=4, num.y=1, num.xi=2, num.eta=1, xi="x1-x2,x3-x4", eta="y1")
fill_matrices(specs)
}

