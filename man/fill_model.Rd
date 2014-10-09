% File src/library/base/man/fill_model.Rd
\name{fill_model}
\alias{fill_model}
\title{Fill a specified structural equation model with parameters.}
\description{
    Fill a structural equation model created with
    \code{\link{specify_model}} and \code{\link{fill_matrices}} with
    parameters. Can be used for data simulations.
}
\usage{
    fill_model(model, parameters)
}
\arguments{
    \item{model}{Structural equation model specified with \code{\link{specify_model}} and \code{\link{fill_matrices}}}
    \item{parameters}{Original parameters which should be filled into model.}
}
\examples{
# specify model
specs <- specify_model(4, 1, 2, 1, xi="x1-x2,x3-x4", eta="y1")
model <- fill_matrices(specs)
# specify parameters
pars.orig <- c(0.6, 0.7,                # Lx
              0.2, 0.4,                 # G
              0.51, 0.64, 0.36, 0.51,   # Td
              0.2,                      # Te
              0.2,                      # Psi
              0.49, 0.235, 0.64,        # A
              0.7, 0.3, 0.3             # Omega
              )
# put parameters in model
mod.filled <- fill_model(model, pars.orig)
}
