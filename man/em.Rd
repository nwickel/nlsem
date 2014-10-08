% File src/library/base/man/lms.Rd
\name{em}
\alias{em}
\title{Maximum likelihood estimation of latent interaction effects.}
\description{
    Fits a structural equation model with latent interaction effects using LMS.
}
\usage{
    em(model, data, start, logger=FALSE, threshold=10e-7, max.iter=40, m=16, ...)
}
\arguments{
    \item{model}{A specified structural equation model of class "lms".}
    \item{data}{The data the model should be fitted to.}
    \item{start}{Starting values for parameters.}
    \item{logger}{If output of EM algorithm should be shown during fitting.}
    \item{threshold}{Convergence threshold.}
    \item{max.iter}{Maximum number of Iterations before EM algorithm stops.}
    \item{m}{Number of nodes for Hermite-Gaussian quadrature.}
    \item{...}{Additional arguments. See Details.}
}
\references{
    Klein, A. &, Moosbrugger, H. (2000). Maximum likelihood estimation of
    latent interaction effects with the LMS method. \emph{Psychometrika, 65},
    457--474.
}
\examples{
# specify model
specs <- specify_model(num.x=6, num.y=3, num.xi=2, num.eta=1, xi="x1-x3,x4-x6", eta="y1-y3")
model <- fill_matrices(specs)

# simulate data
pars.orig <- c(0.6, 0.7, 0.8, 0.7,                                   # Lx
              0.5, 0.3,                                              # Ly
              0.2, 0.4,                                              # G
              0.51, 0.64, 0.36, 0.25, 0.67, 0.23,                    # Td
              0.31, 0.45, 0.15,                                      # Te
              0.2,                                                   # Psi
              0.49, 0.64, 0.54,                                      # A
              0.7                                                    # Omega
              )
mod.filled <- fill_model(model, pars.orig)
dat <- simulate(mod.filled)

# fit model
pars.start <- runif(22, 0.4, 0.6)

\dontrun{
system.time({
    res <- em(model, dat, pars.start, logger=TRUE)
})
summary(res)
}
}


