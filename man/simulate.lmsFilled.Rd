% File src/library/base/man/simulate.lmsFilled.Rd
\name{simulate.lmsFilled}
\alias{simulate.lmsFilled}
\title{Simulate data from a structural equation model.}
\description{
    Simulate data from a structural equation model.
}
\usage{
\method{simulate}{lmsFilled}(object, nsim=1, seed=NULL, n=400, m=16, ...)
}
\arguments{
    \item{object}{Structural equation model of class "lmsFilled".}
    \item{nsim}{Don't know. Do I need it?}
    \item{seed}{Set seed. Default is NULL.}
    \item{n}{Data for how many observations should be simulated.}
    \item{m}{Number of nodes for Hermite-Gaussian quadrature.}
    \item{...}{Additional arguments.}
}
\examples{
# specify model
specs <- specify_model(num.x=6, num.y=3, num.xi=2, num.eta=1,
                       xi="x1-x3,x4-x6", eta="y1-y3",
                       interaction="xi1:xi2")
model <- fill_matrices(specs)

# original parameters
pars.orig <- c(.6, .5, .4, .5, .4, .6, .5, .2, .6, .7, .3, .2, .5, .7, .3, .4,
                .6, .2, .3, .4, .6, .8)

# put parameters into model
mod.filled <- fill_model(model, pars.orig)

# simulate data from model
dat <- simulate(mod.filled)
}

