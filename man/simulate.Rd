% File src/library/base/man/simulate.Rd
\name{simulate}
\alias{simulate.lms}
\alias{simulate.semm}
\alias{simulate.nsemm}
\title{Simulate data from a structural equation model}
\description{
    Simulate data from a structural equation mixture model.
}
\usage{
\method{simulate}{lms}(object, nsim=1, seed=NULL, n=400, m=16, parameters, ...)
\method{simulate}{semm}(object, nsim=1, seed=NULL, n=400, parameters, ...)
\method{simulate}{nsemm}(object, nsim=1, seed=NULL, n=400, m=16, parameters, ...)
}
\arguments{
    \item{object}{structural equation model of class \code{lms},
    \code{semm}, or \code{nsemm}.}
    \item{parameters}{`true' parameters which should be used to simulate
    data.}
    \item{nsim}{number of response vectors to simulate.  Defaults to 1.}
    \item{seed}{set seed. Default is NULL.}
    \item{n}{data for how many observations should be simulated.}
    \item{m}{number of nodes for Hermite-Gaussian quadrature. Only needed
    for \code{lms} and \code{nsemm}.}
    \item{...}{additional arguments.}
}
\value{
    Returns a matrix with n rows and as many columns as indicators are
    entered into the model.
}
\examples{
# specify model
model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                     xi="x1-x3,x4-x6", eta="y1-y3",
                     interaction="xi1:xi2")

# original parameters
pars.orig <- c(.6, .5, .4, .5, .4, .6, .5, .2, .6, .7, .3, .2, .5,
               .7, .3, .4, .6, .2, .3, .4, .6, .8)

# simulate data from model
dat <- simulate(model, parameters=pars.orig)
}

