% File src/library/base/man/specify_model.Rd
\name{specify_model}
\alias{specify_model}
\title{Specify a structural equation model}
\description{
    Specify a structural equation model with constraints.
}
\usage{
specify_model(num.x, num.y, num.xi, num.eta, xi, eta, interaction="all",
              constraints="default", interc_obs=FALSE, interc_lat=FALSE)
}
\arguments{
    \item{num.x}{number of observed variables for xi.}
    \item{num.y}{number of observed variables for eta}
    \item{num.xi}{number of latent exogenous variabeles.}
    \item{num.eta}{number of latent endognous variables.}
    \item{xi}{which observed variables are indicators for which exogenous
    variable. See Details}
    \item{eta}{which observed variables are indicators for which endogenous
    variable. See Details.}
    \item{interaction}{define which interaction terms should be included.
    Default is "all" which means that all interactions between exogenous
    variables including quadratic terms are included. See Details for how
    to enter interaction terms.}
    \item{constraints}{a vector of constraints. If set to "default"
    (default) then parameters are constrained to common practice. See
    Details.}
    \item{interc_obs}{logical. If intercepts for measurement model should be
    estimated (mu's). Default is FALSE (intercepts are fixed at 0).}
    \item{interc_lat}{logical. If intercepts for latent variable model should be
    estimated (alpha's). Default is FALSE (intercepts are fixed at 0).}
}
\examples{
# with default constraints
specs <- specify_model(num.x=6, num.y=3, num.xi=2, num.eta=1, xi="x1-x3,x4-x6", eta="y1-y3")

# add custom constraints
# example 1
constr <- c(1, NA, NA, 0, 0, 0, 0, 0, 0, 1, NA, NA, 1, NA, NA, NA, NA, 0, NA,
    0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0,
    0, 0, 0, NA, 0, 0, 0, 0, 0, 0, NA, NA, 0, 0, 0, NA, 0, 0, 0, NA, NA, NA,
    NA, 0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0, 0, 0, 0, NA, 0)
specs$ustart <- constr

# eaxmple 2 (in one step)
specs <- specify_model(num.x=6, num.y=3, num.xi=2, num.eta=1, xi="x1-x3,x4-x6", eta="y1-y3", constraints=constr)
}

