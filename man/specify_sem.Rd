% File src/library/base/man/specify_sem.Rd
\name{specify_sem}
\alias{specify_sem}
\title{Specify a structural equation model}
\description{
    Specify a structural equation model with constraints.
}
\usage{
specify_sem(num.x, num.y, num.xi, num.eta, xi, eta, num.groups=1,
            interaction="all", interc_obs=FALSE,
            interc_lat=FALSE, relation_lat="default")
}
\arguments{
    \item{num.x}{number of observed variables for xi.}
    \item{num.y}{number of observed variables for eta}
    \item{num.xi}{number of latent exogenous variabeles.}
    \item{num.eta}{number of latent endognous variables.}
    \item{num.groups}{number of latent groups}
    \item{xi}{which observed variables are indicators for which exogenous
    variable. See Details}
    \item{eta}{which observed variables are indicators for which endogenous
    variable. See Details.}
    \item{interaction}{define which interaction terms should be included.
    Default is 'all' which means that all interactions between exogenous
    variables including quadratic terms are included. See Details for how
    to enter interaction terms.}
    \item{interc_obs}{logical. If intercepts for measurement model should be
    estimated (mu's). Default is \code{FALSE} (intercepts are fixed at 0).}
    \item{interc_lat}{logical. If intercepts for latent variable model should be
    estimated (alpha's). Default is \code{FALSE} (intercepts are fixed at 0).}
    \item{relation_lat}{define relations between latent variables.
    Influences Beta and Gamma matrices. For 'defaults' and how to define
    see Details.}
}
\details{
    Which indicators belong to which latent variable is defined by
    \code{xi} and \code{eta}. Must be specified in the following way:
    \code{xi='x1-x2,x3-x4'} which means that variables \code{x1, x2} are
    indicators for \code{xi1} and \code{x3, x4} are indicators for
    \code{xi2}. And accordingly for the endogenous variables \code{eta}.

    Interactions between latent exogenous variables are defined like
    \code{interaction='xi1:xi2,xi1:xi1}. It is important to note, that
    interactions must always start with \code{xi1} and build from there. A
    definitin like \code{interaction='xi1:xi2,xi2:xi3} is not feasible and
    must be changed to \code{interaction='xi1:xi2,xi1:xi3} (by simple
    switching \code{xi1} and \code{xi2} in my definitions).

    \code{relation_lat} defines which latent variables influence each
    other. It must be defined like
    \code{relation_lat='xi1,xi2>eta1;eta1>eta2'}. Free parameters will be
    set accordingly in Beta and Gamma matrices. When nothing is defined
    Gamma defaults to all NAs (which means all \code{xi}'s influence all
    \code{eta}'s) and Beta is an identity matrix.
}
\value{
    An object of class 'lms', 'stemm', or 'nsemm' which can be used to
    estimate parameters using \code{\link{em}} that consists of the
    following components:
    \item{matrices}{list of matrices specifying the structural equation
    model.}
    \item{info}{list of informations about structural equation model.}
}
\examples{
# with default constraints
model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1, xi="x1-x3,x4-x6", eta="y1-y3")

# create data frame
specs <- as.data.frame(model)
# and add custom constraints
constr <- c(1, NA, NA, 0, 0, 0, 0, 0, 0, 1, NA, NA, 1, NA, NA, NA, NA, 0, NA,
    0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0,
    0, 0, 0, NA, 0, 0, 0, 0, 0, 0, NA, NA, 0, 0, 0, NA, 0, 0, 0, NA, NA, NA,
    NA, 0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0, 0, 0, 0, NA, 0)
specs$group1 <- constr
}

