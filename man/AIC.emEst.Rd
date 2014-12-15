% File src/library/base/man/AIC.emEst.Rd
\name{AIC.emEst}
\alias{AIC.emEst}
\alias{BIC.emEst}
\alias{logLik.emEst}
\title{Akaike's Information Criterion}
\description{
    Calculate AIC for one or several structural equation models.
}
\usage{
\method{AIC}{emEst}(object, ..., k=2)
\method{AIC}{emEst}(object, ...)
\method{logLik}{emEst}(object, ...)
}
\arguments{
    \item{object}{estimated structural equation model of class `emEst'.}
    \item{...}{optionally more fitted model objects.}
    \item{k}{numeric, the _penalty_ per parameter to be used; the default
          k = 2 is the classical AIC.}
}
\value{
     A data frame with rows corresponding to the objects and columns
     representing the number of parameters in the model (df) and the
     AIC/BIC.
}
\seealso{\code{\link{anova.emEst}}}

