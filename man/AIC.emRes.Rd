% File src/library/base/man/AIC.emRes.Rd
\name{AIC.emRes}
\alias{AIC.emRes}
\alias{logLik.emRes}
\title{Akaike's Information Criterion}
\description{
    Calculate AIC for one or several structural equation models.
}
\usage{
\method{AIC}{emRes}(object, ..., k=2)
\method{logLik}{emRes}(object, ...)
}
\arguments{
    \item{object}{estimated structural equation model of class "emRes".}
    \item{...}{optionally more fitted model objects.}
    \item{k}{numeric, the _penalty_ per parameter to be used; the default
          k = 2 is the classical AIC.}
}
\value{
     A data.frame with rows corresponding to the objects and columns
     representing the number of parameters in the model (df) and the AIC.
}
\seealso{\code{\link{anova.emRes}}}

