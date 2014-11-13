% File src/library/base/man/anova.emEst.Rd
\name{coef.emEst}
\alias{coef.emEst}
\alias{coefficients.emEst}
\title{Extract Model Coeffcients}
\description{
    Extracts model coefficients from a structural equation model fitted
    with \code{\link{em}}.}
\usage{
\method{coef}{emEst}(object, ...)
\method{coefficients}{emEst}(object, ...)
}
\arguments{
    \item{object}{estimated structural equation model of class 'emEst'.}
    \item{...}{additional arguments.}
}
\value{
     Returns a numeric vector that contains estimated parameters of
     structural equation model.
}
\seealso{\code{\link{em}}}

