% File src/library/base/man/anova.emRes.Rd
\name{coef.emRes}
\alias{coefficients.emRes}
\title{Extract Model Coeffcients}
\description{
    Extracts model coefficients form a structural equation model fitted
    with \code{\link{em}}.}
\usage{
\method{coef}{emRes}(object, ...)
}
\arguments{
    \item{object}{estimated structural equation model of class "emRes".}
    \item{...}{additional arguments.}
}
\value{
     Returns a numeric vector that contains estimated parameters of
     structural equation model.
}
\seealso{\code{\link{em}}}

