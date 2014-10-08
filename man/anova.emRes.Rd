% File src/library/base/man/anova.emRes.Rd
\name{anova.emRes}
\alias{anova.emRes}
\title{Anova Tables}
\description{
    Calculate likelihood ratio statistic for one or several structural equation models.
}
\usage{
\method{anova}{emRes}(object, ..., test=c("Chisq", "none"))
}
\arguments{
    \item{object}{estimated structural equation model of class "emRes".}
    \item{...}{additional objects of the same type.}
    \item{test}{a character string, (partially) matching one of "Chisq" or
    "none". Should the p-values of the chi-square distributions be
    reported?}
}
\value{
     Returns an object of class "anova".  These objects represent
     analysis-of-variance and analysis-of-deviance tables. It is not
     implemented for a single argument. 
}
\seealso{\code{\link{AIC.emRes}}, \code{\link{logLik.emRes}}}

