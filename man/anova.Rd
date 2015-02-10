% File src/library/base/man/anova.Rd
\name{anova}
\alias{anova.emEst}
\title{Anova Tables}
\description{
    Calculate likelihood ratio tests to compare two or more structural equation models.
}
\usage{
\method{anova}{emEst}(object, ..., test = c("Chisq", "none"))
}
\arguments{
    \item{object}{estimated structural equation model of class \code{emEst}.}
    \item{...}{additional objects of the same type.}
    \item{test}{a character string, (partially) matching one of `Chisq' or
    `none'. Should the p-values of the chi-square distributions be
    reported?}
}
\value{
     Returns an object of class code{anova}.  These objects represent
     analysis-of-variance and analysis-of-deviance tables. It is not
     implemented for a single argument. 
}

