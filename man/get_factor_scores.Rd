% File src/library/base/man/get_factor_scores.Rd
\name{get_factor_scores}
\alias{get_factor_scores}
\title{Obtain factor scores from fitted nlsem object}
\description{
  Bartlett factor scores can be obtained from a fitted model.
}
\usage{
  get_factor_scores(object, model, data)
}
\arguments{
  \item{object}{estimated structural equation model of class \code{emEst}
  obtained from \code{\link{em}}.}
  \item{model}{a model created by \code{\link{specify_sem}} or
  \code{\link{create_sem}}.}
  \item{data}{the model was fitted to.}
}
\value{
  Gives back a matrix with as many columns as number of exogenous
  variables and as many rows as number of observations.
}
\examples{
# specify model
model <- specify_sem(num.x = 11, num.y = 4, num.xi = 2, num.eta = 1, 
  xi = "x1-x5,x6-x11", eta = "y1-y4", interaction = "xi1:xi2")

data("jordan")

start <- runif(count_free_parameters(model))
\dontrun{
res <- em(model, jordan, start)

fs <- get_factor_scores(res, model, jordan)
}
}
