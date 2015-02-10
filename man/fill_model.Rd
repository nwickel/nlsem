% File src/library/base/man/fill_model.Rd
\name{fill_model}
\alias{fill_model}
\title{Fills an empty structural equation model with parameters}
\description{
    Creates a model of the same class as \code{model} and puts
    \code{parameters} where \code{model} has \code{NA}'s.
}
\usage{
    fill_model(model, parameters)
}
\arguments{
    \item{model}{a model created by \code{\link{specify_sem}} or
    \code{\link{create_sem}}.}
    \item{parameters}{numeric vector with length of number of free parameters in
    \code{model}.\cr See \code{\link{count_free_parameters}}.}
}
\value{
    Gives back an object of class \code{lms}, \code{semm}, or \code{nsemm}.
}
\seealso{\code{\link{specify_sem}}, \code{\link{create_sem}},
            \code{\link{count_free_parameters}}}
\examples{
# specify model
model <- specify_sem(num.x = 4, num.y = 1, num.xi = 2, num.eta = 1, 
                     xi = "x1-x2,x3-x4", eta = "y1",
                     interaction = "xi1:xi2")
pars <- runif(count_free_parameters(model))
fill_model(model, parameters = pars)
}
