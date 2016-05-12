% File src/library/base/man/lav2nlsem.Rd
\name{lav2nlsem}
\alias{lav2nlsem}
\title{Create a structural equation model from lavaan syntax}
\description{
  Create model matrices from a string specifying a structural equation
  model in lavaan syntax.
}
\usage{
  lav2nlsem(model, constraints=c("indirect", "direct1", "direct2"),
            class.spec="class")
}
\arguments{
  \item{model}{A description of the user-specified model. The model is
  described using the lavaan model syntax. See Details in
  \code{?model.syntax} in lavaan for more information.}
  \item{constraints}{which should be set for a model with more than one
  latent class. See Details in \code{?specify_sem}.}
  \item{class.spec}{String used to specify latent classes. Can be any
  string e.g. `class', `mixture', etc. Default is `class'.}
}
\details{
TODO
}
\value{
  Gives back an object of class \code{singleClass}, \code{semm}, or
  \code{nsemm} which can be fitted using \code{\link{em}}.
}
\seealso{\code{\link{specify_sem}}, \code{\link{create_sem}}}
\examples{
# create model with three latent classes
lav.model <- '
          class: 1
          eta =~ y1 + y2 + y3 + y4
          xi1 =~ x1 + x2 + x3 + x4 + x5
          xi2 =~ x6 + x7 + x8 + x9 + x10 + x11

          eta ~ xi1 + xi2 + xi1:xi1

          class: 2
          eta =~ y1 + y2 + y3 + y4
          xi1 =~ x1 + x2 + x3 + x4 + x5
          xi2 =~ x6 + x7 + x8 + x9 + x10 + x11

          eta ~ xi1 + xi2 + xi1:xi2 + xi1:xi1

          class: 3
          eta =~ y1 + y2 + y3 + y4
          xi1 =~ x1 + x2 + x3 + x4 + x5
          xi2 =~ x6 + x7 + x8 + x9 + x10 + x11

          eta ~ xi1 + xi2 + xi1:xi2'

model <- lav2nlsem(lav.model, constraints = "direct1", class.spec = "class")
}
