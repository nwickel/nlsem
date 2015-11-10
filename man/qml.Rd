% File src/library/base/man/qml.Rd
\name{qml}
\alias{qml}
\title{Quasi-maximum likelihood estimation of a nonlinear structural equation model}
\description{
  Fits a structural equation model with latent interaction effects using
  Quasi-maximum likelihood estimation.
}
\usage{
  qml(model, data, start, max.iter = 150, optimizer = c("nlminb",
      "optim"), neg.hessian = TRUE, ...) }
\arguments{
  \item{model}{a specified structural equation model of class
  \code{singleClass}.}
  \item{data}{the data the model should be fitted to. Data needs to be a
  matrix and variables need to be in the order x1, x2, ..., y1, y2, ...
  as specified in \code{\link{specify_sem}}. Data matrix needs no column
  names (will be ignored anyways).}
  \item{start}{starting values for parameters.}
  \item{max.iter}{maximum number of iterations for optimizer.}
  \item{optimizer}{which optimizer should be used for maximization of
  parameters: \code{\link{nlminb}} or \code{\link{optim}}.}
  \item{neg.hessian}{should negative Hessian be calculated.}
  \item{...}{additional arguments. See Details.}
}
\details{
  Additional arguments can be passed to \dots for these optimizers. See
  documentation for \code{\link{optim}} and \code{\link{nlminb}}.

  Quasi-maximum likelihood (QML) estimation is in principle a faster
  version for LMS, but might be less accurate for normal data. For
  practical purposes differences are negligible, though. For nonnormal
  data QML outperforms LMS.
}
\value{
  An object of class \code{qmlEst} that consists of the following components:
  \item{model.class}{class of model that was fitted. Will always be
  \code{singleClass}.}
  \item{coefficients}{estimated parameters.}
  \item{objective}{final loglikelihood obtained with EM algorithm.}
  \item{convergence}{convergence code for optimizer. See documentation
  for \code{\link{optim}} and \code{\link{nlminb}}.}
  \item{Hessian}{negative Hessian matrix for final parameter estimation.}
  \item{info}{list of number of exogenous (\code{num.xi}) and endogenous
  (\code{num.eta}) variables and of indicators (\code{num.x} and
  \code{num.y}). Corresponds to specifications given to
  \code{\link{specify_sem}} when specifiying structural equation model.}
}
\references{
  Klein, A. &, Muthen, B. O. (2007). Quasi-Maximum Likelihood Estimation
  of Structural Equation Models With Multiple Interaction and Quadratic
  Effects \emph{Multivariate Behavioral Research, 42}, 647--673.
}
\seealso{\code{\link{specify_sem}}}
\examples{
# specify model of class singleClass
sc <- specify_sem(num.x=4, num.y=2, num.xi=2, num.eta=1, xi="x1-x2,x3-x4",
                  eta="y1-y2", interaction="eta1~xi1:xi2") 

# simulate data
pars.orig <- c(0.6, 0.7,                 # Lx
               0.8,                      # Ly
               0.2, 0.4,                 # G
               0.25, 0.25, 0.25, 0.25,   # Td
               0.25, 0.25,               # Te
               0.2,                      # Psi
               0.49, 0.235, 0.64,        # Phi
               0, 0,                     # nu.x
               0,                        # nu.x
               1,                        # alpha
               1, 1,                     # tau
               0.7                       # Omega
              )

dat <- simulate(sc, parameters=pars.orig, seed=81)

# fit model
set.seed(1609)
start <- runif(count_free_parameters(sc))
\dontrun{
qml1 <- qml(sc, dat, start)
summary(qml1)
}
}

