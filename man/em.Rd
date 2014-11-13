% File src/library/base/man/em.Rd
\name{em}
\alias{em}
\title{Maximum likelihood estimation of latent interaction effects and
mixture structural equation models}
\description{
    Fits a structural equation model with latent interaction effects using
    mixture approaches (LMS, STEMM, NSEMM).
}
\usage{
    em(model, data, start, logger=FALSE, threshold=1e-05, max.iter=40,
       m=16, optimizer=c("nlminb", "optim"), ...)
}
\arguments{
    \item{model}{a specified structural equation model of class "lms".}
    \item{data}{the data the model should be fitted to.}
    \item{start}{starting values for parameters.}
    \item{logger}{if output of EM algorithm should be shown during fitting.}
    \item{threshold}{convergence threshold.}
    \item{max.iter}{maximum number of Iterations before EM algorithm stops.}
    \item{m}{number of nodes for Hermite-Gaussian quadrature.}
    \item{optimizer}{which optimizer should be used in maximization step of
    EM algorithm: \code{\link{nlminb}} or \code{\link{optim}}.}
    \item{...}{additional arguments. See Details.}
}
\value{
    An object of class 'emEst' that consists of the following components:
    \item{model.class}{What class model has that was fitted, can be 'lms',
    'stemm', or 'nsemm'.}
    \item{par}{estimated parameters.}
    \item{objective}{final loglikelihood obtained with EM algorithm.}
    \item{convergence_final_step}{convergence message of final step of
    maximization step either done with \code{\link{optim}}} or
    \code{\link{nlminb}.}
    \item{message_final_step}{message of final step. See above.}
    \item{Hessian}{Hesse matrix for final parameter estimation.}
    \item{loglikelihoods}{loglikelihoods obtained during each iteratin of
    EM algorithm}
    \item{info}{list of number of exogenous (\code{num.xi}) and endogenous
    (\code{num.eta}) variables and of indicators (\code{num.x} and
    \code{num.y}). Corresponds to specifications given to
    \code{\link{specify_sem}} when specifiying structural equation model.}
}
\details{
    \code{em} can be used to estimate parameters for structural equation
    mixture models with latent interaction effects with an EM algorithm.
    The maximization step of the EM algorithm can use two different
    optimizers: \code{\link{optim}} or \code{\link{nlminb}}. Default is
    \code{\link{nlminb}}.

    Additional arguments can be passed to \dots for these optimizers, like
    iterations etc. See documentation for \code{\link{optim}} and
    \code{\link{nlminb}}.
}
\references{
    Klein, A. &, Moosbrugger, H. (2000). Maximum likelihood estimation of
    latent interaction effects with the LMS method. \emph{Psychometrika, 65},
    457--474.

    Jedidi, K., Jagpal, H. S., & DeSarbo, W. S. (1997). STEMM: A General
    Finite Mixture Structural Equation Model, \emph{Journal of
    Classification, 14}, 23--50.

    Kelava, A., Nagengast, B., & Brandt, H. (2014). A nonlinear structural
    equation mixture modeling approach for non-normally distributed latent
    predictor variables. \emph{Structural Equation Modeling, 21}, 468-481.
}
\examples{
# specify model
model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                       xi="x1-x3,x4-x6", eta="y1-y3",
                       interaction="xi1:xi2")

# simulate data
pars.orig <- c(0.6, 0.7, 0.8, 0.7,                                   # Lx
              0.5, 0.3,                                              # Ly
              0.2, 0.4,                                              # G
              0.51, 0.64, 0.36, 0.25, 0.67, 0.23,                    # Td
              0.31, 0.45, 0.15,                                      # Te
              0.2,                                                   # Psi
              0.49, 0.64, 0.54,                                      # A
              0.7                                                    # Omega
              )

dat <- simulate(model, parameters=pars.orig)

# fit model
pars.start <- runif(22, 0.4, 0.6)

\dontrun{
system.time({
    res <- em(model, dat, pars.start, logger=TRUE)
})
summary(res)
}
}


