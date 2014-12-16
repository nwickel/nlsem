% File src/library/base/man/em.Rd
\name{em}
\alias{em}
\title{Maximum likelihood estimation of latent interaction effects and
mixture structural equation models}
\description{
    Fits a structural equation model with latent interaction effects using
    mixture approaches (LMS, SEMM, NSEMM).
}
\usage{
    em(model, data, start, logger=FALSE, convergence=1e-02, max.iter=100,
       m=16, optimizer=c("nlminb", "optim"), max.mstep=1, 
       neg.hessian=TRUE, ...) }
\arguments{
    \item{model}{a specified structural equation model of class \code{lms},
    \code{semm}, or \code{nsemm}.}
    \item{data}{the data the model should be fitted to.}
    \item{start}{starting values for parameters.}
    \item{logger}{if output of EM algorithm should be shown during fitting.}
    \item{convergence}{convergence threshold.}
    \item{max.iter}{maximum number of iterations before EM algorithm stops.}
    \item{m}{number of nodes for Hermite-Gaussian quadrature. Defaults to
    16. See Datails.}
    \item{optimizer}{which optimizer should be used in maximization step of
    EM algorithm: \code{\link{nlminb}} or \code{\link{optim}}.}
    \item{max.mstep}{maximum iteration steps the optimizer should use in
    its mstep during one EM iteration. Defaults to 1.}
    \item{neg.hessian}{should negative Hessian be calculated in last step
    of iteration.}
    \item{...}{additional arguments. See Details.}
}
\details{
    \code{em} can be used to estimate parameters for structural equation
    mixture models with latent interaction effects with an EM algorithm.
    The maximization step of the EM algorithm can use two different
    optimizers: \code{\link{optim}} or \code{\link{nlminb}}. Default is
    \code{\link{nlminb}}.

    Additional arguments can be passed to \dots for these optimizers. See
    documentation for \code{\link{optim}} and \code{\link{nlminb}}.

    The LMS approach (Klein & Moosbrugger, 2000) uses Hermite-Gauss
    quadrature for numerical approximation. The nodes used in this
    approximation need to be prespecified by the user. The more nodes are
    used the better the numerical approximation but also the slower the
    calculations.
}
\value{
    An object of class \code{emEst} that consists of the following components:
    \item{model.class}{class of model that was fitted, can be \code{lms},
    \code{semm}, or \code{nsemm}.}
    \item{coefficients}{estimated parameters.}
    \item{objective}{final loglikelihood obtained with EM algorithm.}
    \item{em_convergence}{yes or no. Did EM algorithm converge?}
    \item{Hessian}{Hessian matrix for final parameter estimation.}
    \item{loglikelihoods}{loglikelihoods obtained during each iteration of
    EM algorithm}
    \item{info}{list of number of exogenous (\code{num.xi}) and endogenous
    (\code{num.eta}) variables and of indicators (\code{num.x} and
    \code{num.y}). Corresponds to specifications given to
    \code{\link{specify_sem}} when specifiying structural equation model.}
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
\seealso{\code{\link{specify_sem}}}
\examples{

## LMS
## ===
# specify model
model1 <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                       xi="x1-x3,x4-x6", eta="y1-y3",
                       interaction="xi1:xi2")
# simulate data
pars.orig1 <- c(0.6, 0.7, 0.8, 0.7,                  # Lambda.x
              0.5, 0.3,                              # Lambda.y
              0.2, 0.4,                              # Gamma
              0.51, 0.64, 0.36, 0.25, 0.67, 0.23,    # Theta.d
              0.31, 0.45, 0.15,                      # Theta.e
              0.2,                                   # Psi
              0.49, 0.64, 0.54,                      # Phi
              0.7                                    # Omega
              )
dat1 <- simulate(model1, seed=34, parameters=pars.orig1)

# fit model
set.seed(4)
pars.start1 <- runif(count_free_parameters(model1), 0.4, 0.6)
\dontrun{
    system.time({
        res1 <- em(model1, dat1, pars.start1, logger=TRUE)
    })
    summary(res1)
}

## SEMM
## =====
# specify model
model2 <- specify_sem(num.x=4, num.y=2, num.xi=2, num.eta=1,
                     xi="x1-x2,x3-x4", eta="y1-y2", num.classes=2,
                     interaction="", interc.obs=FALSE, interc.lat=FALSE,
                     relation.lat="eta1~xi1+xi2")
# simulate data
pars.orig2 <- c(
                # class 1
                1, 1,           # Lambda.x
                1,              # Lambda.y
                0.2, 0.3,       # Gamma
                rep(0.5, 4),    # Theta.d
                0.5, 0.5,       # Theta.e
                0.5,            # Psi
                rep(0.8, 3),    # Phi
                # class 2
                1, 1,           # Lambda.x
                1,              # Lambda.y
                0.8, 0.9,       # Gamma
                rep(0.5, 4),    # Theta.d
                0.5, 0.5,       # Theta.e
                0.5,            # Psi
                rep(0.8, 3)     # Phi
)

dat2 <- simulate(model2, seed=1, parameters=pars.orig2)

# fit model
set.seed(7)
pars.start2 <- runif(count_free_parameters(model2), 0.1, 1)

\dontrun{
    system.time({
        res2 <- em(model2, dat2, pars.start2, logger=TRUE, max.iter=1000)
    })
    summary(res2)
}

## NSEMM
## =====
# specify model
model3 <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                     xi="x1-x3,x4-x6", eta="y1-y3",
                     num.classes=2, interaction="xi1:xi2",
                     interc.obs=FALSE, interc.lat=FALSE,
                     relation.lat="eta1~xi1+xi2")
# simulate data
pars.orig3 <- c(
                # class 1
                1, 1, 1, 1,         # Lambda.x
                1, 1,               # Lambda.y
                0.2, 0.8,           # Gamma
                rep(0.5, 6),        # Theta.d
                0.5, 0.5, 0.5,      # Theta.e
                0.5,                # Psi
                0.5, 0.5, 0.5,      # Phi
                0.3,                # Omega
                # class 2
                1, 1, 1, 1,         # Lambda.x
                1, 1,               # Lambda.y
                0.8, 0.2,           # Gamma
                rep(0.5, 6),        # Theta.d
                0.5, 0.5, 0.5,      # Theta.e
                0.5,                # Psi
                0.5, 0.5, 0.5,      # Phi
                0.3                 # Omega
)

dat3 <- simulate(model3, seed=1, parameters=pars.orig3)

# fit model
set.seed(6)
pars.start3 <- runif(count_free_parameters(model3), 0.1, 0.8)
\dontrun{
    system.time({
        res3 <- em(model3, dat3, pars.start3, logger=TRUE)
    })
    summary(res3)
}
}



