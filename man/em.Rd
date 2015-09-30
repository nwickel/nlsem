% File src/library/base/man/em.Rd
\name{em}
\alias{em}
\title{Maximum likelihood estimation of structural equation mixture models}
\description{
  Fits a structural equation model with latent interaction effects using
  mixture approaches (LMS, SEMM, NSEMM).
}
\usage{
  em(model, data, start, qml = FALSE, verbose = FALSE, convergence = 1e-02,
     max.iter = 100, m = 16, optimizer = c("nlminb", "optim"),
     max.mstep = 1, max.singleClass = 1, neg.hessian = TRUE, ...) }
\arguments{
  \item{model}{a specified structural equation model of class
  \code{singleClass}, \code{semm}, or \code{nsemm}.}
  \item{data}{the data the model should be fitted to. Data needs to be a
  matrix and variables need to be in the order x1, x2, ..., y1, y2, ...
  as specified in \code{\link{specify_sem}}. Data matrix needs no column
  names (will be ignored anyways).}
  \item{start}{starting values for parameters.}
  \item{qml}{logical. Indicating if QML estimation should be used instead
  of LMS for estimation of nonlinear effects. Defaults to FALSE. QML is
  much faster, though.}
  \item{verbose}{if output of EM algorithm should be shown during fitting.}
  \item{convergence}{convergence threshold.}
  \item{max.iter}{maximum number of iterations before EM algorithm stops.}
  \item{m}{number of nodes for Hermite-Gaussian quadrature. Defaults to
  16. See Datails.}
  \item{optimizer}{which optimizer should be used in maximization step of
  EM algorithm: \code{\link{nlminb}} or \code{\link{optim}}.}
  \item{max.mstep}{maximum iteration steps the optimizer should use in
  its mstep during one EM iteration. Defaults to 1.}
  \item{max.singleClass}{maximum iteration steps for singleClass model
  inside of NSEMM model. Defaults to 1 (and should only be changed for
  valid reasons).}
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
  \item{model.class}{class of model that was fitted, can be
  \code{singleClass}, \code{semm}, or \code{nsemm}.}
  \item{coefficients}{estimated parameters.}
  \item{objective}{final loglikelihood obtained with EM algorithm.}
  \item{em_convergence}{yes or no. Did EM algorithm converge?}
  \item{Hessian}{Hessian matrix for final parameter estimation.}
  \item{loglikelihoods}{loglikelihoods obtained during each iteration of
  EM algorithm.}
  \item{info}{list of number of exogenous (\code{num.xi}) and endogenous
  (\code{num.eta}) variables and of indicators (\code{num.x} and
  \code{num.y}). Corresponds to specifications given to
  \code{\link{specify_sem}} when specifiying structural equation model.}
}
\references{
  Jedidi, K., Jagpal, H. S., & DeSarbo, W. S. (1997). STEMM: A General
  Finite Mixture Structural Equation Model, \emph{Journal of
  Classification, 14}, 23--50.

  Kelava, A., Nagengast, B., & Brandt, H. (2014). A nonlinear structural
  equation mixture modeling approach for non-normally distributed latent
  predictor variables. \emph{Structural Equation Modeling, 21}, 468-481.

  Klein, A. &, Moosbrugger, H. (2000). Maximum likelihood estimation of
  latent interaction effects with the LMS method. \emph{Psychometrika, 65},
  457--474.
}
\seealso{\code{\link{specify_sem}}}
\examples{

###### Example for SEMM ######
# load data
data("PoliticalDemocracy", package = "lavaan")
dat <- as.matrix(PoliticalDemocracy[ ,c(9:11,1:8)])

# specify model of class SEMM
model <- specify_sem(num.x = 3, num.y = 8, num.xi = 1, num.eta = 2, 
  xi = "x1-x3", eta = "y1-y4,y5-y8", rel.lat =
  "eta1~xi1,eta2~xi1,eta2~eta1", num.classes = 2)

# fit model
set.seed(911)
start <- runif(count_free_parameters(model))
\dontrun{
res <- em(model, dat, start, convergence = 0.1, max.iter = 200)
summary(res)
plot(res)
}

###### Example for LMS ######
# specify model
model <- specify_sem(num.x = 11, num.y = 4, num.xi = 2, num.eta = 1,
  xi = "x1-x5,x6-x11", eta = "y1-y4", interaction = "eta~xi1:xi2")

data("jordan")

set.seed(110)
start <- runif(count_free_parameters(model))
\dontrun{
res <- em(model, jordan, start, convergence=1, verbose=TRUE)
summary(res)
plot(res)
}
}
