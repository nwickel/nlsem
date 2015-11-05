\name{nlsem-package}
\alias{nlsem-package}
\docType{package}
\title{
Fitting structural equation mixture models
}
\description{
Estimation of structural equation models with nonlinear effects and underlying
nonnormal distributions.
}
\details{
This is a package for estimating nonlinear structural equation mixture
models using an expectation-maximization (EM) algorithm. Four different
approaches are implemented. Firstly, the Latent Moderated Structural
Equations (LMS) approach (Klein & Moosbrugger, 2000) and the Quasi-Maximum
Likelihood (QML) approach (Klein & Muthen, 2007), which allow for two-way
interaction and quadratic terms in the structural model. Due to the
nonlinearity, the latent criterion variables cannot be assumed to be
normally distributed. Therefore, the latent criterins's distribution is
approximated with a mixture of normal distributions in LMS. Secondly, the
Structural Equation finite Mixture Model (STEMM or SEMM) approach (Jedidi,
Jagpal & DeSarbo, 1997), which uses mixtures to model latent classes. In
this way it can deal with heterogeneity in the sample or nonlinearity and
nonnormality of the latent variables and their indicators. And thirdly, a
combination of these two approaches, the Nonlinear Structural Equation
Mixture Model (NSEMM) approach (Kelava, Nagengast & Brandt, 2014). Here,
interaction and quadratic terms as well as latent classes can be modeled.

The models can be specified with \code{\link{specify_sem}}. Depending
on the specification of \code{interaction} and the number of latent classes
(\code{num.classes}) the returned object will be of class
\code{singleClass}, \code{semm}, or \code{nsemm}. Each of these can be
estimated using \code{\link{em}} and models of type \code{singleClass}
can additionally be fitted with the function \code{\link{qml}}.
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

  Klein, A. &, Muthen, B. O. (2007). Quasi-Maximum Likelihood Estimation
  of Structural Equation Models With Multiple Interaction and Quadratic
  Effects \emph{Multivariate Behavioral Research, 42}, 647--673.
}
\section{Future Features}{
\itemize{
  \item NSEMM, LMS and QML for more than one latent endogenous variable.
  \item Parameter standardization.
}
}
\keyword{package}

