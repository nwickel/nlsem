# nonnormal.R
# 
# Generating multivariate nonnormal data based on Vale & Maurelli (1983).
# Based on the code from Cengiz Zopluoglu (2011).
#
# created Jun/22/2015, NU

#--------------- main functions ---------------

generate_nonnormal <- function(n, cor, skew, kurt, k, ...) {

  stopifnot(length(skew) == k)
  stopifnot(length(kurt) == k)

  # Inputs:
  # n    = sample size
  # k    = number of variables
  # cor  = desired correlation matrix between bivariate non-normal variables
  # skew = a vector of k elements, skewness for the variables
  # kurt = a vector of k elements, kurtosis for the variables

  # Compute the constants a, b, c, and d for each variable with a desired
  # skewness and kurtosis
  
  constants <- matrix(nrow=k, ncol=4)
  for(i in 1:k) {
    constants[i, ] <- solve_constants(skew[i], kurt[i])
  }
  
  # Compute intermediate intercorrelation matrix required for normal
  # variables. Normal variables are used to construct nonnormal
  # variables.
  
  inter <- matrix(nrow=k, ncol=k)
  for(i in 1:k) {
    for(j in 1:k) {
      inter[i,j] <- solve_rho(cor[i,j], constants[i,1], constants[j,1],
        constants[i,2], constants[j,2], constants[i,3], constants[j,3],
        constants[i,4], constants[j,4])
    }
  }

  diag(inter) <- 1
  
  # Compute multivariate normal variables based on intermediate
  # intercorrelation matrix
  
  # Eigen decomposition of correlation matrix
  U <- eigen(inter)$vectors
  L <- eigen(inter)$values
  b <- U %*% diag(sqrt(L))
  
  # Creating independent multivariate normal variables
  normal <- matrix(nrow=n, ncol=k)
  for(i in 1:k) normal[,i] <- rnorm(n)
  
  # Creating correlated multivariate normal variables
  dat <- as.data.frame(normal %*% t(b))
  
  # Creating correlated nonnormal multivariate variables from correlated
  # multivariate normal variables using constants a, b, c, and d 
  nonnormal <- as.data.frame(matrix(nrow=n, ncol=k))
  
  for(i in 1:k) {
    nonnormal[,i] <- constants[i,1] + constants[i,2]*dat[,i] +
    constants[i,3]*dat[,i]^2 + constants[i,4]*dat[,i]^3 
  }
  
  nonnormal
}

#--------------- helper functions ---------------

# Function to compute a, b, c, d for a variable given skewness and
# kurtosis 
# Use Newton-Raphson Iteration with a Jacobian matrix to solve system
# of nonlinear Equations 2, 3, and 4 in Vale & Maurelli (1983)
  
solve_constants <- function(skew, kurt, start=c(1,0,0), tol=.00001, max.iter=500){

  # sk = desired skewness
  # ku = desired kurtosis
  # start = starting values for the iteration, based on Fleishman (1978)

  fn <- function(x, skew, kurt){
  
    f <- matrix(0, nrow=length(x))
    b <- x[1]
    c <- x[2]
    d <- x[3]
    # Eq. 2-4
    f[1] <- b^2 + 6*b*d + 2*c^2 + 15*d^2 - 1
    f[2] <- 2*c*(b^2 + 24*b*d + 105*d^2 + 2) - skew
    f[3] <- 24*(b*d + c^2*(1 + b^2 + 28*b*d) + d^2*(12 + 48*b*d + 141*c^2 +
      225*d^2)) - kurt
    
    # Jacobian matrix
    j <- matrix(0, ncol=length(x), nrow=length(x))
    j[1,1] <- 2*b + 6*d
    j[1,2] <- 4*c
    j[1,3] <- 6*b + 30*d
    j[2,1] <- 4*b*c + 48*c*d
    j[2,2] <- 2*b^2 + 48*b*d + 210*d^2+4
    j[2,3] <- 48*b * c+420*c*d
    j[3,1] <- 24*d + 48 * c^2*b
    j[3,2] <- 48*c + 48 * c * b^2 + 1344*c*b*d + 6768*c*d^2
    j[3,3] <- 24*b + 672*c^2*b + 576*d + 3456*b*d^2 + 6768*d*c^2 + 21600*d^3

    solve(j) %*% f
  }
  
  x0   <- start
  d    <- fn(x0, skew=skew, kurt=kurt)
  iter <- 0
  
  while((abs(d[1]) > tol) && (iter < max.iter)) {
    x0   <- x0 - fn(x0, skew=skew, kurt=kurt)
    d    <- fn(x0, skew=skew, kurt=kurt)
    iter <- iter + 1
  }
  retval <- c(-x0[2], x0)
  retval
} 

# Function to solve the polynomial function to find the intermediate
# correlation between two normal variables for a given desired correlation
# between two nonnormal variables and constants a, b, c, d 
# Use Newton-Raphson iteration to approximate the root

solve_rho <- function(r12, a1, a2, b1, b2, c1, c2, d1, d2, tol=.00001,
                      max.iter=500, start=.5) {

  ftn <- function(rho) {
    a <-((b1*b2 + 3*b1*d2 + 3*d1*b2 + 9*d1*d2)*rho) + ((2*c1*c2)*rho^2) +
      ((6*d1*d2)*rho^3) - r12
    b <-(b1*b2 + 3*b1*d2 + 3*d1*b2 + 9*d1*d2) + ((4*c1*c2)*rho) +
      ((12*d1*d2)*rho^2) 
    c(a,b)
  }
  
  rho  <- start
  fx   <- ftn(rho)
  iter <- 0
  
  while((abs(fx[1]) > tol) && (iter < max.iter)) {
    rho  <- rho - fx[1]/fx[2]
    fx   <- ftn(rho)
    iter <- iter + 1
  }
  rho
}

