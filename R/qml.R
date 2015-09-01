# qml.R
#
# created: Feb/04/2015, NU
# poluted: Jul/22/2015, HB
# last mod: Sep/01/2015, NU

#--------------- main functions ---------------

qml <- function(model, data, start, max.iter=150, 
                optimizer=c("nlminb", "optim"), neg.hessian=TRUE, ...) {

  if (is.matrix(data)) {
    data <- data
  } else if (is.data.frame(data)) {
    data <- as.matrix(data)
  } else {
    stop("data need to be a matrix or a data frame.")
  }

  suppressWarnings(
  est <- mstep_qml(model=model, data=data, parameters=start,
    neg.hessian=neg.hessian, optimizer=optimizer,
    max.iter=max.iter, ...)
  )
  # TODO Should these warnings really be supressed?

  names(est$par) <- model$info$par.names

  if (sum(est$par - start) == 0) {
    stop("NA/NaN function evaluation. Please try different set of starting parameters.")
  }
  # TODO: Maybe this should be inside mstep_qml under "nlminb"??
  
  out <- list(model.class=class(model), coefficients=est$par,
    objective=-est$objective,
    convergence=est$convergence,
    neg.hessian=est$hessian,
    iterations=est$iterations,
    info=model$info[c("num.xi", "num.eta", "num.x", "num.y",
      "num.classes")])

  class(out) <- "qmlEst"
  out
}

mu_qml <- function(model, data) {

  m <- model$matrices$class1 

  # extract x and y from data frame
  x <- data[, 1:model$info$num.x]
  y <- data[, (model$info$num.x + 1):dim(data)[2]]
  
  ny <- model$info$num.y
  ne <- model$info$num.eta
  
  
  if (model$info$num.y > 1) {
    # transformation of y
    # cf. Eq 7
    beta <- m$Lambda.y[(ne+1):ny,] 
    R <- cbind(-beta, diag(ny-ne))
    u <- y %*% t(R)
  } else {
    u <- 0
  # TODO: What happens to beta and R in this case??
  }

  
  # Eqs 15, 16
  Sigma1 <- m$Phi - m$Phi %*% t(m$Lambda.x) %*% solve(m$Lambda.x %*%
            m$Phi %*% t(m$Lambda.x) + m$Theta.d) %*% m$Lambda.x %*% m$Phi 

  L1 <- m$Phi %*% t(m$Lambda.x) %*% solve(m$Lambda.x %*% 
        m$Phi %*% t(m$Lambda.x) + m$Theta.d)
  # cf. Eq 24
  L2 <- -m$Theta.e[(1:ne),(1:ne)] %*% t(beta) %*% solve(R %*% m$Theta.e %*% t(R))
  
  # m.x is unconditional
  # m.u mean vector for R %*% epsilon (0)
  # Eq 14: but mu.x <- 0, i.e., without means for xi
  mu.x  <- m$nu.x + m$Lambda.x %*% m$tau 
  
  mu.u  <- as.matrix(rep(0, (ny-ne)))
  

  # m.y1 is conditional given x and u
  # Eq 12 
  # Beta must be transformed since it is defined for SEMM approach
  Beta <- m$Beta - diag(model$info$num.eta)
  B <- solve(diag(ne) - Beta)
  # TODO: Was passiert hier generell wenn diese Differenz nicht
  # invertierbar ist?
  alpha.b <- B %*% m$alpha
  ga1.b   <- B %*% m$Gamma
  if (is.matrix(m$Omega)) {
    Gamma2 <- t(vech(m$Omega))
  } else {
    Gamma2 <- t(apply(m$Omega, 3, vech))
  }
  ga2.b   <- B %*% Gamma2       # contains interaction and quadratic effects
  
  N <- nrow(data)

  mt.m   <- matrix(rep(m$tau, N), model$info$num.xi, N, byrow=F)
  ma.m   <- matrix(rep(alpha.b, N), ne, N, byrow=T)
  mux.m  <- matrix(rep(mu.x, N), N, model$info$num.x, byrow=T)
  muu.m  <- matrix(rep(mu.u, N), N, (ny-ne), byrow=T)
  mvy1.m <- matrix(rep(m$nu.y[1:ne], N), ne, N, byrow=T) 
  
  tmp <- matrix(nrow=ne, ncol=N)
  for (i in seq_len(N)) {
    tmp[,i] <- ga2.b %*% vech((mt.m + L1 %*% t(x[i,] - mux.m)) %*% 
      t(mt.m + L1 %*% t(x[i,] - mux.m)) + Sigma1)
  }

  # Eq. 20
  mu.y1 <- 
    mvy1.m +
    ma.m +
    ga1.b %*% (mt.m + L1 %*% t(x - mux.m)) +
    tmp +
    L2 %*% t(u-muu.m)
    
  mu.xu <- c(mu.x, mu.u)    # mean for f2
  mu <- list(mu.xu, mu.y1)  

  mu
}


sigma_qml <- function(model, data) {

  m <- model$matrices$class1

  # extract x and y from data frame
  x <- data[, 1:model$info$num.x]
  y <- data[, (model$info$num.x + 1):dim(data)[2]]

  ny <- model$info$num.y
  ne <- model$info$num.eta
  
  if (model$info$num.y > 1) {
    # transformation of y
    # cf. Eq 7
    beta <- m$Lambda.y[(ne + 1):ny,] 
    R <- cbind(-beta, diag(ny-ne))
    u <- y %*% t(R)
  } else {
      u <- 0
      # TODO: s.o
  }
  
  # Eqs 15, 16
  Sigma1 <- m$Phi - m$Phi %*% t(m$Lambda.x) %*% solve(m$Lambda.x %*% 
    m$Phi %*% t(m$Lambda.x) + m$Theta.d) %*% m$Lambda.x %*% m$Phi 
  L1 <- m$Phi %*% t(m$Lambda.x) %*% solve(m$Lambda.x %*% m$Phi %*%
    t(m$Lambda.x) + m$Theta.d)

  # s. Eq 17,33 
  Beta <- m$Beta - diag(model$info$num.eta)
  B <- solve(diag(ne) - Beta)
  psi.b  <- B %*% m$Psi %*% t(B)

  # s. Eq 34
  Sigma2 <- m$Theta.e[(1:ne),(1:ne)] - m$Theta.e[(1:ne),(1:ne)] %*% 
    t(beta) %*% solve(R %*% m$Theta.e %*% t(R)) %*% beta %*%
    m$Theta.e[(1:ne),(1:ne)]
  # s. Eq 24
  L2 <- -m$Theta.e[(1:ne),(1:ne)] %*% t(beta) %*% solve(R %*% m$Theta.e %*% t(R))
  
  # cf. Eq 42
  # one covariance matrix per person (i in 1:N)
      
  mu.x  <- m$nu.x + m$Lambda.x %*% m$tau 
  N <- nrow(data)

  Sigma3 <- list()
  for(i in 1:N){
    a <- m$tau + L1 %*% (x[i,] - mu.x)
    Sigma3[[i]] <- sigma.xi(a, Sigma1)
  }
  
  # Eq 14
  # Cov(x), Cov(u)
  sigma.x <- m$Lambda.x %*% m$Phi %*% t(m$Lambda.x) + m$Theta.d
  sigma.u <- R %*% m$Theta.e %*% t(R)

  # Cov(x,u)
  p1 <- length(x[1,])
  p2 <- length(u[1,])
  sigma.xu <- matrix(0, p1+p2, p1+p2)
  sigma.xu[1:p1, 1:p1] <- sigma.x
  sigma.xu[(p1+1):(p1+p2), (p1+1):(p1+p2)] <- sigma.u
  
  # sigma.y1 is conditional given x and u
  # Eq 17
  # sigma.y1 <- diag((c(m$Gamma) + 2*x %*% t(L1) %*% m$Omega) %*% Sigma1
  #             %*% t(c(m$Gamma) + 2*x %*% t(L1) %*% m$Omega)) + Sigma2 + 
  #             Sigma3
  # --> original equation from paper: faulty!
  
  ga1.b   <- B %*% m$Gamma
  if (is.matrix(m$Omega)) {
    Gamma2 <- t(vech(m$Omega))
  } else {
    Gamma2 <- t(apply(m$Omega, 3, vech))
  }
  ga2.b   <- B %*% Gamma2       # contains interaction and quadratic effects
  
  sigma.y1 <- list()  # TODO more efficient version? Transform to vector?
  for(i in 1:N){
    sigma.y1[[i]] <- 
      ga1.b %*% Sigma1 %*% t(ga1.b) +
      Sigma2 +
      psi.b +
      ga2.b %*% Sigma3[[i]] %*% t(ga2.b)
      
  }
  
  sigma.xy <- list(sigma.xu, sigma.y1)

  sigma.xy
}

loglikelihood_qml <- function(parameters, model, data) {
    
  mod.filled <- fill_model(model = model, parameters = parameters)

  # extract x and y from data frame
  x <- data[, 1:model$info$num.x]
  y <- data[, (model$info$num.x + 1):dim(data)[2]]
     
  mean.qml <- mu_qml(model = mod.filled, data=data)
  sigma.qml  <- sigma_qml(model = mod.filled, data=data)
  
  if (model$info$num.y > 1) {
    # transformation of y
    beta <- mod.filled$matrices$class1$Lambda.y[(model$info$num.eta + 1):model$info$num.y,] 
    R <- cbind(-beta, diag(model$info$num.y-model$info$num.eta))
    u <- y %*% t(R)
  } else {
    u <- 0
    # TODO: s.o.
  }

  N <- nrow(data)

  # Eq 10: densities
  f2 <- dmvnorm(cbind(x, u), mean=mean.qml[[1]], sigma=sigma.qml[[1]],
    log=TRUE)
  # original implementation: produces NaN when sds are negative
  
  lls <- NULL
  for (i in seq_len(N)) {
    f3 <- dmvnorm(as.matrix(y[,1:model$info$num.eta]), mean=mean.qml[[2]][,i],
      sigma=sigma.qml[[2]][[i]], log=TRUE)
    lls <- c(lls, sum(f2 + f3))
  }
  
  res <- sum(lls)
  
  return(-res)
}

mstep_qml <- function(model, parameters, data, neg.hessian=FALSE,
                      optimizer=c("nlminb", "optim"), max.iter=1,
                      control=list(), ...) {

  # optimizer
  optimizer <- match.arg(optimizer)

  if (optimizer == "nlminb") {

    if (is.null(control$iter.max)) {
        control$iter.max <- max.iter
    } else warning("iter.max is set for nlminb. max.iter will be ignored.")

    est <- nlminb(start=parameters, objective=loglikelihood_qml,
                  data=data, model=model,
                  upper=model$info$bounds$upper,
                  lower=model$info$bounds$lower, control=control, ...)
  } else {

    if (is.null(control$maxit)){
        control$maxit <- max.iter
    } else warning("maxit is set for optim. max.iter will be ignored.")

    est <- optim(par=parameters, fn=loglikelihood_qml, data=data,
                 model=model, upper=model$info$bounds$upper,
                 lower=model$info$bounds$lower, method="L-BFGS-B",
                 control=control, ...) 
    # fit est to nlminb output
    names(est) <- gsub("value", "objective", names(est))
    names(est) <- gsub("counts", "iterations", names(est))
  }

  if (neg.hessian == TRUE){
      est$hessian <- fdHess(pars=est$par, fun=loglikelihood_qml,
                            model=model, data=data)$Hessian
  }
  est
}


#--------------- helper functions ---------------

var.z <- function(Omega, Sigma1){

  ds <- dim(Sigma1)[1]
  varz <- 0
  
  # Eq. 18
  for(i in 1:ds){
    for(j in 1:ds){
      for(k in 1:ds){
        for(s in 1:ds){
          varzij <- Omega[i,j]*Omega[k,s]*(Sigma1[i,j]*Sigma1[k,s] +
                    Sigma1[i,k]*Sigma1[j,s] + Sigma1[i,s]*Sigma1[j,k])
          varz <- varz + varzij
        }
      }
    }
  }
  varz <- varz - sum(diag(Omega %*% Sigma1))^2
  varz
}


vech <- function (x, diag=TRUE) {
  x[upper.tri(x, diag=diag)]
}

sigma.xi <- function(a,sigma){
  # cov(vech(xixi')|x,u)
  # input: covariance matrix of z
  # output: covariance matrix of vech(az+za+zz==xi)
  dimz  <- nrow(sigma)
  dimzz <- length(vech(sigma))
  sigmazz <- matrix(0, dimzz, dimzz)
  asigmaz <- matrix(0, dimzz, dimzz)
  
  p <- q <- 1
  for(i in 1:dimz){
    for(j in i:dimz){
      for(k in 1:dimz){
        for(l in k:dimz){
          sigmazz[p,q] <- sigma[i,k]*sigma[j,l]+sigma[i,l]*sigma[j,k]
          
          asigmaz[p,q] <- 
            a[i]*a[k]*sigma[j,l] +
            a[i]*a[l]*sigma[j,k] +
            a[j]*a[k]*sigma[i,l] +
            a[j]*a[l]*sigma[i,k]
          
          q <- q+1
          if(q>dimzz){p <- p + 1
                      q <- 1}
          
        }
      }
    }
  }
  
  sigma.xixi <- asigmaz + sigmazz
  sigma.xixi
  
}


