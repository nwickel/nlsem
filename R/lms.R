# lms.R
#
# created: Sep/11/2014, NU
# last mod: Oct/07/2014, NU

mu_lms <- function(model, z) {
    
    stopifnot(class(model) == "lmsFilled")

    matrices <- model$matrices$group1
    k    <- get_k(matrices$Omega)     # number of nonzero rows in Omega
    n    <- nrow(matrices$A)          # number of zero rows in Omega
    z.1  <- c(z, rep(0, n - k))       # [z_1 0]'
    A.z  <- matrices$A %*% z.1 
    mu.x <- matrices$nu.x + matrices$Lambda.x %*% A.z 
    mu.y <- matrices$nu.y + matrices$Lambda.y %*% (matrices$alpha +
            matrices$Gamma %*% A.z + t(A.z) %*% matrices$Omega %*% A.z)
    mu   <- c(mu.x, mu.y)

    mu
}

sigma_lms <- function(model, z) {

    stopifnot(class(model) == "lmsFilled")

    matrices <- model$matrices$group1
    k     <- get_k(matrices$Omega)    # number of nonzero rows in Omega
    n     <- nrow(matrices$A)         # number of zero rows in Omega
    z.1   <- c(z, rep(0, n - k))      # [z_1 0]'
    A.z   <- matrices$A %*% z.1 
    d.mat <- get_d(n=n, k=k)
    Lx.A  <- matrices$Lambda.x %*% matrices$A
    temp  <- matrices$Gamma %*% matrices$A + t(A.z) %*% matrices$Omega %*% matrices$A
    s11   <- Lx.A %*% d.mat %*% t(Lx.A) + matrices$Theta.d 
    s12   <- Lx.A %*% d.mat %*% t(temp) %*% t(matrices$Lambda.y)
    s21   <- t(s12)
    s22   <- matrices$Lambda.y %*% temp %*% d.mat %*% t(temp) %*%
             t(matrices$Lambda.y) + matrices$Lambda.y %*% matrices$Psi %*%
             t(matrices$Lambda.y) + matrices$Theta.e
    sigma <- rbind(cbind(s11,s12), cbind(s21,s22))
    
    # check if sigma is symmetric
    if(!isSymmetric(sigma)) stop("Sigma has to be symmetric")

    sigma
}

get_k <- function(Omega) which(rowSums(Omega) == 0)[1] - 1

get_d <- function(n, k) {
    mat <- diag(n)
    mat[1:k, 1:k] <- 0
    mat
}

quadrature <- function(m, k) {

    one.dim         <- hermite.h.quadrature.rules(m)[[m]]
    test            <- as.matrix(expand.grid(lapply(vector("list", k), function(x) {x <- 1:m; x})))
    final.nodes     <- matrix(one.dim$x[test], ncol=k, byrow=FALSE)
    permute.weights <- matrix(one.dim$w[test], ncol=k, byrow=FALSE)
    final.weights   <- apply(permute.weights, 1, prod)

    n               <- final.nodes * sqrt(2)
    w               <- final.weights * pi^(-k/2)
  
    out <- list(n=n, w=w, k=k, m=m)
    out
}


estep_lms <- function(model, parameters, dat, m, ...) {

    stopifnot(count_free_parameters(model) == length(parameters))

    mod.filled <- fill_model(model=model, parameters=parameters)

    k <- get_k(mod.filled$matrices$group1$Omega)
    quad <- quadrature(m, k)

    V <- quad$n       # matrix of node vectors m x k
    w <- quad$w       # weights

    stopifnot(sum(w) - 1 < 1e-5)

    P <- NULL
    for(node.num in 1:m) {
        v.par <- V[node.num, ]
        rho   <- w[node.num]
        p.ij  <- rho * dmvnorm(dat, mean=mu_lms(model=mod.filled, z=v.par), 
                               sigma=sigma_lms(model=mod.filled, z=v.par))
        P     <- cbind(P, p.ij)
    }

    P <- P / rowSums(P)   # divide each rho_j*phi(x_i, y_i) by whole density (row)

    stopifnot(all.equal(rowSums(P), rep(1, nrow(P))))
    # TODO Should get a meaningful error message

    P
}

# log likelihood function which will be optimized
loglikelihood <- function(model, parameters, dat, P, m=16, ...) {
    
    mod.filled <- fill_model(model=model, parameters=parameters)

    k <- get_k(mod.filled$matrices$group1$Omega)
    quad <- quadrature(m, k)
    V <- quad$n

    res <- 0
    for(node.num in 1:m) {
        v.par    <- V[node.num,]
        lls      <- dmvnorm(dat, mean=mu_lms(model=mod.filled, z=v.par),
                            sigma=sigma_lms(model=mod.filled, z=v.par),
                            log=TRUE) * P[,node.num]
        res      <- res + sum(lls)
    }

    -res
}


mstep_lms <- function(model, parameters, dat, P, m, Hessian=FALSE, ...) {


     # # optimizer
     # if (Hessian == FALSE){
     #     est <- nlminb(start=parameters, objective=fun, dat=dat, model=model, P=P, ...)
     #     out <- est
     # } else {
     #     est <- optim(par=parameters, fn=fun, dat=dat, model=model, P=P, hessian=TRUE, method="L-BFGS-B", ...)
     #     out <- list(par=est$par, objective=est$value,
     #                 convergence=est$convergence, evaluations=est$counts,
     #                 message=est$message, hessian=est$hessian) }
     ## --> Alternative calculation of Hessian (does not require nlme)
 

    # constrain variances to +Inf
    upper <- rep(Inf, count_free_parameters(model))
    lower <- rep(-Inf, count_free_parameters(model))
    lower[grep("Theta.[de]", model$info$par.names)] <- 0
    lower[grep("Psi", model$info$par.names)] <- 0
    # TODO What about A? Does that have to be positiv as well??? Since Phi
    # should be...
    # TODO What about if Theta.d and Theta.e are not diagonal matrices?
    # TODO What about Psi, when eta > 1?
    # optimizer
    est <- nlminb(start=parameters, objective=loglikelihood, dat=dat,
                  model=model, P=P, upper=upper, lower=lower, ...)

    # TODO Try out constrOptim and compare results...
    # ??? How to use boundaries in constrOptim?

    if (Hessian == TRUE){
        est$hessian <- nlme::fdHess(pars=est$par, fun=loglikelihood, model=model,
                                    dat=dat, P=P)
    }
    # TODO Maybe add analytical solution for Hessian matrix?

    est
}

