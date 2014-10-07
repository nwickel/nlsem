# lms.R
#
# created: Sep/11/2014, NU
# last mod: Oct/10/2014, NU

mu_lms <- function(model, z) {
    
    stopifnot(class(model) == "lmsFilled")

    matrices <- model$matrices
    k    <- get_k(matrices$O)         # number of nonzero rows in Omega
    n    <- nrow(matrices$A)          # number of zero rows in Omega
    z.1  <- c(z, rep(0, n - k))       # [z_1 0]'
    A.z  <- matrices$A %*% z.1 
    mu.x <- matrices$vx + matrices$Lx %*% A.z 
    mu.y <- matrices$vy + matrices$Ly %*% (matrices$alpha + matrices$G %*% A.z +
            t(A.z) %*% matrices$O %*% A.z)
    mu   <- c(mu.x, mu.y)

    mu
}

sigma_lms <- function(model, z) {

    stopifnot(class(model) == "lmsFilled")

    matrices <- model$matrices
    k     <- get_k(matrices$O)        # number of nonzero rows in Omega
    n     <- nrow(matrices$A)         # number of zero rows in Omega
    z.1   <- c(z, rep(0, n - k))      # [z_1 0]'
    A.z   <- matrices$A %*% z.1 
    d.mat <- get_d(n=n, k=k)
    Lx.A  <- matrices$Lx %*% matrices$A
    temp  <- matrices$G %*% matrices$A + t(A.z) %*% matrices$O %*% matrices$A
    s11   <- Lx.A %*% d.mat %*% t(Lx.A) + matrices$Td 
    s12   <- Lx.A %*% d.mat %*% t(temp) %*% t(matrices$Ly)
    s21   <- t(s12)
    s22   <- matrices$Ly %*% temp %*% d.mat %*% t(temp) %*% t(matrices$Ly) + 
             matrices$Ly %*% matrices$Psi %*% t(matrices$Ly) + matrices$Te
    sigma <- rbind(cbind(s11,s12), cbind(s21,s22))
    
    # check if sigma is symmetric
    if(!isSymmetric(sigma)) stop("Sigma has to be symmetric")

    sigma
}

get_k <- function(omega) which(rowSums(omega) == 0)[1] - 1

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


estep_lms <- function(model, parameters, dat, m, k, ...) {

    stopifnot(free_parameters(model) == length(parameters))

    mod.filled <- fill_model(model=model, parameters=parameters)

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
    # ??? What is this for? What are we trying to check?
    # Should all P's be near 1? Why?
    # Anyway, this should get a meaningful error message TODO

    P
}

loglikelihood <- function(model, parameters, data, m, k){

    quad <- quadrature(m, k)
    V <- quad$n

    P <- estep_lms(model, parameters, data, m, k)

    mod.filled <- fill_model(model=model, parameters=parameters)

    res <- 0
    for(node.num in 1:m) {
        v.par    <- V[node.num,]
        lls      <- dmvnorm(data, mean=mu_lms(model=mod.filled, z=v.par),
                            sigma=sigma_lms(model=mod.filled, z=v.par),
                            log=TRUE) * P[,node.num]
        res      <- res + sum(lls)
    }
    res
}
# FIX ME: I do not give the correct minimum, since I calculate a different
# P... (or because of something else...)

mstep_lms <- function(model, parameters, dat, P, m, k, Hessian=FALSE, ...) {

    # likelihood function which will be optimized
    fun <- function(parameters, dat, model, P) {
        
        quad <- quadrature(m, k)
        V <- quad$n
  
        mod.filled <- fill_model(model=model, parameters=parameters)

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
    upper <- rep(Inf, free_parameters(model))
    lower <- rep(-Inf, free_parameters(model))
    lower[grep("[Td, Te, Psi]", model$info$par.names)] <- 0
    # TODO What about A? Does that have to be positiv as well??? Since Phi
    # should be...
    # TODO What about if Td and Te are not diagonal matrices?
    # TODO What about Psi, when eta > 1?
    # optimizer
    est <- nlminb(start=parameters, objective=fun, dat=dat, model=model,
                  P=P, upper=upper, lower=lower, ...)

    # TODO Try out constrOptim and compare results...
    # ??? How to use boundaries in constrOptim?

    if (Hessian == TRUE){
        est$hessian <- nlme::fdHess(pars=est$par, fun=fun, model=model,
                                    dat=dat, P=P)
    }
    # TODO Maybe add analytical solution for Hessian matrix?

    est
}


em <- function(model, data, start, logger=FALSE, threshold=10e-7,
                max.iter=40, m=16, k=1, ...) {

    cat("-----------------------------------\n")
    cat("Starting EM-algorithm\n")
    cat(paste("Threshold: ", threshold, "\n"))
    cat("-----------------------------------\n")
    cat("-----------------------------------\n")
    
    ll.old   <- 0     # likelihood of the last iteration
    ll.new   <- 1     # likelihood of the current iteration
    ll.ret   <- NULL
    num.iter <- 0     # number of iterations
    
    while(abs(ll.old - ll.new) > threshold) { # as long as no convergence reached
        if(ll.new - ll.old > 0.001 && num.iter > 3) {
            warning("Likelihood should be decreasing")
        }
      
        if(logger == TRUE) {
            cat(paste("Iteration", num.iter+1, "\n"))
            cat("Doing expectation-step \n")
        }
  
        # Update likelihood
        ll.old <- ll.new

        # E-step
        P      <- estep_lms(model=model, parameters=start, dat=data, m=m, k=k, ...)
  
        if(logger == TRUE){
            cat("Doing maximization-step \n")
        }
        
        # M-step
        m.step <- mstep_lms(model=model, P=P, dat=data, parameters=start, m=m, k=k, ...)
  
        if(logger == TRUE) {
            cat("Results of maximization \n")
            cat(paste0("Final likelihood: ", m.step$objective, "\n"))
            cat(paste0("Convergence: ", m.step$convergence, "\n"))
            cat(paste0("Number of iterations: ", m.step$iterations, "\n"))
            cat("----------------------------------- \n")
        }
      
        ll.new     <- m.step$objective
        ll.ret     <- c(ll.ret, ll.new)
        start      <- m.step$par
        num.iter   <- num.iter + 1
  
        if(num.iter == max.iter) break
    }

    cat("-----------------------------------\n")
    cat("EM completed \n")
    cat(paste0("Previous Likelihood: ", ll.old, "\n"))
    cat(paste0("Final Likelihood: ", ll.new,"\n"))
    cat("-----------------------------------\n")
    cat("-----------------------------------\n")
    cat("Computing Hessian \n")
    cat("-----------------------------------\n")
    
    # Compute hessian of final parameters
    final <- mstep_lms(model=model, P=P, dat=data, parameters=start, Hessian=TRUE, m=m, k=k, ...)

    names(final$par) <- model$info$par.names
    
    out <- list(par=final$par, objective=-final$objective,
           convergence=final$convergence, message=final$message,
           Hessian=final$hessian$Hessian, gradient=final$hessian$gradient,
           likelihoods=-ll.ret, info=model$info[1:4])
  
    class(out) <- "emRes"
  
    out
}

simulate.lmsFilled <- function(object, nsim=1, seed=NULL, n=400, m=16, k=1, ...){
    
    # set seed
    set.seed(seed)
    
    # Gauss-Hermite quadrature
    quad <- quadrature(m=m, k=k)
    V <- quad$n
    w <- quad$w

    # calculate n for mixture components
    n.mix <- ceiling(w*n)
    
    # simulate data (compare Equation 15 in Klein & Moosbrugger (2000)
    dat.sim <- sapply(1:m, function(i){
                           rmvnorm(n.mix[i], 
                           mean=mu_lms(model=object, z=V[i,]),
                           sigma=sigma_lms(model=object, z=V[i,]))
                           })

    dat <- Reduce(rbind, dat.sim)[1:n,]         # ceiling "inflates" number of observations 
    dat
}


rel_change <- function(x){
    
    rel.change <- numeric(length(x))
    for (i in 2:length(rel.change)){
        
        rel.change[i] <- abs(x[i-1]-x[i])/max(abs(x[i-1]), abs(x[i]))
    }
    rel.change
}

summary.emRes <- function(object, ...){
    
    # estimates
    est <- object$par
    
    # calculate Phi
    A <- matrix(0, nrow=object$info$num.xi, ncol=object$info$num.xi)
    A[lower.tri(A, diag=TRUE)] <- est[grep("A", names(est))]
    Phi <- A %*% t(A)
    est[grep("A", names(est))] <- Phi[lower.tri(Phi, diag=TRUE)] 
    names(est)[grep("A", names(est))] <- paste0("Phi", 1:sum(lower.tri(Phi, diag=TRUE)))

    # standard errors
    s.error <- sqrt(diag(solve(object$Hessian)))
    # TODO Do I have to take the neg. Hessian? Or do I already have
    # the negativ Hessian? Covariance matrix is the inverse of the negative
    # Hessian, isn't it?
    tvalue <- est / s.error
    pvalue <- 2 * pnorm(-abs(tvalue))
    # TODO Do I have to take it twice here? Is that because we only test if
    # different from 0 and do not care about direction?
    est.table <- cbind(est, s.error, tvalue, pvalue)
    dimnames(est.table)  <- list(names(est), c("Estimate", "Std. Error", "t value", "Pr(>|z|)"))

    # likelihoods
    iterations   <- length(object$likelihoods) 
    abs.change   <- c(0, abs(diff(object$likelihoods)))
    rel.change   <- rel_change(object$likelihoods)
    logLik.table <- cbind(object$likelihoods, abs.change, rel.change)
    dimnames(logLik.table) <- list(1:iterations, c("loglikelihood", "absolute change", "relative change"))

    ans <- list(estimates=est.table,
                iterations=iterations,
                finallogLik=object$objective,
                logLikelihoods=logLik.table)

    class(ans) <- "summary.emRes"

    ans
}

print.summary.emRes <- function(x, digits=max(3, getOption("digits") - 3),
                               cs.ind=2:3, ...){

    cat("\nEstimates:\n")
    printCoefmat(x$estimates, digits=digits, cs.ind=cs.ind, ...)
    cat("\nNumber of iterations:", x$iterations, "\nFinal loglikelihood:",
        format(x$finallogLik, digits=digits), "\n\n")
    cat("\nLikelihoods:\n")
    printCoefmat(x$logLikelihoods, digits=6, cs.ind=2:3, ...)

}

logLik.emRes <- function(object, ...){
    if(length(list(...)))
        warning("extra arguments discarded")

    out <- object$objective
    attr(out, "df") <- length(object$par)
    class(out) <- "logLik"
    out
}

anova.emRes <- function(object, ..., test=c("Chisq", "none")){

    test <- match.arg(test)
    dots <- list(...)
    if (length(dots) == 0)
        stop('anova is not implemented for a single "emRes" object')
    
    mlist <- list(object, ...)
    names(mlist) <- c(deparse(substitute(object)),
                as.character(substitute(...[]))[2:length(mlist)])
    if (any(!sapply(mlist, inherits, "emRes")))
        stop('not all objects are of class "emRes"')
    nt <- length(mlist)

    dflist <- sapply(mlist, function(x) length(x$par))
    
    s <- order(dflist, decreasing=TRUE)
    mlist <- mlist[s]
    dflist <- dflist[s]

    lls <- sapply(mlist, function(x) -2*x$objective)
    tss <- c("", paste(1:(nt - 1), 2:nt, sep = " vs "))
    df <- c(NA, -diff(dflist))
    x2 <- c(NA, -diff(lls))
    pr <- c(NA, 1 - pchisq(x2[-1], df[-1]))
    out <- data.frame(Model=names(mlist), Resid.df=dflist, Deviance=lls,
                      Test=tss, Df=df, LRtest=x2, Prob=pr)
    names(out) <- c("Model", "Resid. df", "Resid. Dev", "Test",
                    "   Df", "LR stat.", "Pr(>Chi)")
    rownames(out) <- 1:nt
    if (test == "none") out <- out[, 1:6]
    class(out) <- c("Anova", "data.frame")
    attr(out, "heading") <- "Analysis of deviance table for LMS models\n"
    # TODO LMS might be misleading here!!!
    # FIX ME Header does not show.
    out
}

AIC.emRes <- function(object, ..., k=2){
    
    dots <- list(...)
    if (length(dots) == 0){
        out <- as.numeric(-2*logLik(object) + k*length(object$par))
    } else {
        mlist <- list(object, ...)
        names(mlist) <- c(deparse(substitute(object)),
                          as.character(substitute(...[]))[2:length(mlist)])
        nt <- length(mlist)

        dflist <- sapply(mlist, function(x) length(x$par))
        
        aic <- sapply(mlist, function(x) -2*logLik(x) + k*length(x$par))
        s <- order(aic, decreasing=TRUE)
        aic <- aic[s]
        dflist <- dflist[s]
        out <- data.frame(df=dflist, AIC=aic)
        rownames(out) <- names(mlist)
    }
    out
}


