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
    
    k <- get_k(mod.filled$matrices$group1$Omega)
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
    lower[grep("[Theta.d, Theta.e, Psi]", model$info$par.names)] <- 0
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

simulate.lmsFilled <- function(object, nsim=1, seed=NULL, n=400, m=16, ...){
    
    # set seed
    set.seed(seed)
      
    # Gauss-Hermite quadrature
    k <- get_k(object$matrices$group1$Omega)
    quad <- quadrature(m=m, k=k)
    V <- quad$n
    w <- quad$w

    # calculate n for mixture components
    n.mix <- ceiling(w*n)
    
    # simulate data (compare Equation 15 in Klein & Moosbrugger (2000))
    dat.sim <- sapply(1:length(w), function(i){
                           rmvnorm(n.mix[i], 
                           mean=mu_lms(model=object, z=V[i,]),
                           sigma=sigma_lms(model=object, z=V[i,]))
                           })

    dat <- Reduce(rbind, dat.sim)[sample(n),]         # ceiling "inflates" number of observations 
    # TODO This seems random --> better solution?
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
    attr(out, "heading") <- "Analysis of deviance table for SEMs\n"
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

coef.emRes <- coefficients.emRes <- function(object, ...){

    coef <- object$par

    # calculate Phi
    A <- matrix(0, nrow=object$info$num.xi, ncol=object$info$num.xi)
    A[lower.tri(A, diag=TRUE)] <- coef[grep("A", names(coef))]
    Phi <- A %*% t(A)
    coef[grep("A", names(coef))] <- Phi[lower.tri(Phi, diag=TRUE)] 
    names(coef)[grep("A", names(coef))] <- paste0("Phi", 1:sum(lower.tri(Phi, diag=TRUE)))

    coef
}


