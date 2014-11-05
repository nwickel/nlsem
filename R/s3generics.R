# s3generics.R
#
# created Nov/03/2014, KN
# last mod Nov/03/2014, KN

as.data.frame.lms <- as.data.frame.stemm <- as.data.frame.nsemm <- function(object, ...) {
    data <- data.frame(
        label = names(unlist(object$matrices$group1)))
    for (g in seq_len(length(object$matrices))) {
        temp <- data.frame(unlist(object$matrices[[g]], use.names=FALSE))
        names(temp) <- paste0("group", g)
        data <- cbind(data, temp)
    }
    data
}

simulate.stemmFilled <- function(object, nsim=1, seed=NULL, n=400, ...) {

    # set seed
    set.seed(seed)

    num.groups <- object$info$num.groups
    w <- object$info$w

    dat.sim <- lapply(1:num.groups, function(g) {
                      rmvnorm(round(n*w[g]),
                              mean=mu_stemm(model=object, g),
                              sigma=sigma_stemm(model=object, g))
                              })
    Reduce(rbind, dat.sim)
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

summary.emEst <- function(object, ...){

    # estimates
    est <- object$par

    # calculate Phi
    if (object$model.class == "lms") {
        A <- matrix(0, nrow=object$info$num.xi, ncol=object$info$num.xi)
        A[lower.tri(A, diag=TRUE)] <- est[grep("A", names(est))]
        Phi <- A %*% t(A)
        est[grep("A", names(est))] <- Phi[lower.tri(Phi, diag=TRUE)] 
        names(est)[grep("A", names(est))] <- paste0("Phi", 1:sum(lower.tri(Phi, diag=TRUE)))
    }

    # standard errors
    s.error <- sqrt(diag(solve(object$Hessian)))
    tvalue <- est / s.error
    pvalue <- 2 * pnorm(-abs(tvalue))
    # TODO Do I have to take it twice here? Is that because we only test if
    # different from 0 and do not care about direction?
    est.table <- cbind(est, s.error, tvalue, pvalue)
    dimnames(est.table)  <- list(names(est), c("Estimate", "Std. Error", "t value", "Pr(>|z|)"))

    # loglikelihoods
    iterations   <- length(object$loglikelihoods) 
    abs.change   <- c(0, diff(object$loglikelihoods))
    rel.change   <- rel_change(object$loglikelihoods)
    logLik.table <- cbind(object$loglikelihoods, abs.change, rel.change)
    dimnames(logLik.table) <- list(1:iterations, c("loglikelihood", "difference", "relative change"))

    ans <- list(estimates=est.table,
                iterations=iterations,
                finallogLik=object$objective,
                logLikelihoods=logLik.table)

    if (object$model.class == "stemm" || object$model.class == "nsemm") {
        ans$group.weights <- object$info$w
    }

    class(ans) <- "summary.emEst"

    ans
}

print.summary.emEst <- function(x, digits=max(3, getOption("digits") - 3),
                               cs.ind=2:3, ...){

    cat("\nEstimates:\n")
    printCoefmat(x$estimates, digits=digits, cs.ind=cs.ind, ...)
    cat("\nNumber of iterations:", x$iterations, "\nFinal loglikelihood:",
        format(x$finallogLik, digits=digits), "\n\n")
    cat("\nLikelihoods:\n")
    printCoefmat(x$logLikelihoods, digits=6, cs.ind=2:3, ...)

}

logLik.emEst <- function(object, ...){
    if(length(list(...)))
        warning("extra arguments discarded")

    out <- object$objective
    attr(out, "df") <- length(object$par)
    class(out) <- "logLik"
    out
}

anova.emEst <- function(object, ..., test=c("Chisq", "none")){

    test <- match.arg(test)
    dots <- list(...)
    if (length(dots) == 0)
        stop('anova is not implemented for a single "emEst" object')

    mlist <- list(object, ...)
    names(mlist) <- c(deparse(substitute(object)),
                as.character(substitute(...[]))[2:length(mlist)])
    if (any(!sapply(mlist, inherits, "emEst")))
        stop('not all objects are of class "emEst"')
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

AIC.emEst <- function(object, ..., k=2){

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

coef.emEst <- coefficients.emEst <- function(object, ...){

    coef <- object$par

    # calculate Phi
    if (object$model.class == "lms") {
        A <- matrix(0, nrow=object$info$num.xi, ncol=object$info$num.xi)
        A[lower.tri(A, diag=TRUE)] <- coef[grep("A", names(coef))]
        Phi <- A %*% t(A)
        coef[grep("A", names(coef))] <- Phi[lower.tri(Phi, diag=TRUE)] 
        names(coef)[grep("A", names(coef))] <- paste0("Phi", 1:sum(lower.tri(Phi, diag=TRUE)))
    }
    coef
}

plot.emEst <- function(x, y, ...){

    plot(x$loglikelihoods, type="l", xlab="Number of iterations", 
         ylab="log likelihood", axes=F)
    axis(1, at=1:length(x$loglikelihoods))
    axis(2)
    box()

}
