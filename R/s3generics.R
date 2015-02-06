# s3generics.R
#
# created Nov/03/2014, KN
# last mod Feb/06/2015, KN

#--------------- main functions ---------------

as.data.frame.lms <- as.data.frame.semm <- as.data.frame.nsemm <- function(x, ...) {
    data <- data.frame(
        label = names(unlist(x$matrices$class1)))
    for (c in seq_len(length(x$matrices))) {
        temp <- data.frame(unlist(x$matrices[[c]], use.names=FALSE))
        names(temp) <- paste0("class", c)
        data <- cbind(data, temp)
    }
    data
}

simulate.nsemm <- function(object, nsim=1, seed=NULL, n=400, m=16, parameters, ...) {

    set.seed(seed)

    mod.filled <- fill_model(object, parameters)

    num.classes <- mod.filled$info$num.classes
    w <- mod.filled$info$w

    # simulate n data points for each mixture as lms
    dat.sim <- lapply(seq_len(num.classes), function(c) {
                      simulate(lms_ify(object, c),
                               nsim=nsim, seed=seed, n=n, m=m,
                               parameters=get_class_parameters(object, parameters)[[c]],
                               ...)

    })

    # see simulate_lms for explanation
    border <- cumsum(w)
    prob <- runif(n)

    dat <- NULL
    for (i in seq_len(n)) {
        ind <- sum(prob[i] > border) + 1
        dat <- rbind(dat, dat.sim[[ind]][i,])
    }

    dat
}

simulate.semm <- function(object, nsim=1, seed=NULL, n=400, parameters, ...) {

    set.seed(seed)

    mod.filled <- fill_model(object, parameters)

    num.classes <- mod.filled$info$num.classes
    w <- mod.filled$info$w

    # simulate n data points from each mixture distribution
    dat.sim <- lapply(1:num.classes, function(c) {
                      rmvnorm(n,
                              mean=mu_semm(matrices=mod.filled$matrices[[c]]),
                              sigma=sigma_semm(matrices=mod.filled$matrices[[c]]))
                              })

    # see simulate.lms for explanation
    border <- cumsum(w)
    prob <- runif(n)

    dat <- NULL
    for (i in seq_len(n)){
        ind <- sum(prob[i] > border) + 1
        dat <- rbind(dat, dat.sim[[ind]][i,])
    }

    dat
}

simulate.lms <- function(object, nsim=1, seed=NULL, n=400, m=16, parameters, ...) {

    # set seed
    set.seed(seed)

    # Gauss-Hermite quadrature
    k <- get_k(object$matrices$class1$Omega)
    quad <- quadrature(m=m, k=k)
    V <- quad$n
    w <- quad$w
    
    parameters <- convert_parameters_lms(object, parameters)
    names(object$matrices$class1)[grep("Phi", names(object$matrices$class1))] <- "A"

    mod.filled <- fill_model(object, parameters)

    # simulate n data points from each mixture distribution
    dat.sim <- lapply(seq_along(w), function(i){
                           rmvnorm(n, 
                           mean=mu_lms(model=mod.filled, z=V[i,]),
                           sigma=sigma_lms(model=mod.filled, z=V[i,]))
                           })

    # decide which data points from each mixture should be included in
    # simulated data set: weights give intervall borders between 0 and 1;
    # we draw random numbers from a uniform distribution and check in what
    # intervall they lie: the ith element from that distribution will be
    # taken and put into the data frame
    border <- cumsum(w)
    prob <- runif(n)

    dat <- NULL
    for (i in seq_len(n)){
        ind <- sum(prob[i] > border) + 1
        dat <- rbind(dat, dat.sim[[ind]][i,])
    }

    dat
}


summary.emEst <- function(object, ...) {

    # estimates
    est <- object$coefficients

    # standard errors
    if (is.numeric(est)) {
        s.error <- calc_standard_error(object$neg.hessian)
        tvalue <- est / s.error
        pvalue <- 2 * pnorm(-abs(tvalue))
        est.table <- cbind(est, s.error, tvalue, pvalue)
        dimnames(est.table)  <- list(names(est), c("Estimate", "Std. Error", "t value", "Pr(>|z|)"))
    } else {
        est.table <- Reduce('rbind', lapply(seq_along(est), function(c) {
                            s.error <- calc_standard_error(object$neg.hessian[[c]])
                            tvalue <- est[[c]] / s.error
                            pvalue <- 2 * pnorm(-abs(tvalue))
                            est.table <- cbind(est[[c]], s.error, tvalue, pvalue)
                            dimnames(est.table)  <- list(paste0("class", c, ".", names(est[[c]])),
                                                         c("Estimate", "Std. Error",
                                                           "t value", "Pr(>|z|)"))
                            est.table
                           }))
    }

    # loglikelihoods
    iterations   <- length(object$loglikelihoods) 
    abs.change   <- c(0, diff(object$loglikelihoods))
    rel.change   <- rel_change(object$loglikelihoods)
    logLik.table <- cbind(object$loglikelihoods, abs.change, rel.change)
    dimnames(logLik.table) <- list(1:iterations, c("loglikelihood", "difference", "relative change"))

    ans <- list(model=object$model.class,
                estimates=est.table,
                iterations=iterations,
                finallogLik=object$objective,
                logLikelihoods=logLik.table)

    if (object$model.class == "semm" || object$model.class == "nsemm") {
        ans$class.weights <- object$info$w
    }

    class(ans) <- "summary.emEst"

    ans
}

print.summary.emEst <- function(x, digits=max(3, getOption("digits") - 3),
                               cs.ind=2:3, ...) {
    
    cat("\nSummary for model of class", x$model, "\n")
    cat("\nEstimates:\n")
    printCoefmat(x$estimates, digits=digits, cs.ind=cs.ind, ...)
    cat("\nNumber of iterations:", x$iterations,
        "\nFinal loglikelihood:", round(x$finallogLik, 3)) 
    if (x$model == "semm" || x$model == "nsemm"){
        cat("\nFinal weights:", round(x$class.weights, 3))
    }
    cat("\n\n", "\nLikelihoods:\n")
    printCoefmat(x$logLikelihoods, digits=6, cs.ind=2:3, ...)

}

logLik.emEst <- function(object, ...){
    if(length(list(...)))
        warning("extra arguments discarded")

    out <- object$objective
    attr(out, "df") <- length(unlist(object$coef))
    class(out) <- "logLik"
    out
}

anova.emEst <- function(object, ..., test=c("Chisq", "none")) {
    # Adapted from anova.polr by Brian Ripley
    
    test <- match.arg(test)
    dots <- list(...)
    if (length(dots) == 0)
        stop('anova is not implemented for a single "emEst" object')

    mlist <- list(object, ...)
    if (any(!sapply(mlist, function(x) x$model.class == "lms"))) {
        stop('Likelihood Ratio Test only meaningful for models of class "lms".')
    }

    nlist <- sapply(mlist, function(x) x$info$n)
    if (!all(nlist == object$info$n)) {
        stop("SEM have not all been fitted to the same data set.")
    }

    names(mlist) <- c(deparse(substitute(object)),
                as.character(substitute(...[]))[2:length(mlist)])
    if (any(!sapply(mlist, inherits, "emEst")))
        stop('not all objects are of class "emEst"')
    nt <- length(mlist)

    dflist <- sapply(mlist, function(x) length(unlist(x$coef)))

    s <- order(dflist)
    mlist <- mlist[s]
    dflist <- dflist[s]

    lls <- sapply(mlist, function(x) -2*x$objective)
    tss <- c("", paste(1:(nt - 1), 2:nt, sep = " vs "))
    df <- c(NA, diff(dflist))
    x2 <- c(NA, -diff(lls))
    pr <- c(NA, 1 - pchisq(x2[-1], df[-1]))
    out <- data.frame(Model=names(mlist), Resid.df=dflist, Deviance=lls,
                      Test=tss, Df=df, LRtest=x2, Prob=pr)
    names(out) <- c("Model", "Numb. coef", "-2logLik", "Test",
                    "   Df", "LR stat.", "Pr(>Chi)")
    rownames(out) <- 1:nt
    if (test == "none") out <- out[, 1:6]
    class(out) <- c("Anova", "data.frame")
    attr(out, "heading") <- "Chi Square test statistic for (nonlinear) SEM\n"
    out
}

AIC.emEst <- function(object, ..., k=2) {

    dots <- list(...)
    if (length(dots) == 0){
        out <- as.numeric(-2*logLik(object) + k*length(object$coef))
    } else {
        mlist <- list(object, ...)
        names(mlist) <- c(deparse(substitute(object)),
                          as.character(substitute(...[]))[2:length(mlist)])
        nt <- length(mlist)

        dflist <- sapply(mlist, function(x) length(unlist(x$coef)))

        aic <- sapply(mlist, function(x) -2*logLik(x) + k*length(unlist(x$coef)))
        s <- order(aic, decreasing=TRUE)
        aic <- aic[s]
        dflist <- dflist[s]
        out <- data.frame(df=dflist, AIC=aic)
        rownames(out) <- names(mlist)
    }
    out
}

BIC.emEst <- function(object, ...) {

    dots <- list(...)
    if (length(dots) == 0){
        out <- as.numeric(-2*logLik(object) +
               log(object$info$n)*length(unlist(object$coef)))
    } else {
        mlist <- list(object, ...)
        names(mlist) <- c(deparse(substitute(object)),
                          as.character(substitute(...[]))[2:length(mlist)])
        nt <- length(mlist)

        dflist <- sapply(mlist, function(x) length(unlist(x$coef)))

        bic <- sapply(mlist, function(x) -2*logLik(x) + log(object$info$n)*length(unlist(x$coef)))
        s <- order(bic, decreasing=TRUE)
        bic <- bic[s]
        dflist <- dflist[s]
        out <- data.frame(df=dflist, BIC=bic)
        rownames(out) <- names(mlist)
    }
    out
}

plot.emEst <- function(x, y, ...) {

    plot(x$loglikelihoods, type="l", xlab="Number of iterations", 
         ylab="log likelihood", axes=F)
    axis(1, at=1:length(x$loglikelihoods))
    axis(2)
    box()

}

#--------------- helper functions ---------------

# calculates relative change defined as absolute difference divided by
# maximum absolute value
rel_change <- function(x) {

    rel.change <- numeric(length(x))
    for (i in 2:length(rel.change)){

        rel.change[i] <- abs(x[i-1]-x[i])/max(abs(x[i-1]), abs(x[i]))
    }
    rel.change
}

calc_standard_error <- function(neg.hessian) {
    warn.msg <- "Standard errors could not be computed, because negative Hessian was either not available or singular"
    s.error <- tryCatch({
                    sqrt(diag(solve(neg.hessian)))
                }, error=function(e) {
                    NA
                }, warning=function(w) {
                     if (grepl("NaN", conditionMessage(w))) {
                       suppressWarnings(sqrt(diag(solve(neg.hessian))))
                    } else{
                        sqrt(diag(solve(neg.hessian)))
                    }
                })
                if (all(is.na(s.error))) warning(warn.msg)
                if (any(is.nan(s.error))) warning("Standard errors for some coefficients could not be computed.") 
    s.error
}
