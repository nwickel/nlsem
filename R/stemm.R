# stemm.R
#
# created: Okt/20/2014, KN
# last mod: Dec/03/2014, NU

#--------------- main functions ---------------

# Calculate mu of multivariate normal distribution for joint vector of
# indicators (See equation 4 in Jedidi, Jagpal & DeSarbo, 1997).
# The order is (x, y) as opposed to the paper.
mu_stemm <- function(matrices) {

    check_filled(matrices)

    mu.y <- matrices$nu.y + matrices$Lambda.y %*% solve(matrices$Beta) %*%
            (matrices$alpha + matrices$Gamma %*% matrices$tau)
    mu.x <- matrices$nu.x + matrices$Lambda.x %*% matrices$tau
    mu <- rbind(mu.x, mu.y) # vertical vector

    mu
}

# Calculate sigma of multivariate normal distribution for joint vector of
# indicators (y, x). (See equation 5 in Jedidi, Jagpal & DeSarbo, 1997)
# The order is (x, y), as opposed to the paper. Therefore the rows and cols in
# sigma are switched.
sigma_stemm <- function(matrices) {

    check_filled(matrices)

    # Lambda.y * B^-1
    Ly.Binv <- matrices$Lambda.y %*% solve(matrices$Beta)

    s22 <- Ly.Binv %*% (matrices$Gamma %*% matrices$Phi %*% t(matrices$Gamma) +
                        matrices$Psi) %*% t(Ly.Binv) + matrices$Theta.e

    s21 <- Ly.Binv %*% matrices$Gamma %*% matrices$Phi %*% t(matrices$Lambda.x)
    s12 <- t(s21)
    s11 <- matrices$Lambda.x %*% matrices$Phi %*% t(matrices$Lambda.x) + matrices$Theta.d
    sigma <- rbind(cbind(s11,s12), cbind(s21, s22))

    # TODO check is this warning is really necessary
    if (!isSymmetric(sigma)) warning("Sigma was not symmetric. This is probably due to
                                     numerical calculation")
    tryCatch(solve(sigma), error = function(e) stop("Sigma is not nonsingular."))

    sigma
}

# Expectation step of the EM-algorithm (see Jedidi, Jagpal & DeSarbo, 1997)
estep_stemm <- function(model, parameters, data) {

    model.filled <- fill_model(model=model, parameters=parameters)

    P <- NULL
    for (c in seq_len(model$info$num.classes)) {
        # class weight
        w.c <- model$info$w[c]

        p.ij <- w.c * dmvnorm(data, mean=mu_stemm(model.filled$matrices[[c]]),
                              sigma=sigma_stemm(model.filled$matrices[[c]]))
        if (sum(p.ij) == 0) stop("Posterior probability could not be calculated
                                 properly. Choose different starting
                                 parameters.")
        P <- cbind(P, p.ij, deparse.level=0)
    }
    P <- P / rowSums(P)
    P
}

# negative log likelihood function which will be optimized in M-step (see below)
loglikelihood_stemm <- function(parameters, matrices, data, p, w) {
    # fill matrices
    for (i in seq_along(matrices)) {
        matrix.i <- matrices[[i]]
        # number of NA's in matrix
        num.na <- length(matrix.i[is.na(matrix.i)])
        if (num.na > 0) {
            matrix.i[is.na(matrix.i)] <- parameters[1:num.na]
            parameters <- parameters[-(1:num.na)]
            matrices[[i]] <- matrix.i
        }
    }
    matrices$Psi <- fill_symmetric(matrices$Psi)
    matrices$Phi <- fill_symmetric(matrices$Phi)

    N <- nrow(data)
    N.c <- sum(p)
    mu <- mu_stemm(matrices)
    sigma <- sigma_stemm(matrices)
    T.c <- 1/N.c * Reduce('+', lapply(1:N, function(i)(
                                       p[i] * (data[i,]-mu) %*% 
                                       t(data[i,]-mu))))
    res <- 1/2 * N.c * (log(det(sigma)) + sum(diag(T.c %*% solve(sigma))) -
                        2*log(w))

    res
}

# negative log likelihood function for maximization of all classes at once
loglikelihood_stemm_constraints <- function(parameters, model, data, P) {
    model.filled <- fill_model(model, parameters)
    N <- nrow(data)
    res <- 0

    for (c in seq_len(model$info$num.classes)) {
        w.c <- model$info$w[c]
        N.c <- sum(P[,c])
        mu.c <- mu_stemm(model.filled$matrices[[c]])
        sigma.c <- sigma_stemm(model.filled$matrices[[c]])
        T.c <- 1/N.c * Reduce('+', lapply(1:N, function(i)(
                                           P[i,c] * (data[i,]-mu.c) %*% 
                                           t(data[i,]-mu.c))))
        res <- res+(1/2 * N.c * (log(det(sigma.c)) + sum(diag(T.c %*%
                                           solve(sigma.c))) - 2*log(w.c)))
    }

    res
}


# Maximization step of the EM-algorithm (see Jedidi, Jagpal & DeSarbo, 1997)
mstep_stemm <- function(model, parameters, data, P, neg.hessian=FALSE,
                        optimizer=c("nlminb", "optim"), constraints=FALSE,
                        max.mstep, control=list(), ...) {
    # --> TODO add constraints argument to doku, if it stays!

    optimizer <- match.arg(optimizer)

    if (constraints == FALSE) {
    ## maximizing each class seperately
        num.classes <- model$info$num.classes
        class.pars <- get_class_parameters(model, parameters)

        est <- lapply(seq_len(num.classes), function(c) {
                    if (optimizer == "nlminb") {
                            if (is.null(control$iter.max)) {
                                control$iter.max <- max.mstep
                            } else {
                                warning("iter.max is set for nlminb. max.mstep will be ignored.")
                            }

                        res <- nlminb(start=class.pars[[c]],
                                      objective=loglikelihood_stemm,
                                      data=data, matrices=model$matrices[[c]],
                                      p=P[,c], w=model$info$w[[c]],
                                      upper=model$info$bounds$upper[[c]],
                                      lower=model$info$bounds$lower[[c]],
                                      control=control, ...)
                    } else {
                            if (is.null(control$maxit)){
                                control$maxit <- max.mstep
                            } else {
                                warning("maxit is set for optim. max.mstep will be ignored.")
                            }

                        res <- optim(par=class.pars[[c]],
                                       fn=loglikelihood_stemm, data=data,
                                       matrices=model$matrices[[c]],
                                       p=P[,c], w=model$info$w[[c]],
                                       upper=model$info$bounds$upper[[c]],
                                       lower=model$info$bounds$lower[[c]],
                                       method="L-BFGS-B", control=control, ...)
                    }
        })
        if (optimizer == "optim") {
            for (c in seq_len(num.classes)) {
                names(est[[c]]) <- gsub("value", "objective", names(est[[c]]))
            }
        }
        res <- list(objective=0)
        for (c in seq_len(num.classes)) {
            res$par[[c]] <- est[[c]]$par
            res$objective <- res$objective + est[[c]]$objective
            res$convergence[[c]] <- est[[c]]$convergence
            res$message[[c]] <- est[[c]]$message
            res$iterations <- est[[c]]$iterations
        }
        names(res$par) <- paste0("class", seq_len(num.classes))

        if (neg.hessian == TRUE) {
            for (c in seq_len(num.classes)) {
                if (optimizer == "nlminb") {
                    res$hessian[[c]] <- fdHess(pars=est[[c]]$par,
                                               fun=loglikelihood_stemm,
                                               matrices=model$matrices[[c]],
                                               data=data, p=P[,c],
                                               w=model$info$w[[c]])$Hessian
                } else {
                    res$hessian[[c]] <- optimHess(par=est[[c]]$par,
                                                fn=loglikelihood_stemm,
                                                matrices=model$matrices[[c]],
                                                data=data, p=P[,c],
                                                w=model$info$w[c])
                }
            }
        names(res$hessian) <- paste0("class", seq_len(num.classes))
        }
        res

    } else {
    # Maximization of all classes together
        if (optimizer == "nlminb") {
            est <- nlminb(start=parameters,
                          objective=loglikelihood_stemm_constraints, data=data,
                          model=model, P=P,
                          upper=unlist(model$info$bounds$upper),
                          lower=unlist(model$info$bounds$lower), ...)
            if (neg.hessian == TRUE){
                est$hessian <- fdHess(pars=est$par,
                                      fun=loglikelihood_stemm_constraints,
                                      model=model, data=data, P=P)$Hessian
            }
        } else {
            est <- optim(par=parameters, fn=loglikelihood_stemm_constraints,
                         model=model, data=data, P=P,
                         upper=unlist(model$info$bounds$upper),
                         lower=unlist(model$info$bounds$lower),
                         method="L-BFGS-B", ...)
            # fit est to nlminb output
            names(est) <- gsub("value", "objective", names(est))
            if (neg.hessian == TRUE) {
                est$hessian <- optimHess(est$par,
                                         fn=loglikelihood_stemm_constraints,
                                         model=model, P=P, data=data)
            }
        }
        est
    }
}

#--------------- helper functions ---------------

# make a list of class specific parameter vectors
get_class_parameters <- function(model, parameters) {
    class.pars <- list()
    for (c in seq_len(model$info$num.classes)) {
        class.pars[[c]] <- parameters[1:length(model$info$par.names[[c]])]
        parameters <- parameters[(length(model$info$par.names[[c]]) + 1):length(parameters)]
    }
    class.pars
}






