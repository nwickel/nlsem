# stemm.R
#
# created: Okt/20/2014, KN
# last mod: Nov/11/2014, KN

# Calculate mu of multivariate normal distribution for joint vector of
# indicators (See equation 4 in Jedidi, Jagpal & DeSarbo, 1997).
# The order is (x, y) as opposed to the paper.
mu_stemm <- function(matrices) {
    # stopifnot(class(model) == "stemmFilled" || class(model) == "nsemmFilled")

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
    # stopifnot(class(model) == "stemmFilled" || class(model) == "nsemmFilled")

    # Lambda.y * B^-1
    Ly.Binv <- matrices$Lambda.y %*% solve(matrices$Beta)

    s22 <- Ly.Binv %*% (matrices$Gamma %*% matrices$Phi %*% t(matrices$Gamma) +
                        matrices$Psi) %*% t(Ly.Binv) + matrices$Theta.e

    s21 <- Ly.Binv %*% matrices$Gamma %*% matrices$Phi %*% t(matrices$Lambda.x)
    s12 <- t(s21)
    s11 <- matrices$Lambda.x %*% matrices$Phi %*% t(matrices$Lambda.x) + matrices$Theta.d
    sigma <- rbind(cbind(s11,s12), cbind(s21, s22))

    if (!isSymmetric(sigma)) stop("Sigma has to be symmetric")
    tryCatch(solve(sigma), error = function(e) stop("Sigma is not nonsingular."))

    sigma
}

# Expectation step of the EM-algorithm (see Jedidi, Jagpal & DeSarbo, 1997)
estep_stemm <- function(model, parameters, data) {

    model.filled <- fill_model(model=model, parameters=parameters)

    P <- NULL
    for (g in seq_len(model$info$num.groups)) {
        # group weight
        w.g <- model$info$w[g]

        p.ij <- w.g * dmvnorm(data, mean=mu_stemm(model.filled$matrices[[g]]),
                              sigma=sigma_stemm(model.filled$matrices[[g]]))
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

    N <- nrow(data)

    N.g <- sum(p)
    mu <- mu_stemm(matrices)
    sigma <- sigma_stemm(matrices)
    T.g <- 1/N.g * Reduce('+', lapply(1:N, function(i)(
                                       p[i] * (data[i,]-mu) %*% 
                                       t(data[i,]-mu))))
    res <- 1/2 * N.g * (log(det(sigma)) + sum(diag(T.g %*% solve(sigma))) -
                        2*log(w))

    res
}

# negative log likelihood function for maximization of all groups at once
loglikelihood_stemm_constraints <- function(parameters, model, data, P) {
    model.filled <- fill_model(model, parameters)
    N <- nrow(data)
    res <- 0

    for (g in seq_len(model$info$num.groups)) {
        w.g <- model$info$w[g]
        N.g <- sum(P[,g])
        mu.g <- mu_stemm(model.filled$matrices[[g]])
        sigma.g <- sigma_stemm(model.filled$matrices[[g]])
        T.g <- 1/N.g * Reduce('+', lapply(1:N, function(i)(
                                           P[i,g] * (data[i,]-mu.g) %*% 
                                           t(data[i,]-mu.g))))
        res <- res+(1/2 * N.g * (log(det(sigma.g)) + sum(diag(T.g %*%
                                           solve(sigma.g))) - 2*log(w.g)))
    }

    res
}


# Maximization step of the EM-algorithm (see Jedidi, Jagpal & DeSarbo, 1997)
mstep_stemm <- function(model, parameters, data, P, Hessian=FALSE,
                        optimizer=c("nlminb", "optim"), constraints=FALSE, ...) {
    # --> TODO add constraints argument to doku, if it stays!

    optimizer <- match.arg(optimizer)

    if (constraints == FALSE) {
    ## maximizing each group seperately
        num.groups <- model$info$num.groups

        # getting parameters for each group
        pars <- parameters
        group.pars <- list()
        for (g in seq_len(num.groups)) {
            group.pars[[g]] <- pars[1:length(model$info$par.names[[g]])]
            pars <- pars[(length(model$info$par.names[[g]]) + 1):length(pars)]
        }

        est <- lapply(seq_len(num.groups), function(g) {
                    if (optimizer == "nlminb") {
                        res <- nlminb(start=group.pars[[g]],
                                      objective=loglikelihood_stemm,
                                      data=data, matrices=model$matrices[[g]],
                                      p=P[,g], w=model$info$w[[g]],
                                      upper=model$info$bounds$upper[[g]],
                                      lower=model$info$bounds$lower[[g]], ...)
                    } else {
                        res <- optim(par=group.pars[[g]],
                                       fn=loglikelihood_stemm, data=data,
                                       matrices=model$matrices[[g]],
                                       p=P[,g], w=model$info$w[[g]],
                                       upper=model$info$bounds$upper[[g]],
                                       lower=model$info$bounds$lower[[g]],
                                       method="L-BFGS-B", ...)
                    }
        })
        if (optimizer == "optim") {
            for (g in seq_len(num.groups)) {
                names(est[[g]]) <- gsub("value", "objective", names(est[[g]]))
            }
        }
        res <- list(objective=0)
        for (g in seq_len(num.groups)) {
            res$par[[g]] <- est[[g]]$par
            res$objective <- res$objective + est[[g]]$objective
            res$convergence[[g]] <- est[[g]]$convergence
            res$message[[g]] <- est[[g]]$message
            # res$iterations <- est[[g]]$iterations
        }
        names(res$par) <- paste0("group", seq_len(num.groups))

        if (Hessian == TRUE) {
            for (g in seq_len(num.groups)) {
                if (optimizer == "nlminb") {
                    res$hessian[[g]] <- fdHess(pars=est[[g]]$par,
                                               fun=loglikelihood_stemm,
                                               matrices=model$matrices[[g]],
                                               data=data, p=P[,g],
                                               w=model$info$w[[g]])$Hessian
                } else {
                    res$hessian[[g]] <- optimHess(par=est[[g]]$par,
                                                fn=loglikelihood_stemm,
                                                matrices=model$matrices[[g]],
                                                data=data, p=P[,g],
                                                w=model$info$w[g])
                }
            }
        names(res$hessian) <- paste0("group", seq_len(num.groups))
        }
        res

    } else {
    # Maximization of all groups together
        if (optimizer == "nlminb") {
            est <- nlminb(start=parameters,
                          objective=loglikelihood_stemm_constraints, data=data,
                          model=model, P=P,
                          upper=unlist(model$info$bounds$upper),
                          lower=unlist(model$info$bounds$lower), ...)
            if (Hessian == TRUE){
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
            if (Hessian == TRUE) {
                est$hessian <- optimHess(est$par,
                                         fn=loglikelihood_stemm_constraints,
                                         model=model, P=P, data=data)
            }
        }
        est
    }
}

