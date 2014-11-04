# stemm.R
#
# created: Okt/20/2014, KN
# last mod: Okt/30/2014, KN

mu_stemm <- function(model, group) {
    stopifnot(class(model) == "stemmFilled")

    matrices <- model$matrices[[group]]

    # TODO case for no y's (as in Davids code). Necessary here?
    # TODO catch error if B is not nonsingular
    mu.y <- matrices$nu.y + matrices$Lambda.y %*% solve(matrices$Beta) %*%
            (matrices$alpha + matrices$Gamma %*% matrices$tau)
    mu.x <- matrices$nu.x + matrices$Lambda.x %*% matrices$tau
    mu <- rbind(mu.y, mu.x) # vertical vector

    mu
}

sigma_stemm <- function(model, group) {
    stopifnot(class(model) == "stemmFilled")

    matrices <- model$matrices[[group]]
    # Lambda.y * B^-1
    Ly.Binv <- matrices$Lambda.y %*% solve(matrices$Beta)

    s11 <- Ly.Binv %*% (matrices$Gamma %*% matrices$Phi %*% t(matrices$Gamma) +
                        matrices$Psi) %*% t(Ly.Binv) + matrices$Theta.e
    s12 <- Ly.Binv %*% matrices$Gamma %*% matrices$Phi %*% t(matrices$Lambda.x)
    s21 <- matrices$Lambda.x %*% matrices$Phi %*% t(matrices$Gamma) %*% t(Ly.Binv)
    s22 <- matrices$Lambda.x %*% matrices$Phi %*% t(matrices$Lambda.x) + matrices$Theta.d
    sigma <- rbind(cbind(s11,s12), cbind(s21, s22))

    if (!isSymmetric(sigma)) stop("Sigma has to be symmetric")

    sigma
}

estep_stemm <- function(model, parameters, data) {
    # TODO do I need the ... here?

    model.filled <- fill_model(model=model, parameters=parameters)

    P <- NULL
    for (g in seq_len(model$info$num.groups)) {
        # group weight
        w.g <- model$info$w[g]

        p.ij <- w.g * dmvnorm(data, mean=mu_stemm(model.filled, g),
                              sigma=sigma_stemm(model.filled, g))
        P <- cbind(P, p.ij, deparse.level=0)
    }
    P <- P / rowSums(P)
    P
}

loglikelihood_stemm <- function(parameters, model, data, P) {
    # TODO model or model.filled as input parameter?
    model.filled <- fill_model(model, parameters)
    N <- nrow(data)
    res <- 0
    for (g in seq_len(model$info$num.groups)) {
        w.g <- model$info$w[g]
        N.g <- sum(P[,g])
        mu.g <- mu_stemm(model.filled, g)
        sigma.g <- sigma_stemm(model.filled, g)
        T.g <- 1/N.g * Reduce('+', lapply(1:N, function(i)(
                                           P[i,g] * (data[i,]-mu.g) %*% 
                                           t(data[i,]-mu.g))))
        res <- res+(1/2 * N.g * (log(det(sigma.g)) + sum(diag(T.g %*%
                                           solve(sigma.g))) - 2*log(w.g)))
    }
    res
}

mstep_stemm <- function(model, parameters, data, P, Hessian=FALSE,
                        optimizer, ...) {

    optimizer <- match.arg(optimizer)

    if (optimizer == "nlminb") {
        est <- nlminb(start=parameters, objective=loglikelihood_stemm, data=data,
                      model=model, P=P, upper=model$info$bounds$upper,
                      lower=model$info$bounds$lower, ...)

        if (Hessian == TRUE){
            est$hessian <- nlme::fdHess(pars=est$par, fun=loglikelihood_stemm,
                                        model=model, data=data, P=P)
        }
    } else {
        est <- optim(par=parameters, fn=loglikelihood_stemm, model=model,
                     data=data, P=P, upper=model$info$bounds$upper,
                     lower=model$info$bounds$lower, ...)
        # fit est to nlminb output
        names(est) <- gsub("value", "objective", names(est))
        if (Hessian == TRUE) {
            est$hessian <- optimHess(est$par, fn=loglikelihood_stemm, model=model,
                                     P=P, data=data)
        }
    }
    est
}

