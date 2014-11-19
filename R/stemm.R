# stemm.R
#
# created: Okt/20/2014, KN
# last mod: Nov/11/2014, KN

# Calculate mu of multivariate normal distribution for joint vector of
# indicators (See equation 4 in Jedidi, Jagpal & DeSarbo, 1997).
# The order is (x, y) as opposed to the paper.
mu_stemm <- function(model, group) {
    stopifnot(class(model) == "stemmFilled" || class(model) == "nsemmFilled")

    matrices <- model$matrices[[group]]

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
sigma_stemm <- function(model, group) {
    stopifnot(class(model) == "stemmFilled" || class(model) == "nsemmFilled")

    matrices <- model$matrices[[group]]
    # Lambda.y * B^-1
    Ly.Binv <- matrices$Lambda.y %*% solve(matrices$Beta)

    s22 <- Ly.Binv %*% (matrices$Gamma %*% matrices$Phi %*% t(matrices$Gamma) +
                        matrices$Psi) %*% t(Ly.Binv) + matrices$Theta.e

    if (!isSymmetric(s22)) stop("S22 has to be symmetric")
    # --> TODO remove this error
    s21 <- Ly.Binv %*% matrices$Gamma %*% matrices$Phi %*% t(matrices$Lambda.x)
    s12 <- t(s21)
    s11 <- matrices$Lambda.x %*% matrices$Phi %*% t(matrices$Lambda.x) + matrices$Theta.d
    if (!isSymmetric(s11)) stop("S11 has to be symmetric")
    sigma <- rbind(cbind(s11,s12), cbind(s21, s22))

    if (!isSymmetric(sigma)) stop("Sigma has to be symmetric")
    tryCatch(solve(sigma), error = function(e) stop("Sigma is not nonsingular."))

    sigma
}

# Expectation step of the EM-algorithm (see Jedidi, Jagpal & DeSarbo, 1997)
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

# negative log likelihood function which will be optimized in M-step (see below)
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

# Maximization step of the EM-algorithm (see Jedidi, Jagpal & DeSarbo, 1997)
mstep_stemm <- function(model, parameters, data, P, Hessian=FALSE,
                        optimizer=c("nlminb", "optim"), ...) {

    optimizer <- match.arg(optimizer)

    if (optimizer == "nlminb") {
        est <- nlminb(start=parameters, objective=loglikelihood_stemm, data=data,
                      model=model, P=P,
                      upper=unlist(model$info$bounds$upper),
                      lower=unlist(model$info$bounds$lower), ...)
        if (Hessian == TRUE){
            est$hessian <- fdHess(pars=est$par, fun=loglikelihood_stemm,
                                        model=model, data=data, P=P)$Hessian
        }
    } else {
        est <- optim(par=parameters, fn=loglikelihood_stemm, model=model,
                     data=data, P=P, upper=unlist(model$info$bounds$upper),
                     lower=unlist(model$info$bounds$lower), method="L-BFGS-B", ...)
        # fit est to nlminb output
        names(est) <- gsub("value", "objective", names(est))
        if (Hessian == TRUE) {
            est$hessian <- optimHess(est$par, fn=loglikelihood_stemm, model=model,
                                     P=P, data=data)
        }
    }
    est
}

