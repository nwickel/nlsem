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

    model.filled <- fill_model(model=model, parameters=parameters)

    P <- NULL
    for (g in seq_len(model$info$num.groups)) {
        # group weight
        w.g <- model$info$w[g]

        p.ij <- w.g * dmvnorm(data, mean=mu_stemm(model.filled, g),
                              sigma=sigma_stemm(model.filled, g))
        P <- cbind(P, p.ij)
    }
    P <- P / rowSums(P)
    P
}

likelihood <- function(model, parameters, data, P) {
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

