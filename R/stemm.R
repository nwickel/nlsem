# stemm.R
#
# created: Okt/20/2014, KN
# last mod: Okt/20/2014, KN

mu_stemm <- function(model, group) {
    stopifnot(class(model) == "stemmFilled")

    matrices <- model$matrices[[group]]

    # TODO case for no y's (as in Davids code). Necessary here?
    # TODO catch error if B is not nonsingular
    mu.y <- matrices$nu.y + matrices$Lambda.y %*% solve(matrices$Beta) %*%
            (matrices$alpha + matrices$Gamma %*% matrices$tau)
    mu.x <- matrices$nu.x + matrices$Lambda.x %*% matrices$tau
    mu <- c(mu.y, mu.x)
    # mu <- rbind(mu.y, mu.x) # vertical vector

    mu
}

sigma_stemm <- function(model, group) {
    stopifnot(class(model) == "stemmFilled")

    matrices <- model$matrices[[group]]
    
    Ly.Binv <- matrices$Lambda.y %*% solve(matrices$Beta)

    s11 <- Ly.Binv %*% (matrices$Gamma %*% matrices$Phi %*% t(matrices$Gamma) +
                        matrices$Psi) %*% t(Ly.Binv) + matrices$Theta.e
    s21 <- matrices$Lambda.x %*% matrices$Phi %*% t(matrices$Gamma) %*% t(Ly.Binv)

}

