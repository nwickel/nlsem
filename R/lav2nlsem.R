lav2nlsem <- function(model, constraints=c("indirect", "direct1",
                      "direct2"), class.spec="class") {

  dat0 <- lavMatrixRepresentation(lavaanify(model, meanstructure=TRUE,
    auto=TRUE))

  if (class.spec %in% names(dat0)) {
    num.classes <- length(unique(dat0[, class.spec]))
  } else {
    num.classes <- 1
  }

  matrices <- list()
  for (c in seq_len(num.classes)) {

    if (class.spec %in% names(dat0)) {
      dat <- dat0[dat0[, class.spec] == c, ]
    } else {
      dat <- dat0
    }

    name.eta <- unique(dat$lhs[grep("^~$", dat$op)])
    #name.eta <- lavNames(dat, "lv.y")
    num.y    <- sum(dat$lhs[grep("^=~$", dat$op)] == name.eta)
    num.x    <- sum(dat$lhs[grep("^=~$", dat$op)] != name.eta)
    name.xi  <- lavNames(dat, "lv.x")[!grepl(":", lavNames(dat, "lv.x"))]
    num.xi   <- length(name.xi)
    num.eta  <- length(name.eta)
    name.x <- lavNames(dat, "ov")[dat$lhs %in% name.xi & dat$mat == "lambda"]
    name.y <- lavNames(dat, "ov")[dat$lhs %in% name.eta & dat$mat == "lambda"]

    # create matrices
    # Lambda.x and Lambda.y
    Lambda <- matrix(0, nrow=num.x + num.y, ncol=num.xi + num.eta)
    d.lambda <- dat[dat$mat == "lambda",]
    Lambda[cbind(d.lambda$row, d.lambda$col)] <- d.lambda$ustart
    colnames(Lambda) <- lavNames(dat, "lv")[!grepl(":", lavNames(dat, "lv"))]
    rownames(Lambda) <- lavNames(dat, "ov")

    Lambda.y <- as.matrix(Lambda[name.y, name.eta])
    dimnames(Lambda.y) <- NULL
    Lambda.x <- as.matrix(Lambda[name.x, name.xi])
    dimnames(Lambda.x) <- NULL

    # Gamma
    d.ga <- dat[dat$rhs %in% name.xi & dat$mat == "beta",]
    Gamma <- matrix(0, nrow=num.eta, ncol=num.xi)
    colnames(Gamma) <- name.xi
    rownames(Gamma) <- name.eta

    for (eta in name.eta) {
      for (xi in name.xi) {
        tryCatch({
        Gamma[eta, xi] <- d.ga$ustart[d.ga$lhs == eta & d.ga$rhs == xi]
        }, error=function(e) e)
      }
    }
    dimnames(Gamma) <- NULL

    # Beta
    d.be <- dat[dat$rhs %in% name.eta & dat$mat == "beta",]
    Beta  <- diag(num.eta)
    colnames(Beta) <- name.eta
    rownames(Beta) <- name.eta

    for (e1 in name.eta) {
      for (e2 in name.eta) {
        tryCatch({
        Beta[e1, e2] <- d.be$ustart[d.be$lhs == e1 & d.be$rhs == e2]
        }, error=function(e) e)
      }
    }
    dimnames(Beta) <- NULL

    # Omega
    nl.effects <- lavNames(dat, "lv.interaction")
    Omega <- matrix(0, nrow=num.xi, ncol=num.xi)
    e.s <- strsplit(nl.effects, ":")
    o.ind <- NULL
    for (i in seq_along(e.s)) {
      ind <- which(name.xi %in% e.s[[i]])
      if (length(ind) == 1) {
        ind <- rep(ind, 2)
      }
      o.ind <- rbind(o.ind, ind)
    }
    Omega[o.ind] <- NA

    # Theta.d and Theta.e
    Theta <- diag(nrow=num.x + num.y)
    d.theta <- dat[dat$mat == "theta",]
    Theta[cbind(d.theta$row, d.theta$col)] <- d.theta$ustart
    colnames(Theta) <- c(name.x, name.y)
    rownames(Theta) <- c(name.x, name.y)

    Theta.e <- as.matrix(Theta[name.y, name.y])
    dimnames(Theta.e) <- NULL
    Theta.d <- as.matrix(Theta[name.x, name.x])
    dimnames(Theta.d) <- NULL

    # Psi
    Psi <- matrix(NA, nrow=num.eta, ncol=num.eta)
    Psi[upper.tri(Psi)] <- 0

    # Phi
    Phi <- matrix(NA, nrow=num.xi, num.xi)
    Phi[upper.tri(Phi)] <- 0
    
    # nu's
    nu <- matrix(0, nrow=num.x + num.y, ncol=1)
    d.nu <- dat[dat$mat == "nu",]
    nu[cbind(d.nu$row, d.nu$col)] <- d.nu$ustart
    rownames(nu) <- c(name.x, name.y)

    nu.y <- as.matrix(nu[name.y, ])
    nu.y[which(Lambda.y == 1, arr.ind=T)[, "row"]] <- 0
    dimnames(nu.y) <- NULL
    nu.x <- as.matrix(nu[name.x, ])
    nu.x[which(Lambda.x == 1, arr.ind=T)[, "row"]] <- 0
    dimnames(nu.x) <- NULL

    # alpha
    alpha <- matrix(NA, nrow=num.eta, ncol=1)
    
    # tau
    tau <- matrix(NA, nrow=num.xi, ncol=1)

    # define model class
    if (num.classes == 1) {
      model.class <- "singleClass"
    } else {
      if (anyNA(Omega)) {
        model.class <- "nsemm"
      } else model.class <- "semm"
    }
    
    if (model.class == "singleClass") {
      matrices[[c]] <- list(Lambda.x=Lambda.x, Lambda.y=Lambda.y,
                            Gamma=Gamma, Beta=Beta, Theta.d=Theta.d,
                            Theta.e=Theta.e, Psi=Psi, Phi=Phi,
                            nu.x=nu.x, nu.y=nu.y, alpha=alpha, tau=tau,
                            Omega=Omega)
    } else if (model.class == "semm") {
      matrices[[c]] <- list(Lambda.x=Lambda.x, Lambda.y=Lambda.y,
                            Gamma=Gamma, Beta=Beta, Theta.d=Theta.d,
                            Theta.e=Theta.e, Psi=Psi, Phi=Phi,
                            nu.x=nu.x, nu.y=nu.y, alpha=alpha, tau=tau)
    } else {
      matrices[[c]] <- list(Lambda.x=Lambda.x, Lambda.y=Lambda.y,
                            Gamma=Gamma, Beta=Beta, Theta.d=Theta.d,
                            Theta.e=Theta.e, Psi=Psi, Phi=Phi,
                            nu.x=nu.x, nu.y=nu.y, alpha=alpha, tau=tau,
                            Omega=Omega)
    }
  }     # end loop for classes

  names(matrices) <- paste0("class", seq_len(num.classes))

  # class weights w
  w <- matrix(1/num.classes, nrow=num.classes, ncol=1)

  constraints <- match.arg(constraints)

  model <- list(matrices=matrices, info=list(num.xi=num.xi, num.eta=num.eta,
                                             num.x=num.x, num.y=num.y,
                                             constraints=constraints,
                                             num.classes=num.classes,
                                             par.names=list(), w=w))

  class(model) <- model.class
  # add parameter names to model
  model$info$par.names <- get_parnames(model=model, constraints=constraints)
  # bounds for parameters (variances > 0)
  model$info$bounds <- bounds(model, constraints=constraints)

  model
}


