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

    name.eta <- lavNames(dat, "lv.y")
    num.y    <- sum(dat$lhs[grep("^=~$", dat$op)] == name.eta)
    num.x    <- sum(dat$lhs[grep("^=~$", dat$op)] != name.eta)
    name.xi  <- lavNames(dat, "lv.x")[!grepl(":", lavNames(dat, "lv.x"))]
    num.xi   <- length(name.xi)
    num.eta  <- length(name.eta)
    name.x <- lavNames(dat, "ov")[dat$lhs %in% name.xi & dat$mat == "lambda"]
    name.y <- lavNames(dat, "ov")[dat$lhs %in% name.eta & dat$mat == "lambda"]

    # create matrices
    d.en <- dat[dat$lhs == name.eta & dat$mat == "lambda",]
    Lambda.y <- matrix(0, nrow=num.y, ncol=num.eta)

    Lambda.y[cbind(d.en$row, d.en$col)] <- d.en$ustart

    # Lambda.x
    d.ex <- dat[dat$lhs %in% name.xi & dat$mat == "lambda",]
    d.ex$row.new <- d.ex$row - num.y
    d.ex$col.new <- d.ex$col - num.eta
    Lambda.x <- matrix(0, nrow=num.x, ncol=num.xi)

    Lambda.x[cbind(d.ex$row.new, d.ex$col.new)] <- d.ex$ustart

    effects <- lavNames(dat, "lv.x")
    d.ga <- dat[dat$rhs %in% name.xi & dat$mat == "beta",]
    # Gamma and Beta
    Gamma <- matrix(nrow=num.eta, ncol=num.xi)
    d.ga$col.new <- d.ga$col - num.eta
    Gamma[cbind(d.ga$row, d.ga$col.new)] <- d.ga$ustart

    # if (num.eta == 1) {
    #   Beta  <- diag(num.eta)
    # } else {
      Beta  <- diag(num.eta)
      d.be <- dat[dat$rhs %in% name.eta & dat$mat == "beta",]
      Beta[cbind(d.be$row, d.be$col)] <- d.be$ustart
    # }
    
    # Omega
    nl.effects <- lavNames(dat, "lv.interaction")
    Omega <- matrix(0, nrow=num.xi, ncol=num.xi)
    e.s <- strsplit(nl.effects, ":")
    o.ind <- NULL
    for (i in seq_along(e.s)) {
      o.ind <- rbind(o.ind, which(name.xi %in% e.s[[i]]))
    }
    Omega[o.ind] <- NA

    # Theta.e
    Theta.e <- diag(nrow=num.y)
    d.Te <- dat[dat$lhs %in% name.y & dat$rhs %in% name.y,]
    Theta.e[cbind(d.Te$row, d.Te$col)] <- d.Te$ustart

    # Theta.d
    Theta.d <- diag(nrow=num.x)
    d.Td <- dat[dat$lhs %in% name.x & dat$rhs %in% name.x,]
    d.Td$row.new <- d.Td$row - num.y
    d.Td$col.new <- d.Td$col - num.y
    Theta.d[cbind(d.Td$row.new, d.Td$col.new)] <- d.Td$ustart

    # Psi
    Psi <- matrix(NA, nrow=num.eta, ncol=num.eta)
    Psi[upper.tri(Psi)] <- 0

    # Phi
    Phi <- matrix(NA, nrow=num.xi, num.xi)
    Phi[upper.tri(Phi)] <- 0
    
    # nu's
    nu.y <- matrix(0, nrow=num.y, ncol=1)
    d.nuy <- dat[dat$lhs %in% name.y & dat$mat == "nu",]
    nu.y[cbind(d.nuy$row, d.nuy$col)] <- d.nuy$ustart

    nu.y[which(d.en$ustart == 1), 1] <- 0

    nu.x <- matrix(0, nrow=num.x, ncol=1)
    d.nux <- dat[dat$lhs %in% name.x & dat$mat == "nu",]
    d.nux$row.new <- d.nux$row - num.y
    nu.x[cbind(d.nux$row.new, d.nux$col)] <- d.nux$ustart

    nu.x[which(d.ex$ustart == 1), 1] <- 0

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


