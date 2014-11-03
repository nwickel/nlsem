# model.R
#
# created Sep/23/2014, NU
# last mod Nov/03/2014, KN

grep_ind <- function(x){
    tryCatch({
        if (length(unlist(strsplit(x, "-"))) > 1){
            as.numeric(gsub("^.*[x,y]([0-9]+).*[x,y]([0-9]+)$", "\\1",
            x)):as.numeric(gsub("^.*[x,y]([0-9]+).*[x,y]([0-9]+)$", "\\2", x))
        } else {
            as.numeric(gsub("^.([0-9]+).*$", "\\1", x))
        }
    }, warning = function(war) {
        # TODO more clear error message
        stop("Wrong input for specifying exogenous or endogonous latent
             variables (xi or etas). See ?specify_sem.")
    })
}

calc_interaction_matrix <- function(x){
    tryCatch({
        rows <- as.numeric(gsub("^.*xi([0-9]+):xi[0-9]+$", "\\1", x))
        cols <- as.numeric(gsub("^.*xi.*:xi([0-9]+)$", "\\1", x))
        mat  <- sort_interaction_effects(rows, cols)
        mat
    }, warning = function(war) {
        stop("Wrong input for interaction. See ?specify_sem.")
    }, error = function(err) { # perhaps error catching is unnecessary
        stop("Wrong input for interaction. See ?specify_sem.")
    })
}

sort_interaction_effects <- function(rows, cols){
    
    for (i in seq_along(rows)){
        if (rows[i] > cols[i]){
            rows_i <- rows[i]
            cols_i <- cols[i]
            rows[i] <- cols_i
            cols[i] <- rows_i
        }
    }
    cbind(rows, cols)
}

test_omega <- function(Omega){
    
    if (any(is.na(diag(Omega)))){
        stop("Can't handle quadratic interaction effects (yet). See
             ?specify_sem for details.")
    }

    if (any(is.na(Omega))){

        ind <- which(is.na(Omega), arr.ind=TRUE)
        dim <- nrow(Omega)
        msg <- "Interactions are not well-defined. Please change order of xi's.
        See ?specify_sem for details."

        if (nrow(ind) == 1){
            if (!is.na(Omega[1, dim]))
                stop(msg)

        } else {
            for (i in 2:nrow(ind)){
                if(ind[i,1] >= ind[i-1,2]){
                    stop(msg)
                }
            }
            for (i in 1:max(ind[,1])){
                if (max(which(is.na(Omega[i,]))) < dim){
                    stop(msg)
                }
            }
            if (max(ind[,1]) > 1){
                for (i in 2:max(ind[,1])){
                    if (min(which(is.na(Omega[i-1,]))) >= min(which(is.na(Omega[i,])))) {
                        stop(msg)
                    }
                }
            }
        }
        invisible(NULL)
    } else {
        invisible(NULL)
    }
}

rel_lat <- function(x, num.eta, num.xi){
    
    x.s       <- unlist(strsplit(x, ";"))
    which.xi  <- which(grepl("xi", x.s))
    which.eta <- which(!grepl("xi", x.s))
    
    G <- matrix(0, nrow=num.eta, ncol=num.xi)
    for (i in which.xi){
        xi.s <- unlist(strsplit(x.s[i], ">"))
            if (length(xi.s) < 2){
                stop("Latent variables misspecified. Must be of the form
                'xi1>eta1'. See ?specify_sem for details.")
            }
        xis  <- unlist(strsplit(xi.s[1], ","))
        etas <- unlist(strsplit(xi.s[2], ","))
        ind.xi <- as.numeric(gsub("^.*xi([0-9]).*$", "\\1", xis))
        ind.eta <- as.numeric(gsub("^.*eta([0-9]).*$", "\\1", etas))
        G[ind.eta, ind.xi] <- NA
    }
    
    B <- diag(1, num.eta)
    for (i in which.eta){
        eta.s <- unlist(strsplit(x.s[i], ">"))
            if (length(eta.s) < 2){
                stop("Latent variables misspecified. Must be of the form
                'xi1>eta1'. See ?specify_sem for details.")
            }
        eta.rows <- unlist(strsplit(eta.s[1], ","))
        eta.cols <- unlist(strsplit(eta.s[2], ","))
        ind.rows <- as.numeric(gsub("^.*eta([0-9]).*$", "\\1", eta.rows))
        ind.cols <- as.numeric(gsub("^.*eta([0-9]).*$", "\\1", eta.cols))
        B[ind.rows, ind.cols] <- NA
    }
    out <- list(Gamma=G, Beta=B)
    out
}

get_model_class <- function(num.groups, interaction) {
    if (num.groups == 1) {
            if (interaction == "") {
                # TODO should still work for stemm
                stop("Model needs either more than one latent group or at least one
                     latent interaction (e.g. 'xi1:xi2'). For other models please
                     use lavaan or the like.")
            } else {
                model_class <- "lms"
            }
        } else if (interaction == "") {
            model_class <- "stemm"
        } else {
            model_class <- "nsemm"
        }
    model_class
}

specify_sem <- function(num.x, num.y, num.xi, num.eta, xi, eta, num.groups=1,
                          interaction="all", interc_obs=FALSE,
                          interc_lat=FALSE, relation_lat="default"){

    # check arguments
    if (!is.numeric(num.x) || !is.numeric(num.y) || !is.numeric(num.xi) 
       || !is.numeric(num.eta) || !is.numeric(num.groups)) {
        stop("Number of variables or groups must be numeric.")
    } else if (num.x < num.xi || num.y < num.eta) {
        stop("The model contains not enough observed variables.")
    }
    stopifnot(num.x > 0, num.y >= 0, num.xi > 0, num.eta >= 0, num.groups > 0)

    # check if only defined xi's are in the interaction
    if (interaction != "all" && interaction != "") {
        interact.matrix <- calc_interaction_matrix(unlist(strsplit(interaction, ",")))
        if (max(interact.matrix) > num.xi) {
            stop("Interaction effects contain more xi's than defined.")
        }
    }

    # class of model
    model_class <- get_model_class(num.groups, interaction)

    xi.s <- unlist(strsplit(xi, ","))
    if (length(xi.s) != num.xi) {
        stop("Number of xi's and assignation of x's to xi's does not match.
             See ?specify_sem.")
    }
    xi.ind <- list()
    for (i in seq_len(num.xi)) xi.ind[[i]] <- grep_ind(xi.s[i])
    # TODO Use sapply instead of loop!

    #eta.ind <- grep_ind(eta)
    eta.s <- unlist(strsplit(eta, ","))
    if (length(eta.s) != num.eta) {
        stop("Number of eta's and assignation of y's to eta's does not match.
             See ?specify_sem.")
    }
    eta.ind <- list()
    for (i in seq_len(num.eta)) eta.ind[[i]] <- grep_ind(eta.s[i])
    # TODO Use sapply instead of loop!

    # create matrices with default constraints
    # Lambda.x
    # TODO catch misspecification
    Lambda.x <- matrix(0, nrow=num.x, ncol=num.xi)
    for (i in seq_len(num.xi)){
        Lambda.x[xi.ind[[i]], i] <- c(1, rep(NA, length(xi.ind[[i]])-1))
    }
    # Lambda.y
    # TODO catch misspecification
    Lambda.y <- matrix(0, nrow=num.y, ncol=num.eta)
    for (i in seq_len(num.eta)){
        Lambda.y[eta.ind[[i]], i] <- c(1, rep(NA, length(eta.ind[[i]])-1))
    }
    # Gamma
    # TODO specification for stemm
    if (relation_lat == "default"){
        Gamma <- matrix(nrow=num.eta, ncol=num.xi)
        # Beta
        Beta <- diag(1, nrow=num.eta)
    }
    else {
        GB <- rel_lat(relation_lat, num.eta=num.eta, num.xi=num.xi)
        tryCatch({
        Gamma <- GB[[grep("G", names(GB))]]}, error=function(e){
        Gamma <- matrix(nrow=num.eta, ncol=num.xi);Gamma})
        tryCatch({
        Beta <- GB[[grep("B", names(GB))]]}, error=function(e){
        Beta <- diag(1, nrow=num.eta);Beta})
    }
    # Theta.d
    Theta.d <- diag(NA, nrow=num.x)
    # Theta.e
    Theta.e <- diag(NA, nrow=num.y)
    # Psi
    Psi <- matrix(NA, nrow=num.eta, ncol=num.eta)
    # Psi must be symmetrical, upper.tri = lower.tri -> in fill_model
    Psi[upper.tri(Psi)] <- 0
    # Phi
    Phi <- matrix(NA, nrow=num.xi, num.xi)
    # Phi must be symmetrical, upper.tri = lower.tri -> in fill_model
    Phi[upper.tri(Phi)] <- 0
    # A
    A <- matrix(NA, nrow=num.xi, ncol=num.xi)
    A[upper.tri(A)] <- 0
    # nu's
    if (interc_obs){
        nu.x <- matrix(NA, nrow=num.x, ncol=1)
        nu.y <- matrix(NA, nrow=num.y, ncol=1)
    } else {
        nu.x <- matrix(0, nrow=num.x, ncol=1)
        nu.y <- matrix(0, nrow=num.y, ncol=1)
    }
    # alpha
    if (interc_lat){
        alpha <- matrix(NA, nrow=num.eta, ncol=1)
    } else {
        alpha <- matrix(0, nrow=num.eta, ncol=1)
    }
    # tau
    # TODO specify default constraints
    tau <- matrix(0, nrow=num.xi, ncol=1)
    # Omega
    Omega <- matrix(0, nrow=num.xi, ncol=num.xi)
    if (interaction == "all"){
        Omega[upper.tri(Omega)] <- NA
    } else if (interaction != "") {
        interaction.s <- unlist(strsplit(interaction, ","))
        ind <- calc_interaction_matrix(interaction.s)
        Omega[ind] <- NA
    }
    # check if Omega has row echelon form
    test_omega(Omega)

    # make a list of the matrices for each group
    matrices <- list()
    for (g in seq_len(num.groups)) {
        if (model_class == "lms") {
            matrices[[g]] <- list(Lambda.x=Lambda.x, Lambda.y=Lambda.y,
                                  Gamma=Gamma, Theta.d=Theta.d,
                                  Theta.e=Theta.e, Psi=Psi, A=A,
                                  nu.x=nu.x, nu.y=nu.y, alpha=alpha, tau=tau,
                                  Omega=Omega)
        } else if (model_class == "stemm") {
            matrices[[g]] <- list(Lambda.x=Lambda.x, Lambda.y=Lambda.y,
                                  Gamma=Gamma, Beta=Beta, Theta.d=Theta.d,
                                  Theta.e=Theta.e, Psi=Psi, Phi=Phi,
                                  nu.x=nu.x, nu.y=nu.y, alpha=alpha, tau=tau)
        } else {
            matrices[[g]] <- list(Lambda.x=Lambda.x, Lambda.y=Lambda.y,
                                  Gamma=Gamma, Beta=Beta, Theta.d=Theta.d,
                                  Theta.e=Theta.e, Psi=Psi, Phi=Phi, A=A,
                                  nu.x=nu.x, nu.y=nu.y, alpha=alpha, tau=tau,
                                  Omega=Omega)
        }
    }
    names(matrices) <- paste0("group",1:num.groups)

    # group weights w
    w <- matrix(1/num.groups, nrow=num.groups, ncol=1)

    # create model with matrices and info
    model <- list(matrices=matrices, info=list(num.xi=num.xi, num.eta=num.eta,
                                               num.x=num.x, num.y=num.y,
                                               num.groups=num.groups,
                                               par.names=list(), w=w))

    # add parameter names to model
    # if (num.groups == 1) {
    if (model_class == "lms") {
        model$info$par.names <- get_parnames(model)$group1
    } else {
        model$info$par.names <- get_parnames(model)
    }

    class(model) <- model_class
    model
}

get_parnames <- function(model) {
    # TODO if this function is not needed in fill_matrices, move into
    # specify_sem
    par.names <- list()
    for (g in seq_len(model$info$num.groups)) {
        lst <- unlist(lapply(model$matrices[[g]], is.na))
        par.names[[g]] <- names(lst[lst])
    }
    names(par.names) <- paste0("group", 1:model$info$num.groups)
    par.names
}

fill_matrices <- function(dat, model){
    # TODO model is input here so "info" is not lost. Better solution required.

    stopifnot(is.data.frame(dat))

    Lambda.x <- as.character(dat$label[grep("Lambda.x", dat$label)])
    Lambda.y <- as.character(dat$label[grep("Lambda.y", dat$label)])
    Gamma    <- as.character(dat$label[grep("Gamma", dat$label)])
    Beta     <- as.character(dat$label[grep("Beta", dat$label)])
    Theta.d  <- as.character(dat$label[grep("Theta.d", dat$label)])
    Theta.e  <- as.character(dat$label[grep("Theta.e", dat$label)])
    Psi      <- as.character(dat$label[grep("Psi", dat$label)])
    Phi      <- as.character(dat$label[grep("Phi", dat$label)])
    A        <- as.character(dat$label[grep("A", dat$label)])
    nu.x     <- as.character(dat$label[grep("nu.x", dat$label)])
    nu.y     <- as.character(dat$label[grep("nu.y", dat$label)])
    alpha    <- as.character(dat$label[grep("alpha", dat$label)])
    tau      <- as.character(dat$label[grep("tau", dat$label)])
    Omega    <- as.character(dat$label[grep("Omega", dat$label)])

    # number of latent and indicator variables and groups
    # num.x      <- length(nu.x)
    # num.y      <- length(nu.y)
    # num.xi     <- length(tau)
    # num.eta    <- length(alpha)
    # num.groups <- ncol(dat) - 1
    num.x <- model$info$num.x
    num.y <- model$info$num.y
    num.xi <- model$info$num.xi
    num.eta <- model$info$num.eta
    num.groups <- model$info$num.groups

    # create matrices
    matrices <- list()
    for (g in seq_len(num.groups)) {
        Lambda.x.matrix <- matrix(dat[dat$label %in% Lambda.x, paste0("group",g)],
                                  nrow=num.x, ncol=num.xi)
        Lambda.y.matrix <- matrix(dat[dat$label %in% Lambda.y, paste0("group",g)],
                                  nrow=num.y, ncol=num.eta)
        Gamma.matrix    <- matrix(dat[dat$label %in% Gamma, paste0("group",g)],
                                  nrow=num.eta, ncol=num.xi)
        Beta.matrix     <- matrix(dat[dat$label %in% Beta, paste0("group",g)],
                                  nrow=num.eta, ncol=num.eta)
        Theta.d.matrix  <- matrix(dat[dat$label %in% Theta.d, paste0("group",g)],
                                  nrow=num.x, ncol=num.x)
        Theta.e.matrix  <- matrix(dat[dat$label %in% Theta.e, paste0("group",g)],
                                  nrow=num.y, ncol=num.y)
        Psi.matrix      <- matrix(dat[dat$label %in% Psi, paste0("group",g)],
                                  nrow=num.eta, ncol=num.eta)
        Phi.matrix      <- matrix(dat[dat$label %in% Phi, paste0("group",g)],
                                  nrow=num.xi, ncol=num.xi)
        A.matrix        <- matrix(dat[dat$label %in% A, paste0("group",g)],
                                  nrow=num.xi, ncol=num.xi)
        nu.x.matrix     <- matrix(dat[dat$label %in% nu.x, paste0("group",g)],
                                  nrow=num.x, ncol=1)
        nu.y.matrix     <- matrix(dat[dat$label %in% nu.y, paste0("group",g)],
                                  nrow=num.y, ncol=1)
        alpha.matrix    <- matrix(dat[dat$label %in% alpha, paste0("group",g)],
                                  nrow=num.eta, ncol=1)
        tau.matrix      <- matrix(dat[dat$label %in% tau, paste0("group",g)],
                                  nrow=num.xi, ncol=1)
        Omega.matrix    <- matrix(dat[dat$label %in% Omega, paste0("group",g)],
                                  nrow=num.xi, ncol=num.xi)

        if (class(model) == "lms") {
            matrices[[g]] <- list(Lambda.x=Lambda.x.matrix,
                                  Lambda.y=Lambda.y.matrix, Gamma=Gamma.matrix,
                                  Theta.d=Theta.d.matrix,
                                  Theta.e=Theta.e.matrix, Psi=Psi.matrix,
                                  A=A.matrix, nu.x=nu.x.matrix,
                                  nu.y=nu.y.matrix, alpha=alpha.matrix,
                                  tau=tau.matrix, Omega=Omega.matrix)
        } else if (class(model) == "stemm") {
            matrices[[g]] <- list(Lambda.x=Lambda.x.matrix,
                                  Lambda.y=Lambda.y.matrix, Gamma=Gamma.matrix,
                                  Beta=Beta.matrix, Theta.d=Theta.d.matrix,
                                  Theta.e=Theta.e.matrix, Psi=Psi.matrix,
                                  Phi=Phi.matrix, nu.x=nu.x.matrix,
                                  nu.y=nu.y.matrix, alpha=alpha.matrix,
                                  tau=tau.matrix)
        } else {
             matrices[[g]] <- list(Lambda.x=Lambda.x.matrix,
                                   Lambda.y=Lambda.y.matrix, Gamma=Gamma.matrix,
                                   Beta=Beta.matrix, Theta.d=Theta.d.matrix,
                                   Theta.e=Theta.e.matrix, Psi=Psi.matrix,
                                   Phi=Phi.matrix, A=A.matrix, nu.x=nu.x.matrix,
                                   nu.y=nu.y.matrix, alpha=alpha.matrix,
                                   tau=tau.matrix, Omega=Omega.matrix)
        }
    }
    names(matrices) <- paste0("group",1:num.groups)

    model.new <- list(matrices=matrices, info=model$info)
    model.new$info$par.names <- get_parnames(model.new)
    # TODO perhaps different for lms

    class(model.new) <- class(model) # TODO this must not be the case!
    model.new
}

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

count_free_parameters <- function(model) {
    # w is not counted
    res <- 0
    for (g in seq_len(model$info$num.groups))
        res <- res + sum(unlist(lapply(model$matrices[[g]], is.na)))
    res
    }

fill_model <- function(model, parameters, version="new") {

    stopifnot(class(model) == "lms" || class(model) == "stemm"
              || class(model) == "nsemm")

    stopifnot(count_free_parameters(model) == length(parameters))

    matrices <- model$matrices

    for (g in seq_len(model$info$num.groups)) {
        if (class(model) == "lms") {
            par.names <- model$info$par.names
        } else {
            par.names <- model$info$par.names[[g]]
        }
        matrices.g <- matrices[[g]] # to avoid multiple [[]]

        for (j in seq_along(matrices.g)) {
            matrix.j <- matrices.g[[j]]
            # number of NA's in matrix
            num.na <- length(matrix.j[is.na(matrix.j)])
            if (num.na > 0) {
                matrix.j[is.na(matrix.j)] <- parameters[1:num.na]
                parameters <- parameters[-(1:num.na)]
                matrices.g[[j]] <- matrix.j
            }
        }
        # make 'symmetric' matrices symmetric
        matrices.g$Psi <- fill_symmetric(matrices.g$Psi)
        if (class(model) == "stemm") matrices.g$Phi <- fill_symmetric(matrices.g$Phi)
        matrices[[g]] <- matrices.g
    }
    out <- list(matrices=matrices, info=model$info)
    switch(class(model),
                 "lms" = {class(out) <- "lmsFilled"},
                 "stemm" = {class(out) <- "stemmFilled"},
                 "nsemm" = {class(out) <- "nsemmFilled"})
    out
}

# fill upper.tri of a (filled) matrix which should be symmetric
fill_symmetric <- function(mat) {
    for (i in seq_len(nrow(mat))) {
                for (j in i:ncol(mat)) {
                        mat[i,j] <- mat[j,i]
                }
            }
    mat
}

## TODO Want fill_matrices to work with a data frame created with lavaanify()...

# # specify model with lavaan
# # needs certain labels to work properly
# my.model <- '# measurement
#             xi1 =~ 1*x1 + Lx1*x1 + Lx2*x2
#             xi2 =~ 1*x3 + Lx3*x3 + Lx4*x4
#             eta =~ 1*y + Ly1*y
# 
#             # structural
#             !interaction := xi1*xi2
#             eta ~ a1*1 + G1*xi1 + G2*xi2 
# 
#             # variances and covariances
#             x1 ~~ Td1*x1
#             x2 ~~ Td2*x2
#             x3 ~~ Td3*x3
#             x4 ~~ Td4*x4
#             y ~~ 0*y + Te1*y
#             eta ~~ Psi1*eta
#             xi1 ~~ Phi1*xi1
#             xi2 ~~ Phi2*xi2'
# 
# #lavaanify(myModel)
# specs <- lavaanify(myModel)[,c("label", "ustart")]

