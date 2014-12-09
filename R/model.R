# model.R
#
# created Sep/23/2014, NU
# last mod Nov/25/2014, NU

#--------------- main functions ---------------

# Define model specification for different SEMs with nonlinear effects;
# possible objects classes are 'lms', 'stemm', 'nsemm'; exported function
specify_sem <- function(num.x, num.y, num.xi, num.eta, xi, eta, num.classes=1,
                          interaction="all", interc.obs=FALSE,
                          interc.lat=FALSE, relation.lat="default"){

    # check arguments
    if (!is.numeric(num.x) || !is.numeric(num.y) || !is.numeric(num.xi) 
       || !is.numeric(num.eta) || !is.numeric(num.classes)) {
        stop("Number of variables or classes must be numeric.")
    } else if (num.x < num.xi || num.y < num.eta) {
        stop("The model contains not enough observed variables.")
    }
    stopifnot(num.x > 0, num.y >= 0, num.xi > 0, num.eta >= 0, num.classes > 0)

    # check if only defined xi's are in the interaction
    if (interaction != "all" && interaction != "") {
        interact.matrix <- calc_interaction_matrix(unlist(strsplit(interaction, ",")))
        if (max(interact.matrix) > num.xi) {
            stop("Interaction effects contain more xi's than defined.")
        }
    }

    # class of model
    model.class <- get_model_class(num.classes, interaction)

    # check latent variables and get indices for Lambda.x and Lambda.y
    xi.s <- unlist(strsplit(xi, ","))
    if (length(xi.s) != num.xi) {
        stop("Number of xi's and assignation of x's to xi's does not match.
             See ?specify_sem.")
    }
    xi.ind <- list()
    for (i in seq_len(num.xi)) {
        xi.ind[[i]] <- grep_ind(xi.s[i])
        if (max(xi.ind[[i]]) > num.x) {
            stop("Number of x's assinged to xi exceeds x's specified. See
                 ?specify_sem")
        }
    }

    eta.s <- unlist(strsplit(eta, ","))
    if (length(eta.s) != num.eta) {
        stop("Number of eta's and assignation of y's to eta's does not match.
             See ?specify_sem.")
    }
    eta.ind <- list()
    for (i in seq_len(num.eta)) {
        eta.ind[[i]] <- grep_ind(eta.s[i])
        if (max(eta.ind[[i]]) > num.y) {
            stop("Number of y's assinged to eta exceeds y's specified. See
                 ?specify_sem.")
        }
    }

    # create matrices with default constraints
    # Lambda.x
    Lambda.x <- matrix(0, nrow=num.x, ncol=num.xi)
    for (i in seq_len(num.xi)){
        Lambda.x[xi.ind[[i]], i] <- c(1, rep(NA, length(xi.ind[[i]])-1))
    }
    # Lambda.y
    Lambda.y <- matrix(0, nrow=num.y, ncol=num.eta)
    for (i in seq_len(num.eta)){
        Lambda.y[eta.ind[[i]], i] <- c(1, rep(NA, length(eta.ind[[i]])-1))
    }
    # Gamma and Beta
    if (relation.lat == "default"){
        Gamma <- matrix(nrow=num.eta, ncol=num.xi)
        Beta  <- diag(1, nrow=num.eta)
    }
    else {
        GB    <- rel_lat(relation.lat, num.eta=num.eta, num.xi=num.xi)
        Gamma <- tryCatch({ GB[[grep("G", names(GB))]] },
                            error=function(e) matrix(nrow=num.eta, ncol=num.xi) )
        Beta  <- tryCatch({ GB[[grep("B", names(GB))]] },
                            error=function(e) diag(1, nrow=num.eta))
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
    # nu's
    if (interc.obs){
        nu.x <- matrix(NA, nrow=num.x, ncol=1)
        nu.y <- matrix(NA, nrow=num.y, ncol=1)
    } else {
        nu.x <- matrix(0, nrow=num.x, ncol=1)
        nu.y <- matrix(0, nrow=num.y, ncol=1)
    }
    # alpha
    if (interc.lat){
        alpha <- matrix(NA, nrow=num.eta, ncol=1)
    } else {
        alpha <- matrix(0, nrow=num.eta, ncol=1)
    }
    # tau
    tau <- matrix(0, nrow=num.xi, ncol=1)
    # Omega
    Omega <- matrix(0, nrow=num.xi, ncol=num.xi)
    if (interaction == "all"){
        Omega[upper.tri(Omega)] <- NA
        # --> TODO add the following so "all" does what description says
        # diag(Omega) <- NA
    } else if (interaction != "") {
        interaction.s <- unlist(strsplit(interaction, ","))
        ind <- calc_interaction_matrix(interaction.s)
        Omega[ind] <- NA
    }
    # check if Omega has row echelon form
    test_omega(Omega)

    # make a list of the matrices for each class
    matrices <- list()
    for (c in seq_len(num.classes)) {
        if (model.class == "lms") {
            matrices[[c]] <- list(Lambda.x=Lambda.x, Lambda.y=Lambda.y,
                                  Gamma=Gamma, Theta.d=Theta.d,
                                  Theta.e=Theta.e, Psi=Psi, Phi=Phi,
                                  nu.x=nu.x, nu.y=nu.y, alpha=alpha, tau=tau,
                                  Omega=Omega)
        } else if (model.class == "stemm") {
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
    }
    names(matrices) <- paste0("class",seq_len(num.classes))

    # class weights w
    w <- matrix(1/num.classes, nrow=num.classes, ncol=1)

    # create model with matrices and info
    model <- list(matrices=matrices, info=list(num.xi=num.xi, num.eta=num.eta,
                                               num.x=num.x, num.y=num.y,
                                               num.classes=num.classes,
                                               par.names=list(), w=w))

    # add parameter names to model
    if (model.class == "lms") {
        model$info$par.names <- get_parnames(model)$class1
    } else {
        model$info$par.names <- get_parnames(model)
    }

    ## --> TODO set parameters that do not vary between classes

    class(model) <- model.class
    
    # bounds for parameters (variances > 0)
    model$info$bounds <- bounds(model)

    model
}

# Create model matrices from a dataframe with columns label (for parameter
# labels) and class 1 to class n; only needed when user wants to have full
# control over constraints, etc.; exported function
create_sem <- function(dat){

    stopifnot(is.data.frame(dat))

    Lambda.x <- as.character(dat$label[grep("Lambda.x", dat$label)])
    Lambda.y <- as.character(dat$label[grep("Lambda.y", dat$label)])
    Gamma    <- as.character(dat$label[grep("Gamma", dat$label)])
    Beta     <- as.character(dat$label[grep("Beta", dat$label)])
    Theta.d  <- as.character(dat$label[grep("Theta.d", dat$label)])
    Theta.e  <- as.character(dat$label[grep("Theta.e", dat$label)])
    Psi      <- as.character(dat$label[grep("Psi", dat$label)])
    Phi      <- as.character(dat$label[grep("Phi", dat$label)])
    nu.x     <- as.character(dat$label[grep("nu.x", dat$label)])
    nu.y     <- as.character(dat$label[grep("nu.y", dat$label)])
    alpha    <- as.character(dat$label[grep("alpha", dat$label)])
    tau      <- as.character(dat$label[grep("tau", dat$label)])
    Omega    <- as.character(dat$label[grep("Omega", dat$label)])

    # number of latent and indicator variables and classes
    num.x      <- length(nu.x)
    num.y      <- length(nu.y)
    num.xi     <- length(tau)
    num.eta    <- length(alpha)
    num.classes <- ncol(dat) - 1

    # create matrices
    matrices <- list()
    for (c in seq_len(num.classes)) {
        Lambda.x.matrix <- matrix(dat[dat$label %in% Lambda.x, paste0("class",c)],
                                  nrow=num.x, ncol=num.xi)
        Lambda.y.matrix <- matrix(dat[dat$label %in% Lambda.y, paste0("class",c)],
                                  nrow=num.y, ncol=num.eta)
        Gamma.matrix    <- matrix(dat[dat$label %in% Gamma, paste0("class",c)],
                                  nrow=num.eta, ncol=num.xi)
        Beta.matrix     <- matrix(dat[dat$label %in% Beta, paste0("class",c)],
                                  nrow=num.eta, ncol=num.eta)
        Theta.d.matrix  <- matrix(dat[dat$label %in% Theta.d, paste0("class",c)],
                                  nrow=num.x, ncol=num.x)
        Theta.e.matrix  <- matrix(dat[dat$label %in% Theta.e, paste0("class",c)],
                                  nrow=num.y, ncol=num.y)
        Psi.matrix      <- matrix(dat[dat$label %in% Psi, paste0("class",c)],
                                  nrow=num.eta, ncol=num.eta)
        Phi.matrix      <- matrix(dat[dat$label %in% Phi, paste0("class",c)],
                                  nrow=num.xi, ncol=num.xi)
        nu.x.matrix     <- matrix(dat[dat$label %in% nu.x, paste0("class",c)],
                                  nrow=num.x, ncol=1)
        nu.y.matrix     <- matrix(dat[dat$label %in% nu.y, paste0("class",c)],
                                  nrow=num.y, ncol=1)
        alpha.matrix    <- matrix(dat[dat$label %in% alpha, paste0("class",c)],
                                  nrow=num.eta, ncol=1)
        tau.matrix      <- matrix(dat[dat$label %in% tau, paste0("class",c)],
                                  nrow=num.xi, ncol=1)
        Omega.matrix    <- matrix(dat[dat$label %in% Omega, paste0("class",c)],
                                  nrow=num.xi, ncol=num.xi)

        if (num.classes == 1) {
            matrices[[c]] <- list(Lambda.x=Lambda.x.matrix,
                                  Lambda.y=Lambda.y.matrix, Gamma=Gamma.matrix,
                                  Theta.d=Theta.d.matrix,
                                  Theta.e=Theta.e.matrix, Psi=Psi.matrix,
                                  Phi=Phi.matrix, nu.x=nu.x.matrix,
                                  nu.y=nu.y.matrix, alpha=alpha.matrix,
                                  tau=tau.matrix, Omega=Omega.matrix)
        } else if (all(is.na(Omega.matrix))) {
            matrices[[c]] <- list(Lambda.x=Lambda.x.matrix,
                                  Lambda.y=Lambda.y.matrix, Gamma=Gamma.matrix,
                                  Beta=Beta.matrix, Theta.d=Theta.d.matrix,
                                  Theta.e=Theta.e.matrix, Psi=Psi.matrix,
                                  Phi=Phi.matrix, nu.x=nu.x.matrix,
                                  nu.y=nu.y.matrix, alpha=alpha.matrix,
                                  tau=tau.matrix)
        } else {
             matrices[[c]] <- list(Lambda.x=Lambda.x.matrix,
                                   Lambda.y=Lambda.y.matrix, Gamma=Gamma.matrix,
                                   Beta=Beta.matrix, Theta.d=Theta.d.matrix,
                                   Theta.e=Theta.e.matrix, Psi=Psi.matrix,
                                   Phi=Phi.matrix, nu.x=nu.x.matrix,
                                   nu.y=nu.y.matrix, alpha=alpha.matrix,
                                   tau=tau.matrix, Omega=Omega.matrix)
        }
    }
    names(matrices) <- paste0("class", 1:num.classes)

    par.names <- list()
    for (c in seq_len(num.classes)) {
        class <- paste0("class",c)
        par.names[[c]] <- as.character(dat$label[is.na(dat[,class])])
    }
    names(par.names) <- paste0("class", 1:num.classes)
    w <- matrix(1/num.classes, nrow=num.classes, ncol=1)

    info <- list(num.xi=num.xi, num.eta=num.eta, num.x=num.x, num.y=num.y,
                 num.classes=num.classes, par.names=par.names, w=w)

    model <- list(matrices=matrices, info=info)

    if (all(is.na(Omega.matrix))) {
        interaction <- ""
    } else interaction <- "not_empty"

    class(model) <- get_model_class(num.classes, interaction)

    # bounds for parameters (variances > 0)
    model$info$bounds <- bounds(model)

    model
}

# Count free parameters of a model created with specify_sem (i.e. NAs in
# the model are counted); exported function
count_free_parameters <- function(model) {
    res <- 0
    for (c in seq_len(model$info$num.classes))
        res <- res + sum(unlist(lapply(model$matrices[[c]], is.na)))
    res
    }

# Fill a model created with specify_sem with parameters given as a vector;
# mostly needed to simulate data from a prespecified model; NOT exported
fill_model <- function(model, parameters) {

    stopifnot(class(model) == "lms" || class(model) == "stemm"
              || class(model) == "nsemm")

    stopifnot(count_free_parameters(model) == length(parameters))

    matrices <- model$matrices

    for (c in seq_len(model$info$num.classes)) {
        if (class(model) == "lms") {
            par.names <- model$info$par.names
        } else {
            par.names <- model$info$par.names[[c]]
        }
        matrices.c <- matrices[[c]] # to avoid multiple [[]]

        for (j in seq_along(matrices.c)) {
            matrix.j <- matrices.c[[j]]
            # number of NA's in matrix
            num.na <- length(matrix.j[is.na(matrix.j)])
            if (num.na > 0) {
                matrix.j[is.na(matrix.j)] <- parameters[1:num.na]
                parameters <- parameters[-(1:num.na)]
                matrices.c[[j]] <- matrix.j
            }
        }
        # make 'symmetric' matrices symmetric
        matrices.c$Psi <- fill_symmetric(matrices.c$Psi)
        tryCatch({ matrices.c$Phi <- fill_symmetric(matrices.c$Phi) },
                                       error=function(e) e,
                                       warning=function(w) w)
        # catching error if there is no Phi (like in "intern" LMS)
        matrices[[c]] <- matrices.c
    }
    out <- list(matrices=matrices, info=model$info)
    class(out) <- class(model)
    out
}

#--------------- helper functions ---------------

# all NOT exported

# Check if model or matrices are filled
check_filled <- function(x) {
    if (anyNA(unlist(x))) stop("model is not filled")
}

# Grep indices for Lambda matrices from input that defines which indicators
# are asociated with which latent variable
grep_ind <- function(x){
    tryCatch({
        if (length(unlist(strsplit(x, "-"))) > 1){
            as.numeric(gsub("^.*[x,y]([0-9]+).*[x,y]([0-9]+)$", "\\1",
            x)):as.numeric(gsub("^.*[x,y]([0-9]+).*[x,y]([0-9]+)$", "\\2", x))
        } else {
            as.numeric(gsub("^.([0-9]+).*$", "\\1", x))
        }
    }, warning = function(war) {
        stop("Wrong input for specifying exogenous or endogonous latent
             variables (xi or etas). See ?specify_sem.")
             # this might never be evaluated, error catched earlier
    })
}

# Returns matrix which specifies which latent variables interact with each
# other
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

# Ensures that interaction effects are in the correct order when passed to
# Omega
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

# Tests if input for Omega is in the correct format; Omega needs to be in
# row echelon form; returns nothing if Omega has the correct form
test_omega <- function(Omega){
    
    if (anyNA(Omega)){

        interactions <- Omega
        diag(interactions) <- 0
        ind <- which(is.na(interactions), arr.ind=TRUE)
        dim <- nrow(interactions)
        msg <- "Interactions are not well-defined. Please change order of xi's.
        See ?specify_sem for details."

        # test if any rows are 0 in between interaction effects
        na.r <- rowSums(interactions)
        for (i in seq_along(na.r)){
            if (!is.na(na.r[i]))
                if (is.na(sum(na.r[-c(1:i)])))
                    stop(msg)
        }

        if (nrow(ind) == 1){
            if (!is.na(interactions[1, dim]))
                stop(msg)

        } else {
            if (max(ind[,1]) > 1){
                for (i in 2:max(ind[,1])){
                    if (min(which(is.na(interactions[i-1,]))) >=
                    min(which(is.na(interactions[i,])))) {
                        stop(msg)
                    }
                }
            }
        }

        # test if quadratic effects are defined for xi's that are not
        # involved in interactions
        diag.o <- diag(Omega)
        na.c <- colSums(interactions)
        for (i in 2:dim){
            if (!is.na(na.c[i]) & is.na(diag.o[i]))
                stop("Quadratic effects defined for xi's which are not
                involved in any interactions.")
            }
        invisible(NULL)
    } else {
        invisible(NULL)
    }
}

# Creates Beta and Gamma matrices according to the input obtained by
# relation.lat; matrices define relationships for latent variables except
# for interaction effects
rel_lat <- function(x, num.eta, num.xi){

    error.msg <- "Latent variables misspecified. Must be of the form
                    'xi1>eta1' or 'eta1>eta2'. See ?specify_sem for details."

    x.s       <- unlist(strsplit(x, ";"))
    which.xi  <- which(grepl("xi", x.s))
    which.eta <- which(!grepl("xi", x.s))

    G <- matrix(NA, nrow=num.eta, ncol=num.xi)
    for (i in which.xi){
        xi.s <- unlist(strsplit(x.s[i], ">"))
        if (length(xi.s) < 2) stop(error.msg)

        xis  <- unlist(strsplit(xi.s[1], ","))
        etas <- unlist(strsplit(xi.s[2], ","))
        tryCatch({
            ind.xi <- as.numeric(gsub("^.*xi([0-9]+).*$", "\\1", xis))
            ind.eta <- as.numeric(gsub("^.*eta([0-9]+).*$", "\\1", etas))
            G[-ind.eta, -ind.xi] <- 0
        }, error = function(e) stop(error.msg)
        , warning = function(w) stop(error.msg)
        )
    }
    
    B <- diag(1, num.eta)
    for (i in which.eta){
        eta.s <- unlist(strsplit(x.s[i], ">"))
        if (length(eta.s) < 2) stop(error.msg)

        eta.cols <- unlist(strsplit(eta.s[1], ","))
        eta.rows <- unlist(strsplit(eta.s[2], ","))
        tryCatch({
            ind.rows <- as.numeric(gsub("^.*eta([0-9]+).*$", "\\1", eta.rows))
            ind.cols <- as.numeric(gsub("^.*eta([0-9]+).*$", "\\1", eta.cols))
            if (any(sapply(ind.rows, function(x) {is.element(x, ind.cols)})))
                stop(error.msg)
            B[ind.rows, ind.cols] <- NA
        }, error = function(e) stop(error.msg)
        , warning = function(w) stop(error.msg)
        )
    }
    out <- list(Gamma=G, Beta=B)
    out
}

# Defines model class of a given specification; possible output: 'lms',
# 'stemm', 'nsemm'
get_model_class <- function(num.classes, interaction) {
    if (interaction != "") {
        if (num.classes == 1) model.class <- "lms"
        else model.class <- "nsemm"
    } else model.class <- "stemm"
}

# Obtains parameter names from a given model; used in specify_sem
get_parnames <- function(model) {

    par.names <- list()
    for (c in seq_len(model$info$num.classes)) {
        lst <- unlist(lapply(model$matrices[[c]], is.na))
        par.names[[c]] <- names(lst[lst])
    }
    names(par.names) <- paste0("class", 1:model$info$num.classes)
    par.names
}

# Fills upper.tri of a (filled) matrix which should be symmetric
fill_symmetric <- function(mat) {
    for (i in seq_len(nrow(mat))) {
                for (j in i:ncol(mat)) {
                        mat[i,j] <- mat[j,i]
                }
            }
    mat
}

# Vector for diagonal indices
diag_ind <- function(num) diag(matrix(seq_len(num^2), num))

# Set bounds for parameters: variances > 0
bounds <- function(model) {

    # variances to (0, Inf)
    if (class(model) == "lms") {

        lower <- rep(-Inf, count_free_parameters(model))
        upper <- rep(Inf, count_free_parameters(model))

        if (model$info$num.x > 1){
            t.d <- paste0("Theta.d", diag_ind(model$info$num.x))
        } else t.d <- "Theta.d"
        if (model$info$num.y > 1){
            t.e <- paste0("Theta.e", diag_ind(model$info$num.y))
        } else t.e <- "Theta.e"
        if (model$info$num.eta > 1){
            psi <- paste0("Psi", diag_ind(model$info$num.eta))
        } else psi <- "Psi"
        lower[model$info$par.names %in% c(t.d, t.e, psi)] <- 0
        out <- list(upper=upper, lower=lower)
    } else if (class(model) == "stemm" || class(model) == "nsemm") {

        lower <- rep(-Inf, length(model$info$par.names$class1))
        upper <- rep(Inf, length(model$info$par.names$class1))

        lower.class <- list()
        upper.class <- list()

        for (c in seq_len(model$info$num.classes)) {
            if (model$info$num.x > 1){
                t.d <- paste0("Theta.d", diag_ind(model$info$num.x))
            } else t.d <- "Theta.d"
            if (model$info$num.y > 1){
                t.e <- paste0("Theta.e", diag_ind(model$info$num.y))
            } else t.e <- "Theta.e"
            if (model$info$num.eta > 1){
                psi <- paste0("Psi", diag_ind(model$info$num.eta))
            } else psi <- "Psi"
            if (model$info$num.xi > 1){
                phi <- paste0("Phi", diag_ind(model$info$num.eta))
            } else phi <- "Phi"

            lower[model$info$par.names[[c]] %in% c(t.d, t.e, psi, phi)] <- 0
            lower.class[[c]] <- lower
            upper.class[[c]] <- upper
        }
        names(lower.class) <- paste0("class",seq_len(model$info$num.classes))
        names(upper.class) <- paste0("class",seq_len(model$info$num.classes))
        out <- list(upper=upper.class, lower=lower.class)
    }
    out
}

