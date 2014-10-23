# model.R
#
# created Sep/23/2014, NU
# last mod Oct/23/2014, KN

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

specify_sem <- function(num.x, num.y, num.xi, num.eta, xi, eta, num.groups=1,
                          interaction="all", interc_obs=FALSE,
                          interc_lat=FALSE){

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

    # create data frame with variable names and one column for each group
    # specs <- data.frame(
    #     label = c(paste0("Lambda.x", 1:(num.x*num.xi)), paste0("Lambda.y",
    #         1:(num.y*num.eta)), paste0("Gamma", 1:(num.xi*num.eta)),
    #         paste0("Beta", 1:(num.eta*num.eta)), paste0("Theta.d",
    #         1:(num.x*num.x)), paste0("Theta.e", 1:(num.y*num.y)),
    #         paste0("Psi", 1:(num.eta*num.eta)), paste0("Phi", 1:(num.xi*num.xi)),
    #         paste0("A", 1:(num.xi*num.xi)), paste0("nu.x", 1:num.x),
    #         paste0("nu.y", 1:num.y), paste0("alpha", 1:num.eta),
    #         paste0("tau", 1:num.xi), paste0("Omega", 1:(num.xi*num.xi)))
    # )
    # for (g in seq_len(num.groups)) {
    #     ustart.temp <- data.frame(ustart = 0)
    #     names(ustart.temp) <- paste0("group", g)
    #     specs <- cbind(specs, ustart.temp)
    # }

    # class of model
    if (num.groups == 1) {
        if (interaction == "") {
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
    Gamma    <- matrix(NA, nrow=num.eta, ncol=num.xi)
    # Beta
    Beta     <- diag(1, nrow=num.eta)
    # Theta.d
    Theta.d  <- diag(NA, nrow=num.x)
    # Theta.e
    Theta.e <- diag(NA, nrow=num.y)
    # Psi
    Psi     <- matrix(NA, nrow=num.eta, ncol=num.eta)
    # Phi
    Phi <- matrix(NA, nrow=num.xi, num.xi)
    # A
    A        <- matrix(NA, nrow=num.xi, ncol=num.xi)
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
    tau      <- matrix(0, nrow=num.xi, ncol=1)
    # Omega
    Omega    <- matrix(0, nrow=num.xi, ncol=num.xi)
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
    matrices[[g]] <- list(Lambda.x=Lambda.x, Lambda.y=Lambda.y,
                          Gamma=Gamma, Beta=Beta, Theta.d=Theta.d,
                          Theta.e=Theta.e, Psi=Psi, Phi=Phi, A=A,
                          nu.x=nu.x, nu.y=nu.y, alpha=alpha, tau=tau,
                          Omega=Omega)
    }
    names(matrices) <- paste0("group",1:num.groups)

    # group weights w
    w <- matrix(NA, nrow=num.groups, ncol=1)

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

fill_matrices <- function(dat){
    
    stopifnot(is.data.frame(dat))

    Lambda.x <- as.character(dat$label[grep("Lambda.x", dat$label)])
    Lambda.y <- as.character(dat$label[grep("Lambda.y", dat$label)])
    Gamma    <- as.character(dat$label[grep("Gamma", dat$label)])
    Beta     <- as.character(dat$label[grep("Beta", dat$label)])
    Theta.d  <- as.character(dat$label[grep("Theta.d", dat$label)])
    Theta.e  <- as.character(dat$label[grep("Theta.e", dat$label)])
    Psi      <- as.character(dat$label[grep("Psi", dat$label)])
    A        <- as.character(dat$label[grep("A", dat$label)])
    nu.x     <- as.character(dat$label[grep("nu.x", dat$label)])
    nu.y     <- as.character(dat$label[grep("nu.y", dat$label)])
    alpha    <- as.character(dat$label[grep("alpha", dat$label)])
    tau      <- as.character(dat$label[grep("tau", dat$label)])
    Omega    <- as.character(dat$label[grep("Omega", dat$label)])

    # number of latent and indicator variables and groups
    num.x      <- length(nu.x)
    num.y      <- length(nu.y)
    num.xi     <- length(tau)
    num.eta    <- length(alpha)
    num.groups <- ncol(dat) - 1

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

        matrices[[g]] <- list(Lambda.x=Lambda.x.matrix,
                              Lambda.y=Lambda.y.matrix, Gamma=Gamma.matrix,
                              Beta=Beta.matrix, Theta.d=Theta.d.matrix,
                              Theta.e=Theta.e.matrix, Psi=Psi.matrix,
                              Phi=Phi.matrix, A=A.matrix, nu.x=nu.x.matrix,
                              nu.y=nu.y.matrix, alpha=alpha.matrix,
                              tau=tau.matrix, Omega=Omega.matrix)
    }
    names(matrices) <- paste0("group",1:num.groups)
    matrices
}

as.data.frame.lms <- as.data.frame.stemm <- as.data.frame.nsemm <- function(object, ...) {
    specs <- data.frame(
        label = names(unlist(object$matrices$group1)))
    for (g in seq_len(length(object$matrices))) {
        temp <- data.frame(unlist(object$matrices[[g]], use.names=FALSE))
        names(temp) <- paste0("group", g)
        specs <- cbind(specs, temp)
    }
    specs
}

count_free_parameters <- function(model) {
    # w is not counted
    res <- 0
    for (g in seq_len(model$info$num.groups))
        res <- res + sum(unlist(lapply(model$matrices[[g]], is.na)))
    res
    }

fill_model <- function(model, parameters) {

    stopifnot(class(model) == "lms" || class(model) == "stemm"
              || class(model) == "nsemm")

    stopifnot(count_free_parameters(model) == length(parameters))


    for (g in seq_len(model$info$num.groups)) {
        par.names <- model$info$par.names[[g]]
        matrix.names <- unlist(lapply(par.names, function(x){
                                      gsub("([0-9]*$)", "", x)}))
        index <- unlist(lapply(par.names, function(x){
                               res <- gsub("(^.*[A-Za-z])", "", x)
                               if (res == "") res <- 1
                               as.numeric(res) }))
        # data <- data.frame(matrix.names=matrix.names, index=index,
        #                    parameters=parameters)
        for (i in seq_along(matrix.names)) {
            for (j in seq_along(matrices)) {
                if (names(matrices[j]) == matrix.names[i])
                    matrices[[j]][index[i]] <- parameters[i]
            }
        }
        # TODO cut parameters or adjust index for parameters to group
    }


    # old fill_model
    # specs <- as.data.frame(model)
    # specs[is.na(specs)] <- parameters

    # out_matrices <- fill_matrices(specs, num.x=model$info$num.x,
    #                               num.y=model$info$num.y,
    #                               num.xi=model$info$num.xi,
    #                               num.eta=model$info$num.eta,
    #                               num.groups=model$info$num.groups)
    # names(out_matrices) <- names(model$matrices)

    # out <- list(matrices = out_matrices, info = model$info)

    # switch(class(model),
    #              "lms" = {class(out) <- "lmsFilled"},
    #              "stemm" = {class(out) <- "stemmFilled"},
    #              "nsemm" = {class(out) <- "nsemmFilled"})
    # out
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

