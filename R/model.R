# model.R
#
# created Sep/23/2014, NU
# last mod Oct/14/2014, KN

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

test_omega <- function(omega){

    if (any(is.na(omega))){

        ind <- which(is.na(omega), arr.ind=TRUE)
        dim <- nrow(omega)
        
        if (nrow(ind) == 1){
            if (!is.na(omega[1, dim]))
                stop("Interactions are not well-defined. Please change order of xi's. See ?specify_sem for details.")

        } else {
            for (i in 2:nrow(ind)){
                if(ind[i,1] >= ind[i-1,2]){
                    stop("Interactions are not well-defined. Please change order of xi's. See ?specify_sem for details.")
                } 
            }
            for (i in 1:max(ind[,1])){
                if (max(which(is.na(omega[i,]))) < dim){
                    stop("Interactions are not well-defined. Please change order of xi's. See ?specify_sem for details.")
                }
            }
            if (max(ind[,1]) > 1){
                for (i in 2:max(ind[,1])){
                    if (min(which(is.na(omega[i-1,]))) >= min(which(is.na(omega[i,])))) {
                        stop("Interactions are not well-defined. Please change order of xi's. See ?specify_sem for details.")
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
                          interaction="all", constraints="default",
                          interc_obs=FALSE, interc_lat=FALSE, dataframe=FALSE){

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

    # create list of matrices for each group (therefore index g)
    matrices <- list()
    for (g in seq_len(num.groups)) {
        Lx.matrix    <- matrix(0, nrow=num.x, ncol=num.xi)
        Ly.matrix    <- matrix(0, nrow=num.y, ncol=num.eta)
        G.matrix     <- matrix(0, nrow=num.eta, ncol=num.xi)
        B.matrix     <- matrix(0, nrow=num.eta, ncol=num.eta)
        Td.matrix    <- matrix(0, nrow=num.x, ncol=num.x)
        Te.matrix    <- matrix(0, nrow=num.y, ncol=num.y)
        Psi.matrix   <- matrix(0, nrow=num.eta, ncol=num.eta)
        A.matrix     <- matrix(0, nrow=num.xi, ncol=num.xi)
        vx.matrix    <- matrix(0, nrow=num.x, ncol=1)
        vy.matrix    <- matrix(0, nrow=num.y, ncol=1)
        alpha.matrix <- matrix(0, nrow=num.eta, ncol=1)
        t.matrix     <- matrix(0, nrow=num.xi, ncol=1)
        O.matrix     <- matrix(0, nrow=num.xi, ncol=num.xi)
        temp <- list(Lx=Lx.matrix, Ly=Ly.matrix, G=G.matrix, B=B.matrix,
                     Td=Td.matrix, Te=Te.matrix, Psi=Psi.matrix, A=A.matrix,
                     vx=vx.matrix, vy=vy.matrix, alpha=alpha.matrix,
                     t=t.matrix, O=O.matrix)
        matrices[[g]] <- temp
    }
    names(matrices) <- paste0("group",1:num.groups)

    # create data frame with variable names and one column for each group
    specs <- data.frame(
        label = c(paste0("Lx", 1:(num.x*num.xi)), paste0("Ly",
            1:(num.y*num.eta)), paste0("G", 1:(num.xi*num.eta)), paste0("B",
            1:(num.eta*num.eta)), paste0("Td", 1:(num.x*num.x)), paste0("Te",
            1:(num.y*num.y)), paste0("Psi", 1:(num.eta*num.eta)), paste0("A",
            1:(num.xi*num.xi)), paste0("vx", 1:num.x), paste0("vy", 1:num.y),
            paste0("alpha", 1:num.eta), paste0("t", 1:num.xi), paste0("O",
            1:(num.xi*num.xi)))
    )
    for (g in seq_len(num.groups)) {
        ustart.temp <- data.frame(ustart = 0)
        names(ustart.temp) <- paste0("group", g)
        specs <- cbind(specs, ustart.temp)
    }

    # fill matrices and data frame with constraints or default
    if (is.data.frame(constraints)){

        # check if number of columns in constraints is equal to num.groups
        if (ncol(constraints) != num.groups) {
            stop("Data frame for constraints does not match number of
                 latent groups. See ?specify_sem.")
        }
        # check if number of rows in constraints is correct
        if (nrow(constraints) != nrow(specs)) {
            print(specs)
            stop("Data frame for constraints does not have correct number of
                 rows. See correct number above or see ?specify_sem.")
        }
        # fill specs with constraints
        specs[,2:(num.groups+1)] <- constraints[,1:(num.groups)]

        # fill matrices
        matrices <- fill_matrices(specs)

    } else if (constraints == "default"){
    # default constraints

        xi.s <- unlist(strsplit(xi, ","))
        xi.ind <- list()
        # TODO Test if length(xi.s) == num.xi
        for (i in seq_len(num.xi)) xi.ind[[i]] <- grep_ind(xi.s[i])
        # TODO Use sapply instead of loop!

        #eta.ind <- grep_ind(eta)
        eta.s <- unlist(strsplit(eta, ","))
        eta.ind <- list()
        for (i in seq_len(num.eta)) eta.ind[[i]] <- grep_ind(eta.s[i])
        # TODO Use sapply instead of loop!

        # fill in default constraints
        for (g in seq_len(num.groups)) {

            # Lx
            for (i in seq_len(num.xi)){
                matrices[[g]]$Lx[xi.ind[[i]], i] <- c(1, rep(NA, length(xi.ind[[i]])-1))
            }
            # Ly
            for (i in seq_len(num.eta)){
                matrices[[g]]$Ly[eta.ind[[i]], i] <- c(1, rep(NA, length(eta.ind[[i]])-1))
            }
            # G
            # TODO why not: matrices[[g]]$G[1:num.eta,1:num.xi] <- NA
            matrices[[g]]$G[1:dim(matrices[[g]]$G)[1],1:dim(matrices[[g]]$G)[2]] <- NA
            # Td
            matrices[[g]]$Td <- diag(NA, num.x)
            # Te
            matrices[[g]]$Te <- diag(NA, num.y)
            # Psi
            # TODO why not: matrices[[g]]$Psi[1:num.xi,1:num.xi]
            matrices[[g]]$Psi[1:dim(matrices[[g]]$Psi)[1],1:dim(matrices[[g]]$Psi)[2]] <- NA
            # A
            matrices[[g]]$A[1:dim(matrices[[g]]$A)[1],1:dim(matrices[[g]]$A)[2]] <- NA
            matrices[[g]]$A[upper.tri(matrices[[g]]$A)] <- 0
            # Omega
            if (interaction == "all"){
                matrices[[g]]$O[upper.tri(matrices[[g]]$O)] <- NA
            } else if (interaction != "") {
                interaction.s <- unlist(strsplit(interaction, ","))
                ind <- calc_interaction_matrix(interaction.s)
                matrices[[g]]$O[ind] <- NA
            }
            test_omega(matrices[[g]]$O)
            # check if O has row echelon form
            # nu's
            if (interc_obs){
                matrices[[g]]$vx[1:num.x] <- NA
                matrices[[g]]$vy[1:num.y] <- NA
            }
            # alpha
            if (interc_lat){
                # TODO why not: empty.model[[g]]$alpha[1:num.eta] <- NA
                matrices[[g]]$alpha[1:dim(matrices[[g]]$alpha)[1],1:dim(matrices[[g]]$alpha)[2]] <- NA
            }
            # put constraints in data frame
            specs[,g+1] <- unlist(matrices[[g]])
        }
    } else {
        stop("constraints need to be a data.frame or set to 'default'.")
    }

    # group weights w
    w.matrix <- matrix(NA, nrow=num.groups, ncol=1)

    # create model with matrices and info
    model <- list(matrices=matrices, info=list(num.xi=num.xi, num.eta=num.eta,
                                               num.x=num.x, num.y=num.y,
                                               num.groups=num.groups,
                                               w=w.matrix))
    # TODO add par.names to info

    # class of model
    # TODO testing if there is no interaction should probably be different
    if (num.groups == 1) {
        if (interaction == "") {
            stop("Model needs either more than one latent group or at least one
                 latent interaction (e.g. 'xi1:xi2'). For other models please
                 use lavaan or the like.")
        } else {
            class(model) <- "lms"
        }
    } else if (interaction == "") {
        class(model) <- "stemm"
    } else {
        class(model) <- "nsemm"
    }

    if (dataframe) {
        specs
    } else {
        model
    }
}

fill_matrices <- function(dat){
    
    stopifnot(is.data.frame(dat))

    # grep names
    Lx    <- as.character(dat$label[grep("Lx", dat$label)])
    Ly    <- as.character(dat$label[grep("Ly", dat$label)])
    G     <- as.character(dat$label[grep("G", dat$label)])
    B     <- as.character(dat$label[grep("B", dat$label)])
    Td    <- as.character(dat$label[grep("Td", dat$label)])
    Te    <- as.character(dat$label[grep("Te", dat$label)])
    Psi   <- as.character(dat$label[grep("Psi", dat$label)])
    A     <- as.character(dat$label[grep("A", dat$label)])
    vx    <- as.character(dat$label[grep("vx", dat$label)])
    vy    <- as.character(dat$label[grep("vy", dat$label)])
    alpha <- as.character(dat$label[grep("alpha", dat$label)])
    t     <- as.character(dat$label[grep("t", dat$label)])
    O     <- as.character(dat$label[grep("O", dat$label)])

    # number of latent and indicator variables
    # TODO could also be passed
    num.x      <- length(vx)
    num.y      <- length(vy)
    num.xi     <- length(t)
    num.eta    <- length(alpha)
    num.groups <- ncol(dat) - 1

    # create matrices
    matrices <- list()
    for (g in seq_len(num.groups)) {
        Lx.matrix    <- matrix(dat[dat$label %in% Lx, paste0("group",g)],
                               nrow=num.x, ncol=num.xi)
        Ly.matrix    <- matrix(dat[dat$label %in% Ly, paste0("group",g)],
                               nrow=num.y, ncol=num.eta)
        G.matrix     <- matrix(dat[dat$label %in% G, paste0("group",g)],
                               nrow=num.eta, ncol=num.xi)
        B.matrix     <- matrix(dat[dat$label %in% B, paste0("group",g)],
                               nrow=num.eta, ncol=num.eta)
        Td.matrix    <- matrix(dat[dat$label %in% Td, paste0("group",g)],
                               nrow=num.x, ncol=num.x)
        Te.matrix    <- matrix(dat[dat$label %in% Te, paste0("group",g)],
                               nrow=num.y, ncol=num.y)
        Psi.matrix   <- matrix(dat[dat$label %in% Psi, paste0("group",g)],
                               nrow=num.eta, ncol=num.eta)
        A.matrix     <- matrix(dat[dat$label %in% A, paste0("group",g)],
                               nrow=num.xi, ncol=num.xi)
        vx.matrix    <- matrix(dat[dat$label %in% vx, paste0("group",g)],
                               nrow=num.x, ncol=1)
        vy.matrix    <- matrix(dat[dat$label %in% vy, paste0("group",g)],
                               nrow=num.y, ncol=1)
        alpha.matrix <- matrix(dat[dat$label %in% alpha, paste0("group",g)],
                               nrow=num.eta, ncol=1)
        t.matrix     <- matrix(dat[dat$label %in% t, paste0("group",g)],
                               nrow=num.xi, ncol=1)
        O.matrix     <- matrix(dat[dat$label %in% O, paste0("group",g)],
                               nrow=num.xi, ncol=num.xi)

        matrices[[g]] <- list(Lx=Lx.matrix, Ly=Ly.matrix, G=G.matrix,
                              B=B.matrix, Td=Td.matrix, Te=Te.matrix,
                              Psi=Psi.matrix, A=A.matrix, vx=vx.matrix,
                              vy=vy.matrix, alpha=alpha.matrix, t=t.matrix,
                              O=O.matrix)
    }
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

free_parameters <- function(model) sum(unlist(lapply(model$matrices, is.na)))

fill_model <- function(model, parameters) {
    
    stopifnot(class(model) == "lms")

    stopifnot(free_parameters(model) == length(parameters))

    mod.unlist <- unlist(model$matrices)
    mod.unlist[is.na(mod.unlist)] <- parameters
    
    specs <- data.frame(label=names(mod.unlist), ustart=mod.unlist)

    out <- fill_matrices(specs)
    out$info$par.names <- NULL

    class(out) <- "lmsFilled"

    out
  
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

