# model.R
#
# created Sep/23/2014, NU
# last mod Oct/08/2014, NU

grep_ind <- function(x){
    
    if (length(strsplit(x, "-")[[1]]) > 1){
        as.numeric(gsub("^x([0-9]+).*x([0-9]+)$", "\\1",
        x)):as.numeric(gsub("^x([0-9]+).*x([0-9]+)$", "\\2", x))
    } else {
        as.numeric(gsub("^.([0-9]).*$", "\\1", x))
    }
}

interaction_matrix <- function(x){
    
    rows <- as.numeric(gsub("^xi([0-9]+).*$", "\\1", x))
    cols <- as.numeric(gsub("^xi.*:xi([0-9]+)$", "\\1", x))
    mat  <- cbind(rows, cols)
    mat

}

specify_model <- function(num.x, num.y, num.xi, num.eta, xi, eta,
                          interaction="all", constraints="default",
                          interc_obs=FALSE, interc_lat=FALSE){

    # check arguments
    if (!is.numeric(num.x) || !is.numeric(num.y) || !is.numeric(num.xi) 
       || !is.numeric(num.eta)) { 
        stop("Number of variables must be numeric.")
    } else if (num.x < num.xi || num.y < num.eta) {
        stop("The model contains not enough observed variables.")
    }
    stopifnot(num.x > 0, num.y >= 0, num.xi > 0, num.eta >= 0)
    
    if (is.numeric(constraints)){

        # create data frame with variable names and constraints
        specs <- data.frame(
        label = c(paste0("Lx", 1:(num.x*num.xi)), paste0("Ly",
            1:(num.y*num.eta)), paste0("G", 1:(num.xi*num.eta)), paste0("B",
            1:(num.eta*num.eta)), paste0("Td", 1:(num.x*num.x)), paste0("Te",
            1:(num.y*num.y)), paste0("Psi", 1:(num.eta*num.eta)), paste0("A",
            1:(num.xi*num.xi)), paste0("vx", 1:num.x), paste0("vy", 1:num.y),
            paste0("alpha", 1:num.eta), paste0("t", 1:num.xi), paste0("O",
            1:(num.xi*num.xi))),
        ustart = constraints)
    } else if (constraints == "default"){
    # default constraints
        
        specs <- data.frame(
        label = c(paste0("Lx", 1:(num.x*num.xi)), paste0("Ly",
            1:(num.y*num.eta)), paste0("G", 1:(num.xi*num.eta)), paste0("B",
            1:(num.eta*num.eta)), paste0("Td", 1:(num.x*num.x)), paste0("Te",
            1:(num.y*num.y)), paste0("Psi", 1:(num.eta*num.eta)), paste0("A",
            1:(num.xi*num.xi)), paste0("vx", 1:num.x), paste0("vy", 1:num.y),
            paste0("alpha", 1:num.eta), paste0("t", 1:num.xi), paste0("O",
            1:(num.xi*num.xi))),
        ustart = 0)

        xi.s <- strsplit(xi, ",")[[1]]
        xi.ind <- list()
        for (i in seq_len(num.xi)) xi.ind[[i]] <- grep_ind(xi.s[i])
        # TODO Use sapply instead of loop!

        #eta.ind <- grep_ind(eta)
        eta.s <- strsplit(eta, ",")[[1]]
        eta.ind <- list()
        for (i in seq_len(num.eta)) eta.ind[[i]] <- grep_ind(eta.s[i])
        # TODO Use sapply instead of loop!
        
        # create empty model matrices
        empty.model <- fill_matrices(specs)$matrices

        # fill in default constraints
        # Lx
        for (i in seq_len(num.xi)){
            empty.model$Lx[xi.ind[[i]], i] <- c(1, rep(NA, length(xi.ind[[i]])-1))
        }
        # Ly
        #empty.model$Ly[eta.ind] <- c(1, rep(NA, length(eta.ind)-1))
        for (i in seq_len(num.eta)){
            empty.model$Ly[eta.ind[[i]], i] <- c(1, rep(NA, length(eta.ind[[i]])-1))
        }
        # G
        empty.model$G[1:dim(empty.model$G)[1],1:dim(empty.model$G)[2]] <- NA
        # Td
        empty.model$Td <- diag(NA, num.x)
        # Te
        empty.model$Te <- diag(NA, num.y)
        # Psi
        empty.model$Psi[1:dim(empty.model$Psi)[1],1:dim(empty.model$Psi)[2]] <- NA
        # A
        empty.model$A[1:dim(empty.model$A)[1],1:dim(empty.model$A)[2]] <- NA
        empty.model$A[upper.tri(empty.model$A)] <- 0
        # Omega
        if (interaction == "all"){
            empty.model$O[upper.tri(empty.model$O, diag=TRUE)] <- NA
        } else {
            interaction.s <- strsplit(interaction, ",")[[1]]
            ind <- interaction_matrix(interaction.s)
            empty.model$O[ind] <- NA
        }
        
        # nu's
        if (interc_obs == TRUE){
            empty.model$vx[1:num.x] <- NA
            empty.model$vy[1:num.y] <- NA
        }
        # alpha
        if (interc_lat == TRUE){
            empty.model$alpha[1:dim(empty.model$alpha)[1],1:dim(empty.model$alpha)[2]] <- NA
        }

        # put constraints in data frame
        specs$ustart <- unlist(empty.model)
    } else {
        stop("constraints need to be a numeric vector or set to 'default'.")
    }
    specs
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
    num.x    <- length(vx)
    num.y    <- length(vy)
    num.xi   <- length(t)
    num.eta  <- length(alpha)

    # create matrices
    Lx    <- matrix(dat[dat$label %in% Lx, "ustart"], nrow=num.x, ncol=num.xi)
    Ly    <- matrix(dat[dat$label %in% Ly, "ustart"], nrow=num.y, ncol=num.eta)
    G     <- matrix(dat[dat$label %in% G, "ustart"], nrow=num.eta, ncol=num.xi)
    B     <- matrix(dat[dat$label %in% B, "ustart"], nrow=num.eta, ncol=num.eta)
    Td    <- matrix(dat[dat$label %in% Td, "ustart"], nrow=num.x, ncol=num.x)
    Te    <- matrix(dat[dat$label %in% Te, "ustart"], nrow=num.y, ncol=num.y)
    Psi   <- matrix(dat[dat$label %in% Psi, "ustart"], nrow=num.eta, ncol=num.eta)
    A     <- matrix(dat[dat$label %in% A, "ustart"], nrow=num.xi, ncol=num.xi)
    vx    <- matrix(dat[dat$label %in% vx, "ustart"], nrow=num.x, ncol=1)
    vy    <- matrix(dat[dat$label %in% vy, "ustart"], nrow=num.y, ncol=1)
    alpha <- matrix(dat[dat$label %in% alpha, "ustart"], nrow=num.eta, ncol=1)
    t     <- matrix(dat[dat$label %in% t, "ustart"], nrow=num.xi, ncol=1)
    O     <- matrix(dat[dat$label %in% O, "ustart"], nrow=num.xi, ncol=num.xi)
    
    out <- list(matrices=list(Lx=Lx, Ly=Ly, G=G, B=B, Td=Td, Te=Te, Psi=Psi, A=A,
                vx=vx, vy=vy, alpha=alpha, t=t, O=O), info=list(num.xi=num.xi,
                num.eta=num.eta, num.x=num.x, num.y=num.y,
                par.names=as.character(dat$label[is.na(dat$ustart)])))

    class(out) <- "lms"

    out
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

