# em.R
#
# last mod: Nov/25/2014, NU

# Performs EM-algorithm for different models of class 'lms', 'stemm', and
# soon 'nsemm'
em <- function(model, data, start, logger=FALSE, threshold=1e-03,
                max.iter=40, m=16, optimizer=c("nlminb", "optim"),
                neg.hessian=TRUE, ...) {

    stopifnot(class(model) == "lms" || class(model) == "stemm" ||
              class(model) == "nsemm")

    if (!count_free_parameters(model) == length(start)){
        stop("Number of starting parameters is not equal to number of free
        parameters in model.")
    }

    if (ncol(data) != (model$info$num.x + model$info$num.y)) {
        stop("Number of columns in data does not match number of x's and y's")
    }

    if (class(model) == "lms" || class(model) == "nsemm"){
        n.na <- length(which(is.na(model$matrices$group1$Omega)))
        if (any(start[-c(1:(length(start) - 3))] == 0)){
            stop("Starting parameters for Omega should not be 0.")
        }
    }

    cat("-----------------------------------\n")
    cat("Starting EM-algorithm for", class(model), "\n")
    cat(paste("Threshold: ", threshold, "\n"))
    cat("-----------------------------------\n")
    cat("-----------------------------------\n")
    
    ll.old   <- 0     # loglikelihood of the last iteration
    ll.new   <- 1     # loglikelihood of the current iteration
    ll.ret   <- NULL
    num.iter <- 0     # number of iterations
    if (class(model) == "stemm" || class(model) == "nsemm") {
        par.new <- start
    } else {
        par.new <- convert_parameters_lms(model, start)
    }
    par.old  <- 0
    
    # when function aborts
    on.exit({
        cat("-----------------------------------\n")
        if (neg.hessian == TRUE) {
            cat("Computing Hessian \n")
        } else {
            cat("Computing final model \n")
        }
        cat("-----------------------------------\n")

        switch(class(model),
           "lms" = {
                final <- mstep_lms(model=model, P=P, dat=data,
                                   parameters=par.new, neg.hessian=neg.hessian, m=m,
                                   optimizer=optimizer, ...)
                names(final$par) <- model$info$par.names
                # Transform parameters back to Phi
                A <- matrix(0, nrow=model$info$num.xi, ncol=model$info$num.xi)
                A[lower.tri(A, diag=TRUE)] <- final$par[grep("Phi", names(final$par))]
                Phi <- A %*% t(A)
                final$par[grep("Phi", names(final$par))] <- Phi[lower.tri(Phi, diag=TRUE)] 
            },
            "stemm" = {
                final <- mstep_stemm(model=model, parameters=par.old, P=P,
                                     data=data, neg.hessian=neg.hessian,
                                     optimizer=optimizer, ...)
                if (is.numeric(final$par)) {
                    par.names <- NULL
                    for (g in seq_len(model$info$num.groups)) {
                        par.names <- c(par.names,
                                       paste0("group", g, ".", model$info$par.names[[g]]))
                    }
                    names(final$par) <- par.names
                } else {
                    for (g in seq_len(model$info$num.groups)) {
                        names(final$par[[g]]) <- model$info$par.names[[g]]
                    }
                }
            },
            "nsemm" = {
                final <- mstep_nsemm(model=model, parameters=par.old, P=P,
                                     data=data, neg.hessian=neg.hessian,
                                     optimizer=optimizer, ...)
                if (is.numeric(final$par)) {
                    par.names <- NULL
                    for (g in seq_len(model$info$num.groups)) {
                        par.names <- c(par.names,
                                       paste0("group", g, ".", model$info$par.names[[g]]))
                    }
                    names(final$par) <- par.names
                } else {
                    for (g in seq_len(model$info$num.groups)) {
                        names(final$par[[g]]) <- model$info$par.names[[g]]
                    }
                }
                # -> TODO check if par is not already named!
            }
        )

        # in case break happens before first m-step
        if (is.null(ll.ret)) {ll.ret <- final$objective}

        out <- list(model.class=class(model), coefficients=final$par,
                    objective=-final$objective,
                    convergence_final_step=final$convergence,
                    message_final_step=final$message,
                    negHessian=final$hessian,
                    loglikelihoods=-ll.ret,
                    info=model$info[1:4])

        # attach w for stemm and nsemm
        if (class(model) == "stemm" || class(model) == "nsemm") {
            out$info <- model$info[c(1:4,7)]
        }

        class(out) <- "emEst"
        return(out)
    })

    while(abs(ll.old - ll.new) > threshold) { # as long as no convergence reached
    #while(sum((par.old - par.new)^2) > threshold) { # as long as no convergence reached
        if(ll.new - ll.old > 0.001 && num.iter > 3) {
            warning("Likelihood should be decreasing")
        }

        if(logger == TRUE) {
            cat(paste("Iteration", num.iter+1, "\n"))
            cat("Doing expectation-step \n")
        }

        # Update loglikelihood
        ll.old <- ll.new
        par.old <- par.new

        # E-step
        switch(class(model),
           "lms" = {
                names(model$matrices$group1)[grep("Phi", names(model$matrices$group1))] <- "A"
                P <- estep_lms(model=model, parameters=par.old, dat=data, m=m, ...)
            },
           "stemm" = {
                P <- estep_stemm(model=model, parameters=par.old, data=data)
                model$info$w <- colSums(P) / nrow(data)
                if (logger == TRUE) {
                    cat("Group weights: ", round(model$info$w, digits=4), "\n")
                }
            },
            "nsemm" = {
                res <- estep_nsemm(model=model, parameters=par.old, data=data,
                                   logger=logger, ...)
                P            <- res$P
                model$info$w <- res$w.g
                par.old      <- res$par.old
                if (logger == TRUE) {
                    cat("Group weights: ", round(model$info$w, digits=4), "\n")
                }
            }
        )
  
        if(logger == TRUE){
            cat("Doing maximization-step \n")
        }
        
        # M-step
        switch(class(model),
            "lms" = {
                m.step <- mstep_lms(model=model, P=P, dat=data, parameters=par.old,
                                m=m, optimizer=optimizer, ...)
            },
            "stemm" = {
                m.step <- mstep_stemm(model=model, parameters=par.old, P=P,
                                  data=data, optimizer=optimizer, ...)
            },
            "nsemm" = {
                m.step <- mstep_nsemm(model=model, parameters=par.old, P=P,
                                  data=data, optimizer=optimizer, ...)
            }
        )

        if(logger == TRUE) {
            cat("Results of maximization \n")
            cat(paste0("Final loglikelihood: ", round(-m.step$objective, 3), "\n"))
            cat(paste0("Convergence: ", m.step$convergence, "\n"))
            cat(paste0("Number of iterations: ", m.step$iterations, "\n"))
            cat("----------------------------------- \n")
        }
      
        ll.new     <- m.step$objective
        ll.ret     <- c(ll.ret, ll.new)
        par.new    <- unlist(m.step$par)
        num.iter   <- num.iter + 1
  
        if(num.iter == max.iter){
            warning("Maximum number of iterations was reached. EM algorithm might not have
            converged.")
            break
        }
    }
    cat("-----------------------------------\n")
    cat("EM completed \n")
    cat(paste0("Previous loglikelihood: ", round(-ll.old, 3), "\n"))
    cat(paste0("Final loglikelihood: ", round(-ll.new, 3),"\n"))
    cat("-----------------------------------\n")

    # When EM is completed, on.exit above is called
    # TOTHINK should official version have the on.exit functionality? If not,
    # move code from on.exit to here

}


