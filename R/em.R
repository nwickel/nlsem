# em.R
#
# last mod: Nov/11/2014, NU

# Performs EM-algorithm for different models of class 'lms', 'stemm', and
# soon 'nsemm'
em <- function(model, data, start, logger=FALSE, threshold=1e-05,
                max.iter=40, m=16, optimizer=c("nlminb", "optim"), ...) {

    if (!count_free_parameters(model) == length(start)){
        stop("Number of starting parameters is not equal to number of free
        parameters in model.")
    }

    if (class(model) == "lms"){
        n.na <- length(which(is.na(model$matrices$group1$Omega)))
        if (any(start[-c(1:(length(start) - 3))] == 0)){
            stop("Starting parameters for Omega should not be 0.")
        }
    }

    cat("-----------------------------------\n")
    cat("Starting EM-algorithm\n")
    cat(paste("Threshold: ", threshold, "\n"))
    cat("-----------------------------------\n")
    cat("-----------------------------------\n")
    
    ll.old   <- 0     # loglikelihood of the last iteration
    ll.new   <- 1     # loglikelihood of the current iteration
    ll.ret   <- NULL
    num.iter <- 0     # number of iterations
    par.new  <- start
    par.old  <- 0
    
    # when function aborts
    on.exit({
        cat("-----------------------------------\n")
        cat("Computing Hessian \n")
        cat("-----------------------------------\n")


        if (class(model) == "lms") {
            final <- mstep_lms(model=model, P=P, dat=data, parameters=par.new,
                               Hessian=TRUE, m=m, optimizer=optimizer, ...)
            names(final$par) <- model$info$par.names
        } else if (class(model) == "stemm") {
            final <- mstep_stemm(model=model, parameters=par.old, P=P, data=data,
                                 Hessian=TRUE, optimizer=optimizer, ...)
            par.names <- NULL
            for (g in seq_len(model$info$num.groups)) {
                par.names <- c(par.names, paste0("group", g, ".",
                                                 model$info$par.names[[g]]))
            }
            names(final$par) <- par.names
        }

        # in case break happens before first m-step
        if (is.null(ll.ret)) {ll.ret <- final$objective}

        out <- list(model.class=class(model), par=final$par,
                    objective=-final$objective,
                    convergence_final_step=final$convergence,
                    message_final_step=final$message,
                    Hessian=final$hessian,
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
        if (class(model) == "lms") {
            P <- estep_lms(model=model, parameters=par.old, dat=data, m=m, ...)
        } else if (class(model) == "stemm") {
            P <- estep_stemm(model=model, parameters=par.old, data=data)
            w.g <- colSums(P) / nrow(data)
            model$info$w <- w.g
            # TODO this changes w from "matrix" to "numeric". Problematic?
            if (logger == TRUE) {
                cat("Group weights: ", round(w.g, digits=4), "\n")
            }
        } else {
            stop("E-step not implemented for other classes than lms or stemm")
        }
  
        if(logger == TRUE){
            cat("Doing maximization-step \n")
        }
        
        # M-step
        if (class(model) == "lms") {
            m.step <- mstep_lms(model=model, P=P, dat=data, parameters=par.old,
                                m=m, optimizer=optimizer, ...)
        } else if (class(model) == "stemm") {
            m.step <- mstep_stemm(model=model, parameters=par.old, P=P,
                                  data=data, optimizer=optimizer, ...)
        } else {
            stop("M-step not implemented for other classes than lms or stemm")
        }


        if(logger == TRUE) {
            cat("Results of maximization \n")
            cat(paste0("Final loglikelihood: ", round(-m.step$objective, 3), "\n"))
            cat(paste0("Convergence: ", m.step$convergence, "\n"))
            cat(paste0("Number of iterations: ", m.step$iterations, "\n"))
            cat("----------------------------------- \n")
        }
      
        ll.new     <- m.step$objective
        ll.ret     <- c(ll.ret, ll.new)
        par.new    <- m.step$par
        num.iter   <- num.iter + 1
  
        if(num.iter == max.iter) break
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


