
em <- function(model, data, start, logger=FALSE, threshold=10e-7,
                max.iter=40, m=16, ...) {

    cat("-----------------------------------\n")
    cat("Starting EM-algorithm\n")
    cat(paste("Threshold: ", threshold, "\n"))
    cat("-----------------------------------\n")
    cat("-----------------------------------\n")
    
    ll.old   <- 0     # likelihood of the last iteration
    ll.new   <- 1     # likelihood of the current iteration
    ll.ret   <- NULL
    num.iter <- 0     # number of iterations
    
    while(abs(ll.old - ll.new) > threshold) { # as long as no convergence reached
        if(ll.new - ll.old > 0.001 && num.iter > 3) {
            warning("Likelihood should be decreasing")
        }
      
        if(logger == TRUE) {
            cat(paste("Iteration", num.iter+1, "\n"))
            cat("Doing expectation-step \n")
        }
  
        # Update likelihood
        ll.old <- ll.new

        # E-step
        P      <- estep_lms(model=model, parameters=start, dat=data, m=m, ...)
  
        if(logger == TRUE){
            cat("Doing maximization-step \n")
        }
        
        # M-step
        m.step <- mstep_lms(model=model, P=P, dat=data, parameters=start, m=m, ...)
  
        if(logger == TRUE) {
            cat("Results of maximization \n")
            cat(paste0("Final likelihood: ", m.step$objective, "\n"))
            cat(paste0("Convergence: ", m.step$convergence, "\n"))
            cat(paste0("Number of iterations: ", m.step$iterations, "\n"))
            cat("----------------------------------- \n")
        }
      
        ll.new     <- m.step$objective
        ll.ret     <- c(ll.ret, ll.new)
        start      <- m.step$par
        num.iter   <- num.iter + 1
  
        if(num.iter == max.iter) break
    }

    cat("-----------------------------------\n")
    cat("EM completed \n")
    cat(paste0("Previous Likelihood: ", ll.old, "\n"))
    cat(paste0("Final Likelihood: ", ll.new,"\n"))
    cat("-----------------------------------\n")
    cat("-----------------------------------\n")
    cat("Computing Hessian \n")
    cat("-----------------------------------\n")
    
    # Compute hessian of final parameters
    final <- mstep_lms(model=model, P=P, dat=data, parameters=start, Hessian=TRUE, m=m, ...)

    names(final$par) <- model$info$par.names
    
    out <- list(par=final$par, objective=-final$objective,
           convergence=final$convergence, message=final$message,
           Hessian=final$hessian$Hessian, gradient=final$hessian$gradient,
           likelihoods=-ll.ret, info=model$info[1:4])
  
    class(out) <- "emRes"
  
    out
}


