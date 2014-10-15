
em <- function(model, data, start, logger=FALSE, threshold=1e-05,
                max.iter=40, m=16, ...) {

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
    
    tryCatch({
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
        P      <- estep_lms(model=model, parameters=par.old, dat=data, m=m, ...)
  
        if(logger == TRUE){
            cat("Doing maximization-step \n")
        }
        
        # M-step
        m.step <- mstep_lms(model=model, P=P, dat=data, parameters=par.old, m=m, ...)
  
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
    cat("-----------------------------------\n")
    cat("Computing Hessian \n")
    cat("-----------------------------------\n")
    
    # Compute hessian of final parameters
    final <- mstep_lms(model=model, P=P, dat=data, parameters=par.new, Hessian=TRUE, m=m, ...)

    names(final$par) <- model$info$par.names
    
    out <- list(par=final$par, objective=-final$objective,
           convergence_final_step=final$convergence,
           message_final_step=final$message, Hessian=final$hessian$Hessian,
           gradient=final$hessian$gradient, loglikelihoods=-ll.ret,
           info=model$info[1:4])
  
    class(out) <- "emRes"
    }, error=function(e){out <- list(par=par.new,
                         objective=m.step$objective, loglikelihoods=-ll.ret)})
                            
  
    out
}


