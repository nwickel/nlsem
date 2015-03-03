# nsemm.R
#
# created: Nov/14/2014, KN
# last mod: Dec/03/2014, NU

#--------------- main functions ---------------

estep_nsemm <- function(model, parameters, data, max.lms, convergence,
                        logger=FALSE, ...) {
    num.classes <- model$info$num.classes

    class.parameters <- get_class_parameters(model, parameters)
    par.new <- NULL

    # lms for each class
    # Note that B is not estimated
    for (c in seq_len(num.classes)) {
        lms.model <- lms_ify(model, c)

        # em for lms
        est <- em(model=lms.model, data=data, start=class.parameters[[c]],
                  logger=logger, neg.hessian=FALSE, max.iter=max.lms,
                  convergence=convergence, ...)

        par.new <- c(par.new, est$coefficients)
    }

    # e-step for semm
    # Note that Omega and A are not estimated
    P <- estep_semm(model=model, parameters=par.new, data=data)
    w.c <- colSums(P) / nrow(data)

    res <- list(P=P, w.c=w.c, par.old=par.new)
    res
}


mstep_nsemm <- function(model, parameters, P, data, optimizer, max.mstep,
                        control=list(), ...) {

    est <- mstep_semm(model=model, parameters=parameters, P=P,
                                data=data, optimizer=optimizer,
                                max.mstep=max.mstep, control=control, ...)

    est
}

#--------------- helper functions ---------------

# create lms model for a specific class of an nsemm model
lms_ify <- function(model, c) {
    lms.model <- list(matrices=list(class1=model$matrices[[c]]),
                      info=model$info)
    lms.model$info$par.names <- model$info$par.names[[c]]
    lms.model$info$bounds$upper <- model$info$bounds$upper[[c]]
    lms.model$info$bounds$lower <- model$info$bounds$lower[[c]]
    lms.model$info$num.classes <- 1
    class(lms.model) = "lms"
    lms.model
}
