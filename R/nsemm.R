# nsemm.R
#
# created: Nov/14/2014, KN
# last mod: Nov/27/2014, KN

#--------------- main functions ---------------

estep_nsemm <- function(model, parameters, data, logger, ...) {
    num.groups <- model$info$num.groups

    group.parameters <- get_group_parameters(model, parameters)
    par.new <- NULL

    # lms for each group
    # Note that B is not estimated
    for (g in seq_len(num.groups)) {
        if (logger == TRUE) {
            cat("LMS for group ", g, "\n")
        }
        lms.model <- lms_ify(model, g)

        # em for lms
        est <- em(model=lms.model, data=data, start=group.parameters[[g]],
                  logger=logger, Hessian=FALSE, ...)

        par.new <- c(par.new, est$coefficients)
    }

    # e-step for stemm
    # Note that Omega and A are not estimated
    P <- estep_stemm(model=model, parameters=par.new, data=data)
    w.g <- colSums(P) / nrow(data)

    res <- list(P=P, w.g=w.g, par.old=par.new)
    res
}


mstep_nsemm <- function(model, parameters, P, data, optimizer, ...) {
    est <- mstep_stemm(model=model, parameters=parameters, P=P,
                                data=data, optimizer=optimizer, ...)

    cat("===================================\n")
    est
}

#--------------- helper functions ---------------

# create lms model for a specific group of an nsemm model
lms_ify <- function(model, g) {
    lms.model <- list(matrices=list(group1=model$matrices[[g]]),
                      info=model$info)
    lms.model$info$par.names <- model$info$par.names[[g]]
    lms.model$info$bounds$upper <- model$info$bounds$upper[[g]]
    lms.model$info$bounds$lower <- model$info$bounds$lower[[g]]
    lms.model$info$num.groups <- 1
    class(lms.model) = "lms"
    lms.model
}
