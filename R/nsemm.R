# nsemm.R
#
# created: Nov/14/2014, KN
# last mod: Nov/14/2014, KN

estep_nsemm <- function(model, parameters, data, ...) {
    num.groups <- model$info$num.groups
    # TODO more args needed here?

    par.old <- parameters
    par.new <- NULL
    # lms for each group
    # Note that B is not estimated
    for (g in seq_len(num.groups)) {
        if (logger == TRUE) {
            cat("LMS for group ", g, "\n")
        }
        lms.model <- list(matrices=list(group1=model$matrices[[g]]),
                          info=model$info)
        lms.model$info$par.names <- model$info$par.names[[g]]
        lms.model$info$bounds <- model$info$bounds[[g]]
        class(lms.model) = "lms"
        group.parameters <- par.old[0:length(lms.model$info$par.names)]
        par.old <- par.old[(length(lms.model$info$par.names) + 1):length(par.old)]

        # em for lms
        est <- em(model=lms.model, data=data, start=group.parameters, ...)

        # get new values for Phi
        par.new[grep("Phi", par.new)] <- coef(est)[grep("Phi", names(coef(est)))]

        par.new <- c(par.new, est$par)
    }

    # e-step for stemm
    # Note that Omega and A are not estimated









    # TODO calculate A from Phi
}
