# Testing functions for stemm
# ===========================

# create model
# ------------
model <- specify_sem(num.x=4, num.y=4, num.xi=2, num.eta=2,
                     xi="x1-x2,x3-x4", eta="y1-y2,y3-y4", num.groups=3,
                     interaction="", interc_obs=FALSE, interc_lat=FALSE,
                     relation_lat="xi1>eta1; xi2>eta2")

# parameters for testing so sigma is positive definite for all groups
parameters <- c(0.5158426, 0.2671909, 1.1394707, 0.9567445, 1.4950897,
                1.2707014, 0.6446239, 0.9593543, 0.3433933, 1.3785417,
                0.9340849, 1.1918894, 1.4115227, 0.7500661, 0.2658292,
                0.3880807, 0.9806666, 0.7433136, 0.2214245, 1.1364833,
                0.5087990, 1.1211797, 0.3396354, 0.6745136, 1.3445510,
                0.4380809, 0.1933387, 0.7142123, 0.5302632, 0.5673774,
                0.2291646, 0.9650112, 0.5051136, 0.4465713, 1.1487853,
                0.5149691, 0.5705360, 0.3755698, 0.8203366, 1.4476701,
                0.9771912, 0.5439453, 0.8588292, 1.0953211, 1.3972953,
                1.4726688, 0.4516402, 1.0838227, 0.2123353, 0.5227697,
                0.7848908, 1.1258676, 0.6165608, 0.4482860, 1.3358925,
                1.4823573, 1.3351624, 0.1592073, 0.2160227, 0.7526599)
# parameters <- runif(count_free_parameters(model), 0.1, 1.5)

data <- simulate(model, parameters)
# P <- estep_stemm(model, parameters, data)

system.time({
    res.nlminb <- em(model, data, parameters, logger=TRUE, optimizer="nlminb")
})
system.time({
res.optim <- em(model, data, parameters, logger=TRUE, optimizer="optim")
})

# EMSEM example: STEMM model for structural equation models
# =========================================================
mod <- specify_sem(num.x=4, num.y=4, num.xi=2, num.eta=2,
                     xi="x1-x2,x3-x4", eta="y1-y2,y3-y4", num.groups=2,
                     interaction="",
                     interc_obs=FALSE, interc_lat=FALSE)
dat <- as.data.frame(mod)

dat[c(2,8,10,16),2:3] <- 1    # Lambda.x.21 & 42 and Lambda.y.21 & 42
dat[18:19,2:3]        <- 0    # gamma.12 and gamma.21
dat[22:23,2:3]        <- NA   # beta.12 and beta.21
dat[c(58,62),2:3]     <- 0    # psi.21 and phi.21
dat[65:72,2:3]        <- 1    # nu.x and nu.y
dat[75:76,2:3]        <- NA   # tau

model <- create_sem(dat)
parameters <- c(
                # group 1
                rep(0.5, 2),    # Gamma
                c(-0.3, 0.7),   # Beta
                rep(0.5, 10),   # Theta.d, Theta.e, Psi
                rep(1, 4),      # Phi, tau
                # group 2
                rep(-0.5, 2),   # Gamma
                c(0.7, -0.3),   # Beta
                rep(0.5, 10),   # Theta.d, Theta.e, Psi
                rep(1, 2),      # Phi
                rep(4, 2)       # tau
)
# add noise
# parameters <- parameters + rnorm(count_free_parameters(model), 0, 0.3)

# constrain lower bounds for variances to 0
lower <- rep(-Inf, count_free_parameters(model))
    start.index <- 0
    for (g in seq_len(model$info$num.groups)) {
        if (g == 1) {
            indices <- c(grep("Theta", model$info$par.names[[g]]),
                         grep("Psi", model$info$par.names[[g]]),
                         grep("Phi", model$info$par.names[[g]]))
        } else {
            start.index <- start.index + length(model$info$par.names[[g-1]])
            new.indices <- start.index + c(grep("Theta", model$info$par.names[[g]]),
                                           grep("Psi", model$info$par.names[[g]]),
                                           grep("Phi", model$info$par.names[[g]]))
            indices <- c(indices, new.indices)
        }
    }
    lower[indices] <- 0
model$info$bounds$lower <- lower
# constrain upper bounds for variances to 1
upper <- rep(Inf, count_free_parameters(model))
# upper[indices] <- 1
# this leads to: Error in solve.default(matrices$Beta) : system is
# computationally singular: reciprocal condition number = 0
model$info$bounds$upper <- upper

# model.filled <- fill_model(model, parameters)
# data <- simulate(model.filled)
data <- as.matrix(read.table("stemm_data", header=TRUE))
P <- estep_stemm(model, parameters, data)
LL <- loglikelihood_stemm(parameters, model, data, P)

system.time({
res.nlminb <- em(model, data, parameters, logger=TRUE, optimizer="nlminb")
})
system.time({
res.optim <- em(model, data, parameters, logger=TRUE, optimizer="optim")
})

# Estimating the model above
# on HPC:
#    user  system elapsed 
# 441.125   0.426 438.858 
# on macserver
#    User      System verstrichen 
# 260.083       4.278     264.511 


# small model
# ------------
model <- specify_sem(num.x=4, num.y=2, num.xi=2, num.eta=1,
                     xi="x1-x2,x3-x4", eta="y1-y2", num.groups=2,
                     interaction="", interc_obs=FALSE, interc_lat=FALSE,
                     relation_lat="xi1,xi2>eta1")

parameters <- runif(count_free_parameters(model), 0.1, 1.5)

data <- simulate(model, parameters)
P <- estep_stemm(model, parameters, data)

system.time({
    res.nlminb <- em(model, data, parameters, logger=TRUE, optimizer="nlminb")
})
system.time({
res.optim <- em(model, data, parameters, logger=TRUE, optimizer="optim")
})


