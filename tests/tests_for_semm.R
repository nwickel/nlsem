# Testing functions for semm
# ===========================

library(nlsem)

# --> TODO fix tests

# # create model
# # ------------
# model <- specify_sem(num.x=4, num.y=4, num.xi=2, num.eta=2,
#                      xi="x1-x2,x3-x4", eta="y1-y2,y3-y4", num.classes=3,
#                      interaction="", interc.obs=FALSE, interc.lat=FALSE,
#                      relation.lat="xi1>eta1; xi2>eta2")
# # set.seed(127)
# pars.orig <- runif(count_free_parameters(model), 0.1, 1.5)
# data <- simulate(model, parameters=pars.orig)
# 
# parameters <- runif(count_free_parameters(model), 0.1, 1.5)
# system.time(
#     res.nlminb <- em(model, data, parameters, logger=TRUE, optimizer="nlminb")
# )
# system.time(
#     res.optim <- em(model, data, parameters, logger=TRUE, optimizer="optim")
# )
# 
# EMSEM example: STEMM model for structural equation models
# =========================================================
mod <- specify_sem(num.x=4, num.y=4, num.xi=2, num.eta=2,
                     xi="x1-x2,x3-x4", eta="y1-y2,y3-y4", num.classes=2,
                     interaction="",
                     interc.obs=TRUE, interc.lat=FALSE,
                     relation.lat="eta1~xi1;eta2~xi2;eta2~eta2;eta1~eta2")

dat <- as.data.frame(mod)

dat[c(2,8,10,16),2:3] <- 1    # Lambda.x.21 & 42 and Lambda.y.21 & 42
dat[c(58,62),2:3]     <- 0    # psi.21 and phi.21
dat[75:76,2:3]        <- NA   # tau

model <- create_sem(dat)
parameters <- c(
                # class 1
                rep(0.5, 2),    # Gamma
                c(-0.3, 0.7),   # Beta
                rep(0.5, 10),   # Theta.d, Theta.e, Psi
                rep(1, 4),      # Phi, tau
                # class 2
                rep(-0.5, 2),   # Gamma
                c(0.7, -0.3),   # Beta
                rep(0.5, 10),   # Theta.d, Theta.e, Psi
                rep(1, 2),      # Phi
                rep(4, 2)       # tau
)
parameters <- runif(count_free_parameters(model), 0.1, 1.5)
data <- simulate(model, parameters=parameters)

# system.time({
#     res <- em(model, data, parameters, logger=TRUE, optimizer="nlminb")
# })

# small model
# ------------
model_semm <- specify_sem(num.x=4, num.y=2, num.xi=2, num.eta=1,
                     xi="x1-x2,x3-x4", eta="y1-y2", num.classes=2,
                     interaction="", interc.obs=FALSE, interc.lat=FALSE,
                     relation.lat="eta1~xi1+xi2")

pars.orig_semm <- c(
                     # class 1
                     1.5, 1,        # Lambda.x
                     0.7,           # Lambda.y
                     0.2, 0.3,      # Gamma
                     rep(0.5, 4),   # Theta.d
                     0.5, 0.5,      # Theta.e
                     0.5,           # Psi
                     rep(0.8, 3),   # Phi
                     # class 2
                     1, 1.2,        # Lambda.x
                     1.4,           # Lambda.y
                     0.8, 0.9,      # Gamma
                     rep(0.5, 4),   # Theta.d
                     0.5, 0.5,      # Theta.e
                     0.5,           # Psi
                     rep(0.8, 3)    # Phi
)

data_semm <- simulate(model_semm, parameters=pars.orig_semm)

set.seed(3)
parameters_semm <- runif(count_free_parameters(model_semm), 0.1, 1.5)


# system.time({
#     res_semm <- em(model_semm, data_semm, parameters_semm, logger=TRUE)
# })

