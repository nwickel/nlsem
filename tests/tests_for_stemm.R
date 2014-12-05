# Testing functions for stemm
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
# # EMSEM example: STEMM model for structural equation models
# # =========================================================
# mod <- specify_sem(num.x=4, num.y=4, num.xi=2, num.eta=2,
#                      xi="x1-x2,x3-x4", eta="y1-y2,y3-y4", num.classes=2,
#                      interaction="",
#                      interc.obs=FALSE, interc.lat=FALSE)
# dat <- as.data.frame(mod)
# 
# dat[c(2,8,10,16),2:3] <- 1    # Lambda.x.21 & 42 and Lambda.y.21 & 42
# dat[18:19,2:3]        <- 0    # gamma.12 and gamma.21
# dat[22:23,2:3]        <- NA   # beta.12 and beta.21
# dat[c(58,62),2:3]     <- 0    # psi.21 and phi.21
# dat[65:72,2:3]        <- 1    # nu.x and nu.y
# dat[75:76,2:3]        <- NA   # tau
# 
# model <- create_sem(dat)
# parameters <- c(
#                 # class 1
#                 rep(0.5, 2),    # Gamma
#                 c(-0.3, 0.7),   # Beta
#                 rep(0.5, 10),   # Theta.d, Theta.e, Psi
#                 rep(1, 4),      # Phi, tau
#                 # class 2
#                 rep(-0.5, 2),   # Gamma
#                 c(0.7, -0.3),   # Beta
#                 rep(0.5, 10),   # Theta.d, Theta.e, Psi
#                 rep(1, 2),      # Phi
#                 rep(4, 2)       # tau
# )
# parameters <- runif(count_free_parameters(model), 0.1, 1.5)
# data <- simulate(model, parameters=parameters)
# 
# system.time({
#     res <- em(model, data, parameters, logger=TRUE, optimizer="nlminb",
#                  control=list(iter.max=1))
# })
# 
# # small model
# # ------------
# model_stemm <- specify_sem(num.x=4, num.y=2, num.xi=2, num.eta=1,
#                      xi="x1-x2,x3-x4", eta="y1-y2", num.classes=2,
#                      interaction="", interc.obs=FALSE, interc.lat=FALSE,
#                      relation.lat="xi1,xi2>eta1")
# 
# set.seed(3249)
# parameters_stemm <- runif(count_free_parameters(model_stemm), 0.1, 1.5)
# 
# data_stemm <- simulate(model_stemm, parameters=parameters_stemm)
# 
# system.time({
#     res_stemm <- em(model_stemm, data_stemm, parameters_stemm, logger=TRUE,
#                     optimizer="nlminb", control=list(iter.max=1))
# })
# 
