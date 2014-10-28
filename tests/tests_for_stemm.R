# Testing functions for stemm
# ===========================

# create model
# ------------
mod <- specify_sem(num.x=4, num.y=4, num.xi=2, num.eta=2,
                     xi="x1-x2,x3-x4", eta="y1-y2,y3-y4", num.groups=3,
                     interaction="",
                     interc_obs=FALSE, interc_lat=FALSE)

dat <- as.data.frame(mod)

# Gamma -> diagonal
dat[18:19, 2:4] <- 0

model <- fill_matrices(dat, mod)

parameters <- runif(count_free_parameters(model), 0.1, 1.5)
model.filled <- fill_model(model, parameters)

data <- simulate(model.filled, parameters)
