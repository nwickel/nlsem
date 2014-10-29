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
model.filled <- fill_model(model, parameters)

data <- simulate(model.filled)

estep_stemm(model, parameters, data)
