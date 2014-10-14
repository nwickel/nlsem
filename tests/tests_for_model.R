# ordinary lms
# ============
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.groups=1,
                         interaction="xi1:xi2", constraints="default",
                         interc_obs=FALSE, interc_lat=FALSE)
# runs without problem

# ordinary lms with wrong input for observed variables
# ====================================================
# too many x's in xi / too few num.x
# ----------------------------------
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x4,x5-x8", eta="y1-y3", num.groups=1,
                         interaction="xi1:xi2", constraints="default",
                         interc_obs=FALSE, interc_lat=FALSE)
# TODO subscript out of bounds

# too few x's in xi / too many num.x
# ----------------------------------
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x2,x3-x4", eta="y1-y3", num.groups=1,
                         interaction="xi1:xi2", constraints="default",
                         interc_obs=FALSE, interc_lat=FALSE)
# TODO this runs, but should throw an error (?)

# too many y's in eta / too few num.y
# -----------------------------------
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y5", num.groups=1,
                         interaction="xi1:xi2", constraints="default",
                         interc_obs=FALSE, interc_lat=FALSE)
# TODO subscript out of bounds

# too few y's in eta / too many num.y
# -----------------------------------
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y2", num.groups=1,
                         interaction="xi1:xi2", constraints="default",
                         interc_obs=FALSE, interc_lat=FALSE)
# TODO this runs, but should throw an error (?)

# ordinary lms with wrong input for latent variables
# ==================================================
# too few num.xi
# --------------
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=1, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.groups=1,
                         interaction="xi1:xi2", constraints="default",
                         interc_obs=FALSE, interc_lat=FALSE)
# error: Interaction effects contain more xi's than defined.

# nonsense input for xi
# ---------------------
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x-x4", eta="y1-y3", num.groups=1,
                         interaction="xi1:xi2", constraints="default",
                         interc_obs=FALSE, interc_lat=FALSE)
# Wrong input for specifying exogenous or endogonous latent variables (xi or
# etas). See ?specify_sem.
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1x2-x3x4", eta="y1-y3", num.groups=1,
                         interaction="xi1:xi2", constraints="default",
                         interc_obs=FALSE, interc_lat=FALSE)
# TODO this runs, but should throw an error

# too few num.eta
# ---------------
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1, xi="x1-x3,x4-x6",
            eta="y1-y3,y4-y6", num.groups=1, interaction="xi1:xi2",
            constraints="default", interc_obs=FALSE, interc_lat=FALSE)
# TODO this should throw an error

# nonsense input for eta
# ----------------------
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y-y3", num.groups=1,
                         interaction="xi1:xi2", constraints="default",
                         interc_obs=FALSE, interc_lat=FALSE)
# Wrong input for specifying exogenous or endogonous latent variables (xi or
# etas). See ?specify_sem.
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1y3-y4", num.groups=1,
                         interaction="xi1:xi2", constraints="default",
                         interc_obs=FALSE, interc_lat=FALSE)
# TODO subscript out of bounds

# ordinary lms with different input for interaction
# =================================================
# "all"
# -----
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.groups=1,
                         interaction="all", constraints="default",
                         interc_obs=FALSE, interc_lat=FALSE)
# TODO this runs, but is it really a lms?

# ""
# --
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.groups=1,
                         interaction="", constraints="default",
                         interc_obs=FALSE, interc_lat=FALSE)
#   Model needs either more than one latent group or at least one latent
#   interaction (e.g. 'xi1:xi2'). For other models please use lavaan or the
#   like

# more interaction terms than xi's
# --------------------------------
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.groups=1,
                         interaction="xi2:xi3", constraints="default",
                         interc_obs=FALSE, interc_lat=FALSE)
# Interaction effects contain more xi's than defined.

# nonsense
# --------
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.groups=1,
                         interaction="x3:x4", constraints="default",
                         interc_obs=FALSE, interc_lat=FALSE)
# Wrong input for interaction. See ?specify_sem.








# Parameters for specify_sem (nonlin)
num.x <- 6
num.xi <- 2
num.eta <- 2
num.y <- 4
num.groups <- 3
xi <- "x1-x3,x4-x6"
eta <- "y1-y2,y3-y4"
interc_obs <- TRUE
interc_lat <- TRUE
# interaction <- "xi1:xi1, xi2:xi1"
interaction <- "all"

specify_sem(num.x, num.y, num.xi, num.eta, xi, eta, num.groups, interaction,
            constraints="default", interc_obs, interc_lat)
