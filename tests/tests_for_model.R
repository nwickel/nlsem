# =====================
# Tests for specify_sem
# =====================

# ordinary lms
# ============
num.groups <- 1
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.groups,
                         interaction="xi1:xi2",
                         interc_obs=FALSE, interc_lat=FALSE)
# OK runs without problem
as.data.frame(lms_model)
# OK runs without problem

# ordinary lms with wrong input for observed variables
# ====================================================
# too many x's in xi / too few num.x
# ----------------------------------
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x4,x5-x8", eta="y1-y3", num.groups,
                         interaction="xi1:xi2",
                         interc_obs=FALSE, interc_lat=FALSE)
# TODO subscript out of bounds

# too few x's in xi / too many num.x
# ----------------------------------
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x2,x3-x4", eta="y1-y3", num.groups,
                         interaction="xi1:xi2",
                         interc_obs=FALSE, interc_lat=FALSE)
# TODO this runs, but should throw an error (?)

# too many y's in eta / too few num.y
# -----------------------------------
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y5", num.groups,
                         interaction="xi1:xi2",
                         interc_obs=FALSE, interc_lat=FALSE)
# TODO subscript out of bounds

# too few y's in eta / too many num.y
# -----------------------------------
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y2", num.groups,
                         interaction="xi1:xi2",
                         interc_obs=FALSE, interc_lat=FALSE)
# TODO this runs, but should throw an error (?)

# ordinary lms with wrong input for latent variables
# ==================================================
# too few num.xi
# --------------
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=1, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.groups,
                         interaction="xi1:xi2",
                         interc_obs=FALSE, interc_lat=FALSE)
# OK error: Interaction effects contain more xi's than defined.

# nonsense input for xi
# ---------------------
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x-x4", eta="y1-y3", num.groups,
                         interaction="xi1:xi2",
                         interc_obs=FALSE, interc_lat=FALSE)
# OK Number of xi's and assignation of x's to xi's does not match.
# See ?specify_sem.
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1x2-x3x4", eta="y1-y3", num.groups,
                         interaction="xi1:xi2",
                         interc_obs=FALSE, interc_lat=FALSE)
# OK Number of xi's and assignation of x's to xi's does not match. See
# ?specify_sem.

# too few num.eta / too many eta's
# --------------------------------
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1, xi="x1-x3,x4-x6",
            eta="y1-y3,y4-y6", num.groups, interaction="xi1:xi2",
            interc_obs=FALSE, interc_lat=FALSE)
# OK Number of eta's and assignation of y's to eta's does not match. See
# ?specify_sem.

# nonsense input for eta
# ----------------------
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y-y3", num.groups,
                         interaction="xi1:xi2",
                         interc_obs=FALSE, interc_lat=FALSE)
# OK Wrong input for specifying exogenous or endogonous latent variables (xi or
# etas). See ?specify_sem.
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1y3-y4", num.groups,
                         interaction="xi1:xi2",
                         interc_obs=FALSE, interc_lat=FALSE)
# TODO subscript out of bounds

# ordinary lms with different inputs for interaction
# =================================================
# "all"
# -----
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.groups,
                         interaction="all",
                         interc_obs=FALSE, interc_lat=FALSE)
# OK this runs

# ""
# --
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.groups,
                         interaction="",
                         interc_obs=FALSE, interc_lat=FALSE)
# OK Model needs either more than one latent group or at least one latent
# interaction (e.g. 'xi1:xi2'). For other models please use lavaan or the
# like

# more interaction terms than xi's
# --------------------------------
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.groups,
                         interaction="xi2:xi3",
                         interc_obs=FALSE, interc_lat=FALSE)
# OK Interaction effects contain more xi's than defined.

# nonsense
# --------
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.groups,
                         interaction="x3:x4",
                         interc_obs=FALSE, interc_lat=FALSE)
# OK Wrong input for interaction. See ?specify_sem.

# ordinary lms with nonsense input for interc_obs, interc_obs
# ===========================================================
# (should be boolean)
# ===================
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1, xi="x1-x3,x4-x6",
            eta="y1-y3", num.groups, interaction="xi1:xi2",
            interc_obs="bla", interc_lat=FALSE)
# Error in if (interc_obs) { : argument is not interpretable as logical

lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.groups,
                         interaction="xi1:xi2",
                         interc_obs=TRUE, interc_lat="bla")
# Error in if (interc_lat) { : argument is not interpretable as logical
# These error messages are quite informative, so just leave them like that

# stemm model
# ===========
stemm_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.groups=3,
                         interaction="xi1:xi2",
                         interc_obs=FALSE, interc_lat=FALSE)
# OK

# nonsense input for num.groups
# -----------------------------
stemm_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.groups="three",
                         interaction="xi1:xi2",
                         interc_obs=FALSE, interc_lat=FALSE)
# OK Number of variables or groups must be numeric.


# ==============================================
# Tests for count_free_parameters and fill_model
# ==============================================

# lms
# ---
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.groups=1,
                         interaction="xi1:xi2",
                         interc_obs=FALSE, interc_lat=FALSE)
class(lms_model)
count_free_parameters(lms_model)
parameters <- 1:count_free_parameters(lms_model)
fill_model(lms_model, parameters)
rm(parameters, lms_model)
# OK

# stemm
# -----
stemm_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.groups=3,
                         interaction="",
                         interc_obs=FALSE, interc_lat=FALSE)
class(stemm_model)
count_free_parameters(stemm_model)
parameters <- 1:count_free_parameters(stemm_model)
fill_model(stemm_model, parameters)
rm(parameters, stemm_model)
# OK

# nsemm
# -----
nsemm_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.groups=3,
                         interaction="xi1:xi2",
                         interc_obs=FALSE, interc_lat=FALSE)
class(nsemm_model)
count_free_parameters(nsemm_model)
parameters <- 1:count_free_parameters(nsemm_model)
fill_model(nsemm_model, parameters)
rm(parameters, nsemm_model)
# OK


# General specification
# =====================

# Parameters for specify_sem (nonlin)
# -----------------------------------
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
            interc_obs, interc_lat)
