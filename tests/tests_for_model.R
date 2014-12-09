# =====================
# Tests for specify_sem
# =====================

library(nlsem)

# ordinary lms
# ============
num.classes <- 1
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.classes,
                         interaction="xi1:xi2",
                         interc.obs=FALSE, interc.lat=FALSE)
as.data.frame(lms_model)

# ordinary lms with different inputs for interaction
# =================================================
# "all"
# -----
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.classes,
                         interaction="all",
                         interc.obs=FALSE, interc.lat=FALSE)

# ""
# --
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.classes,
                         interaction="",
                         interc.obs=FALSE, interc.lat=FALSE)

# not all interactions with xi defined
# ------------------------------------
lms_model <- specify_sem(num.x=8, num.y=6, num.xi=4, num.eta=3,
                         xi="x1-x2,x3-x4,x5-x6,x7-x8", eta="y1-y2,y3-y4,y5-y6",
                         num.classes=1, interaction="xi1:xi2,xi1:xi1",
                         interc.obs=FALSE, interc.lat=FALSE)

# stemm model
# ===========
stemm_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.classes=3,
                         interaction="",
                         interc.obs=FALSE, interc.lat=FALSE)

# ===============================
# Tests for count_free_parameters
# ===============================

# lms
# ---
lms_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.classes=1,
                         interaction="xi1:xi2",
                         interc.obs=FALSE, interc.lat=FALSE)
class(lms_model)
count_free_parameters(lms_model)
parameters <- 1:count_free_parameters(lms_model)
rm(parameters, lms_model)

# stemm
# -----
stemm_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.classes=3,
                         interaction="",
                         interc.obs=FALSE, interc.lat=FALSE)
class(stemm_model)
count_free_parameters(stemm_model)
parameters <- 1:count_free_parameters(stemm_model)
rm(parameters, stemm_model)

# nsemm
# -----
nsemm_model <- specify_sem(num.x=6, num.y=3, num.xi=2, num.eta=1,
                         xi="x1-x3,x4-x6", eta="y1-y3", num.classes=3,
                         interaction="xi1:xi2",
                         interc.obs=FALSE, interc.lat=FALSE)
class(nsemm_model)
count_free_parameters(nsemm_model)
parameters <- 1:count_free_parameters(nsemm_model)
rm(parameters, nsemm_model)

# nsemm with relation.lat
# =======================
nsemm_model <- specify_sem(num.x=8, num.y=6, num.xi=4, num.eta=3,
                         xi="x1-x2,x3-x4,x5-x6,x7-x8", eta="y1-y2,y3-y4,y5-y6",
                         num.classes=3, interaction="", interc.obs=FALSE,
                         interc.lat=FALSE,
                         relation.lat="xi1>eta1;xi2>eta2;xi3,xi4>eta3;eta1,eta3>eta2")
