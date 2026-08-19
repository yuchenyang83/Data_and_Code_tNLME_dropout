################################################################################
#
#   Filename    :    fit_actg398_sensitivity.R
#   Project     :    Biometrics article "Modeling HIV Viral Dynamics Using a
#                    Nonlinear Mixed-Effects Framework for Heavy-Tailed Data
#                    with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#   Date        :    18.08.2026
#   Purpose     :    perform the Appendix F SAEM tuning sensitivity analysis
#                    for the Scenario (III) AR(1) tNLME model under MCAR, MAR,
#                    and MNAR by varying K = 5, 10, 15 and c = 1.8, 2.4, 3.0
#
#   Input data files  :  Data_and_Code/Data/source/actg398.txt
#
#   Output data files :  Data_and_Code/Data/sensitivity_SAEM/
#                         fit.t.III.ARp.MCAR_k*_c*.RData
#                         fit.t.III.ARp.MAR_k*_c*.RData
#                         fit.t.III.ARp.MNAR_k*_c*.RData
#
#
#   R Version   :    R-4.6.0
#
################################################################################
rm(list = ls())
set.seed(20260518)
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
actg398 <- read.table(paste(PATH, "/Data/source/actg398.txt", sep = ""), header = T)
table(actg398$calwk)

setdiff(1:481, actg398$patid[actg398$calwk == 0])
setdiff(1:481, actg398$patid[actg398$txday == 0])
actg398$patid[actg398$calwk == 0]

actg398[actg398$patid == 56, ]
actg398[actg398$patid == 229, ]

# actg398[actg398$patid == 56, ]$calwk <- c(0, 2, 8)
# actg398[actg398$patid == 56, ]$txday <- actg398[actg398$patid == 56, ]$txday - 13
# actg398[actg398$patid == 229, ]$calwk <- c(0, 2, 8, 16, 24)
# actg398[actg398$patid == 229, ]$txday <- actg398[actg398$patid == 229, ]$txday - 15


# actg398[actg398$patid == setdiff(1:481, actg398$patid[actg398$txday == 0]), ]

setdiff(1:481, actg398$patid[actg398$txday == 0])
actg398[actg398$patid == 1, ]

table(actg398$trtarm[actg398$calwk == 0])
table(actg398$trtarm[actg398$calwk == 0], actg398$nnrti[actg398$calwk == 0])
actg398$trtarm[actg398$trtarm == 1] <- 0
actg398$trtarm[actg398$trtarm == 2] <- 0
actg398$trtarm[actg398$trtarm == 3] <- 0
actg398$trtarm[actg398$trtarm == 4] <- 1
table(actg398$trtarm[actg398$calwk == 0])

# iid <- which(actg398$txday[actg398$calwk == 0] != 0)
# for (i in 1:length(iid))
# {
#   actg398$txday[actg398$patid == iid[i]] <- actg398$txday[actg398$patid == iid[i]] - min(actg398$txday[actg398$patid == iid[i]])
# }
table(actg398$calwk)
table(actg398$txday)


ll <- c(0, 2, 4, 8, 16, 24)
actg398$Subject <- actg398$patid
actg398$Time <- actg398$calwk

Subject <- unique(actg398$Subject)
N <- length(Subject)
ni <- numeric(N)
for (i in 1:N) ni[i] <- length(actg398$Time[actg398$Subject == Subject[i]])

actg398
actg398$D <- NA
for (i in 1:N)
{
  actg398$D[actg398$Subject == Subject[i]] <- max(actg398$Time[actg398$Subject == Subject[i]])
}

actg398[actg398$D != 24, ]
actg398[actg398$patid == 19, ]
actg398[actg398$patid == 56, ]


D.max <- 24
actg398.miss <- NULL
ll <- c(0, 2, 4, 8, 16, 24)
for (i in 1:N)
{
  if (unique(actg398[actg398$Subject == i, ]$D) < D.max) {
    actg398.i <- actg398[actg398$Subject == i, ]
    kk <- which(ll == actg398.i$D[1])
    cc <- ll[(kk + 1):length(ll)]
    yy <- which(actg398.i$Time == max(actg398.i$D))
    
    for (ii in 1:length(cc)) actg398.i <- rbind(actg398.i, actg398.i[dim(actg398.i)[1], ])
    actg398.i$Time[(yy + 1):length(actg398.i$Time)] <- ll[(kk + 1):length(ll)]
    actg398.i$logrna[(yy + 1):length(actg398.i$Time)] <- NA
    actg398.i$D <- actg398.i$Time[yy + 1]
    actg398.i <- cbind(actg398.i, miss = 1)
    actg398.miss <- rbind(actg398.miss, actg398.i)
  } else {
    actg398.i <- actg398[actg398$Subject == i, ]
    actg398.i <- cbind(actg398.i, miss = 0)
    actg398.miss <- rbind(actg398.miss, actg398.i)
  }
}
actg398.miss[actg398.miss$patid == 19, ]
actg398.miss[actg398.miss$patid == 169, ]
actg398.miss[actg398.miss$patid == 175, ]

actg398.miss[actg398.miss$patid == 56, ]
actg398.miss[actg398.miss$patid == 229, ]

actg398.miss[which(is.na(actg398.miss$logrna)), ]
actg398.miss[which(is.na(actg398.miss$logrna)), ]$txday <- actg398.miss[which(is.na(actg398.miss$logrna)), ]$Time * 7

## change to matrix
n <- max(ni)
N <- length(unique(actg398$Subject))
Subject <- unique(actg398.miss$Subject)
n <- length(Subject)
nj <- numeric(n)
for (i in 1:n) nj[i] <- length(actg398.miss$Dayt[actg398$Subject == Subject[i]])

actg398.miss$patid[actg398.miss$miss == 1]
actg398.miss[actg398.miss$miss == 1, ]

########  fit model
library(nlme)
actg398.miss <- data.frame(
  actg398.miss$Subject, actg398.miss$txday / 7, actg398.miss$logrna,
  actg398.miss$D, actg398.miss$miss, actg398.miss$trtarm, actg398.miss$nnrti
)
colnames(actg398.miss) <- c("Subject", "Time", "Var1", "D", "miss", "trtarm", "nnrti")
actg398.miss <- groupedData(Var1 ~ Time | Subject, data = actg398.miss)
actg398.miss$R <- 0
actg398.miss$R[which(is.na(actg398.miss$Var1))] <- 1
Data <- actg398.miss
Data$Subject <- as.numeric(as.character(Data$Subject))

Data[Data$Subject==56,]

#################### ####################
#################### ####################
## ===== initial values =====
library(nlme)

Data.obs <- Data[which(!is.na(Data$Var1)), ]

Data.obs$Subject <- factor(Data.obs$Subject)
Data.obs$Time    <- as.numeric(Data.obs$Time)
Data.obs$Var1    <- as.numeric(Data.obs$Var1)
Data.obs$nnrti   <- as.numeric(Data.obs$nnrti)
Data.obs$trtarm  <- as.numeric(Data.obs$trtarm)

Data.obs$Yexp <- 10^Data.obs$Var1

n.obs <- table(Data.obs$Subject)
keep.id <- names(n.obs[n.obs >= 4])
Data.obs4 <- Data.obs[Data.obs$Subject %in% keep.id, ]

fm1.list <- nlsList(
  Yexp ~ SSbiexp(Time, A1s, lrc1, A2s, lrc2) | Subject,
  data = Data.obs4,
  control = nls.control(
    maxiter = 500,
    warnOnly = TRUE,
    minFactor = 1 / 10000
  ),
  na.action = na.omit,
  pool = FALSE
)

coef.ss <- as.data.frame(coef(fm1.list))
coef.ss$Subject <- rownames(coef.ss)

head(coef.ss)
summary(coef.ss)

coef.ss <- coef.ss[complete.cases(coef.ss), ]

coef.ss <- coef.ss[
  is.finite(coef.ss$A1s) &
    is.finite(coef.ss$A2s) &
    is.finite(coef.ss$lrc1) &
    is.finite(coef.ss$lrc2) &
    coef.ss$A1s > 0 &
    coef.ss$A2s > 0,
]

coef.ss$logA1 <- log(coef.ss$A1s)
coef.ss$logA2 <- log(coef.ss$A2s)
coef.ss$L1    <- exp(coef.ss$lrc1)
coef.ss$L2    <- exp(coef.ss$lrc2)

idx.swap <- which(coef.ss$L1 < coef.ss$L2)

if (length(idx.swap) > 0) {
  tmp.logA <- coef.ss$logA1[idx.swap]
  tmp.L    <- coef.ss$L1[idx.swap]
  
  coef.ss$logA1[idx.swap] <- coef.ss$logA2[idx.swap]
  coef.ss$L1[idx.swap]    <- coef.ss$L2[idx.swap]
  
  coef.ss$logA2[idx.swap] <- tmp.logA
  coef.ss$L2[idx.swap]    <- tmp.L
}

Data.obs.tmp <- Data.obs[order(Data.obs$Subject, Data.obs$Time), ]

subj.cov <- Data.obs.tmp[!duplicated(Data.obs.tmp$Subject), c("Subject", "nnrti", "trtarm")]
subj.cov$Subject <- as.character(subj.cov$Subject)
rownames(subj.cov) <- NULL

coef.ss$Subject <- as.character(coef.ss$Subject)

coef.ss <- merge(coef.ss, subj.cov, by = "Subject")
coef.ss <- coef.ss[complete.cases(coef.ss), ]

## ===== regress subject-level rates on covariates =====
lm.L1 <- lm(L1 ~ nnrti + trtarm, data = coef.ss)
lm.L2 <- lm(L2 ~ nnrti + trtarm, data = coef.ss)

get.coef <- function(object, name, default = 0) {
  cc <- coef(object)
  if (name %in% names(cc)) {
    val <- unname(cc[name])
    ifelse(is.finite(val), val, default)
  } else {
    default
  }
}

start.selfStart <- c(
  phi1 = median(coef.ss$logA1, na.rm = TRUE),
  phi2 = get.coef(lm.L1, "(Intercept)"),
  phi3 = get.coef(lm.L1, "nnrti"),
  phi4 = get.coef(lm.L1, "trtarm"),
  phi5 = median(coef.ss$logA2, na.rm = TRUE),
  phi6 = get.coef(lm.L2, "(Intercept)"),
  phi7 = get.coef(lm.L2, "nnrti"),
  phi8 = get.coef(lm.L2, "trtarm")
)

round(start.selfStart, 3)


start.nlme <- start.selfStart

f1 <- Var1 ~ log10(
  exp(phi1 - phi2 * Time - phi3 * nnrti * Time - phi4 * trtarm * Time) +
    exp(phi5 - phi6 * Time - phi7 * nnrti * Time - phi8 * trtarm * Time)
)

fm1.nlme <- nlme(
  f1, fixed = phi1 + phi2 + phi3 + phi4 + phi5 + phi6 + phi7 + phi8 ~ 1,
  random = phi1 + phi5 ~ 1 | Subject,
  data = Data.obs, start = start.nlme, verbose = TRUE,
  control = nlmeControl(
    pnlsMaxIter = 50,
    msMaxIter = 200
  )
)

summary(fm1.nlme)
fixed.effects(fm1.nlme)


Beta <- fixed.effects(fm1.nlme)
q <- 2
DD <- matrix(NA, q, q)
DD[1, 1] <- as.numeric(VarCorr(fm1.nlme)[1, 1])
DD[2, 2] <- as.numeric(VarCorr(fm1.nlme)[2, 1])
DD[1, 2] <- DD[2, 1] <- sqrt(as.numeric(VarCorr(fm1.nlme)[1, 2]) * as.numeric(VarCorr(fm1.nlme)[2, 2])) * as.numeric(VarCorr(fm1.nlme)[2, 3])
sigma <- sigma(fm1.nlme)^2

## prediction of missing response
source(paste0(PATH, "/function/analyze_realdata_AIDS.R"))
alpha <- as.numeric(prediction_ym(Data, fm1.nlme, "UNC")$alpha.hat)
Data$yc <- prediction_ym(Data, fm1.nlme, "UNC")$yc
prediction_ym(Data, fm1.nlme, "UNC")$fm.result

############################################################
## Source functions
############################################################

source(paste0(PATH, "/function/analyze_realdata_AIDS.R"))
source(paste0(PATH, "/function/NLMMmissingSAEM.R"))
source(paste0(PATH, "/function/tNLMMmissingSAEM.R"))


############################################################
## Structure (I): random effects on phi1 and phi5
############################################################

fm1.I <- nlme(
  f1,
  fixed = phi1 + phi2 + phi3 + phi4 + phi5 + phi6 + phi7 + phi8 ~ 1,
  random = phi1 + phi5 ~ 1 | Subject,
  data = Data.obs,
  start = start.nlme,
  verbose = TRUE,
  control = nlmeControl(pnlsMaxIter = 50, msMaxIter = 200)
)

Beta.I <- fixed.effects(fm1.I)
vc.I <- VarCorr(fm1.I)
var.I <- suppressWarnings(as.numeric(vc.I[, "Variance"]))
var.I <- var.I[is.finite(var.I)]
DD.I <- diag(pmax(var.I[1:2], 1e-4), 2)
sigma.I <- sigma(fm1.I)^2
q <- 2
DD <- matrix(NA, q, q)
DD[1, 1] <- as.numeric(VarCorr(fm1.nlme)[1, 1])
DD[2, 2] <- as.numeric(VarCorr(fm1.nlme)[2, 1])
DD[1, 2] <- DD[2, 1] <- sqrt(as.numeric(VarCorr(fm1.nlme)[1, 2]) * as.numeric(VarCorr(fm1.nlme)[2, 2])) * as.numeric(VarCorr(fm1.nlme)[2, 3])


pred.I <- prediction_ym(Data, fm1.I, "UNC")
Data.I <- Data
Data.I$yc <- pred.I$yc

alpha.I.MNAR <- as.numeric(pred.I$alpha.hat)
alpha.I.MAR <- alpha.I.MNAR[-5]
alpha.I.MCAR <- alpha.I.MNAR[1]

init.N.I.MNAR <- list(Beta = Beta.I, DD = DD.I, sigma = sigma.I, Phi = 0.5, ga = 1, alpha = alpha.I.MNAR)
init.N.I.MAR <- list(Beta = Beta.I, DD = DD.I, sigma = sigma.I, Phi = 0.5, ga = 1, alpha = alpha.I.MAR)
init.N.I.MCAR <- list(Beta = Beta.I, DD = DD.I, sigma = sigma.I, Phi = 0.5, ga = 1, alpha = alpha.I.MCAR)

init.t.I.MNAR <- list(Beta = Beta.I, DD = DD, sigma = sigma.I, Phi = 0.5, ga = 1, nu = 50, alpha = alpha.I.MNAR)
init.t.I.MAR <- list(Beta = Beta.I, DD = DD, sigma = sigma.I, Phi = 0.5, ga = 1, nu = 50, alpha = alpha.I.MAR)
init.t.I.MCAR <- list(Beta = Beta.I, DD = DD, sigma = sigma.I, Phi = 0.5, ga = 1, nu = 50, alpha = alpha.I.MCAR)

init.t.I.MNAR.new <- list(Beta = Beta.I, DD = DD.I, sigma = sigma.I, Phi = 0.5, ga = 1, nu = 50, alpha = alpha.I.MNAR)
init.t.I.MAR.new <- list(Beta = Beta.I, DD = DD.I, sigma = sigma.I, Phi = 0.5, ga = 1, nu = 50, alpha = alpha.I.MAR)
init.t.I.MCAR.new <- list(Beta = Beta.I, DD = DD.I, sigma = sigma.I, Phi = 0.5, ga = 1, nu = 50, alpha = alpha.I.MCAR)


############################################################
## Structure (II): random effects on phi2 and phi6
############################################################

fm1.II <- nlme(
  f1,
  fixed = phi1 + phi2 + phi3 + phi4 + phi5 + phi6 + phi7 + phi8 ~ 1,
  random = list(Subject = pdDiag(phi2 + phi6 ~ 1)),
  data = Data.obs,
  start = fixed.effects(fm1.I),
  verbose = TRUE,
  control = nlmeControl(
    pnlsMaxIter = 30,
    msMaxIter = 100,
    msMaxEval = 200,
    pnlsTol = 1e-3,
    tolerance = 1e-4,
    returnObject = TRUE
  )
)

Beta.II <- fixed.effects(fm1.II)
vc.II <- VarCorr(fm1.II)
var.II <- suppressWarnings(as.numeric(vc.II[, "Variance"]))
var.II <- var.II[is.finite(var.II)]
DD.II <- diag(pmax(var.II[1:2], 1e-4), 2)
sigma.II <- sigma(fm1.II)^2

init.N.II.MNAR <- list(Beta = Beta.II, DD = DD.II, sigma = sigma.II, Phi = 0.5, ga = 1, alpha = alpha.I.MNAR)
init.N.II.MAR <- list(Beta = Beta.II, DD = DD.II, sigma = sigma.II, Phi = 0.5, ga = 1, alpha = alpha.I.MAR)
init.N.II.MCAR <- list(Beta = Beta.II, DD = DD.II, sigma = sigma.II, Phi = 0.5, ga = 1, alpha = alpha.I.MCAR)

init.t.II.MNAR <- list(Beta = Beta.II, DD = DD.II, sigma = sigma.II, Phi = 0.5, ga = 1, nu = 50, alpha = alpha.I.MNAR)
init.t.II.MAR <- list(Beta = Beta.II, DD = DD.II, sigma = sigma.II, Phi = 0.5, ga = 1, nu = 50, alpha = alpha.I.MAR)
init.t.II.MCAR <- list(Beta = Beta.II, DD = DD.II, sigma = sigma.II, Phi = 0.5, ga = 1, nu = 50, alpha = alpha.I.MCAR)

init.t.II.MNAR.new <- list(Beta = Beta.II, DD = diag(1,2), sigma = sigma.II, Phi = 0.5, ga = 1, nu = 50, alpha = alpha.I.MNAR)
init.t.II.MAR.new <- list(Beta = Beta.II, DD = diag(1,2), sigma = sigma.II, Phi = 0.5, ga = 1, nu = 50, alpha = alpha.I.MAR)
init.t.II.MCAR.new <- list(Beta = Beta.II, DD = diag(1,2), sigma = sigma.II, Phi = 0.5, ga = 1, nu = 50, alpha = alpha.I.MCAR)

############################################################
## Structure (III): random effects on phi1, phi2, phi5, and phi6
############################################################

fm1.III <- nlme(
  f1,
  fixed = phi1 + phi2 + phi3 + phi4 + phi5 + phi6 + phi7 + phi8 ~ 1,
  random = list(Subject = pdDiag(phi1 + phi2 + phi5 + phi6 ~ 1)),
  data = Data.obs,
  start = fixed.effects(fm1.II),
  verbose = TRUE,
  control = nlmeControl(
    pnlsMaxIter = 30,
    msMaxIter = 150,
    msMaxEval = 300,
    pnlsTol = 1e-3,
    tolerance = 1e-4,
    returnObject = TRUE
  )
)


Beta.III <- fixed.effects(fm1.III)
vc.III <- VarCorr(fm1.III)
var.III <- suppressWarnings(as.numeric(vc.III[, "Variance"]))
var.III <- var.III[is.finite(var.III)]
DD.III <- diag(pmax(var.III[1:4], 1e-4), 4)
sigma.III <- sigma(fm1.III)^2

init.N.III.MNAR <- list(Beta = Beta.III, DD = DD.III, sigma = sigma.III, Phi = 0.5, ga = 1, alpha = alpha.I.MNAR)
init.N.III.MAR <- list(Beta = Beta.III, DD = DD.III, sigma = sigma.III, Phi = 0.5, ga = 1, alpha = alpha.I.MAR)
init.N.III.MCAR <- list(Beta = Beta.III, DD = DD.III, sigma = sigma.III, Phi = 0.5, ga = 1, alpha = alpha.I.MCAR)

init.t.III.MNAR <- list(Beta = Beta.III, DD = DD.III, sigma = sigma.III, Phi = 0.5, ga = 1, nu = 50, alpha = alpha.I.MNAR)
init.t.III.MAR <- list(Beta = Beta.III, DD = DD.III, sigma = sigma.III, Phi = 0.5, ga = 1, nu = 50, alpha = alpha.I.MAR)
init.t.III.MCAR <- list(Beta = Beta.III, DD = DD.III, sigma = sigma.III, Phi = 0.5, ga = 1, nu = 50, alpha = alpha.I.MCAR)


############################################################
## ===== Alternative initial values based on pooled nls =====
############################################################

## Clean data for pooled nls
Data.nls <- Data.obs[
  is.finite(Data.obs$Time) &
    is.finite(Data.obs$Var1) &
    is.finite(Data.obs$Yexp) &
    Data.obs$Yexp > 0,
]

Data.nls$Time <- as.numeric(Data.nls$Time)
Data.nls$Var1 <- as.numeric(Data.nls$Var1)
Data.nls$Yexp <- as.numeric(Data.nls$Yexp)

## Collapse repeated observations by time to obtain a smoother pooled trajectory.
## Median is more robust than mean for heavy-tailed viral-load data.
Data.pool <- aggregate(
  Yexp ~ Time,
  data = Data.nls,
  FUN = median,
  na.rm = TRUE
)

Data.pool <- Data.pool[
  is.finite(Data.pool$Time) &
    is.finite(Data.pool$Yexp) &
    Data.pool$Yexp > 0,
]

Data.pool <- Data.pool[order(Data.pool$Time), ]

## Data-driven crude initial values
y0.pool <- max(Data.pool$Yexp, na.rm = TRUE)
yT.pool <- min(Data.pool$Yexp, na.rm = TRUE)

start.pool <- c(
  logA1 = log(0.7 * y0.pool),
  L1    = 0.5,
  logA2 = log(0.3 * y0.pool),
  L2    = 0.05
)

## Fit pooled biexponential model on original RNA scale
fm0.nls <- nls(
  Yexp ~ exp(logA1 - L1 * Time) + exp(logA2 - L2 * Time),
  data = Data.pool,
  start = start.pool,
  algorithm = "port",
  lower = c(
    logA1 = -20,
    L1    = 1e-6,
    logA2 = -20,
    L2    = 1e-6
  ),
  upper = c(
    logA1 = 30,
    L1    = 10,
    logA2 = 30,
    L2    = 10
  ),
  control = nls.control(
    maxiter = 1000,
    warnOnly = TRUE,
    minFactor = 1 / 10000
  )
)

coef.nls0 <- coef(fm0.nls)

logA1.nls <- coef.nls0["logA1"]
logA2.nls <- coef.nls0["logA2"]
L1.nls    <- coef.nls0["L1"]
L2.nls    <- coef.nls0["L2"]

## Order the two exponential phases so that phase 1 has the faster decay rate
if (L1.nls < L2.nls) {
  tmp.logA <- logA1.nls
  tmp.L    <- L1.nls
  
  logA1.nls <- logA2.nls
  L1.nls    <- L2.nls
  
  logA2.nls <- tmp.logA
  L2.nls    <- tmp.L
}

start.nls0 <- c(
  phi1 = as.numeric(logA1.nls),
  phi2 = as.numeric(L1.nls),
  phi3 = 0,
  phi4 = 0,
  phi5 = as.numeric(logA2.nls),
  phi6 = as.numeric(L2.nls),
  phi7 = 0,
  phi8 = 0
)

round(start.nls0, 3)

############################################################
## ===== Full pooled nls for f1 =====
############################################################

fm1.nls <- nls(
  f1,
  data = Data.nls,
  start = start.nls0,
  algorithm = "port",
  lower = c(
    phi1 = -20, phi2 = 1e-6, phi3 = -10, phi4 = -10,
    phi5 = -20, phi6 = 1e-6, phi7 = -10, phi8 = -10
  ),
  upper = c(
    phi1 = 30, phi2 = 10, phi3 = 10, phi4 = 10,
    phi5 = 30, phi6 = 10, phi7 = 10, phi8 = 10
  ),
  control = nls.control(
    maxiter = 1000,
    warnOnly = TRUE,
    minFactor = 1 / 10000
  ),
  na.action = na.omit
)

start.nls <- coef(fm1.nls)

## Again order the two phases if necessary
if (start.nls["phi2"] < start.nls["phi6"]) {
  tmp <- start.nls[c("phi1", "phi2", "phi3", "phi4")]
  start.nls[c("phi1", "phi2", "phi3", "phi4")] <- start.nls[c("phi5", "phi6", "phi7", "phi8")]
  start.nls[c("phi5", "phi6", "phi7", "phi8")] <- tmp
}

round(start.nls, 3)

init.N.III.MNAR.nls <- list(Beta = start.nls, DD = diag(1, 4), sigma = 1, Phi = 0.5, ga = 1, alpha = alpha.I.MNAR)
init.N.III.MAR.nls <- list(Beta = start.nls, DD = diag(1, 4), sigma = 1, Phi = 0.5, ga = 1, alpha = alpha.I.MAR)
init.N.III.MCAR.nls <- list(Beta = start.nls, DD = diag(1, 4), sigma = 1, Phi = 0.5, ga = 1, alpha = alpha.I.MCAR)

init.t.III.MNAR.nls <- list(Beta = start.nls, DD = diag(1, 4), sigma = 1, Phi = 0.5, ga = 1, nu = 50, alpha = alpha.I.MNAR)
init.t.III.MAR.nls <- list(Beta = start.nls, DD = diag(1, 4), sigma = 1, Phi = 0.5, ga = 1, nu = 50, alpha = alpha.I.MAR)
init.t.III.MCAR.nls <- list(Beta = start.nls, DD = diag(1, 4), sigma = 1, Phi = 0.5, ga = 1, nu = 50, alpha = alpha.I.MCAR)



M.LL <- 1000
tol <- 1e-5
max.iter <- 1000
per <- 500
M <- 10
SAVE_PATH <- file.path(PATH, "Data", "sensitivity_SAEM")

############################################################
## Source functions sensitivity
############################################################

source(paste0(PATH, "/function/analyze_realdata_AIDS_sensitivity.R"))
source(paste0(PATH, "/function/NLMMmissingSAEM_sensitivity.r"))
source(paste0(PATH, "/function/tNLMMmissingSAEM_sensitivity.r"))

K.list <- c(5, 10, 15)
c.list <- c(1.8, 2.4, 3.0)

for (K.label in K.list) {
  for (c.label in c.list) {
    
    K.run <- K.label
    
    if (K.label == 5 && c.label == 2.4) K.run <- 10
    if (K.label == 10 && c.label == 2.4) K.run <- 5
    
    seed.run <- 20260527
    
    if (K.label == 10 && c.label == 2.4) seed.run <- 20250502
    
    k.step <- K.run
    ccc <- c.label
    
    kc.tag <- paste0("_k", K.label, "_c", c.label)
    
    set.seed(seed.run)
    
    fit.t.III.ARp.MNAR <- tNLMM.miss.SAEM(
      Data,
      g = 1,
      init.para = init.t.III.MNAR,
      cor.type = "ARp",
      k.step = k.step,
      ccc = ccc,
      M = 10,
      M.LL = M.LL,
      P = 1,
      tol = tol,
      max.iter = max.iter,
      per = per,
      mechanism = "MNAR",
      random.structure = "III"
    )
    
    save(
      fit.t.III.ARp.MNAR,
      file = file.path(
        SAVE_PATH,
        paste0("fit.t.III.ARp.MNAR", kc.tag, ".RData")
      )
    )
    
    rm(fit.t.III.ARp.MNAR)
    
    set.seed(seed.run)
    
    fit.t.III.ARp.MAR <- tNLMM.miss.SAEM(
      Data,
      g = 1,
      init.para = init.t.III.MAR,
      cor.type = "ARp",
      k.step = k.step,
      ccc = ccc,
      M = 10,
      M.LL = M.LL,
      P = 1,
      tol = tol,
      max.iter = max.iter,
      per = per,
      mechanism = "MAR",
      random.structure = "III"
    )
    
    save(
      fit.t.III.ARp.MAR,
      file = file.path(
        SAVE_PATH,
        paste0("fit.t.III.ARp.MAR", kc.tag, ".RData")
      )
    )
    
    rm(fit.t.III.ARp.MAR)
    
    set.seed(seed.run)
    
    fit.t.III.ARp.MCAR <- tNLMM.miss.SAEM(
      Data,
      g = 1,
      init.para = init.t.III.MCAR,
      cor.type = "ARp",
      k.step = k.step,
      ccc = ccc,
      M = 10,
      M.LL = M.LL,
      P = 1,
      tol = tol,
      max.iter = max.iter,
      per = per,
      mechanism = "MCAR",
      random.structure = "III"
    )
    
    save(
      fit.t.III.ARp.MCAR,
      file = file.path(
        SAVE_PATH,
        paste0("fit.t.III.ARp.MCAR", kc.tag, ".RData")
      )
    )
    
    rm(fit.t.III.ARp.MCAR)
  }
}
