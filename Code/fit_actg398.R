#################################################################################
#
#   Project     :    "Robust HIV Viral Dynamics: A Nonlinear Mixed-Effects Framework for
#                    Heavy-Tailed Data with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#
# Read the ACTG 398 dataset from 'Data_and_Code/Data/source/actg398.txt'.
# You may replace this with your own data if desired.
#
# This R script allows you to perform model fitting for both tNLME and NLME models
# under various missing data mechanisms: MNAR, MAR, and MCAR.
#
# The models can be fitted with one of the following within-subject error structures: UNC, CS, MA(1), or AR(1).
#
# The fitted model results will be saved to 'Data_and_Code/Data/fit_result.RData'.
#
#################################################################################
# rm(list = ls())
PATH <- getwd()
actg398 <- read.table(paste(PATH, "/Data/source/actg398.txt", sep = ""), header = T)
table(actg398$calwk)

setdiff(1:481, actg398$patid[actg398$calwk == 0])
setdiff(1:481, actg398$patid[actg398$txday == 0])
actg398$patid[actg398$calwk == 0]

actg398[actg398$patid == 56, ]
actg398[actg398$patid == 229, ]

actg398[actg398$patid == 56, ]$calwk <- c(0, 2, 8)
actg398[actg398$patid == 56, ]$txday <- actg398[actg398$patid == 56, ]$txday - 13
actg398[actg398$patid == 229, ]$calwk <- c(0, 2, 8, 16, 24)
actg398[actg398$patid == 229, ]$txday <- actg398[actg398$patid == 229, ]$txday - 15


# actg398[actg398$patid == setdiff(1:481, actg398$patid[actg398$txday == 0]), ]

setdiff(1:481, actg398$patid[actg398$txday == 0])
actg398[actg398$patid == 123, ]

table(actg398$trtarm[actg398$calwk == 0])
table(actg398$trtarm[actg398$calwk == 0], actg398$nnrti[actg398$calwk == 0])
actg398$trtarm[actg398$trtarm == 1] <- 0
actg398$trtarm[actg398$trtarm == 2] <- 0
actg398$trtarm[actg398$trtarm == 3] <- 0
actg398$trtarm[actg398$trtarm == 4] <- 1
table(actg398$trtarm[actg398$calwk == 0])

iid <- which(actg398$txday[actg398$calwk == 0] != 0)
for (i in 1:length(iid))
{
  actg398$txday[actg398$patid == iid[i]] <- actg398$txday[actg398$patid == iid[i]] - min(actg398$txday[actg398$patid == iid[i]])
}
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

Data[Data$miss == 1, ]

f1 <- Var1 ~ log10(exp(phi1 - phi2 * Time - phi3 * nnrti * Time - phi4 * trtarm * Time) + exp(phi5 - phi6 * Time - phi7 * nnrti * Time - phi8 * trtarm * Time))
fm1.list <- nlsList(f1, data = Data[which(!is.na(Data$Var1)), ], start = list(
  phi1 = 10, phi2 = 1.8, phi3 = 0.2, phi4 = 0.1,
  phi5 = 6, phi6 = -0.5, phi7 = -0.2, phi8 = -0.1
))
coef(fm1.list)
fm1.nlme <- nlme(f1,
  fixed = phi1 + phi2 + phi3 + phi4 + phi5 + phi6 + phi7 + phi8 ~ 1,
  random = phi1 + phi5 ~ 1 | Subject,
  data = Data[which(!is.na(Data$Var1)), ], start = c(
    phi1 = 10, phi2 = 1.8, phi3 = 0.2, phi4 = 0.1,
    phi5 = 6, phi6 = -0.5, phi7 = -0.2, phi8 = -0.1
  ), verbose = T
)
summary(fm1.nlme)

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

init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha)

M.LL <- 1000
tol <- 1e-6
max.iter <- 1000
per <- 500

M <- 10
cor.type <- "UNC"
mechanism <- "MNAR"

source(paste(PATH, "/function/tNLMMmissingSAEM.r", sep = ""))
source(paste(PATH, "/function/NLMMmissingSAEM.r", sep = ""))
###############
###############
DD1 <- diag(1, 2)
set.seed(202501002)
init.para <- list(Beta = Beta, DD = DD1, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha)

fit.MNAR <- NLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("UNC"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MNAR"
)

set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD1, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha[-1])

fit.MAR <- NLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("UNC"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MAR"
)

set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD1, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha[1])

fit.MCAR <- NLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("UNC"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MCAR"
)

source(paste(PATH, "tLMMmissingSAEM.r", sep = ""))
###############
###############
### UNC
set.seed(202501005)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha)

est.MNAR <- tNLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("UNC"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MNAR"
)

set.seed(20250202)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha[-1])
est.MAR <- tNLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("UNC"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MAR"
)
set.seed(20250202)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha[1])
est.MCAR <- tNLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("UNC"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MCAR"
)


###############
###############
### ARp
source(paste(PATH, "LMMmissingSAEM.r", sep = ""))
cor.type <- c("ARp")
set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha)

fit.MNAR.ARp <- NLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("ARp"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MNAR"
)

set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha[-1])

fit.MAR.ARp <- NLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("ARp"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MAR"
)

set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha[1])

fit.MCAR.ARp <- NLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("ARp"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MCAR"
)


source(paste(PATH, "tLMMmissingSAEM.r", sep = ""))
###############
###############
### ARp
cor.type <- c("ARp")
set.seed(20250206)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha)

est.MNAR.ARp <- tNLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("ARp"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MNAR"
)

set.seed(202501002)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha[-1])
est.MAR.ARp <- tNLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("ARp"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MAR"
)
set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha[1])
est.MCAR.ARp <- tNLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("ARp"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MCAR"
)

###############
###############
### CS
source(paste(PATH, "LMMmissingSAEM.r", sep = ""))
cor.type <- c("CS")
set.seed(20250203)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha)

fit.MNAR.CS <- NLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("CS"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MNAR"
)

set.seed(202502004)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha[-1])

fit.MAR.CS <- NLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("CS"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MAR"
)

set.seed(20250202)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha[1])

fit.MCAR.CS <- NLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("CS"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MCAR"
)


source(paste(PATH, "tLMMmissingSAEM.r", sep = ""))
###############
###############
### CS
cor.type <- c("CS")
set.seed(20250204)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha)

est.MNAR.CS <- tNLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("CS"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MNAR"
)

set.seed(20250202)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha[-1])
est.MAR.CS <- tNLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("CS"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MAR"
)
set.seed(20250202)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha[1])
est.MCAR.CS <- tNLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("CS"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MCAR"
)



###############
###############
### BAND1
source(paste(PATH, "LMMmissingSAEM.r", sep = ""))
cor.type <- c("BAND1")
set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha)

fit.MNAR.BAND1 <- NLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("BAND1"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MNAR"
)

set.seed(202501002)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha[-1])

fit.MAR.BAND1 <- NLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("BAND1"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MAR"
)

set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha[1])

fit.MCAR.BAND1 <- NLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("BAND1"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MCAR"
)


source(paste(PATH, "tLMMmissingSAEM.r", sep = ""))
###############
###############
### BAND1
cor.type <- c("BAND1")
set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha)

est.MNAR.BAND1 <- tNLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("BAND1"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MNAR"
)

set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha[-1])
est.MAR.BAND1 <- tNLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("BAND1"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MAR"
)
set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha[1])
est.MCAR.BAND1 <- tNLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("BAND1"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MCAR"
)

###############
###############
### CAR1
source(paste(PATH, "LMMmissingSAEM.r", sep = ""))
cor.type <- c("CAR1")
set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha)

fit.MNAR.CAR1 <- NLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("CAR1"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MNAR"
)

set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha[-1])

fit.MAR.CAR1 <- NLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("CAR1"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MAR"
)

set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha[1])

fit.MCAR.CAR1 <- NLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("CAR1"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MCAR"
)


source(paste(PATH, "tLMMmissingSAEM.r", sep = ""))
###############
###############
### CAR1
cor.type <- c("CAR1")
set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha)

est.MNAR.CAR1 <- tNLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("CAR1"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MNAR"
)

set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha[-1])
est.MAR.CAR1 <- tNLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("CAR1"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MAR"
)
set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha[1])
est.MCAR.CAR1 <- tNLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("CAR1"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MCAR"
)


###############
###############
### DEC
source(paste(PATH, "LMMmissingSAEM.r", sep = ""))
cor.type <- c("DEC")
set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, ga = 1, nu = 5, alpha = alpha)

fit.MNAR.DEC <- NLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("DEC"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MNAR"
)

set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, ga = 1, nu = 5, alpha = alpha[-1])

fit.MAR.DEC <- NLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("DEC"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MAR"
)

set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, ga = 1, nu = 5, alpha = alpha[1])

fit.MCAR.DEC <- NLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("DEC"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MCAR"
)


source(paste(PATH, "tLMMmissingSAEM.r", sep = ""))
###############
###############
### DEC
cor.type <- c("DEC")
set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, ga = 1, nu = 5, alpha = alpha)

est.MNAR.DEC <- tNLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("DEC"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MNAR"
)

set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, ga = 1, nu = 5, alpha = alpha[-1])
est.MAR.DEC <- tNLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("DEC"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MAR"
)
set.seed(20250201)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, ga = 1, nu = 5, alpha = alpha[1])
est.MCAR.DEC <- tNLMM.miss.SAEM(Data,
  g = 1, init.para, cor.type = c("DEC"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MCAR"
)
