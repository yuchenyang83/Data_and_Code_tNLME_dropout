#################################################################################
#
#   Project     :    "Robust HIV Viral Dynamics: A Nonlinear Mixed-Effects Framework for
#                    Heavy-Tailed Data with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#
# Read the ACTG 398 dataset from 'Data_and_Code/Data/source/actg398.txt'.
# You may replace this with your own data if desired.
#
# This R script conducts a sensitivity analysis by fixing the parameter alpha_2,
# which corresponds to the coefficient of the current response in the
# dropout model. This approach allows users to assess the robustness of the
# model estimates under different fixed values of alpha_2.
#
# Model fitting is performed using the tNLME model under the MNAR mechanism
# with an AR(1) within-subject error structure.
#
# Due to the potentially lengthy iterative process, intermediate steps and
# iteration details are not printed. Only the final results are saved to:
# 'Data_and_Code/Data/fixed_alpha.txt'.
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
source(paste0(PATH, "/function/fix_alpha/run_sensitivity_analysis.r"))
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

source(paste(PATH, "/function/fix_alpha/tNLMMmissingSAEM_fixalpha.r", sep = ""))
###############
###############

cor.type <- "ARp"
alpha[5] <- 0.0001

ll <- c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 4, 6)
ll = c(ll, -c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 4, 6))
result <- NULL
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha)
per <- 500
for (iii in 1:length(ll))
{
  kkk <- 20250101 + iii
  set.seed(kkk)
  print(kkk)
  cor.type <- "ARp"
  alpha[5] <- ll[iii]
  init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 1e-6, nu = 5, alpha = alpha)
  est.MNAR.ARp <- try(tNLMM.miss.SAEM(Data,
    g = 1, init.para, cor.type = c("ARp"), M = 10, M.LL = M.LL, P = 1, tol = tol, max.iter = max.iter, per = per,
    mechanism = "MNAR"
  ), silent = F)
  if (class(est.MNAR.ARp) == "try-error") next
  result <- rbind(result, c(kkk, est.MNAR.ARp$model.inf$bic, sum(is.na(est.MNAR.ARp$IM$out[2, ])), est.MNAR.ARp$IM$out[, 6]))
  save.image(paste0(PATH, "est.MNAR.ARp", kkk, ".RData"))
}

