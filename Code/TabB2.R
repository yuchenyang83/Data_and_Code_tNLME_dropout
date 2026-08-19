################################################################################
#
#   Filename    :    TabB2.R
#   Project     :    Biometrics article "Modeling HIV Viral Dynamics Using a
#                    Nonlinear Mixed-Effects Framework for Heavy-Tailed Data
#                    with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#   Date        :    18.08.2026
#   Purpose     :    produce Table B.2 for the simulation study by comparing
#                    AIC, BIC, and MSPE between the NLME and tNLME models under
#                    the MNAR mechanism across sample sizes and dropout levels
#
#   Input data files  :  Data_and_Code/Data/simulation/SS-simulation-t25/SIM1-SIM5/
#                         Data_and_Code/Data/simulation/SS-simulation-t50/SIM1-SIM5/
#                         Data_and_Code/Data/simulation/SS-simulation-t75/SIM1-SIM5/
#
#   Output data files :  Data_and_Code/Result/TableB2.csv
#
#   R Version   :    R-4.6.0
#   Required R packages : none
#
################################################################################
#################
#################  25 %; nu=5
#################
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
PATH0 <- paste(PATH, "/Data/simulation/SS-simulation-t25/", sep = "")
PATH1 <- paste(PATH0, "SIM1/", sep = "")
PATH2 <- paste(PATH0, "SIM2/", sep = "")
PATH3 <- paste(PATH0, "SIM3/", sep = "")
PATH4 <- paste(PATH0, "SIM4/", sep = "")
PATH5 <- paste(PATH0, "SIM5/", sep = "")


loglik1 <- as.matrix(read.table(paste(PATH1, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE1 <- as.matrix(read.table(paste(PATH1, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara1 <- colMeans(as.matrix(read.table(paste(PATH1, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik2 <- as.matrix(read.table(paste(PATH2, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE2 <- as.matrix(read.table(paste(PATH2, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara2 <- colMeans(as.matrix(read.table(paste(PATH2, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik3 <- as.matrix(read.table(paste(PATH3, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE3 <- as.matrix(read.table(paste(PATH3, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara3 <- colMeans(as.matrix(read.table(paste(PATH3, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik4 <- as.matrix(read.table(paste(PATH4, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE4 <- as.matrix(read.table(paste(PATH4, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara4 <- colMeans(as.matrix(read.table(paste(PATH4, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik5 <- as.matrix(read.table(paste(PATH5, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE5 <- as.matrix(read.table(paste(PATH5, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara5 <- colMeans(as.matrix(read.table(paste(PATH5, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

####################
##########
##########   table1
nn <- c(25, 50, 100, 200, 400)
aic1 <- aic2 <- bic1 <- bic2 <- NULL
mMNAR <- length(realpara1)
mMAR <- length(realpara1) - 1
mMCAR <- length(realpara1) - 2

mtMNAR <- length(realpara1) - 1
mtMAR <- length(realpara1) - 2
mtMCAR <- length(realpara1) - 3
mmm <- c(mMNAR, mMAR, mMCAR, mtMNAR, mtMAR, mtMCAR)

n <- nn[1]
AIc <- 2 * (mmm) - 2 * t(loglik1[, 2:7])
BIc <- mmm * log(n) - 2 * t(loglik1[, 2:7])

MSE1 <- t(MSE1)
sum.tabel <- rbind(
  colMeans(loglik1)[-1][c(1, 4)],
  colMeans(t(AIc))[c(1, 4)],
  apply(AIc, 1, sd)[c(1, 4)],
  c(sum(apply(AIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(AIc[c(1, 4), ], 2, order)[2, ] == 1)),
  colMeans(t(BIc))[c(1, 4)],
  apply(BIc, 1, sd)[c(1, 4)],
  c(sum(apply(BIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(BIc[c(1, 4), ], 2, order)[2, ] == 1))
)


n <- nn[2]
AIc <- 2 * (mmm) - 2 * t(loglik2[, 2:7])
BIc <- mmm * log(n) - 2 * t(loglik2[, 2:7])

MSE2 <- t(MSE2)
sum.tabel <- rbind(
  sum.tabel,
  rbind(
    colMeans(loglik2)[-1][c(1, 4)],
    colMeans(t(AIc))[c(1, 4)],
    apply(AIc, 1, sd)[c(1, 4)],
    c(sum(apply(AIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(AIc[c(1, 4), ], 2, order)[2, ] == 1)),
    colMeans(t(BIc))[c(1, 4)],
    apply(BIc, 1, sd)[c(1, 4)],
    c(sum(apply(BIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(BIc[c(1, 4), ], 2, order)[2, ] == 1))
  )
)


n <- nn[3]
AIc <- 2 * (mmm) - 2 * t(loglik3[, 2:7])
BIc <- mmm * log(n) - 2 * t(loglik3[, 2:7])

MSE3 <- t(MSE3)
sum.tabel <- rbind(
  sum.tabel,
  rbind(
    colMeans(loglik3)[-1][c(1, 4)],
    colMeans(t(AIc))[c(1, 4)],
    apply(AIc, 1, sd)[c(1, 4)],
    c(sum(apply(AIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(AIc[c(1, 4), ], 2, order)[2, ] == 1)),
    colMeans(t(BIc))[c(1, 4)],
    apply(BIc, 1, sd)[c(1, 4)],
    c(sum(apply(BIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(BIc[c(1, 4), ], 2, order)[2, ] == 1))
  )
)
n <- nn[4]
AIc <- 2 * (mmm) - 2 * t(loglik4[, 2:7])
BIc <- mmm * log(n) - 2 * t(loglik4[, 2:7])

MSE4 <- t(MSE4)
sum.tabel <- rbind(
  sum.tabel,
  rbind(
    colMeans(loglik4)[-1][c(1, 4)],
    colMeans(t(AIc))[c(1, 4)],
    apply(AIc, 1, sd)[c(1, 4)],
    c(sum(apply(AIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(AIc[c(1, 4), ], 2, order)[2, ] == 1)),
    colMeans(t(BIc))[c(1, 4)],
    apply(BIc, 1, sd)[c(1, 4)],
    c(sum(apply(BIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(BIc[c(1, 4), ], 2, order)[2, ] == 1))
  )
)
n <- nn[2]
AIc <- 2 * (mmm) - 2 * t(loglik5[, 2:7])
BIc <- mmm * log(n) - 2 * t(loglik5[, 2:7])

MSE5 <- t(MSE5)
sum.tabel <- rbind(
  sum.tabel,
  rbind(
    colMeans(loglik5)[-1][c(1, 4)],
    colMeans(t(AIc))[c(1, 4)],
    apply(AIc, 1, sd)[c(1, 4)],
    c(sum(apply(AIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(AIc[c(1, 4), ], 2, order)[2, ] == 1)),
    colMeans(t(BIc))[c(1, 4)],
    apply(BIc, 1, sd)[c(1, 4)],
    c(sum(apply(BIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(BIc[c(1, 4), ], 2, order)[2, ] == 1))
  )
)

round(sum.tabel, 3)

for (i in (c(4, 7, 11, 14, 18, 21, 25, 28, 32, 35)))
{
  if (sum(sum.tabel[i, ][1:2]) > 100) sum.tabel[i, ] <- round((sum.tabel[i, ] / sum(sum.tabel[i, ][1:2])) * 100)
}

#################
#################  50 %; nu=5
#################
PATH0 <- paste(PATH, "/Data/simulation/SS-simulation-t50/", sep = "")
PATH1 <- paste(PATH0, "SIM1/", sep = "")
PATH2 <- paste(PATH0, "SIM2/", sep = "")
PATH3 <- paste(PATH0, "SIM3/", sep = "")
PATH4 <- paste(PATH0, "SIM4/", sep = "")
PATH5 <- paste(PATH0, "SIM5/", sep = "")

loglik1 <- as.matrix(read.table(paste(PATH1, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE1 <- as.matrix(read.table(paste(PATH1, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara1 <- colMeans(as.matrix(read.table(paste(PATH1, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik2 <- as.matrix(read.table(paste(PATH2, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE2 <- as.matrix(read.table(paste(PATH2, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara2 <- colMeans(as.matrix(read.table(paste(PATH2, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik3 <- as.matrix(read.table(paste(PATH3, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE3 <- as.matrix(read.table(paste(PATH3, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara3 <- colMeans(as.matrix(read.table(paste(PATH3, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik4 <- as.matrix(read.table(paste(PATH4, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE4 <- as.matrix(read.table(paste(PATH4, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara4 <- colMeans(as.matrix(read.table(paste(PATH4, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik5 <- as.matrix(read.table(paste(PATH5, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE5 <- as.matrix(read.table(paste(PATH5, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara5 <- colMeans(as.matrix(read.table(paste(PATH5, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

####################
##########
##########   table1
nn <- c(25, 50, 100, 200, 400)
aic1 <- aic2 <- bic1 <- bic2 <- NULL
mMNAR <- length(realpara1)
mMAR <- length(realpara1) - 1
mMCAR <- length(realpara1) - 2

mtMNAR <- length(realpara1) - 1
mtMAR <- length(realpara1) - 2
mtMCAR <- length(realpara1) - 3
mmm <- c(mMNAR, mMAR, mMCAR, mtMNAR, mtMAR, mtMCAR)

n <- nn[1]
AIc <- 2 * (mmm) - 2 * t(loglik1[, 2:7])
BIc <- mmm * log(n) - 2 * t(loglik1[, 2:7])

MSE1 <- t(MSE1)
sum.tabel1 <- rbind(
  colMeans(loglik1)[-1][c(1, 4)],
  colMeans(t(AIc))[c(1, 4)],
  apply(AIc, 1, sd)[c(1, 4)],
  c(sum(apply(AIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(AIc[c(1, 4), ], 2, order)[2, ] == 1)),
  colMeans(t(BIc))[c(1, 4)],
  apply(BIc, 1, sd)[c(1, 4)],
  c(sum(apply(BIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(BIc[c(1, 4), ], 2, order)[2, ] == 1))
)


n <- nn[2]
AIc <- 2 * (mmm) - 2 * t(loglik2[, 2:7])
BIc <- mmm * log(n) - 2 * t(loglik2[, 2:7])

MSE2 <- t(MSE2)
sum.tabel1 <- rbind(
  sum.tabel1,
  rbind(
    colMeans(loglik2)[-1][c(1, 4)],
    colMeans(t(AIc))[c(1, 4)],
    apply(AIc, 1, sd)[c(1, 4)],
    c(sum(apply(AIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(AIc[c(1, 4), ], 2, order)[2, ] == 1)),
    colMeans(t(BIc))[c(1, 4)],
    apply(BIc, 1, sd)[c(1, 4)],
    c(sum(apply(BIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(BIc[c(1, 4), ], 2, order)[2, ] == 1))
  )
)


n <- nn[3]
AIc <- 2 * (mmm) - 2 * t(loglik3[, 2:7])
BIc <- mmm * log(n) - 2 * t(loglik3[, 2:7])

MSE3 <- t(MSE3)
sum.tabel1 <- rbind(
  sum.tabel1,
  rbind(
    colMeans(loglik3)[-1][c(1, 4)],
    colMeans(t(AIc))[c(1, 4)],
    apply(AIc, 1, sd)[c(1, 4)],
    c(sum(apply(AIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(AIc[c(1, 4), ], 2, order)[2, ] == 1)),
    colMeans(t(BIc))[c(1, 4)],
    apply(BIc, 1, sd)[c(1, 4)],
    c(sum(apply(BIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(BIc[c(1, 4), ], 2, order)[2, ] == 1))
  )
)
n <- nn[4]
AIc <- 2 * (mmm) - 2 * t(loglik4[, 2:7])
BIc <- mmm * log(n) - 2 * t(loglik4[, 2:7])

MSE4 <- t(MSE4)
sum.tabel1 <- rbind(
  sum.tabel1,
  rbind(
    colMeans(loglik4)[-1][c(1, 4)],
    colMeans(t(AIc))[c(1, 4)],
    apply(AIc, 1, sd)[c(1, 4)],
    c(sum(apply(AIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(AIc[c(1, 4), ], 2, order)[2, ] == 1)),
    colMeans(t(BIc))[c(1, 4)],
    apply(BIc, 1, sd)[c(1, 4)],
    c(sum(apply(BIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(BIc[c(1, 4), ], 2, order)[2, ] == 1))
  )
)
n <- nn[2]
AIc <- 2 * (mmm) - 2 * t(loglik5[, 2:7])
BIc <- mmm * log(n) - 2 * t(loglik5[, 2:7])

MSE5 <- t(MSE5)
sum.tabel1 <- rbind(
  sum.tabel1,
  rbind(
    colMeans(loglik5)[-1][c(1, 4)],
    colMeans(t(AIc))[c(1, 4)],
    apply(AIc, 1, sd)[c(1, 4)],
    c(sum(apply(AIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(AIc[c(1, 4), ], 2, order)[2, ] == 1)),
    colMeans(t(BIc))[c(1, 4)],
    apply(BIc, 1, sd)[c(1, 4)],
    c(sum(apply(BIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(BIc[c(1, 4), ], 2, order)[2, ] == 1))
  )
)

round(sum.tabel1, 3)

for (i in (c(4, 7, 11, 14, 18, 21, 25, 28, 32, 35)))
{
  if (sum(sum.tabel1[i, ][1:2]) > 100) sum.tabel1[i, ] <- round((sum.tabel1[i, ] / sum(sum.tabel1[i, ][1:2])) * 100)
}

#################
#################  75%; nu=5
#################
PATH0 <- paste(PATH, "/Data/simulation/SS-simulation-t75/", sep = "")
PATH1 <- paste(PATH0, "SIM1/", sep = "")
PATH2 <- paste(PATH0, "SIM2/", sep = "")
PATH3 <- paste(PATH0, "SIM3/", sep = "")
PATH4 <- paste(PATH0, "SIM4/", sep = "")
PATH5 <- paste(PATH0, "SIM5/", sep = "")

loglik1 <- as.matrix(read.table(paste(PATH1, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE1 <- as.matrix(read.table(paste(PATH1, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara1 <- colMeans(as.matrix(read.table(paste(PATH1, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik2 <- as.matrix(read.table(paste(PATH2, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE2 <- as.matrix(read.table(paste(PATH2, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara2 <- colMeans(as.matrix(read.table(paste(PATH2, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik3 <- as.matrix(read.table(paste(PATH3, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE3 <- as.matrix(read.table(paste(PATH3, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara3 <- colMeans(as.matrix(read.table(paste(PATH3, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik4 <- as.matrix(read.table(paste(PATH4, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE4 <- as.matrix(read.table(paste(PATH4, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara4 <- colMeans(as.matrix(read.table(paste(PATH4, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik5 <- as.matrix(read.table(paste(PATH5, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE5 <- as.matrix(read.table(paste(PATH5, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara5 <- colMeans(as.matrix(read.table(paste(PATH5, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

####################
##########
##########   table1
nn <- c(25, 50, 100, 200, 400)
aic1 <- aic2 <- bic1 <- bic2 <- NULL
mMNAR <- length(realpara1)
mMAR <- length(realpara1) - 1
mMCAR <- length(realpara1) - 2

mtMNAR <- length(realpara1) - 1
mtMAR <- length(realpara1) - 2
mtMCAR <- length(realpara1) - 3
mmm <- c(mMNAR, mMAR, mMCAR, mtMNAR, mtMAR, mtMCAR)

n <- nn[1]
AIc <- 2 * (mmm) - 2 * t(loglik1[, 2:7])
BIc <- mmm * log(n) - 2 * t(loglik1[, 2:7])

MSE1 <- t(MSE1)
sum.tabel2 <- rbind(
  colMeans(loglik1)[-1][c(1, 4)],
  colMeans(t(AIc))[c(1, 4)],
  apply(AIc, 1, sd)[c(1, 4)],
  c(sum(apply(AIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(AIc[c(1, 4), ], 2, order)[2, ] == 1)),
  colMeans(t(BIc))[c(1, 4)],
  apply(BIc, 1, sd)[c(1, 4)],
  c(sum(apply(BIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(BIc[c(1, 4), ], 2, order)[2, ] == 1))
)


n <- nn[2]
AIc <- 2 * (mmm) - 2 * t(loglik2[, 2:7])
BIc <- mmm * log(n) - 2 * t(loglik2[, 2:7])

MSE2 <- t(MSE2)
sum.tabel2 <- rbind(
  sum.tabel2,
  rbind(
    colMeans(loglik2)[-1][c(1, 4)],
    colMeans(t(AIc))[c(1, 4)],
    apply(AIc, 1, sd)[c(1, 4)],
    c(sum(apply(AIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(AIc[c(1, 4), ], 2, order)[2, ] == 1)),
    colMeans(t(BIc))[c(1, 4)],
    apply(BIc, 1, sd)[c(1, 4)],
    c(sum(apply(BIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(BIc[c(1, 4), ], 2, order)[2, ] == 1))
  )
)


n <- nn[3]
AIc <- 2 * (mmm) - 2 * t(loglik3[, 2:7])
BIc <- mmm * log(n) - 2 * t(loglik3[, 2:7])

MSE3 <- t(MSE3)
sum.tabel2 <- rbind(
  sum.tabel2,
  rbind(
    colMeans(loglik3)[-1][c(1, 4)],
    colMeans(t(AIc))[c(1, 4)],
    apply(AIc, 1, sd)[c(1, 4)],
    c(sum(apply(AIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(AIc[c(1, 4), ], 2, order)[2, ] == 1)),
    colMeans(t(BIc))[c(1, 4)],
    apply(BIc, 1, sd)[c(1, 4)],
    c(sum(apply(BIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(BIc[c(1, 4), ], 2, order)[2, ] == 1))
  )
)
n <- nn[4]
AIc <- 2 * (mmm) - 2 * t(loglik4[, 2:7])
BIc <- mmm * log(n) - 2 * t(loglik4[, 2:7])

MSE4 <- t(MSE4)
sum.tabel2 <- rbind(
  sum.tabel2,
  rbind(
    colMeans(loglik4)[-1][c(1, 4)],
    colMeans(t(AIc))[c(1, 4)],
    apply(AIc, 1, sd)[c(1, 4)],
    c(sum(apply(AIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(AIc[c(1, 4), ], 2, order)[2, ] == 1)),
    colMeans(t(BIc))[c(1, 4)],
    apply(BIc, 1, sd)[c(1, 4)],
    c(sum(apply(BIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(BIc[c(1, 4), ], 2, order)[2, ] == 1))
  )
)
n <- nn[2]
AIc <- 2 * (mmm) - 2 * t(loglik5[, 2:7])
BIc <- mmm * log(n) - 2 * t(loglik5[, 2:7])

MSE5 <- t(MSE5)
sum.tabel2 <- rbind(
  sum.tabel2,
  rbind(
    colMeans(loglik5)[-1][c(1, 4)],
    colMeans(t(AIc))[c(1, 4)],
    apply(AIc, 1, sd)[c(1, 4)],
    c(sum(apply(AIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(AIc[c(1, 4), ], 2, order)[2, ] == 1)),
    colMeans(t(BIc))[c(1, 4)],
    apply(BIc, 1, sd)[c(1, 4)],
    c(sum(apply(BIc[c(1, 4), ], 2, order)[1, ] == 1), sum(apply(BIc[c(1, 4), ], 2, order)[2, ] == 1))
  )
)

for (i in (c(4, 7, 11, 14, 18, 21, 25, 28, 32, 35)))
{
  if (sum(sum.tabel2[i, ][1:2]) > 100) sum.tabel2[i, ] <- round((sum.tabel2[i, ] / sum(sum.tabel2[i, ][1:2])) * 100)
}

sum.tabel.log.c <- round(cbind(sum.tabel, sum.tabel1, sum.tabel2), 3)[-c(
  1, 3, 6,
  8, 10, 13,
  15, 17, 20,
  22, 24, 27,
  29, 31, 34
), c(2, 1, 4, 3, 6, 5)]


################################################################################
#################
#################  25 %; nu=5
#################
PATH0 <- paste(PATH, "/Data/simulation/SS-simulation-t25/", sep = "")
PATH1 <- paste(PATH0, "SIM1/", sep = "")
PATH2 <- paste(PATH0, "SIM2/", sep = "")
PATH3 <- paste(PATH0, "SIM3/", sep = "")
PATH4 <- paste(PATH0, "SIM4/", sep = "")
PATH5 <- paste(PATH0, "SIM5/", sep = "")

loglik1 <- as.matrix(read.table(paste(PATH1, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE1 <- as.matrix(read.table(paste(PATH1, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara1 <- colMeans(as.matrix(read.table(paste(PATH1, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik2 <- as.matrix(read.table(paste(PATH2, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE2 <- as.matrix(read.table(paste(PATH2, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara2 <- colMeans(as.matrix(read.table(paste(PATH2, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik3 <- as.matrix(read.table(paste(PATH3, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE3 <- as.matrix(read.table(paste(PATH3, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara3 <- colMeans(as.matrix(read.table(paste(PATH3, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik4 <- as.matrix(read.table(paste(PATH4, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE4 <- as.matrix(read.table(paste(PATH4, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara4 <- colMeans(as.matrix(read.table(paste(PATH4, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik5 <- as.matrix(read.table(paste(PATH5, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE5 <- as.matrix(read.table(paste(PATH5, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara5 <- colMeans(as.matrix(read.table(paste(PATH5, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))


MSE1 <- t(MSE1)
sum.tabel <- rbind(
  colMeans(t(MSE1))[-1][c(1, 4)],
  apply(MSE1, 1, sd)[-1][c(1, 4)],
  c(sum(apply(MSE1[-1, ][c(1, 4), ], 2, order)[1, ] == 1), sum(apply(MSE1[-1, ][c(1, 4), ], 2, order)[2, ] == 1))
)

MSE2 <- t(MSE2)
sum.tabel <- rbind(
  sum.tabel,
  rbind(
    colMeans(t(MSE2))[-1][c(1, 4)],
    apply(MSE2, 1, sd)[-1][c(1, 4)],
    c(sum(apply(MSE2[-1, ][c(1, 4), ], 2, order)[1, ] == 1), sum(apply(MSE2[-1, ][c(1, 4), ], 2, order)[2, ] == 1))
  )
)

MSE3 <- t(MSE3)
sum.tabel <- rbind(
  sum.tabel,
  rbind(
    colMeans(t(MSE3))[-1][c(1, 4)],
    apply(MSE3, 1, sd)[-1][c(1, 4)],
    c(sum(apply(MSE3[-1, ][c(1, 4), ], 2, order)[1, ] == 1), sum(apply(MSE3[-1, ][c(1, 4), ], 2, order)[2, ] == 1))
  )
)

MSE4 <- t(MSE4)
sum.tabel <- rbind(
  sum.tabel,
  rbind(
    colMeans(t(MSE4))[-1][c(1, 4)],
    apply(MSE4, 1, sd)[-1][c(1, 4)],
    c(sum(apply(MSE4[-1, ][c(1, 4), ], 2, order)[1, ] == 1), sum(apply(MSE4[-1, ][c(1, 4), ], 2, order)[2, ] == 1))
  )
)

MSE5 <- t(MSE5)
sum.tabel <- rbind(
  sum.tabel,
  rbind(
    colMeans(t(MSE5))[-1][c(1, 4)],
    apply(MSE5, 1, sd)[-1][c(1, 4)],
    c(sum(apply(MSE5[-1, ][c(1, 4), ], 2, order)[1, ] == 1), sum(apply(MSE5[-1, ][c(1, 4), ], 2, order)[2, ] == 1))
  )
)

round(sum.tabel, 3)



for (i in c(3, 6, 9, 12, 15))
{
  if (sum(sum.tabel[i, ][1:2]) > 100) sum.tabel[i, ] <- round((sum.tabel[i, ] / sum(sum.tabel[i, ][1:2])) * 100)
}

sum.tabel.TT <- NULL
sum.tabel.TT <- cbind(sum.tabel.TT, round(sum.tabel, 3))

#################
#################  50 %; nu=5
#################
PATH0 <- paste(PATH, "/Data/simulation/SS-simulation-t50/", sep = "")
PATH1 <- paste(PATH0, "SIM1/", sep = "")
PATH2 <- paste(PATH0, "SIM2/", sep = "")
PATH3 <- paste(PATH0, "SIM3/", sep = "")
PATH4 <- paste(PATH0, "SIM4/", sep = "")
PATH5 <- paste(PATH0, "SIM5/", sep = "")

loglik1 <- as.matrix(read.table(paste(PATH1, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE1 <- as.matrix(read.table(paste(PATH1, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara1 <- colMeans(as.matrix(read.table(paste(PATH1, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik2 <- as.matrix(read.table(paste(PATH2, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE2 <- as.matrix(read.table(paste(PATH2, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara2 <- colMeans(as.matrix(read.table(paste(PATH2, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik3 <- as.matrix(read.table(paste(PATH3, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE3 <- as.matrix(read.table(paste(PATH3, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara3 <- colMeans(as.matrix(read.table(paste(PATH3, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik4 <- as.matrix(read.table(paste(PATH4, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE4 <- as.matrix(read.table(paste(PATH4, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara4 <- colMeans(as.matrix(read.table(paste(PATH4, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik5 <- as.matrix(read.table(paste(PATH5, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE5 <- as.matrix(read.table(paste(PATH5, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara5 <- colMeans(as.matrix(read.table(paste(PATH5, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))


MSE1 <- t(MSE1)
sum.tabel <- rbind(
  colMeans(t(MSE1))[-1][c(1, 4)],
  apply(MSE1, 1, sd)[-1][c(1, 4)],
  c(sum(apply(MSE1[-1, ][c(1, 4), ], 2, order)[1, ] == 1), sum(apply(MSE1[-1, ][c(1, 4), ], 2, order)[2, ] == 1))
)

MSE2 <- t(MSE2)
sum.tabel <- rbind(
  sum.tabel,
  rbind(
    colMeans(t(MSE2))[-1][c(1, 4)],
    apply(MSE2, 1, sd)[-1][c(1, 4)],
    c(sum(apply(MSE2[-1, ][c(1, 4), ], 2, order)[1, ] == 1), sum(apply(MSE2[-1, ][c(1, 4), ], 2, order)[2, ] == 1))
  )
)

MSE3 <- t(MSE3)
sum.tabel <- rbind(
  sum.tabel,
  rbind(
    colMeans(t(MSE3))[-1][c(1, 4)],
    apply(MSE3, 1, sd)[-1][c(1, 4)],
    c(sum(apply(MSE3[-1, ][c(1, 4), ], 2, order)[1, ] == 1), sum(apply(MSE3[-1, ][c(1, 4), ], 2, order)[2, ] == 1))
  )
)

MSE4 <- t(MSE4)
sum.tabel <- rbind(
  sum.tabel,
  rbind(
    colMeans(t(MSE4))[-1][c(1, 4)],
    apply(MSE4, 1, sd)[-1][c(1, 4)],
    c(sum(apply(MSE4[-1, ][c(1, 4), ], 2, order)[1, ] == 1), sum(apply(MSE4[-1, ][c(1, 4), ], 2, order)[2, ] == 1))
  )
)

MSE5 <- t(MSE5)
sum.tabel <- rbind(
  sum.tabel,
  rbind(
    colMeans(t(MSE5))[-1][c(1, 4)],
    apply(MSE5, 1, sd)[-1][c(1, 4)],
    c(sum(apply(MSE5[-1, ][c(1, 4), ], 2, order)[1, ] == 1), sum(apply(MSE5[-1, ][c(1, 4), ], 2, order)[2, ] == 1))
  )
)

round(sum.tabel, 3)



for (i in c(3, 6, 9, 12, 15))
{
  if (sum(sum.tabel[i, ][1:2]) > 100) sum.tabel[i, ] <- round((sum.tabel[i, ] / sum(sum.tabel[i, ][1:2])) * 100)
}


sum.tabel.TT <- cbind(sum.tabel.TT, round(sum.tabel, 3))


#################
#################  75%; nu=5
#################
PATH0 <- paste(PATH, "/Data/simulation/SS-simulation-t75/", sep = "")
PATH1 <- paste(PATH0, "SIM1/", sep = "")
PATH2 <- paste(PATH0, "SIM2/", sep = "")
PATH3 <- paste(PATH0, "SIM3/", sep = "")
PATH4 <- paste(PATH0, "SIM4/", sep = "")
PATH5 <- paste(PATH0, "SIM5/", sep = "")

loglik1 <- as.matrix(read.table(paste(PATH1, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE1 <- as.matrix(read.table(paste(PATH1, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara1 <- colMeans(as.matrix(read.table(paste(PATH1, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik2 <- as.matrix(read.table(paste(PATH2, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE2 <- as.matrix(read.table(paste(PATH2, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara2 <- colMeans(as.matrix(read.table(paste(PATH2, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik3 <- as.matrix(read.table(paste(PATH3, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE3 <- as.matrix(read.table(paste(PATH3, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara3 <- colMeans(as.matrix(read.table(paste(PATH3, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik4 <- as.matrix(read.table(paste(PATH4, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE4 <- as.matrix(read.table(paste(PATH4, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara4 <- colMeans(as.matrix(read.table(paste(PATH4, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))

loglik5 <- as.matrix(read.table(paste(PATH5, "loglik.txt", sep = ""), na.strings = "NA", sep = ""))
MSE5 <- as.matrix(read.table(paste(PATH5, "MSE.txt", sep = ""), na.strings = "NA", sep = ""))
realpara5 <- colMeans(as.matrix(read.table(paste(PATH5, "realpara.txt", sep = ""), na.strings = "NA", sep = "")))


MSE1 <- t(MSE1)
sum.tabel <- rbind(
  colMeans(t(MSE1))[-1][c(1, 4)],
  apply(MSE1, 1, sd)[-1][c(1, 4)],
  c(sum(apply(MSE1[-1, ][c(1, 4), ], 2, order)[1, ] == 1), sum(apply(MSE1[-1, ][c(1, 4), ], 2, order)[2, ] == 1))
)

MSE2 <- t(MSE2)
sum.tabel <- rbind(
  sum.tabel,
  rbind(
    colMeans(t(MSE2))[-1][c(1, 4)],
    apply(MSE2, 1, sd)[-1][c(1, 4)],
    c(sum(apply(MSE2[-1, ][c(1, 4), ], 2, order)[1, ] == 1), sum(apply(MSE2[-1, ][c(1, 4), ], 2, order)[2, ] == 1))
  )
)

MSE3 <- t(MSE3)
sum.tabel <- rbind(
  sum.tabel,
  rbind(
    colMeans(t(MSE3))[-1][c(1, 4)],
    apply(MSE3, 1, sd)[-1][c(1, 4)],
    c(sum(apply(MSE3[-1, ][c(1, 4), ], 2, order)[1, ] == 1), sum(apply(MSE3[-1, ][c(1, 4), ], 2, order)[2, ] == 1))
  )
)

MSE4 <- t(MSE4)
sum.tabel <- rbind(
  sum.tabel,
  rbind(
    colMeans(t(MSE4))[-1][c(1, 4)],
    apply(MSE4, 1, sd)[-1][c(1, 4)],
    c(sum(apply(MSE4[-1, ][c(1, 4), ], 2, order)[1, ] == 1), sum(apply(MSE4[-1, ][c(1, 4), ], 2, order)[2, ] == 1))
  )
)

MSE5 <- t(MSE5)
sum.tabel <- rbind(
  sum.tabel,
  rbind(
    colMeans(t(MSE5))[-1][c(1, 4)],
    apply(MSE5, 1, sd)[-1][c(1, 4)],
    c(sum(apply(MSE5[-1, ][c(1, 4), ], 2, order)[1, ] == 1), sum(apply(MSE5[-1, ][c(1, 4), ], 2, order)[2, ] == 1))
  )
)

round(sum.tabel, 3)

sum.tabel.TT <- cbind(sum.tabel.TT, round(sum.tabel, 3))

sum.table <- sum.tabel.TT[-c(2, 5, 8, 11, 14), c(2, 1, 4, 3, 6, 5)]

sum.table <- rbind(
  sum.tabel.log.c[1:4, ], sum.table[1:2, ],
  sum.tabel.log.c[5:8, ], sum.table[3:4, ],
  sum.tabel.log.c[9:12, ], sum.table[5:6, ],
  sum.tabel.log.c[13:16, ], sum.table[7:8, ],
  sum.tabel.log.c[17:20, ], sum.table[9:10, ]
)

colnames(sum.table) <- c("NLME", "tNLME", "NLME", "tNLME", "NLME", "tNLME")
row.names(sum.table) <- c(rep(c("AIC", "(Freq)", "BIC", "(Freq)", "MSPE", "(Freq)"), time = 5))


write.csv(sum.table, paste0(PATH, "/Result/TableB2.csv"))
