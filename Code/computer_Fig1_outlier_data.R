################################################################################
#
#   Filename    :    computer_Fig1_outlier_data.R
#   Project     :    Biometrics article "Modeling HIV Viral Dynamics Using a
#                    Nonlinear Mixed-Effects Framework for Heavy-Tailed Data
#                    with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#   Purpose     :    compute the subject-level outlier identifiers used in
#                    Figure 1 from the selected Scenario (III) AR(1) tNLME
#                    model under MNAR dropout
#
#   Input data files  :  Data/fit_III_result.RData
#                        function/analyze_realdata_AIDS.R
#
#   Output data file  :  Data/Figure1_outlier_subjects.txt
#   R Version   :    R-4.6.0
#
################################################################################
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
load(file.path(PATH, "Data", "fit_III_result.RData"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(PATH, "function", "analyze_realdata_AIDS.R"))

Beta  <- as.matrix(fit.t.III.ARp.MNAR$para.est$Beta)
sigma <- fit.t.III.ARp.MNAR$para.est$sigma
DD    <- fit.t.III.ARp.MNAR$para.est$D
Phi   <- fit.t.III.ARp.MNAR$para.est$Phi
b.hat <- as.matrix(fit.t.III.ARp.MNAR$para.est$b)
nu    <- fit.t.III.ARp.MNAR$para.est$nu
ytilde <- as.numeric(fit.t.III.ARp.MNAR$y.c)

p <- length(Beta)
q <- nrow(DD)

N <- length(unique(Data$Subject))

ni <- numeric(N)
for (i in 1:N) {
  ni[i] <- length(Data$Subject[Data$Subject == i])
}

n <- sum(ni)

na.ind <- which(is.na(as.vector(t(Data$Var1))))
cumsum.ni <- cumsum(ni)
cumsum.q <- cumsum(rep(q, N))

ni.o <- numeric(N)
for (i in 1:N) {
  ni.o[i] <- sum(!is.na(Data$Var1[Data$Subject == i]))
}
cumsum.ni.o <- cumsum(ni.o)

A <- diag(p)

B <- matrix(c(
  1, 0, 0, 0, 0, 0, 0, 0,
  0, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 0, 0, 0,
  0, 0, 0, 0, 0, 1, 0, 0
), ncol = q)

MU <- NULL
TXtilde <- 0
TZtilde <- matrix(0, ncol = N * q, nrow = n)
TLam <- TLam.inv <- matrix(0, n, n)

cor.type <- "ARp"

for (i in 1:N) {
  
  if (i == 1) {
    idx1 <- 1:cumsum.ni[1]
    idx2 <- 1:cumsum.q[1]
  } else {
    idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
    idx2 <- (cumsum.q[i - 1] + 1):cumsum.q[i]
  }
  
  eta.ij <- A %*% Beta + B %*% b.hat[idx2, ]
  
  MU <- c(
    MU,
    mu.fn(
      eta.ij,
      Data$Time[Data$Subject == i],
      Data$nnrti[Data$Subject == i],
      Data$trtarm[Data$Subject == i]
    )
  )
  
  dmu.ij <- dmu(
    eta.ij,
    Data$Time[Data$Subject == i],
    Data$nnrti[Data$Subject == i],
    Data$trtarm[Data$Subject == i]
  )
  
  TXtilde <- rbind(TXtilde, dmu.ij %*% A)
  
  Ztilde <- dmu.ij %*% B
  
  TZtilde[idx1, ((i - 1) * q + 1):(i * q)] <- Ztilde
  
  Lam <- Ztilde %*% DD %*% t(Ztilde) +
    sigma * cor.fn(
      Phi,
      dim = ni[i],
      type = cor.type,
      Ti = Data$Time[Data$Subject == i],
      ga = NULL
    )
  
  TLam[idx1, idx1] <- Lam
  TLam.inv[idx1, idx1] <- solve(Lam)
}

TXtilde <- TXtilde[-1, ]

ytilde.o <- ytilde[-na.ind]
Xtilde.o <- TXtilde[-na.ind, ]

TLam.oo <- TLam[-na.ind, -na.ind]

TLam.oo.inv <- matrix(
  0,
  ncol = sum(ni.o),
  nrow = sum(ni.o)
)

for (i in 1:N) {
  
  if (i == 1) {
    idx1 <- 1:cumsum.ni.o[1]
  } else {
    idx1 <- (cumsum.ni.o[i - 1] + 1):cumsum.ni.o[i]
  }
  
  TLam.oo.inv[idx1, idx1] <-
    solve(TLam.oo[idx1, idx1])
}

yo.cent <- ytilde.o - Xtilde.o %*% Beta

Delta <- numeric(N)

for (i in 1:N) {
  
  if (i == 1) {
    idx1 <- 1:cumsum.ni.o[1]
  } else {
    idx1 <- (cumsum.ni.o[i - 1] + 1):cumsum.ni.o[i]
  }
  
  Delta[i] <-
    t(yo.cent[idx1]) %*%
    TLam.oo.inv[idx1, idx1] %*%
    yo.cent[idx1]
}

tau.hat <- (nu + ni) / (nu + Delta)

beta.cut <- (1 + ni / nu) *
  qbeta(
    0.025,
    shape1 = nu / 2,
    shape2 = ni / 2
  )

out.flag <- tau.hat < beta.cut

subject.id <- sort(unique(Data$Subject))
out_subj <- subject.id[out.flag]

outlier.data <- data.frame(
  Subject = out_subj
)

write.table(
  outlier.data,
  file = file.path(
    PATH,
    "Data",
    "Figure1_outlier_subjects.txt"
  ),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

outlier.data
nrow(outlier.data)