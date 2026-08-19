################################################################################
#
#   Filename    :    Tab1.R
#   Project     :    Biometrics article "Modeling HIV Viral Dynamics Using a
#                    Nonlinear Mixed-Effects Framework for Heavy-Tailed Data
#                    with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#   Date        :    18.08.2026
#   Purpose     :    produce Table 1 for the ACTG 398 data by comparing the
#                    72 candidate NLME and tNLME models across random-effects
#                    Scenarios (I)-(III), four within-subject correlation
#                    structures, and three missingness mechanisms
#
#   Input data files  :  Data_and_Code/Data/fit_I_result.RData
#                         Data_and_Code/Data/fit_II_result.RData
#                         Data_and_Code/Data/fit_III_result.RData
#
#   Intermediate file :  Data_and_Code/Data/Table1_raw.txt
#
#   Output data files :  Data_and_Code/Result/Table1.csv
#
#   R Version   :    R-4.6.0
#   Required R packages : none
#
################################################################################

load(paste0(PATH, "/Data/fit_I_result.RData"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
load(paste0(PATH, "/Data/fit_II_result.RData"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
load(paste0(PATH, "/Data/fit_III_result.RData"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

N <- length(unique(Data$Subject))

if (!dir.exists(file.path(PATH, "Result"))) {
  dir.create(file.path(PATH, "Result"), recursive = TRUE)
}


############################################################
## Structure (I)
############################################################

Tab.I.MCAR <- rbind(
  c(dim(fit.N.I.UNC.MCAR$IM$out)[2], dim(fit.t.I.UNC.MCAR$IM$out)[2], dim(fit.N.I.CS.MCAR$IM$out)[2], dim(fit.t.I.CS.MCAR$IM$out)[2], dim(fit.N.I.ARp.MCAR$IM$out)[2], dim(fit.t.I.ARp.MCAR$IM$out)[2], dim(fit.N.I.BAND1.MCAR$IM$out)[2], dim(fit.t.I.BAND1.MCAR$IM$out)[2]),
  c(fit.N.I.UNC.MCAR$model.inf$loglik, fit.t.I.UNC.MCAR$model.inf$loglik, fit.N.I.CS.MCAR$model.inf$loglik, fit.t.I.CS.MCAR$model.inf$loglik, fit.N.I.ARp.MCAR$model.inf$loglik, fit.t.I.ARp.MCAR$model.inf$loglik, fit.N.I.BAND1.MCAR$model.inf$loglik, fit.t.I.BAND1.MCAR$model.inf$loglik),
  c(2 * dim(fit.N.I.UNC.MCAR$IM$out)[2] - 2 * fit.N.I.UNC.MCAR$model.inf$loglik, 2 * dim(fit.t.I.UNC.MCAR$IM$out)[2] - 2 * fit.t.I.UNC.MCAR$model.inf$loglik, 2 * dim(fit.N.I.CS.MCAR$IM$out)[2] - 2 * fit.N.I.CS.MCAR$model.inf$loglik, 2 * dim(fit.t.I.CS.MCAR$IM$out)[2] - 2 * fit.t.I.CS.MCAR$model.inf$loglik, 2 * dim(fit.N.I.ARp.MCAR$IM$out)[2] - 2 * fit.N.I.ARp.MCAR$model.inf$loglik, 2 * dim(fit.t.I.ARp.MCAR$IM$out)[2] - 2 * fit.t.I.ARp.MCAR$model.inf$loglik, 2 * dim(fit.N.I.BAND1.MCAR$IM$out)[2] - 2 * fit.N.I.BAND1.MCAR$model.inf$loglik, 2 * dim(fit.t.I.BAND1.MCAR$IM$out)[2] - 2 * fit.t.I.BAND1.MCAR$model.inf$loglik),
  c(dim(fit.N.I.UNC.MCAR$IM$out)[2] * log(N) - 2 * fit.N.I.UNC.MCAR$model.inf$loglik, dim(fit.t.I.UNC.MCAR$IM$out)[2] * log(N) - 2 * fit.t.I.UNC.MCAR$model.inf$loglik, dim(fit.N.I.CS.MCAR$IM$out)[2] * log(N) - 2 * fit.N.I.CS.MCAR$model.inf$loglik, dim(fit.t.I.CS.MCAR$IM$out)[2] * log(N) - 2 * fit.t.I.CS.MCAR$model.inf$loglik, dim(fit.N.I.ARp.MCAR$IM$out)[2] * log(N) - 2 * fit.N.I.ARp.MCAR$model.inf$loglik, dim(fit.t.I.ARp.MCAR$IM$out)[2] * log(N) - 2 * fit.t.I.ARp.MCAR$model.inf$loglik, dim(fit.N.I.BAND1.MCAR$IM$out)[2] * log(N) - 2 * fit.N.I.BAND1.MCAR$model.inf$loglik, dim(fit.t.I.BAND1.MCAR$IM$out)[2] * log(N) - 2 * fit.t.I.BAND1.MCAR$model.inf$loglik),
  c(fit.N.I.UNC.MCAR$iter, fit.t.I.UNC.MCAR$iter, fit.N.I.CS.MCAR$iter, fit.t.I.CS.MCAR$iter, fit.N.I.ARp.MCAR$iter, fit.t.I.ARp.MCAR$iter, fit.N.I.BAND1.MCAR$iter, fit.t.I.BAND1.MCAR$iter),
  c(fit.N.I.UNC.MCAR$model.inf$time / 3600, fit.t.I.UNC.MCAR$model.inf$time / 3600, fit.N.I.CS.MCAR$model.inf$time / 3600, fit.t.I.CS.MCAR$model.inf$time / 3600, fit.N.I.ARp.MCAR$model.inf$time / 3600, fit.t.I.ARp.MCAR$model.inf$time / 3600, fit.N.I.BAND1.MCAR$model.inf$time / 3600, fit.t.I.BAND1.MCAR$model.inf$time / 3600),
  c(
    mean(colMeans(fit.N.I.UNC.MCAR$Taccept.rate)[colMeans(fit.N.I.UNC.MCAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.I.UNC.MCAR$Taccept.rate)[colMeans(fit.t.I.UNC.MCAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.I.CS.MCAR$Taccept.rate)[colMeans(fit.N.I.CS.MCAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.I.CS.MCAR$Taccept.rate)[colMeans(fit.t.I.CS.MCAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.I.ARp.MCAR$Taccept.rate)[colMeans(fit.N.I.ARp.MCAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.I.ARp.MCAR$Taccept.rate)[colMeans(fit.t.I.ARp.MCAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.I.BAND1.MCAR$Taccept.rate)[colMeans(fit.N.I.BAND1.MCAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.I.BAND1.MCAR$Taccept.rate)[colMeans(fit.t.I.BAND1.MCAR$Taccept.rate) != 0])
  )
)

Tab.I.MAR <- rbind(
  c(dim(fit.N.I.UNC.MAR$IM$out)[2], dim(fit.t.I.UNC.MAR$IM$out)[2], dim(fit.N.I.CS.MAR$IM$out)[2], dim(fit.t.I.CS.MAR$IM$out)[2], dim(fit.N.I.ARp.MAR$IM$out)[2], dim(fit.t.I.ARp.MAR$IM$out)[2], dim(fit.N.I.BAND1.MAR$IM$out)[2], dim(fit.t.I.BAND1.MAR$IM$out)[2]),
  c(fit.N.I.UNC.MAR$model.inf$loglik, fit.t.I.UNC.MAR$model.inf$loglik, fit.N.I.CS.MAR$model.inf$loglik, fit.t.I.CS.MAR$model.inf$loglik, fit.N.I.ARp.MAR$model.inf$loglik, fit.t.I.ARp.MAR$model.inf$loglik, fit.N.I.BAND1.MAR$model.inf$loglik, fit.t.I.BAND1.MAR$model.inf$loglik),
  c(2 * dim(fit.N.I.UNC.MAR$IM$out)[2] - 2 * fit.N.I.UNC.MAR$model.inf$loglik, 2 * dim(fit.t.I.UNC.MAR$IM$out)[2] - 2 * fit.t.I.UNC.MAR$model.inf$loglik, 2 * dim(fit.N.I.CS.MAR$IM$out)[2] - 2 * fit.N.I.CS.MAR$model.inf$loglik, 2 * dim(fit.t.I.CS.MAR$IM$out)[2] - 2 * fit.t.I.CS.MAR$model.inf$loglik, 2 * dim(fit.N.I.ARp.MAR$IM$out)[2] - 2 * fit.N.I.ARp.MAR$model.inf$loglik, 2 * dim(fit.t.I.ARp.MAR$IM$out)[2] - 2 * fit.t.I.ARp.MAR$model.inf$loglik, 2 * dim(fit.N.I.BAND1.MAR$IM$out)[2] - 2 * fit.N.I.BAND1.MAR$model.inf$loglik, 2 * dim(fit.t.I.BAND1.MAR$IM$out)[2] - 2 * fit.t.I.BAND1.MAR$model.inf$loglik),
  c(dim(fit.N.I.UNC.MAR$IM$out)[2] * log(N) - 2 * fit.N.I.UNC.MAR$model.inf$loglik, dim(fit.t.I.UNC.MAR$IM$out)[2] * log(N) - 2 * fit.t.I.UNC.MAR$model.inf$loglik, dim(fit.N.I.CS.MAR$IM$out)[2] * log(N) - 2 * fit.N.I.CS.MAR$model.inf$loglik, dim(fit.t.I.CS.MAR$IM$out)[2] * log(N) - 2 * fit.t.I.CS.MAR$model.inf$loglik, dim(fit.N.I.ARp.MAR$IM$out)[2] * log(N) - 2 * fit.N.I.ARp.MAR$model.inf$loglik, dim(fit.t.I.ARp.MAR$IM$out)[2] * log(N) - 2 * fit.t.I.ARp.MAR$model.inf$loglik, dim(fit.N.I.BAND1.MAR$IM$out)[2] * log(N) - 2 * fit.N.I.BAND1.MAR$model.inf$loglik, dim(fit.t.I.BAND1.MAR$IM$out)[2] * log(N) - 2 * fit.t.I.BAND1.MAR$model.inf$loglik),
  c(fit.N.I.UNC.MAR$iter, fit.t.I.UNC.MAR$iter, fit.N.I.CS.MAR$iter, fit.t.I.CS.MAR$iter, fit.N.I.ARp.MAR$iter, fit.t.I.ARp.MAR$iter, fit.N.I.BAND1.MAR$iter, fit.t.I.BAND1.MAR$iter),
  c(fit.N.I.UNC.MAR$model.inf$time / 3600, fit.t.I.UNC.MAR$model.inf$time / 3600, fit.N.I.CS.MAR$model.inf$time / 3600, fit.t.I.CS.MAR$model.inf$time / 3600, fit.N.I.ARp.MAR$model.inf$time / 3600, fit.t.I.ARp.MAR$model.inf$time / 3600, fit.N.I.BAND1.MAR$model.inf$time / 3600, fit.t.I.BAND1.MAR$model.inf$time / 3600),
  c(
    mean(colMeans(fit.N.I.UNC.MAR$Taccept.rate)[colMeans(fit.N.I.UNC.MAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.I.UNC.MAR$Taccept.rate)[colMeans(fit.t.I.UNC.MAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.I.CS.MAR$Taccept.rate)[colMeans(fit.N.I.CS.MAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.I.CS.MAR$Taccept.rate)[colMeans(fit.t.I.CS.MAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.I.ARp.MAR$Taccept.rate)[colMeans(fit.N.I.ARp.MAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.I.ARp.MAR$Taccept.rate)[colMeans(fit.t.I.ARp.MAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.I.BAND1.MAR$Taccept.rate)[colMeans(fit.N.I.BAND1.MAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.I.BAND1.MAR$Taccept.rate)[colMeans(fit.t.I.BAND1.MAR$Taccept.rate) != 0])
  )
)

Tab.I.MNAR <- rbind(
  c(dim(fit.N.I.UNC.MNAR$IM$out)[2], dim(fit.t.I.UNC.MNAR$IM$out)[2], dim(fit.N.I.CS.MNAR$IM$out)[2], dim(fit.t.I.CS.MNAR$IM$out)[2], dim(fit.N.I.ARp.MNAR$IM$out)[2], dim(fit.t.I.ARp.MNAR$IM$out)[2], dim(fit.N.I.BAND1.MNAR$IM$out)[2], dim(fit.t.I.BAND1.MNAR$IM$out)[2]),
  c(fit.N.I.UNC.MNAR$model.inf$loglik, fit.t.I.UNC.MNAR$model.inf$loglik, fit.N.I.CS.MNAR$model.inf$loglik, fit.t.I.CS.MNAR$model.inf$loglik, fit.N.I.ARp.MNAR$model.inf$loglik, fit.t.I.ARp.MNAR$model.inf$loglik, fit.N.I.BAND1.MNAR$model.inf$loglik, fit.t.I.BAND1.MNAR$model.inf$loglik),
  c(2 * dim(fit.N.I.UNC.MNAR$IM$out)[2] - 2 * fit.N.I.UNC.MNAR$model.inf$loglik, 2 * dim(fit.t.I.UNC.MNAR$IM$out)[2] - 2 * fit.t.I.UNC.MNAR$model.inf$loglik, 2 * dim(fit.N.I.CS.MNAR$IM$out)[2] - 2 * fit.N.I.CS.MNAR$model.inf$loglik, 2 * dim(fit.t.I.CS.MNAR$IM$out)[2] - 2 * fit.t.I.CS.MNAR$model.inf$loglik, 2 * dim(fit.N.I.ARp.MNAR$IM$out)[2] - 2 * fit.N.I.ARp.MNAR$model.inf$loglik, 2 * dim(fit.t.I.ARp.MNAR$IM$out)[2] - 2 * fit.t.I.ARp.MNAR$model.inf$loglik, 2 * dim(fit.N.I.BAND1.MNAR$IM$out)[2] - 2 * fit.N.I.BAND1.MNAR$model.inf$loglik, 2 * dim(fit.t.I.BAND1.MNAR$IM$out)[2] - 2 * fit.t.I.BAND1.MNAR$model.inf$loglik),
  c(dim(fit.N.I.UNC.MNAR$IM$out)[2] * log(N) - 2 * fit.N.I.UNC.MNAR$model.inf$loglik, dim(fit.t.I.UNC.MNAR$IM$out)[2] * log(N) - 2 * fit.t.I.UNC.MNAR$model.inf$loglik, dim(fit.N.I.CS.MNAR$IM$out)[2] * log(N) - 2 * fit.N.I.CS.MNAR$model.inf$loglik, dim(fit.t.I.CS.MNAR$IM$out)[2] * log(N) - 2 * fit.t.I.CS.MNAR$model.inf$loglik, dim(fit.N.I.ARp.MNAR$IM$out)[2] * log(N) - 2 * fit.N.I.ARp.MNAR$model.inf$loglik, dim(fit.t.I.ARp.MNAR$IM$out)[2] * log(N) - 2 * fit.t.I.ARp.MNAR$model.inf$loglik, dim(fit.N.I.BAND1.MNAR$IM$out)[2] * log(N) - 2 * fit.N.I.BAND1.MNAR$model.inf$loglik, dim(fit.t.I.BAND1.MNAR$IM$out)[2] * log(N) - 2 * fit.t.I.BAND1.MNAR$model.inf$loglik),
  c(fit.N.I.UNC.MNAR$iter, fit.t.I.UNC.MNAR$iter, fit.N.I.CS.MNAR$iter, fit.t.I.CS.MNAR$iter, fit.N.I.ARp.MNAR$iter, fit.t.I.ARp.MNAR$iter, fit.N.I.BAND1.MNAR$iter, fit.t.I.BAND1.MNAR$iter),
  c(fit.N.I.UNC.MNAR$model.inf$time / 3600, fit.t.I.UNC.MNAR$model.inf$time / 3600, fit.N.I.CS.MNAR$model.inf$time / 3600, fit.t.I.CS.MNAR$model.inf$time / 3600, fit.N.I.ARp.MNAR$model.inf$time / 3600, fit.t.I.ARp.MNAR$model.inf$time / 3600, fit.N.I.BAND1.MNAR$model.inf$time / 3600, fit.t.I.BAND1.MNAR$model.inf$time / 3600),
  c(
    mean(colMeans(fit.N.I.UNC.MNAR$Taccept.rate)[colMeans(fit.N.I.UNC.MNAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.I.UNC.MNAR$Taccept.rate)[colMeans(fit.t.I.UNC.MNAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.I.CS.MNAR$Taccept.rate)[colMeans(fit.N.I.CS.MNAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.I.CS.MNAR$Taccept.rate)[colMeans(fit.t.I.CS.MNAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.I.ARp.MNAR$Taccept.rate)[colMeans(fit.N.I.ARp.MNAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.I.ARp.MNAR$Taccept.rate)[colMeans(fit.t.I.ARp.MNAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.I.BAND1.MNAR$Taccept.rate)[colMeans(fit.N.I.BAND1.MNAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.I.BAND1.MNAR$Taccept.rate)[colMeans(fit.t.I.BAND1.MNAR$Taccept.rate) != 0])
  )
)

Tab.I <- rbind(
  Tab.I.MCAR,
  Tab.I.MAR,
  Tab.I.MNAR
)


############################################################
## Structure (II)
############################################################

Tab.II.MCAR <- rbind(
  c(dim(fit.N.II.UNC.MCAR$IM$out)[2], dim(fit.t.II.UNC.MCAR$IM$out)[2], dim(fit.N.II.CS.MCAR$IM$out)[2], dim(fit.t.II.CS.MCAR$IM$out)[2], dim(fit.N.II.ARp.MCAR$IM$out)[2], dim(fit.t.II.ARp.MCAR$IM$out)[2], dim(fit.N.II.BAND1.MCAR$IM$out)[2], dim(fit.t.II.BAND1.MCAR$IM$out)[2]),
  c(fit.N.II.UNC.MCAR$model.inf$loglik, fit.t.II.UNC.MCAR$model.inf$loglik, fit.N.II.CS.MCAR$model.inf$loglik, fit.t.II.CS.MCAR$model.inf$loglik, fit.N.II.ARp.MCAR$model.inf$loglik, fit.t.II.ARp.MCAR$model.inf$loglik, fit.N.II.BAND1.MCAR$model.inf$loglik, fit.t.II.BAND1.MCAR$model.inf$loglik),
  c(2 * dim(fit.N.II.UNC.MCAR$IM$out)[2] - 2 * fit.N.II.UNC.MCAR$model.inf$loglik, 2 * dim(fit.t.II.UNC.MCAR$IM$out)[2] - 2 * fit.t.II.UNC.MCAR$model.inf$loglik, 2 * dim(fit.N.II.CS.MCAR$IM$out)[2] - 2 * fit.N.II.CS.MCAR$model.inf$loglik, 2 * dim(fit.t.II.CS.MCAR$IM$out)[2] - 2 * fit.t.II.CS.MCAR$model.inf$loglik, 2 * dim(fit.N.II.ARp.MCAR$IM$out)[2] - 2 * fit.N.II.ARp.MCAR$model.inf$loglik, 2 * dim(fit.t.II.ARp.MCAR$IM$out)[2] - 2 * fit.t.II.ARp.MCAR$model.inf$loglik, 2 * dim(fit.N.II.BAND1.MCAR$IM$out)[2] - 2 * fit.N.II.BAND1.MCAR$model.inf$loglik, 2 * dim(fit.t.II.BAND1.MCAR$IM$out)[2] - 2 * fit.t.II.BAND1.MCAR$model.inf$loglik),
  c(dim(fit.N.II.UNC.MCAR$IM$out)[2] * log(N) - 2 * fit.N.II.UNC.MCAR$model.inf$loglik, dim(fit.t.II.UNC.MCAR$IM$out)[2] * log(N) - 2 * fit.t.II.UNC.MCAR$model.inf$loglik, dim(fit.N.II.CS.MCAR$IM$out)[2] * log(N) - 2 * fit.N.II.CS.MCAR$model.inf$loglik, dim(fit.t.II.CS.MCAR$IM$out)[2] * log(N) - 2 * fit.t.II.CS.MCAR$model.inf$loglik, dim(fit.N.II.ARp.MCAR$IM$out)[2] * log(N) - 2 * fit.N.II.ARp.MCAR$model.inf$loglik, dim(fit.t.II.ARp.MCAR$IM$out)[2] * log(N) - 2 * fit.t.II.ARp.MCAR$model.inf$loglik, dim(fit.N.II.BAND1.MCAR$IM$out)[2] * log(N) - 2 * fit.N.II.BAND1.MCAR$model.inf$loglik, dim(fit.t.II.BAND1.MCAR$IM$out)[2] * log(N) - 2 * fit.t.II.BAND1.MCAR$model.inf$loglik),
  c(fit.N.II.UNC.MCAR$iter, fit.t.II.UNC.MCAR$iter, fit.N.II.CS.MCAR$iter, fit.t.II.CS.MCAR$iter, fit.N.II.ARp.MCAR$iter, fit.t.II.ARp.MCAR$iter, fit.N.II.BAND1.MCAR$iter, fit.t.II.BAND1.MCAR$iter),
  c(fit.N.II.UNC.MCAR$model.inf$time / 3600, fit.t.II.UNC.MCAR$model.inf$time / 3600, fit.N.II.CS.MCAR$model.inf$time / 3600, fit.t.II.CS.MCAR$model.inf$time / 3600, fit.N.II.ARp.MCAR$model.inf$time / 3600, fit.t.II.ARp.MCAR$model.inf$time / 3600, fit.N.II.BAND1.MCAR$model.inf$time / 3600, fit.t.II.BAND1.MCAR$model.inf$time / 3600),
  c(
    mean(colMeans(fit.N.II.UNC.MCAR$Taccept.rate)[colMeans(fit.N.II.UNC.MCAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.II.UNC.MCAR$Taccept.rate)[colMeans(fit.t.II.UNC.MCAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.II.CS.MCAR$Taccept.rate)[colMeans(fit.N.II.CS.MCAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.II.CS.MCAR$Taccept.rate)[colMeans(fit.t.II.CS.MCAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.II.ARp.MCAR$Taccept.rate)[colMeans(fit.N.II.ARp.MCAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.II.ARp.MCAR$Taccept.rate)[colMeans(fit.t.II.ARp.MCAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.II.BAND1.MCAR$Taccept.rate)[colMeans(fit.N.II.BAND1.MCAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.II.BAND1.MCAR$Taccept.rate)[colMeans(fit.t.II.BAND1.MCAR$Taccept.rate) != 0])
  )
)

Tab.II.MAR <- rbind(
  c(dim(fit.N.II.UNC.MAR$IM$out)[2], dim(fit.t.II.UNC.MAR$IM$out)[2], dim(fit.N.II.CS.MAR$IM$out)[2], dim(fit.t.II.CS.MAR$IM$out)[2], dim(fit.N.II.ARp.MAR$IM$out)[2], dim(fit.t.II.ARp.MAR$IM$out)[2], dim(fit.N.II.BAND1.MAR$IM$out)[2], dim(fit.t.II.BAND1.MAR$IM$out)[2]),
  c(fit.N.II.UNC.MAR$model.inf$loglik, fit.t.II.UNC.MAR$model.inf$loglik, fit.N.II.CS.MAR$model.inf$loglik, fit.t.II.CS.MAR$model.inf$loglik, fit.N.II.ARp.MAR$model.inf$loglik, fit.t.II.ARp.MAR$model.inf$loglik, fit.N.II.BAND1.MAR$model.inf$loglik, fit.t.II.BAND1.MAR$model.inf$loglik),
  c(2 * dim(fit.N.II.UNC.MAR$IM$out)[2] - 2 * fit.N.II.UNC.MAR$model.inf$loglik, 2 * dim(fit.t.II.UNC.MAR$IM$out)[2] - 2 * fit.t.II.UNC.MAR$model.inf$loglik, 2 * dim(fit.N.II.CS.MAR$IM$out)[2] - 2 * fit.N.II.CS.MAR$model.inf$loglik, 2 * dim(fit.t.II.CS.MAR$IM$out)[2] - 2 * fit.t.II.CS.MAR$model.inf$loglik, 2 * dim(fit.N.II.ARp.MAR$IM$out)[2] - 2 * fit.N.II.ARp.MAR$model.inf$loglik, 2 * dim(fit.t.II.ARp.MAR$IM$out)[2] - 2 * fit.t.II.ARp.MAR$model.inf$loglik, 2 * dim(fit.N.II.BAND1.MAR$IM$out)[2] - 2 * fit.N.II.BAND1.MAR$model.inf$loglik, 2 * dim(fit.t.II.BAND1.MAR$IM$out)[2] - 2 * fit.t.II.BAND1.MAR$model.inf$loglik),
  c(dim(fit.N.II.UNC.MAR$IM$out)[2] * log(N) - 2 * fit.N.II.UNC.MAR$model.inf$loglik, dim(fit.t.II.UNC.MAR$IM$out)[2] * log(N) - 2 * fit.t.II.UNC.MAR$model.inf$loglik, dim(fit.N.II.CS.MAR$IM$out)[2] * log(N) - 2 * fit.N.II.CS.MAR$model.inf$loglik, dim(fit.t.II.CS.MAR$IM$out)[2] * log(N) - 2 * fit.t.II.CS.MAR$model.inf$loglik, dim(fit.N.II.ARp.MAR$IM$out)[2] * log(N) - 2 * fit.N.II.ARp.MAR$model.inf$loglik, dim(fit.t.II.ARp.MAR$IM$out)[2] * log(N) - 2 * fit.t.II.ARp.MAR$model.inf$loglik, dim(fit.N.II.BAND1.MAR$IM$out)[2] * log(N) - 2 * fit.N.II.BAND1.MAR$model.inf$loglik, dim(fit.t.II.BAND1.MAR$IM$out)[2] * log(N) - 2 * fit.t.II.BAND1.MAR$model.inf$loglik),
  c(fit.N.II.UNC.MAR$iter, fit.t.II.UNC.MAR$iter, fit.N.II.CS.MAR$iter, fit.t.II.CS.MAR$iter, fit.N.II.ARp.MAR$iter, fit.t.II.ARp.MAR$iter, fit.N.II.BAND1.MAR$iter, fit.t.II.BAND1.MAR$iter),
  c(fit.N.II.UNC.MAR$model.inf$time / 3600, fit.t.II.UNC.MAR$model.inf$time / 3600, fit.N.II.CS.MAR$model.inf$time / 3600, fit.t.II.CS.MAR$model.inf$time / 3600, fit.N.II.ARp.MAR$model.inf$time / 3600, fit.t.II.ARp.MAR$model.inf$time / 3600, fit.N.II.BAND1.MAR$model.inf$time / 3600, fit.t.II.BAND1.MAR$model.inf$time / 3600),
  c(
    mean(colMeans(fit.N.II.UNC.MAR$Taccept.rate)[colMeans(fit.N.II.UNC.MAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.II.UNC.MAR$Taccept.rate)[colMeans(fit.t.II.UNC.MAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.II.CS.MAR$Taccept.rate)[colMeans(fit.N.II.CS.MAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.II.CS.MAR$Taccept.rate)[colMeans(fit.t.II.CS.MAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.II.ARp.MAR$Taccept.rate)[colMeans(fit.N.II.ARp.MAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.II.ARp.MAR$Taccept.rate)[colMeans(fit.t.II.ARp.MAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.II.BAND1.MAR$Taccept.rate)[colMeans(fit.N.II.BAND1.MAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.II.BAND1.MAR$Taccept.rate)[colMeans(fit.t.II.BAND1.MAR$Taccept.rate) != 0])
  )
)

Tab.II.MNAR <- rbind(
  c(dim(fit.N.II.UNC.MNAR$IM$out)[2], dim(fit.t.II.UNC.MNAR$IM$out)[2], dim(fit.N.II.CS.MNAR$IM$out)[2], dim(fit.t.II.CS.MNAR$IM$out)[2], dim(fit.N.II.ARp.MNAR$IM$out)[2], dim(fit.t.II.ARp.MNAR$IM$out)[2], dim(fit.N.II.BAND1.MNAR$IM$out)[2], dim(fit.t.II.BAND1.MNAR$IM$out)[2]),
  c(fit.N.II.UNC.MNAR$model.inf$loglik, fit.t.II.UNC.MNAR$model.inf$loglik, fit.N.II.CS.MNAR$model.inf$loglik, fit.t.II.CS.MNAR$model.inf$loglik, fit.N.II.ARp.MNAR$model.inf$loglik, fit.t.II.ARp.MNAR$model.inf$loglik, fit.N.II.BAND1.MNAR$model.inf$loglik, fit.t.II.BAND1.MNAR$model.inf$loglik),
  c(2 * dim(fit.N.II.UNC.MNAR$IM$out)[2] - 2 * fit.N.II.UNC.MNAR$model.inf$loglik, 2 * dim(fit.t.II.UNC.MNAR$IM$out)[2] - 2 * fit.t.II.UNC.MNAR$model.inf$loglik, 2 * dim(fit.N.II.CS.MNAR$IM$out)[2] - 2 * fit.N.II.CS.MNAR$model.inf$loglik, 2 * dim(fit.t.II.CS.MNAR$IM$out)[2] - 2 * fit.t.II.CS.MNAR$model.inf$loglik, 2 * dim(fit.N.II.ARp.MNAR$IM$out)[2] - 2 * fit.N.II.ARp.MNAR$model.inf$loglik, 2 * dim(fit.t.II.ARp.MNAR$IM$out)[2] - 2 * fit.t.II.ARp.MNAR$model.inf$loglik, 2 * dim(fit.N.II.BAND1.MNAR$IM$out)[2] - 2 * fit.N.II.BAND1.MNAR$model.inf$loglik, 2 * dim(fit.t.II.BAND1.MNAR$IM$out)[2] - 2 * fit.t.II.BAND1.MNAR$model.inf$loglik),
  c(dim(fit.N.II.UNC.MNAR$IM$out)[2] * log(N) - 2 * fit.N.II.UNC.MNAR$model.inf$loglik, dim(fit.t.II.UNC.MNAR$IM$out)[2] * log(N) - 2 * fit.t.II.UNC.MNAR$model.inf$loglik, dim(fit.N.II.CS.MNAR$IM$out)[2] * log(N) - 2 * fit.N.II.CS.MNAR$model.inf$loglik, dim(fit.t.II.CS.MNAR$IM$out)[2] * log(N) - 2 * fit.t.II.CS.MNAR$model.inf$loglik, dim(fit.N.II.ARp.MNAR$IM$out)[2] * log(N) - 2 * fit.N.II.ARp.MNAR$model.inf$loglik, dim(fit.t.II.ARp.MNAR$IM$out)[2] * log(N) - 2 * fit.t.II.ARp.MNAR$model.inf$loglik, dim(fit.N.II.BAND1.MNAR$IM$out)[2] * log(N) - 2 * fit.N.II.BAND1.MNAR$model.inf$loglik, dim(fit.t.II.BAND1.MNAR$IM$out)[2] * log(N) - 2 * fit.t.II.BAND1.MNAR$model.inf$loglik),
  c(fit.N.II.UNC.MNAR$iter, fit.t.II.UNC.MNAR$iter, fit.N.II.CS.MNAR$iter, fit.t.II.CS.MNAR$iter, fit.N.II.ARp.MNAR$iter, fit.t.II.ARp.MNAR$iter, fit.N.II.BAND1.MNAR$iter, fit.t.II.BAND1.MNAR$iter),
  c(fit.N.II.UNC.MNAR$model.inf$time / 3600, fit.t.II.UNC.MNAR$model.inf$time / 3600, fit.N.II.CS.MNAR$model.inf$time / 3600, fit.t.II.CS.MNAR$model.inf$time / 3600, fit.N.II.ARp.MNAR$model.inf$time / 3600, fit.t.II.ARp.MNAR$model.inf$time / 3600, fit.N.II.BAND1.MNAR$model.inf$time / 3600, fit.t.II.BAND1.MNAR$model.inf$time / 3600),
  c(
    mean(colMeans(fit.N.II.UNC.MNAR$Taccept.rate)[colMeans(fit.N.II.UNC.MNAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.II.UNC.MNAR$Taccept.rate)[colMeans(fit.t.II.UNC.MNAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.II.CS.MNAR$Taccept.rate)[colMeans(fit.N.II.CS.MNAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.II.CS.MNAR$Taccept.rate)[colMeans(fit.t.II.CS.MNAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.II.ARp.MNAR$Taccept.rate)[colMeans(fit.N.II.ARp.MNAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.II.ARp.MNAR$Taccept.rate)[colMeans(fit.t.II.ARp.MNAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.II.BAND1.MNAR$Taccept.rate)[colMeans(fit.N.II.BAND1.MNAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.II.BAND1.MNAR$Taccept.rate)[colMeans(fit.t.II.BAND1.MNAR$Taccept.rate) != 0])
  )
)

Tab.II <- rbind(
  Tab.II.MCAR,
  Tab.II.MAR,
  Tab.II.MNAR
)


############################################################
## Structure (III)
############################################################

Tab.III.MCAR <- rbind(
  c(dim(fit.N.III.UNC.MCAR$IM$out)[2], dim(fit.t.III.UNC.MCAR$IM$out)[2], dim(fit.N.III.CS.MCAR$IM$out)[2], dim(fit.t.III.CS.MCAR$IM$out)[2], dim(fit.N.III.ARp.MCAR$IM$out)[2], dim(fit.t.III.ARp.MCAR$IM$out)[2], dim(fit.N.III.BAND1.MCAR$IM$out)[2], dim(fit.t.III.BAND1.MCAR$IM$out)[2]),
  c(fit.N.III.UNC.MCAR$model.inf$loglik, fit.t.III.UNC.MCAR$model.inf$loglik, fit.N.III.CS.MCAR$model.inf$loglik, fit.t.III.CS.MCAR$model.inf$loglik, fit.N.III.ARp.MCAR$model.inf$loglik, fit.t.III.ARp.MCAR$model.inf$loglik, fit.N.III.BAND1.MCAR$model.inf$loglik, fit.t.III.BAND1.MCAR$model.inf$loglik),
  c(2 * dim(fit.N.III.UNC.MCAR$IM$out)[2] - 2 * fit.N.III.UNC.MCAR$model.inf$loglik, 2 * dim(fit.t.III.UNC.MCAR$IM$out)[2] - 2 * fit.t.III.UNC.MCAR$model.inf$loglik, 2 * dim(fit.N.III.CS.MCAR$IM$out)[2] - 2 * fit.N.III.CS.MCAR$model.inf$loglik, 2 * dim(fit.t.III.CS.MCAR$IM$out)[2] - 2 * fit.t.III.CS.MCAR$model.inf$loglik, 2 * dim(fit.N.III.ARp.MCAR$IM$out)[2] - 2 * fit.N.III.ARp.MCAR$model.inf$loglik, 2 * dim(fit.t.III.ARp.MCAR$IM$out)[2] - 2 * fit.t.III.ARp.MCAR$model.inf$loglik, 2 * dim(fit.N.III.BAND1.MCAR$IM$out)[2] - 2 * fit.N.III.BAND1.MCAR$model.inf$loglik, 2 * dim(fit.t.III.BAND1.MCAR$IM$out)[2] - 2 * fit.t.III.BAND1.MCAR$model.inf$loglik),
  c(dim(fit.N.III.UNC.MCAR$IM$out)[2] * log(N) - 2 * fit.N.III.UNC.MCAR$model.inf$loglik, dim(fit.t.III.UNC.MCAR$IM$out)[2] * log(N) - 2 * fit.t.III.UNC.MCAR$model.inf$loglik, dim(fit.N.III.CS.MCAR$IM$out)[2] * log(N) - 2 * fit.N.III.CS.MCAR$model.inf$loglik, dim(fit.t.III.CS.MCAR$IM$out)[2] * log(N) - 2 * fit.t.III.CS.MCAR$model.inf$loglik, dim(fit.N.III.ARp.MCAR$IM$out)[2] * log(N) - 2 * fit.N.III.ARp.MCAR$model.inf$loglik, dim(fit.t.III.ARp.MCAR$IM$out)[2] * log(N) - 2 * fit.t.III.ARp.MCAR$model.inf$loglik, dim(fit.N.III.BAND1.MCAR$IM$out)[2] * log(N) - 2 * fit.N.III.BAND1.MCAR$model.inf$loglik, dim(fit.t.III.BAND1.MCAR$IM$out)[2] * log(N) - 2 * fit.t.III.BAND1.MCAR$model.inf$loglik),
  c(fit.N.III.UNC.MCAR$iter, fit.t.III.UNC.MCAR$iter, fit.N.III.CS.MCAR$iter, fit.t.III.CS.MCAR$iter, fit.N.III.ARp.MCAR$iter, fit.t.III.ARp.MCAR$iter, fit.N.III.BAND1.MCAR$iter, fit.t.III.BAND1.MCAR$iter),
  c(fit.N.III.UNC.MCAR$model.inf$time / 3600, fit.t.III.UNC.MCAR$model.inf$time / 3600, fit.N.III.CS.MCAR$model.inf$time / 3600, fit.t.III.CS.MCAR$model.inf$time / 3600, fit.N.III.ARp.MCAR$model.inf$time / 3600, fit.t.III.ARp.MCAR$model.inf$time / 3600, fit.N.III.BAND1.MCAR$model.inf$time / 3600, fit.t.III.BAND1.MCAR$model.inf$time / 3600),
  c(
    mean(colMeans(fit.N.III.UNC.MCAR$Taccept.rate)[colMeans(fit.N.III.UNC.MCAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.III.UNC.MCAR$Taccept.rate)[colMeans(fit.t.III.UNC.MCAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.III.CS.MCAR$Taccept.rate)[colMeans(fit.N.III.CS.MCAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.III.CS.MCAR$Taccept.rate)[colMeans(fit.t.III.CS.MCAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.III.ARp.MCAR$Taccept.rate)[colMeans(fit.N.III.ARp.MCAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.III.ARp.MCAR$Taccept.rate)[colMeans(fit.t.III.ARp.MCAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.III.BAND1.MCAR$Taccept.rate)[colMeans(fit.N.III.BAND1.MCAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.III.BAND1.MCAR$Taccept.rate)[colMeans(fit.t.III.BAND1.MCAR$Taccept.rate) != 0])
  )
)

Tab.III.MAR <- rbind(
  c(dim(fit.N.III.UNC.MAR$IM$out)[2], dim(fit.t.III.UNC.MAR$IM$out)[2], dim(fit.N.III.CS.MAR$IM$out)[2], dim(fit.t.III.CS.MAR$IM$out)[2], dim(fit.N.III.ARp.MAR$IM$out)[2], dim(fit.t.III.ARp.MAR$IM$out)[2], dim(fit.N.III.BAND1.MAR$IM$out)[2], dim(fit.t.III.BAND1.MAR$IM$out)[2]),
  c(fit.N.III.UNC.MAR$model.inf$loglik, fit.t.III.UNC.MAR$model.inf$loglik, fit.N.III.CS.MAR$model.inf$loglik, fit.t.III.CS.MAR$model.inf$loglik, fit.N.III.ARp.MAR$model.inf$loglik, fit.t.III.ARp.MAR$model.inf$loglik, fit.N.III.BAND1.MAR$model.inf$loglik, fit.t.III.BAND1.MAR$model.inf$loglik),
  c(2 * dim(fit.N.III.UNC.MAR$IM$out)[2] - 2 * fit.N.III.UNC.MAR$model.inf$loglik, 2 * dim(fit.t.III.UNC.MAR$IM$out)[2] - 2 * fit.t.III.UNC.MAR$model.inf$loglik, 2 * dim(fit.N.III.CS.MAR$IM$out)[2] - 2 * fit.N.III.CS.MAR$model.inf$loglik, 2 * dim(fit.t.III.CS.MAR$IM$out)[2] - 2 * fit.t.III.CS.MAR$model.inf$loglik, 2 * dim(fit.N.III.ARp.MAR$IM$out)[2] - 2 * fit.N.III.ARp.MAR$model.inf$loglik, 2 * dim(fit.t.III.ARp.MAR$IM$out)[2] - 2 * fit.t.III.ARp.MAR$model.inf$loglik, 2 * dim(fit.N.III.BAND1.MAR$IM$out)[2] - 2 * fit.N.III.BAND1.MAR$model.inf$loglik, 2 * dim(fit.t.III.BAND1.MAR$IM$out)[2] - 2 * fit.t.III.BAND1.MAR$model.inf$loglik),
  c(dim(fit.N.III.UNC.MAR$IM$out)[2] * log(N) - 2 * fit.N.III.UNC.MAR$model.inf$loglik, dim(fit.t.III.UNC.MAR$IM$out)[2] * log(N) - 2 * fit.t.III.UNC.MAR$model.inf$loglik, dim(fit.N.III.CS.MAR$IM$out)[2] * log(N) - 2 * fit.N.III.CS.MAR$model.inf$loglik, dim(fit.t.III.CS.MAR$IM$out)[2] * log(N) - 2 * fit.t.III.CS.MAR$model.inf$loglik, dim(fit.N.III.ARp.MAR$IM$out)[2] * log(N) - 2 * fit.N.III.ARp.MAR$model.inf$loglik, dim(fit.t.III.ARp.MAR$IM$out)[2] * log(N) - 2 * fit.t.III.ARp.MAR$model.inf$loglik, dim(fit.N.III.BAND1.MAR$IM$out)[2] * log(N) - 2 * fit.N.III.BAND1.MAR$model.inf$loglik, dim(fit.t.III.BAND1.MAR$IM$out)[2] * log(N) - 2 * fit.t.III.BAND1.MAR$model.inf$loglik),
  c(fit.N.III.UNC.MAR$iter, fit.t.III.UNC.MAR$iter, fit.N.III.CS.MAR$iter, fit.t.III.CS.MAR$iter, fit.N.III.ARp.MAR$iter, fit.t.III.ARp.MAR$iter, fit.N.III.BAND1.MAR$iter, fit.t.III.BAND1.MAR$iter),
  c(fit.N.III.UNC.MAR$model.inf$time / 3600, fit.t.III.UNC.MAR$model.inf$time / 3600, fit.N.III.CS.MAR$model.inf$time / 3600, fit.t.III.CS.MAR$model.inf$time / 3600, fit.N.III.ARp.MAR$model.inf$time / 3600, fit.t.III.ARp.MAR$model.inf$time / 3600, fit.N.III.BAND1.MAR$model.inf$time / 3600, fit.t.III.BAND1.MAR$model.inf$time / 3600),
  c(
    mean(colMeans(fit.N.III.UNC.MAR$Taccept.rate)[colMeans(fit.N.III.UNC.MAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.III.UNC.MAR$Taccept.rate)[colMeans(fit.t.III.UNC.MAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.III.CS.MAR$Taccept.rate)[colMeans(fit.N.III.CS.MAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.III.CS.MAR$Taccept.rate)[colMeans(fit.t.III.CS.MAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.III.ARp.MAR$Taccept.rate)[colMeans(fit.N.III.ARp.MAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.III.ARp.MAR$Taccept.rate)[colMeans(fit.t.III.ARp.MAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.III.BAND1.MAR$Taccept.rate)[colMeans(fit.N.III.BAND1.MAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.III.BAND1.MAR$Taccept.rate)[colMeans(fit.t.III.BAND1.MAR$Taccept.rate) != 0])
  )
)

Tab.III.MNAR <- rbind(
  c(dim(fit.N.III.UNC.MNAR$IM$out)[2], dim(fit.t.III.UNC.MNAR$IM$out)[2], dim(fit.N.III.CS.MNAR$IM$out)[2], dim(fit.t.III.CS.MNAR$IM$out)[2], dim(fit.N.III.ARp.MNAR$IM$out)[2], dim(fit.t.III.ARp.MNAR$IM$out)[2], dim(fit.N.III.BAND1.MNAR$IM$out)[2], dim(fit.t.III.BAND1.MNAR$IM$out)[2]),
  c(fit.N.III.UNC.MNAR$model.inf$loglik, fit.t.III.UNC.MNAR$model.inf$loglik, fit.N.III.CS.MNAR$model.inf$loglik, fit.t.III.CS.MNAR$model.inf$loglik, fit.N.III.ARp.MNAR$model.inf$loglik, fit.t.III.ARp.MNAR$model.inf$loglik, fit.N.III.BAND1.MNAR$model.inf$loglik, fit.t.III.BAND1.MNAR$model.inf$loglik),
  c(2 * dim(fit.N.III.UNC.MNAR$IM$out)[2] - 2 * fit.N.III.UNC.MNAR$model.inf$loglik, 2 * dim(fit.t.III.UNC.MNAR$IM$out)[2] - 2 * fit.t.III.UNC.MNAR$model.inf$loglik, 2 * dim(fit.N.III.CS.MNAR$IM$out)[2] - 2 * fit.N.III.CS.MNAR$model.inf$loglik, 2 * dim(fit.t.III.CS.MNAR$IM$out)[2] - 2 * fit.t.III.CS.MNAR$model.inf$loglik, 2 * dim(fit.N.III.ARp.MNAR$IM$out)[2] - 2 * fit.N.III.ARp.MNAR$model.inf$loglik, 2 * dim(fit.t.III.ARp.MNAR$IM$out)[2] - 2 * fit.t.III.ARp.MNAR$model.inf$loglik, 2 * dim(fit.N.III.BAND1.MNAR$IM$out)[2] - 2 * fit.N.III.BAND1.MNAR$model.inf$loglik, 2 * dim(fit.t.III.BAND1.MNAR$IM$out)[2] - 2 * fit.t.III.BAND1.MNAR$model.inf$loglik),
  c(dim(fit.N.III.UNC.MNAR$IM$out)[2] * log(N) - 2 * fit.N.III.UNC.MNAR$model.inf$loglik, dim(fit.t.III.UNC.MNAR$IM$out)[2] * log(N) - 2 * fit.t.III.UNC.MNAR$model.inf$loglik, dim(fit.N.III.CS.MNAR$IM$out)[2] * log(N) - 2 * fit.N.III.CS.MNAR$model.inf$loglik, dim(fit.t.III.CS.MNAR$IM$out)[2] * log(N) - 2 * fit.t.III.CS.MNAR$model.inf$loglik, dim(fit.N.III.ARp.MNAR$IM$out)[2] * log(N) - 2 * fit.N.III.ARp.MNAR$model.inf$loglik, dim(fit.t.III.ARp.MNAR$IM$out)[2] * log(N) - 2 * fit.t.III.ARp.MNAR$model.inf$loglik, dim(fit.N.III.BAND1.MNAR$IM$out)[2] * log(N) - 2 * fit.N.III.BAND1.MNAR$model.inf$loglik, dim(fit.t.III.BAND1.MNAR$IM$out)[2] * log(N) - 2 * fit.t.III.BAND1.MNAR$model.inf$loglik),
  c(fit.N.III.UNC.MNAR$iter, fit.t.III.UNC.MNAR$iter, fit.N.III.CS.MNAR$iter, fit.t.III.CS.MNAR$iter, fit.N.III.ARp.MNAR$iter, fit.t.III.ARp.MNAR$iter, fit.N.III.BAND1.MNAR$iter, fit.t.III.BAND1.MNAR$iter),
  c(fit.N.III.UNC.MNAR$model.inf$time / 3600, fit.t.III.UNC.MNAR$model.inf$time / 3600, fit.N.III.CS.MNAR$model.inf$time / 3600, fit.t.III.CS.MNAR$model.inf$time / 3600, fit.N.III.ARp.MNAR$model.inf$time / 3600, fit.t.III.ARp.MNAR$model.inf$time / 3600, fit.N.III.BAND1.MNAR$model.inf$time / 3600, fit.t.III.BAND1.MNAR$model.inf$time / 3600),
  c(
    mean(colMeans(fit.N.III.UNC.MNAR$Taccept.rate)[colMeans(fit.N.III.UNC.MNAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.III.UNC.MNAR$Taccept.rate)[colMeans(fit.t.III.UNC.MNAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.III.CS.MNAR$Taccept.rate)[colMeans(fit.N.III.CS.MNAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.III.CS.MNAR$Taccept.rate)[colMeans(fit.t.III.CS.MNAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.III.ARp.MNAR$Taccept.rate)[colMeans(fit.N.III.ARp.MNAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.III.ARp.MNAR$Taccept.rate)[colMeans(fit.t.III.ARp.MNAR$Taccept.rate) != 0]),
    mean(colMeans(fit.N.III.BAND1.MNAR$Taccept.rate)[colMeans(fit.N.III.BAND1.MNAR$Taccept.rate) != 0]),
    mean(colMeans(fit.t.III.BAND1.MNAR$Taccept.rate)[colMeans(fit.t.III.BAND1.MNAR$Taccept.rate) != 0])
  )
)

Tab.III <- rbind(
  Tab.III.MCAR,
  Tab.III.MAR,
  Tab.III.MNAR
)



############################################################
## Stack Structures I, II, and III from top to bottom
############################################################

sum.table <- rbind(
  Tab.I,
  Tab.II,
  Tab.III
)

sum.table <- round(sum.table, 3)

colnames(sum.table) <- c(
  "UNC_NLME", "UNC_tNLME",
  "CS_NLME",  "CS_tNLME",
  "AR1_NLME", "AR1_tNLME",
  "MA1_NLME", "MA1_tNLME"
)

structure <- c(
  "I", rep("", 20),
  "II", rep("", 20),
  "III", rep("", 20)
)

mechanism.one <- rep(
  c("MCAR", "", "", "", "", "", "",
    "MAR",  "", "", "", "", "", "",
    "MNAR", "", "", "", "", "", ""),
  1
)

mechanism <- rep(mechanism.one, 3)

criterion <- rep(
  c("m", "loglik", "AIC", "BIC", "iter", "time_hours", "acceptance_rate"),
  9
)

Table1 <- data.frame(
  structure = structure,
  mechanism = mechanism,
  criterion = criterion,
  sum.table,
  check.names = FALSE
)

Table1

################################################################################
# 4. Produce Table 1 as reported in the manuscript
################################################################################

Table1.paper <- Table1[Table1$criterion %in% c("m", "loglik", "AIC", "BIC"), ]

Table1.paper$structure <- c("Scenario I", rep("", 11), "Scenario II", rep("", 11), "Scenario III", rep("", 11))
Table1.paper$mechanism <- rep(c("MCAR", "", "", "", "MAR", "", "", "", "MNAR", "", "", ""), 3)
Table1.paper$criterion <- rep(c("m", "l_max", "AIC", "BIC"), 9)

################################################################################
# 5. Format values as reported in Table 1
################################################################################

Table1.output <- Table1.paper

Table1.output[Table1.output$criterion == "m", 4:11] <- format(round(Table1.paper[Table1.paper$criterion == "m", 4:11], 0), nsmall = 0, trim = TRUE)
Table1.output[Table1.output$criterion == "l_max", 4:11] <- format(round(Table1.paper[Table1.paper$criterion == "l_max", 4:11], 3), nsmall = 3, trim = TRUE)
Table1.output[Table1.output$criterion == "AIC", 4:11] <- format(round(Table1.paper[Table1.paper$criterion == "AIC", 4:11], 3), nsmall = 3, trim = TRUE)
Table1.output[Table1.output$criterion == "BIC", 4:11] <- format(round(Table1.paper[Table1.paper$criterion == "BIC", 4:11], 3), nsmall = 3, trim = TRUE)

################################################################################
# 6. Save Table 1
################################################################################

RESULT_DIR <- file.path(PATH, "Result")

write.csv(Table1.output, file = file.path(RESULT_DIR, "Table1.csv"), row.names = FALSE, quote = FALSE, na = "")


