################################################################################
#
#   Filename    :    TabD1.R
#   Project     :    Biometrics article "Modeling HIV Viral Dynamics Using a
#                    Nonlinear Mixed-Effects Framework for Heavy-Tailed Data
#                    with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#   Date        :    18.08.2026
#   Purpose     :    produce Table D.1 for the ACTG 398 data by summarizing the
#                    model-selection results from the Scenario (III) NLME and
#                    tNLME models fitted under the alternative nls-based
#                    initialization
#
#   Input data files  :  Data_and_Code/Data/fit_III_result_nls.RData
#
#   Output data files :  Data_and_Code/Result/TableD1.csv
#
#   R Version   :    R-4.6.0
#   Required R packages : none
#
################################################################################

PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
load(paste0(PATH, "/Data/fit_III_result_nls.RData"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

N <- length(unique(Data$Subject))

if (N != 481) stop("Unexpected number of ACTG 398 subjects.")

################################################################################
# 1. MCAR
################################################################################

TabD1.MCAR <- rbind(
  c(dim(fit.N.III.UNC.MCAR.nls$IM$out)[2], dim(fit.t.III.UNC.MCAR.nls$IM$out)[2], dim(fit.N.III.CS.MCAR.nls$IM$out)[2], dim(fit.t.III.CS.MCAR.nls$IM$out)[2], dim(fit.N.III.ARp.MCAR.nls$IM$out)[2], dim(fit.t.III.ARp.MCAR.nls$IM$out)[2], dim(fit.N.III.BAND1.MCAR.nls$IM$out)[2], dim(fit.t.III.BAND1.MCAR.nls$IM$out)[2]),
  c(fit.N.III.UNC.MCAR.nls$model.inf$loglik, fit.t.III.UNC.MCAR.nls$model.inf$loglik, fit.N.III.CS.MCAR.nls$model.inf$loglik, fit.t.III.CS.MCAR.nls$model.inf$loglik, fit.N.III.ARp.MCAR.nls$model.inf$loglik, fit.t.III.ARp.MCAR.nls$model.inf$loglik, fit.N.III.BAND1.MCAR.nls$model.inf$loglik, fit.t.III.BAND1.MCAR.nls$model.inf$loglik)
)

TabD1.MCAR <- rbind(TabD1.MCAR, 2 * TabD1.MCAR[1, ] - 2 * TabD1.MCAR[2, ], TabD1.MCAR[1, ] * log(N) - 2 * TabD1.MCAR[2, ])

################################################################################
# 2. MAR
################################################################################

TabD1.MAR <- rbind(
  c(dim(fit.N.III.UNC.MAR.nls$IM$out)[2], dim(fit.t.III.UNC.MAR.nls$IM$out)[2], dim(fit.N.III.CS.MAR.nls$IM$out)[2], dim(fit.t.III.CS.MAR.nls$IM$out)[2], dim(fit.N.III.ARp.MAR.nls$IM$out)[2], dim(fit.t.III.ARp.MAR.nls$IM$out)[2], dim(fit.N.III.BAND1.MAR.nls$IM$out)[2], dim(fit.t.III.BAND1.MAR.nls$IM$out)[2]),
  c(fit.N.III.UNC.MAR.nls$model.inf$loglik, fit.t.III.UNC.MAR.nls$model.inf$loglik, fit.N.III.CS.MAR.nls$model.inf$loglik, fit.t.III.CS.MAR.nls$model.inf$loglik, fit.N.III.ARp.MAR.nls$model.inf$loglik, fit.t.III.ARp.MAR.nls$model.inf$loglik, fit.N.III.BAND1.MAR.nls$model.inf$loglik, fit.t.III.BAND1.MAR.nls$model.inf$loglik)
)

TabD1.MAR <- rbind(TabD1.MAR, 2 * TabD1.MAR[1, ] - 2 * TabD1.MAR[2, ], TabD1.MAR[1, ] * log(N) - 2 * TabD1.MAR[2, ])

################################################################################
# 3. MNAR
################################################################################

TabD1.MNAR <- rbind(
  c(dim(fit.N.III.UNC.MNAR.nls$IM$out)[2], dim(fit.t.III.UNC.MNAR.nls$IM$out)[2], dim(fit.N.III.CS.MNAR.nls$IM$out)[2], dim(fit.t.III.CS.MNAR.nls$IM$out)[2], dim(fit.N.III.ARp.MNAR.nls$IM$out)[2], dim(fit.t.III.ARp.MNAR.nls$IM$out)[2], dim(fit.N.III.BAND1.MNAR.nls$IM$out)[2], dim(fit.t.III.BAND1.MNAR.nls$IM$out)[2]),
  c(fit.N.III.UNC.MNAR.nls$model.inf$loglik, fit.t.III.UNC.MNAR.nls$model.inf$loglik, fit.N.III.CS.MNAR.nls$model.inf$loglik, fit.t.III.CS.MNAR.nls$model.inf$loglik, fit.N.III.ARp.MNAR.nls$model.inf$loglik, fit.t.III.ARp.MNAR.nls$model.inf$loglik, fit.N.III.BAND1.MNAR.nls$model.inf$loglik, fit.t.III.BAND1.MNAR.nls$model.inf$loglik)
)

TabD1.MNAR <- rbind(TabD1.MNAR, 2 * TabD1.MNAR[1, ] - 2 * TabD1.MNAR[2, ], TabD1.MNAR[1, ] * log(N) - 2 * TabD1.MNAR[2, ])

################################################################################
# 4. Construct Table D.1
################################################################################

sum.table <- round(rbind(TabD1.MCAR, TabD1.MAR, TabD1.MNAR), 3)

TableD1 <- data.frame(
  Mechanism = c("MCAR", "", "", "", "MAR", "", "", "", "MNAR", "", "", ""),
  Criterion = rep(c("m", "ell_max", "AIC", "BIC"), 3),
  sum.table,
  check.names = FALSE
)

colnames(TableD1) <- c("Mechanism", "Criterion", "UNC_NLME", "UNC_tNLME", "CS_NLME", "CS_tNLME", "AR1_NLME", "AR1_tNLME", "MA1_NLME", "MA1_tNLME")

################################################################################
# 5. Save Table D.1
################################################################################
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
write.csv(TableD1, paste0(PATH, "/Result/TableD1.csv"), row.names = FALSE, na = "")
