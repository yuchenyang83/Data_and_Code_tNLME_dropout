################################################################################
#
#   Filename    :    TabD2.R
#   Project     :    Biometrics article "Modeling HIV Viral Dynamics Using a
#                    Nonlinear Mixed-Effects Framework for Heavy-Tailed Data
#                    with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#   Date        :    18.08.2026
#   Purpose     :    produce Table D.2 for the ACTG 398 data by reporting the ML
#                    estimates, standard errors, and absolute estimate-to-SE
#                    ratios from the Scenario (III) AR(1) tNLME models under
#                    MCAR, MAR, and MNAR using the alternative nls-based
#                    initialization
#
#   Input data files  :  Data_and_Code/Data/fit_III_result_nls.RData
#
#   Output data files :  Data_and_Code/Result/TableD2.csv
#
#   R Version   :    R-4.6.0
#   Required R packages : none
#
################################################################################

PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
load(paste0(PATH, "/Data/fit_III_result_nls.RData"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

################################################################################
# 1. Extract estimates, SEs, and absolute estimate-to-SE ratios
################################################################################

if (ncol(fit.t.III.ARp.MCAR.nls$IM$out) != 22) stop("Unexpected number of parameters in fit.t.III.ARp.MCAR.nls$IM$out.")
if (ncol(fit.t.III.ARp.MAR.nls$IM$out) != 25) stop("Unexpected number of parameters in fit.t.III.ARp.MAR.nls$IM$out.")
if (ncol(fit.t.III.ARp.MNAR.nls$IM$out) != 26) stop("Unexpected number of parameters in fit.t.III.ARp.MNAR.nls$IM$out.")

MCAR.r <- round(t(rbind(fit.t.III.ARp.MCAR.nls$IM$out, abs(fit.t.III.ARp.MCAR.nls$IM$out[1, ] / fit.t.III.ARp.MCAR.nls$IM$out[2, ]))), 3)
MAR.r <- round(t(rbind(fit.t.III.ARp.MAR.nls$IM$out, abs(fit.t.III.ARp.MAR.nls$IM$out[1, ] / fit.t.III.ARp.MAR.nls$IM$out[2, ]))), 3)
MNAR.r <- round(t(rbind(fit.t.III.ARp.MNAR.nls$IM$out, abs(fit.t.III.ARp.MNAR.nls$IM$out[1, ] / fit.t.III.ARp.MNAR.nls$IM$out[2, ]))), 3)

MCAR.r <- rbind(MCAR.r, matrix(NA, nrow = 4, ncol = 3))
MAR.r <- rbind(MAR.r, matrix(NA, nrow = 1, ncol = 3))

################################################################################
# 2. Construct Table D.2 in the same parameter order as the manuscript
################################################################################

sum.table <- cbind(MCAR.r, MAR.r, MNAR.r)

colnames(sum.table) <- c("MCAR_Est", "MCAR_SE", "MCAR_abs_Est_SE", "MAR_Est", "MAR_SE", "MAR_abs_Est_SE", "MNAR_Est", "MNAR_SE", "MNAR_abs_Est_SE")

Parameter <- c(
  "beta1 (intercept)",
  "beta2 (week)",
  "beta3 (NNRTI)",
  "beta4 (Trtarm)",
  "beta5 (intercept)",
  "beta6 (week)",
  "beta7 (NNRTI)",
  "beta8 (Trtarm)",
  "d11",
  "d12",
  "d22",
  "d13",
  "d23",
  "d33",
  "d14",
  "d24",
  "d34",
  "d44",
  "sigma^2",
  "phi",
  "nu",
  "alpha00 (intercept)",
  "alpha01 (NNRTI)",
  "alpha02 (Trtarm)",
  "alpha1 (y_i,j-1)",
  "alpha2 (y_ij)"
)

TableD2 <- data.frame(Parameter = Parameter, sum.table, check.names = FALSE)

################################################################################
# 3. Save Table D.2
################################################################################

write.csv(TableD2, paste0(PATH, "/Result/TableD2.csv"), row.names = FALSE, na = "")


