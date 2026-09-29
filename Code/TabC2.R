################################################################################
#
#   Filename    :    TabC2.R
#   Project     :    Biometrics article "Modeling HIV Viral Dynamics Using a
#                    Nonlinear Mixed-Effects Framework for Heavy-Tailed Data
#                    with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#   Date        :    18.09.2026
#   Purpose     :    produce Table C.2 for the ACTG 398 data by reporting the ML
#                    estimates, standard errors, and absolute estimate-to-SE
#                    ratios from the Scenario (III) AR(1) tNLME models under
#                    MCAR, MAR, and MNAR using the alternative nls-based
#                    initialization, after reparameterizing Trtarm from the old
#                    coding (therapy = 0, placebo = 1) to the new coding
#                    (therapy = 1, placebo = 0)
#
#   Input data files  :  Data_and_Code/Data/fit_III_result_nls.RData
#
#   Output data files :  Data_and_Code/Result/TableC2.csv
#
#   R Version   :    R-4.6.0
#   Required R packages : none
#
################################################################################

PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
load(file.path(PATH, "Data", "fit_III_result_nls.RData"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

################################################################################
# 1. Exact reparameterization for the revised Trtarm coding
################################################################################

reparameterize_IM <- function(fit, mechanism) {
  im.out <- as.matrix(fit$IM$out)
  p <- ncol(im.out)

  if (mechanism == "MCAR" && p != 22) stop("Unexpected number of parameters for MCAR.")
  if (mechanism == "MAR" && p != 25) stop("Unexpected number of parameters for MAR.")
  if (mechanism == "MNAR" && p != 26) stop("Unexpected number of parameters for MNAR.")

  if (!is.null(fit$IM$V.theta)) {
    V.old <- as.matrix(fit$IM$V.theta)
  } else if (!is.null(fit$IM$I.theta)) {
    V.old <- solve(as.matrix(fit$IM$I.theta))
  } else {
    stop("Neither fit$IM$V.theta nor fit$IM$I.theta is available.")
  }

  if (nrow(V.old) != p || ncol(V.old) != p) stop("Dimension mismatch between fit$IM$out and the covariance matrix.")

  se.old <- as.numeric(im.out[2, ])
  if (!isTRUE(all.equal(se.old, sqrt(diag(V.old)), tolerance = 1e-5, check.attributes = FALSE))) {
    warning("The SE row in fit$IM$out is not numerically identical to sqrt(diag(V.theta)); transformed SEs are computed from V.theta.")
  }

  T <- diag(p)

  T[2, ] <- 0
  T[2, 2] <- 1
  T[2, 4] <- 1

  T[4, ] <- 0
  T[4, 4] <- -1

  T[6, ] <- 0
  T[6, 6] <- 1
  T[6, 8] <- 1

  T[8, ] <- 0
  T[8, 8] <- -1

  if (mechanism %in% c("MAR", "MNAR")) {
    T[22, 24] <- 1
    T[24, ] <- 0
    T[24, 24] <- -1
  }

  est.old <- as.numeric(im.out[1, ])
  est.new <- drop(T %*% est.old)
  V.new <- T %*% V.old %*% t(T)
  se.new <- sqrt(diag(V.new))

  out.new <- im.out
  out.new[1, ] <- est.new
  out.new[2, ] <- se.new

  out.new
}

MCAR.out <- reparameterize_IM(fit.t.III.ARp.MCAR.nls, "MCAR")
MAR.out <- reparameterize_IM(fit.t.III.ARp.MAR.nls, "MAR")
MNAR.out <- reparameterize_IM(fit.t.III.ARp.MNAR.nls, "MNAR")

################################################################################
# 2. Extract estimates, SEs, and absolute estimate-to-SE ratios
################################################################################

MCAR.r <- round(t(rbind(MCAR.out, abs(MCAR.out[1, ] / MCAR.out[2, ]))), 3)
MAR.r <- round(t(rbind(MAR.out, abs(MAR.out[1, ] / MAR.out[2, ]))), 3)
MNAR.r <- round(t(rbind(MNAR.out, abs(MNAR.out[1, ] / MNAR.out[2, ]))), 3)

MCAR.r <- rbind(MCAR.r, matrix(NA, nrow = 4, ncol = 3))
MAR.r <- rbind(MAR.r, matrix(NA, nrow = 1, ncol = 3))

################################################################################
# 3. Construct Table C.2 in the same parameter order as the manuscript
################################################################################

sum.table <- cbind(MCAR.r, MAR.r, MNAR.r)

colnames(sum.table) <- c("MCAR_Est", "MCAR_SE", "MCAR_abs_Est_SE",
                         "MAR_Est", "MAR_SE", "MAR_abs_Est_SE",
                         "MNAR_Est", "MNAR_SE", "MNAR_abs_Est_SE")

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

TableC2 <- data.frame(Parameter = Parameter, sum.table, check.names = FALSE)

################################################################################
# 4. Save Table C.2
################################################################################

write.csv(TableC2, file.path(PATH, "Result", "TableC2.csv"), row.names = FALSE, na = "")
