################################################################################
#
#   Filename    :    Tab2.R
#   Project     :    "Robust HIV Viral Dynamics: A Nonlinear Mixed-Effects Framework for
#                    Heavy-Tailed Data with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#   Date        :    14.03.2026
#   Purpose     :    produce Table 2 for AIDS
#
#   Input data files  :  Data_and_Code/Data/fit_result.RData
#   Output data files :  Data_and_Code/results/Table2.csv
#
#   R Version   :    R-4.3.1
#   Required R packages : None
#
################################################################################
load(paste0(PATH, "/Data/fit_result.RData"))
PATH <- getwd()

MCAR.r <- round(t(cbind(rbind(est.MCAR.ARp$IM$out, abs(est.MCAR.ARp$IM$out[1, ] / est.MCAR.ARp$IM$out[2, ])), NA, NA, NA, NA)), 3)
MAR.r <- round(t(cbind(rbind(est.MAR.ARp$IM$out, abs(est.MAR.ARp$IM$out[1, ] / est.MAR.ARp$IM$out[2, ])), NA)), 3)
MNAR.r <- round(t(cbind(rbind(est.MNAR.ARp$IM$out, abs(est.MNAR.ARp$IM$out[1, ] / est.MNAR.ARp$IM$out[2, ])))), 3)

sum.table <- cbind(MCAR.r, MAR.r, MNAR.r)

colnames(sum.table) <- c("Est", "SE", "Est/SE", "Est", "SE", "Est/SE", "Est", "SE", "Est/SE")
row.names(sum.table) <- c("beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "beta8", "d11", "d12", "d22", "sigma^2", "phi", "nu", "alpha00", "alpha01", "alpha02", "alpha1", "alpha2")


write.csv(sum.table, paste0(PATH, "/Result/Table2.csv"))
