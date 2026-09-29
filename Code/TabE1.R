################################################################################
#
#   Filename    :    TabE1.R
#   Project     :    Biometrics article "Modeling HIV Viral Dynamics Using a
#                    Nonlinear Mixed-Effects Framework for Heavy-Tailed Data
#                    with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#   Date        :    18.09.2026
#   Purpose     :    produce Table E.1 for the Metropolis-Hastings tuning
#                    sensitivity analysis by extracting the number of model
#                    parameters, approximated observed-data log-likelihood, and
#                    acceptance rate from the Scenario (III) AR(1) tNLME fits,
#                    and then calculating AIC and BIC
#
#   Input data files  :  Data_and_Code/Data/sensitivity_SAEM/
#                         fit.t.III.ARp.MCAR_k5_c1.8.RData, ..., 
#                         fit.t.III.ARp.MNAR_k15_c3.RData
#
#   Intermediate file :  Data_and_Code/Data/sensitivity_SAEM/TableE1_raw.txt
#
#   Output data files :  Data_and_Code/Result/TableE1.csv
#
#   R Version   :    R-4.6.0
#   Required R packages : none
#
################################################################################

PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
PATH1 <- paste0(PATH, "/Data/sensitivity_SAEM")

options(digits = 16)

TableE1.raw <- data.frame(K = numeric(0), Mechanism = character(0), c = numeric(0), m = numeric(0), loglik = numeric(0), acceptance_rate = numeric(0), stringsAsFactors = FALSE)

################################################################################
# 1. K = 5, MCAR
################################################################################

load(paste0(PATH1, "/fit.t.III.ARp.MCAR_k5_c1.8.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 5, Mechanism = "MCAR", c = 1.8, m = dim(fit.t.III.ARp.MCAR$IM$out)[2], loglik = fit.t.III.ARp.MCAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MCAR$Taccept.rate)[colMeans(fit.t.III.ARp.MCAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MCAR)

load(paste0(PATH1, "/fit.t.III.ARp.MCAR_k5_c2.4.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 5, Mechanism = "MCAR", c = 2.4, m = dim(fit.t.III.ARp.MCAR$IM$out)[2], loglik = fit.t.III.ARp.MCAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MCAR$Taccept.rate)[colMeans(fit.t.III.ARp.MCAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MCAR)

load(paste0(PATH1, "/fit.t.III.ARp.MCAR_k5_c3.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 5, Mechanism = "MCAR", c = 3.0, m = dim(fit.t.III.ARp.MCAR$IM$out)[2], loglik = fit.t.III.ARp.MCAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MCAR$Taccept.rate)[colMeans(fit.t.III.ARp.MCAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MCAR)

################################################################################
# 2. K = 5, MAR
################################################################################

load(paste0(PATH1, "/fit.t.III.ARp.MAR_k5_c1.8.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 5, Mechanism = "MAR", c = 1.8, m = dim(fit.t.III.ARp.MAR$IM$out)[2], loglik = fit.t.III.ARp.MAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MAR$Taccept.rate)[colMeans(fit.t.III.ARp.MAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MAR)

load(paste0(PATH1, "/fit.t.III.ARp.MAR_k5_c2.4.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 5, Mechanism = "MAR", c = 2.4, m = dim(fit.t.III.ARp.MAR$IM$out)[2], loglik = fit.t.III.ARp.MAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MAR$Taccept.rate)[colMeans(fit.t.III.ARp.MAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MAR)

load(paste0(PATH1, "/fit.t.III.ARp.MAR_k5_c3.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 5, Mechanism = "MAR", c = 3.0, m = dim(fit.t.III.ARp.MAR$IM$out)[2], loglik = fit.t.III.ARp.MAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MAR$Taccept.rate)[colMeans(fit.t.III.ARp.MAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MAR)

################################################################################
# 3. K = 5, MNAR
################################################################################

load(paste0(PATH1, "/fit.t.III.ARp.MNAR_k5_c1.8.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 5, Mechanism = "MNAR", c = 1.8, m = dim(fit.t.III.ARp.MNAR$IM$out)[2], loglik = fit.t.III.ARp.MNAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MNAR$Taccept.rate)[colMeans(fit.t.III.ARp.MNAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MNAR)

load(paste0(PATH1, "/fit.t.III.ARp.MNAR_k5_c2.4.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 5, Mechanism = "MNAR", c = 2.4, m = dim(fit.t.III.ARp.MNAR$IM$out)[2], loglik = fit.t.III.ARp.MNAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MNAR$Taccept.rate)[colMeans(fit.t.III.ARp.MNAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MNAR)

load(paste0(PATH1, "/fit.t.III.ARp.MNAR_k5_c3.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 5, Mechanism = "MNAR", c = 3.0, m = dim(fit.t.III.ARp.MNAR$IM$out)[2], loglik = fit.t.III.ARp.MNAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MNAR$Taccept.rate)[colMeans(fit.t.III.ARp.MNAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MNAR)

################################################################################
# 4. K = 10, MCAR
################################################################################

load(paste0(PATH1, "/fit.t.III.ARp.MCAR_k10_c1.8.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 10, Mechanism = "MCAR", c = 1.8, m = dim(fit.t.III.ARp.MCAR$IM$out)[2], loglik = fit.t.III.ARp.MCAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MCAR$Taccept.rate)[colMeans(fit.t.III.ARp.MCAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MCAR)

load(paste0(PATH1, "/fit.t.III.ARp.MCAR_k10_c2.4.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 10, Mechanism = "MCAR", c = 2.4, m = dim(fit.t.III.ARp.MCAR$IM$out)[2], loglik = fit.t.III.ARp.MCAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MCAR$Taccept.rate)[colMeans(fit.t.III.ARp.MCAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MCAR)

load(paste0(PATH1, "/fit.t.III.ARp.MCAR_k10_c3.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 10, Mechanism = "MCAR", c = 3.0, m = dim(fit.t.III.ARp.MCAR$IM$out)[2], loglik = fit.t.III.ARp.MCAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MCAR$Taccept.rate)[colMeans(fit.t.III.ARp.MCAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MCAR)

################################################################################
# 5. K = 10, MAR
################################################################################

load(paste0(PATH1, "/fit.t.III.ARp.MAR_k10_c1.8.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 10, Mechanism = "MAR", c = 1.8, m = dim(fit.t.III.ARp.MAR$IM$out)[2], loglik = fit.t.III.ARp.MAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MAR$Taccept.rate)[colMeans(fit.t.III.ARp.MAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MAR)

load(paste0(PATH1, "/fit.t.III.ARp.MAR_k10_c2.4.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 10, Mechanism = "MAR", c = 2.4, m = dim(fit.t.III.ARp.MAR$IM$out)[2], loglik = fit.t.III.ARp.MAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MAR$Taccept.rate)[colMeans(fit.t.III.ARp.MAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MAR)

load(paste0(PATH1, "/fit.t.III.ARp.MAR_k10_c3.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 10, Mechanism = "MAR", c = 3.0, m = dim(fit.t.III.ARp.MAR$IM$out)[2], loglik = fit.t.III.ARp.MAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MAR$Taccept.rate)[colMeans(fit.t.III.ARp.MAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MAR)

################################################################################
# 6. K = 10, MNAR
################################################################################

load(paste0(PATH1, "/fit.t.III.ARp.MNAR_k10_c1.8.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 10, Mechanism = "MNAR", c = 1.8, m = dim(fit.t.III.ARp.MNAR$IM$out)[2], loglik = fit.t.III.ARp.MNAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MNAR$Taccept.rate)[colMeans(fit.t.III.ARp.MNAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MNAR)

load(paste0(PATH1, "/fit.t.III.ARp.MNAR_k10_c2.4.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 10, Mechanism = "MNAR", c = 2.4, m = dim(fit.t.III.ARp.MNAR$IM$out)[2], loglik = fit.t.III.ARp.MNAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MNAR$Taccept.rate)[colMeans(fit.t.III.ARp.MNAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MNAR)

load(paste0(PATH1, "/fit.t.III.ARp.MNAR_k10_c3.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 10, Mechanism = "MNAR", c = 3.0, m = dim(fit.t.III.ARp.MNAR$IM$out)[2], loglik = fit.t.III.ARp.MNAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MNAR$Taccept.rate)[colMeans(fit.t.III.ARp.MNAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MNAR)

################################################################################
# 7. K = 15, MCAR
################################################################################

load(paste0(PATH1, "/fit.t.III.ARp.MCAR_k15_c1.8.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 15, Mechanism = "MCAR", c = 1.8, m = dim(fit.t.III.ARp.MCAR$IM$out)[2], loglik = fit.t.III.ARp.MCAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MCAR$Taccept.rate)[colMeans(fit.t.III.ARp.MCAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MCAR)

load(paste0(PATH1, "/fit.t.III.ARp.MCAR_k15_c2.4.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 15, Mechanism = "MCAR", c = 2.4, m = dim(fit.t.III.ARp.MCAR$IM$out)[2], loglik = fit.t.III.ARp.MCAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MCAR$Taccept.rate)[colMeans(fit.t.III.ARp.MCAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MCAR)

load(paste0(PATH1, "/fit.t.III.ARp.MCAR_k15_c3.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 15, Mechanism = "MCAR", c = 3.0, m = dim(fit.t.III.ARp.MCAR$IM$out)[2], loglik = fit.t.III.ARp.MCAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MCAR$Taccept.rate)[colMeans(fit.t.III.ARp.MCAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MCAR)

################################################################################
# 8. K = 15, MAR
################################################################################

load(paste0(PATH1, "/fit.t.III.ARp.MAR_k15_c1.8.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 15, Mechanism = "MAR", c = 1.8, m = dim(fit.t.III.ARp.MAR$IM$out)[2], loglik = fit.t.III.ARp.MAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MAR$Taccept.rate)[colMeans(fit.t.III.ARp.MAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MAR)

load(paste0(PATH1, "/fit.t.III.ARp.MAR_k15_c2.4.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 15, Mechanism = "MAR", c = 2.4, m = dim(fit.t.III.ARp.MAR$IM$out)[2], loglik = fit.t.III.ARp.MAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MAR$Taccept.rate)[colMeans(fit.t.III.ARp.MAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MAR)

load(paste0(PATH1, "/fit.t.III.ARp.MAR_k15_c3.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 15, Mechanism = "MAR", c = 3.0, m = dim(fit.t.III.ARp.MAR$IM$out)[2], loglik = fit.t.III.ARp.MAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MAR$Taccept.rate)[colMeans(fit.t.III.ARp.MAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MAR)

################################################################################
# 9. K = 15, MNAR
################################################################################

load(paste0(PATH1, "/fit.t.III.ARp.MNAR_k15_c1.8.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 15, Mechanism = "MNAR", c = 1.8, m = dim(fit.t.III.ARp.MNAR$IM$out)[2], loglik = fit.t.III.ARp.MNAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MNAR$Taccept.rate)[colMeans(fit.t.III.ARp.MNAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MNAR)

load(paste0(PATH1, "/fit.t.III.ARp.MNAR_k15_c2.4.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 15, Mechanism = "MNAR", c = 2.4, m = dim(fit.t.III.ARp.MNAR$IM$out)[2], loglik = fit.t.III.ARp.MNAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MNAR$Taccept.rate)[colMeans(fit.t.III.ARp.MNAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MNAR)

load(paste0(PATH1, "/fit.t.III.ARp.MNAR_k15_c3.RData"))
TableE1.raw <- rbind(TableE1.raw, data.frame(K = 15, Mechanism = "MNAR", c = 3.0, m = dim(fit.t.III.ARp.MNAR$IM$out)[2], loglik = fit.t.III.ARp.MNAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MNAR$Taccept.rate)[colMeans(fit.t.III.ARp.MNAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MNAR)

################################################################################
# 10. Check and save the extracted quantities
################################################################################

if (nrow(TableE1.raw) != 27) stop("TableE1_raw should contain exactly 27 tuning-model fits.")

if (any(TableE1.raw$m[TableE1.raw$Mechanism == "MCAR"] != 22)) stop("Unexpected number of MCAR parameters.")
if (any(TableE1.raw$m[TableE1.raw$Mechanism == "MAR"] != 25)) stop("Unexpected number of MAR parameters.")
if (any(TableE1.raw$m[TableE1.raw$Mechanism == "MNAR"] != 26)) stop("Unexpected number of MNAR parameters.")

write.table(TableE1.raw, paste0(PATH1, "/TableE1_raw.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

################################################################################
# 11. Read TableE1_raw.txt and calculate AIC and BIC
################################################################################

rm(TableE1.raw)

TableE1.raw <- read.table(paste0(PATH1, "/TableE1_raw.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

N <- 481

TableE1.raw$AIC <- 2 * TableE1.raw$m - 2 * TableE1.raw$loglik
TableE1.raw$BIC <- TableE1.raw$m * log(N) - 2 * TableE1.raw$loglik

################################################################################
# 12. Convert to the layout of Table E.1
################################################################################

TableE1.long <- rbind(
  data.frame(K = TableE1.raw$K, Mechanism = TableE1.raw$Mechanism, Criterion = "ell_max", c = TableE1.raw$c, value = TableE1.raw$loglik),
  data.frame(K = TableE1.raw$K, Mechanism = TableE1.raw$Mechanism, Criterion = "AIC", c = TableE1.raw$c, value = TableE1.raw$AIC),
  data.frame(K = TableE1.raw$K, Mechanism = TableE1.raw$Mechanism, Criterion = "BIC", c = TableE1.raw$c, value = TableE1.raw$BIC),
  data.frame(K = TableE1.raw$K, Mechanism = TableE1.raw$Mechanism, Criterion = "Acceptance rate", c = TableE1.raw$c, value = TableE1.raw$acceptance_rate)
)

TableE1.long$Mechanism <- factor(TableE1.long$Mechanism, levels = c("MCAR", "MAR", "MNAR"))
TableE1.long$Criterion <- factor(TableE1.long$Criterion, levels = c("ell_max", "AIC", "BIC", "Acceptance rate"))
TableE1.long <- TableE1.long[order(TableE1.long$K, TableE1.long$Mechanism, TableE1.long$Criterion, TableE1.long$c), ]

TableE1 <- reshape(TableE1.long, idvar = c("K", "Mechanism", "Criterion"), timevar = "c", direction = "wide")

TableE1 <- TableE1[order(TableE1$K, TableE1$Mechanism, TableE1$Criterion), ]
row.names(TableE1) <- NULL

TableE1$Mechanism <- as.character(TableE1$Mechanism)
TableE1$Criterion <- as.character(TableE1$Criterion)

colnames(TableE1) <- c("K", "Mechanism", "Criterion", "c = 1.8", "c = 2.4", "c = 3.0")

TableE1[, c("c = 1.8", "c = 2.4", "c = 3.0")] <- round(TableE1[, c("c = 1.8", "c = 2.4", "c = 3.0")], 3)

################################################################################
# 13. Match the display format of Table E.1
################################################################################

TableE1$K <- as.character(TableE1$K)

TableE1$K[c(2:12, 14:24, 26:36)] <- ""
TableE1$Mechanism[c(2:4, 6:8, 10:12, 14:16, 18:20, 22:24, 26:28, 30:32, 34:36)] <- ""

################################################################################
# 14. Save Table E.1
################################################################################
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

write.csv(TableE1, paste0(PATH, "/Result/TableE1.csv"), row.names = FALSE, na = "")


