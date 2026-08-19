################################################################################
#
#   Filename    :    TabF1.R
#   Project     :    Biometrics article "Modeling HIV Viral Dynamics Using a
#                    Nonlinear Mixed-Effects Framework for Heavy-Tailed Data
#                    with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#   Date        :    18.08.2026
#   Purpose     :    produce Table F.1 for the Metropolis-Hastings tuning
#                    sensitivity analysis by extracting the number of model
#                    parameters, approximated observed-data log-likelihood, and
#                    acceptance rate from the Scenario (III) AR(1) tNLME fits,
#                    and then calculating AIC and BIC
#
#   Input data files  :  Data_and_Code/Data/sensitivity_SAEM/
#                         fit.t.III.ARp.MCAR_k5_c1.8.RData, ..., 
#                         fit.t.III.ARp.MNAR_k15_c3.RData
#
#   Intermediate file :  Data_and_Code/Data/sensitivity_SAEM/TableF1_raw.txt
#
#   Output data files :  Data_and_Code/Result/TableF1.csv
#
#   R Version   :    R-4.6.0
#   Required R packages : none
#
################################################################################

PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
PATH1 <- paste0(PATH, "/Data/sensitivity_SAEM")

options(digits = 16)

TableF1.raw <- data.frame(K = numeric(0), Mechanism = character(0), c = numeric(0), m = numeric(0), loglik = numeric(0), acceptance_rate = numeric(0), stringsAsFactors = FALSE)

################################################################################
# 1. K = 5, MCAR
################################################################################

load(paste0(PATH1, "/fit.t.III.ARp.MCAR_k5_c1.8.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 5, Mechanism = "MCAR", c = 1.8, m = dim(fit.t.III.ARp.MCAR$IM$out)[2], loglik = fit.t.III.ARp.MCAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MCAR$Taccept.rate)[colMeans(fit.t.III.ARp.MCAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MCAR)

load(paste0(PATH1, "/fit.t.III.ARp.MCAR_k5_c2.4.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 5, Mechanism = "MCAR", c = 2.4, m = dim(fit.t.III.ARp.MCAR$IM$out)[2], loglik = fit.t.III.ARp.MCAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MCAR$Taccept.rate)[colMeans(fit.t.III.ARp.MCAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MCAR)

load(paste0(PATH1, "/fit.t.III.ARp.MCAR_k5_c3.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 5, Mechanism = "MCAR", c = 3.0, m = dim(fit.t.III.ARp.MCAR$IM$out)[2], loglik = fit.t.III.ARp.MCAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MCAR$Taccept.rate)[colMeans(fit.t.III.ARp.MCAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MCAR)

################################################################################
# 2. K = 5, MAR
################################################################################

load(paste0(PATH1, "/fit.t.III.ARp.MAR_k5_c1.8.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 5, Mechanism = "MAR", c = 1.8, m = dim(fit.t.III.ARp.MAR$IM$out)[2], loglik = fit.t.III.ARp.MAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MAR$Taccept.rate)[colMeans(fit.t.III.ARp.MAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MAR)

load(paste0(PATH1, "/fit.t.III.ARp.MAR_k5_c2.4.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 5, Mechanism = "MAR", c = 2.4, m = dim(fit.t.III.ARp.MAR$IM$out)[2], loglik = fit.t.III.ARp.MAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MAR$Taccept.rate)[colMeans(fit.t.III.ARp.MAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MAR)

load(paste0(PATH1, "/fit.t.III.ARp.MAR_k5_c3.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 5, Mechanism = "MAR", c = 3.0, m = dim(fit.t.III.ARp.MAR$IM$out)[2], loglik = fit.t.III.ARp.MAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MAR$Taccept.rate)[colMeans(fit.t.III.ARp.MAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MAR)

################################################################################
# 3. K = 5, MNAR
################################################################################

load(paste0(PATH1, "/fit.t.III.ARp.MNAR_k5_c1.8.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 5, Mechanism = "MNAR", c = 1.8, m = dim(fit.t.III.ARp.MNAR$IM$out)[2], loglik = fit.t.III.ARp.MNAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MNAR$Taccept.rate)[colMeans(fit.t.III.ARp.MNAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MNAR)

load(paste0(PATH1, "/fit.t.III.ARp.MNAR_k5_c2.4.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 5, Mechanism = "MNAR", c = 2.4, m = dim(fit.t.III.ARp.MNAR$IM$out)[2], loglik = fit.t.III.ARp.MNAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MNAR$Taccept.rate)[colMeans(fit.t.III.ARp.MNAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MNAR)

load(paste0(PATH1, "/fit.t.III.ARp.MNAR_k5_c3.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 5, Mechanism = "MNAR", c = 3.0, m = dim(fit.t.III.ARp.MNAR$IM$out)[2], loglik = fit.t.III.ARp.MNAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MNAR$Taccept.rate)[colMeans(fit.t.III.ARp.MNAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MNAR)

################################################################################
# 4. K = 10, MCAR
################################################################################

load(paste0(PATH1, "/fit.t.III.ARp.MCAR_k10_c1.8.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 10, Mechanism = "MCAR", c = 1.8, m = dim(fit.t.III.ARp.MCAR$IM$out)[2], loglik = fit.t.III.ARp.MCAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MCAR$Taccept.rate)[colMeans(fit.t.III.ARp.MCAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MCAR)

load(paste0(PATH1, "/fit.t.III.ARp.MCAR_k10_c2.4.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 10, Mechanism = "MCAR", c = 2.4, m = dim(fit.t.III.ARp.MCAR$IM$out)[2], loglik = fit.t.III.ARp.MCAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MCAR$Taccept.rate)[colMeans(fit.t.III.ARp.MCAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MCAR)

load(paste0(PATH1, "/fit.t.III.ARp.MCAR_k10_c3.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 10, Mechanism = "MCAR", c = 3.0, m = dim(fit.t.III.ARp.MCAR$IM$out)[2], loglik = fit.t.III.ARp.MCAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MCAR$Taccept.rate)[colMeans(fit.t.III.ARp.MCAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MCAR)

################################################################################
# 5. K = 10, MAR
################################################################################

load(paste0(PATH1, "/fit.t.III.ARp.MAR_k10_c1.8.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 10, Mechanism = "MAR", c = 1.8, m = dim(fit.t.III.ARp.MAR$IM$out)[2], loglik = fit.t.III.ARp.MAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MAR$Taccept.rate)[colMeans(fit.t.III.ARp.MAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MAR)

load(paste0(PATH1, "/fit.t.III.ARp.MAR_k10_c2.4.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 10, Mechanism = "MAR", c = 2.4, m = dim(fit.t.III.ARp.MAR$IM$out)[2], loglik = fit.t.III.ARp.MAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MAR$Taccept.rate)[colMeans(fit.t.III.ARp.MAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MAR)

load(paste0(PATH1, "/fit.t.III.ARp.MAR_k10_c3.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 10, Mechanism = "MAR", c = 3.0, m = dim(fit.t.III.ARp.MAR$IM$out)[2], loglik = fit.t.III.ARp.MAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MAR$Taccept.rate)[colMeans(fit.t.III.ARp.MAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MAR)

################################################################################
# 6. K = 10, MNAR
################################################################################

load(paste0(PATH1, "/fit.t.III.ARp.MNAR_k10_c1.8.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 10, Mechanism = "MNAR", c = 1.8, m = dim(fit.t.III.ARp.MNAR$IM$out)[2], loglik = fit.t.III.ARp.MNAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MNAR$Taccept.rate)[colMeans(fit.t.III.ARp.MNAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MNAR)

load(paste0(PATH1, "/fit.t.III.ARp.MNAR_k10_c2.4.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 10, Mechanism = "MNAR", c = 2.4, m = dim(fit.t.III.ARp.MNAR$IM$out)[2], loglik = fit.t.III.ARp.MNAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MNAR$Taccept.rate)[colMeans(fit.t.III.ARp.MNAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MNAR)

load(paste0(PATH1, "/fit.t.III.ARp.MNAR_k10_c3.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 10, Mechanism = "MNAR", c = 3.0, m = dim(fit.t.III.ARp.MNAR$IM$out)[2], loglik = fit.t.III.ARp.MNAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MNAR$Taccept.rate)[colMeans(fit.t.III.ARp.MNAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MNAR)

################################################################################
# 7. K = 15, MCAR
################################################################################

load(paste0(PATH1, "/fit.t.III.ARp.MCAR_k15_c1.8.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 15, Mechanism = "MCAR", c = 1.8, m = dim(fit.t.III.ARp.MCAR$IM$out)[2], loglik = fit.t.III.ARp.MCAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MCAR$Taccept.rate)[colMeans(fit.t.III.ARp.MCAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MCAR)

load(paste0(PATH1, "/fit.t.III.ARp.MCAR_k15_c2.4.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 15, Mechanism = "MCAR", c = 2.4, m = dim(fit.t.III.ARp.MCAR$IM$out)[2], loglik = fit.t.III.ARp.MCAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MCAR$Taccept.rate)[colMeans(fit.t.III.ARp.MCAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MCAR)

load(paste0(PATH1, "/fit.t.III.ARp.MCAR_k15_c3.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 15, Mechanism = "MCAR", c = 3.0, m = dim(fit.t.III.ARp.MCAR$IM$out)[2], loglik = fit.t.III.ARp.MCAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MCAR$Taccept.rate)[colMeans(fit.t.III.ARp.MCAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MCAR)

################################################################################
# 8. K = 15, MAR
################################################################################

load(paste0(PATH1, "/fit.t.III.ARp.MAR_k15_c1.8.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 15, Mechanism = "MAR", c = 1.8, m = dim(fit.t.III.ARp.MAR$IM$out)[2], loglik = fit.t.III.ARp.MAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MAR$Taccept.rate)[colMeans(fit.t.III.ARp.MAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MAR)

load(paste0(PATH1, "/fit.t.III.ARp.MAR_k15_c2.4.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 15, Mechanism = "MAR", c = 2.4, m = dim(fit.t.III.ARp.MAR$IM$out)[2], loglik = fit.t.III.ARp.MAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MAR$Taccept.rate)[colMeans(fit.t.III.ARp.MAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MAR)

load(paste0(PATH1, "/fit.t.III.ARp.MAR_k15_c3.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 15, Mechanism = "MAR", c = 3.0, m = dim(fit.t.III.ARp.MAR$IM$out)[2], loglik = fit.t.III.ARp.MAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MAR$Taccept.rate)[colMeans(fit.t.III.ARp.MAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MAR)

################################################################################
# 9. K = 15, MNAR
################################################################################

load(paste0(PATH1, "/fit.t.III.ARp.MNAR_k15_c1.8.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 15, Mechanism = "MNAR", c = 1.8, m = dim(fit.t.III.ARp.MNAR$IM$out)[2], loglik = fit.t.III.ARp.MNAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MNAR$Taccept.rate)[colMeans(fit.t.III.ARp.MNAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MNAR)

load(paste0(PATH1, "/fit.t.III.ARp.MNAR_k15_c2.4.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 15, Mechanism = "MNAR", c = 2.4, m = dim(fit.t.III.ARp.MNAR$IM$out)[2], loglik = fit.t.III.ARp.MNAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MNAR$Taccept.rate)[colMeans(fit.t.III.ARp.MNAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MNAR)

load(paste0(PATH1, "/fit.t.III.ARp.MNAR_k15_c3.RData"))
TableF1.raw <- rbind(TableF1.raw, data.frame(K = 15, Mechanism = "MNAR", c = 3.0, m = dim(fit.t.III.ARp.MNAR$IM$out)[2], loglik = fit.t.III.ARp.MNAR$model.inf$loglik, acceptance_rate = mean(colMeans(fit.t.III.ARp.MNAR$Taccept.rate)[colMeans(fit.t.III.ARp.MNAR$Taccept.rate) != 0])))
rm(fit.t.III.ARp.MNAR)

################################################################################
# 10. Check and save the extracted quantities
################################################################################

if (nrow(TableF1.raw) != 27) stop("TableF1_raw should contain exactly 27 tuning-model fits.")

if (any(TableF1.raw$m[TableF1.raw$Mechanism == "MCAR"] != 22)) stop("Unexpected number of MCAR parameters.")
if (any(TableF1.raw$m[TableF1.raw$Mechanism == "MAR"] != 25)) stop("Unexpected number of MAR parameters.")
if (any(TableF1.raw$m[TableF1.raw$Mechanism == "MNAR"] != 26)) stop("Unexpected number of MNAR parameters.")

write.table(TableF1.raw, paste0(PATH1, "/TableF1_raw.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

################################################################################
# 11. Read TableF1_raw.txt and calculate AIC and BIC
################################################################################

rm(TableF1.raw)

TableF1.raw <- read.table(paste0(PATH1, "/TableF1_raw.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

N <- 481

TableF1.raw$AIC <- 2 * TableF1.raw$m - 2 * TableF1.raw$loglik
TableF1.raw$BIC <- TableF1.raw$m * log(N) - 2 * TableF1.raw$loglik

################################################################################
# 12. Convert to the layout of Table F.1
################################################################################

TableF1.long <- rbind(
  data.frame(K = TableF1.raw$K, Mechanism = TableF1.raw$Mechanism, Criterion = "ell_max", c = TableF1.raw$c, value = TableF1.raw$loglik),
  data.frame(K = TableF1.raw$K, Mechanism = TableF1.raw$Mechanism, Criterion = "AIC", c = TableF1.raw$c, value = TableF1.raw$AIC),
  data.frame(K = TableF1.raw$K, Mechanism = TableF1.raw$Mechanism, Criterion = "BIC", c = TableF1.raw$c, value = TableF1.raw$BIC),
  data.frame(K = TableF1.raw$K, Mechanism = TableF1.raw$Mechanism, Criterion = "Acceptance rate", c = TableF1.raw$c, value = TableF1.raw$acceptance_rate)
)

TableF1.long$Mechanism <- factor(TableF1.long$Mechanism, levels = c("MCAR", "MAR", "MNAR"))
TableF1.long$Criterion <- factor(TableF1.long$Criterion, levels = c("ell_max", "AIC", "BIC", "Acceptance rate"))
TableF1.long <- TableF1.long[order(TableF1.long$K, TableF1.long$Mechanism, TableF1.long$Criterion, TableF1.long$c), ]

TableF1 <- reshape(TableF1.long, idvar = c("K", "Mechanism", "Criterion"), timevar = "c", direction = "wide")

TableF1 <- TableF1[order(TableF1$K, TableF1$Mechanism, TableF1$Criterion), ]
row.names(TableF1) <- NULL

TableF1$Mechanism <- as.character(TableF1$Mechanism)
TableF1$Criterion <- as.character(TableF1$Criterion)

colnames(TableF1) <- c("K", "Mechanism", "Criterion", "c = 1.8", "c = 2.4", "c = 3.0")

TableF1[, c("c = 1.8", "c = 2.4", "c = 3.0")] <- round(TableF1[, c("c = 1.8", "c = 2.4", "c = 3.0")], 3)

################################################################################
# 13. Match the display format of Table F.1
################################################################################

TableF1$K <- as.character(TableF1$K)

TableF1$K[c(2:12, 14:24, 26:36)] <- ""
TableF1$Mechanism[c(2:4, 6:8, 10:12, 14:16, 18:20, 22:24, 26:28, 30:32, 34:36)] <- ""

################################################################################
# 14. Save Table F.1
################################################################################
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

write.csv(TableF1, paste0(PATH, "/Result/TableF1.csv"), row.names = FALSE, na = "")


