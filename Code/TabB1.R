################################################################################
#
#   Filename    :    TabB1.R
#   Project     :    Biometrics article "Modeling HIV Viral Dynamics Using a
#                    Nonlinear Mixed-Effects Framework for Heavy-Tailed Data
#                    with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#   Date        :    18.08.2026
#   Purpose     :    produce Table B.1 for the simulation study by comparing
#                    information-matrix standard errors (IM SE) with empirical
#                    Monte Carlo standard deviations (MC SD) from the fitted
#                    tNLME model under the MNAR mechanism across sample sizes
#                    and dropout levels
#
#   Input data files  :  Data_and_Code/Data/simulation/SS-simulation-t25/SIM1-SIM5/
#                         Data_and_Code/Data/simulation/SS-simulation-t50/SIM1-SIM5/
#                         Data_and_Code/Data/simulation/SS-simulation-t75/SIM1-SIM5/
#
#   Output data files :  Data_and_Code/Result/TableB1.csv
#
#   R Version   :    R-4.6.0
#   Required R packages : none
#
################################################################################
#################
#################  25 %; nu=5
#################
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
PATH1 <- paste(PATH, "/Data/simulation/SS-simulation-t25/", sep = "")

realpara <- colMeans(as.matrix(read.table(paste(PATH1, "SIM1/realpara.txt", sep = ""), na.strings = "NA", sep = "")))
names(realpara) <- c(rep("Beta", 3), "sigma", rep("DD", 1), "nu", rep("alpha", 4))

MAR.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]

MAR.se1 <- as.matrix(read.table(paste(PATH1, "SIM1/MAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.se1 <- as.matrix(read.table(paste(PATH1, "SIM1/MCAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.se1 <- as.matrix(read.table(paste(PATH1, "SIM1/MNAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.se2 <- as.matrix(read.table(paste(PATH1, "SIM2/MAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.se2 <- as.matrix(read.table(paste(PATH1, "SIM2/MCAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.se2 <- as.matrix(read.table(paste(PATH1, "SIM2/MNAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.se3 <- as.matrix(read.table(paste(PATH1, "SIM3/MAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.se3 <- as.matrix(read.table(paste(PATH1, "SIM3/MCAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.se3 <- as.matrix(read.table(paste(PATH1, "SIM3/MNAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.se4 <- as.matrix(read.table(paste(PATH1, "SIM4/MAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.se4 <- as.matrix(read.table(paste(PATH1, "SIM4/MCAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.se4 <- as.matrix(read.table(paste(PATH1, "SIM4/MNAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.se5 <- as.matrix(read.table(paste(PATH1, "SIM5/MAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.se5 <- as.matrix(read.table(paste(PATH1, "SIM5/MCAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.se5 <- as.matrix(read.table(paste(PATH1, "SIM5/MNAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]


sum.table <- rbind(
  round(colMeans(MNAR.se1, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est1, 2, sd))), nrow = 1), 3),
  round(colMeans(MNAR.se2, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est2, 2, sd))), nrow = 1), 3),
  round(colMeans(MNAR.se3, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est3, 2, sd))), nrow = 1), 3),
  round(colMeans(MNAR.se4, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est4, 2, sd))), nrow = 1), 3),
  round(colMeans(MNAR.se5, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est5, 2, sd))), nrow = 1), 3)
)

##########################################################################################
########################################################################
########################################################################  50%
PATH1 <- paste(PATH, "/Data/simulation/SS-simulation-t50/", sep = "")

realpara <- colMeans(as.matrix(read.table(paste(PATH1, "SIM1/realpara.txt", sep = ""), na.strings = "NA", sep = "")))
names(realpara) <- c(rep("Beta", 3), "sigma", rep("DD", 1), "nu", rep("alpha", 4))

MAR.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]

MAR.se1 <- as.matrix(read.table(paste(PATH1, "SIM1/MAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.se1 <- as.matrix(read.table(paste(PATH1, "SIM1/MCAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.se1 <- as.matrix(read.table(paste(PATH1, "SIM1/MNAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.se2 <- as.matrix(read.table(paste(PATH1, "SIM2/MAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.se2 <- as.matrix(read.table(paste(PATH1, "SIM2/MCAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.se2 <- as.matrix(read.table(paste(PATH1, "SIM2/MNAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.se3 <- as.matrix(read.table(paste(PATH1, "SIM3/MAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.se3 <- as.matrix(read.table(paste(PATH1, "SIM3/MCAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.se3 <- as.matrix(read.table(paste(PATH1, "SIM3/MNAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.se4 <- as.matrix(read.table(paste(PATH1, "SIM4/MAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.se4 <- as.matrix(read.table(paste(PATH1, "SIM4/MCAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.se4 <- as.matrix(read.table(paste(PATH1, "SIM4/MNAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.se5 <- as.matrix(read.table(paste(PATH1, "SIM5/MAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.se5 <- as.matrix(read.table(paste(PATH1, "SIM5/MCAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.se5 <- as.matrix(read.table(paste(PATH1, "SIM5/MNAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]

sum.table <- rbind(
  sum.table,
  rbind(
    round(colMeans(MNAR.se1, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est1, 2, sd))), nrow = 1), 3),
    round(colMeans(MNAR.se2, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est2, 2, sd))), nrow = 1), 3),
    round(colMeans(MNAR.se3, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est3, 2, sd))), nrow = 1), 3),
    round(colMeans(MNAR.se4, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est4, 2, sd))), nrow = 1), 3),
    round(colMeans(MNAR.se5, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est5, 2, sd))), nrow = 1), 3)
  )
)


##########################################################################################
########################################################################
########################################################################  75%
PATH1 <- paste(PATH, "/Data/simulation/SS-simulation-t75/", sep = "")

realpara <- colMeans(as.matrix(read.table(paste(PATH1, "SIM1/realpara.txt", sep = ""), na.strings = "NA", sep = "")))
names(realpara) <- c(rep("Beta", 3), "sigma", rep("DD", 1), "nu", rep("alpha", 4))

MAR.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]

MAR.se1 <- as.matrix(read.table(paste(PATH1, "SIM1/MAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.se1 <- as.matrix(read.table(paste(PATH1, "SIM1/MCAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.se1 <- as.matrix(read.table(paste(PATH1, "SIM1/MNAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.se2 <- as.matrix(read.table(paste(PATH1, "SIM2/MAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.se2 <- as.matrix(read.table(paste(PATH1, "SIM2/MCAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.se2 <- as.matrix(read.table(paste(PATH1, "SIM2/MNAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.se3 <- as.matrix(read.table(paste(PATH1, "SIM3/MAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.se3 <- as.matrix(read.table(paste(PATH1, "SIM3/MCAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.se3 <- as.matrix(read.table(paste(PATH1, "SIM3/MNAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.se4 <- as.matrix(read.table(paste(PATH1, "SIM4/MAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.se4 <- as.matrix(read.table(paste(PATH1, "SIM4/MCAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.se4 <- as.matrix(read.table(paste(PATH1, "SIM4/MNAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.se5 <- as.matrix(read.table(paste(PATH1, "SIM5/MAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.se5 <- as.matrix(read.table(paste(PATH1, "SIM5/MCAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.se5 <- as.matrix(read.table(paste(PATH1, "SIM5/MNAR.se.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]


sum.table <- rbind(
  sum.table,
  rbind(
    round(colMeans(MNAR.se1, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est1, 2, sd))), nrow = 1), 3),
    round(colMeans(MNAR.se2, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est2, 2, sd))), nrow = 1), 3),
    round(colMeans(MNAR.se3, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est3, 2, sd))), nrow = 1), 3),
    round(colMeans(MNAR.se4, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est4, 2, sd))), nrow = 1), 3),
    round(colMeans(MNAR.se5, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est5, 2, sd))), nrow = 1), 3)
  )
)


colnames(sum.table) <- c("beta1", "beta2", "beta3", "sigma^2", "d", "nu", "alpha00", "alpha01", "alpha1", "alpha2")
row.names(sum.table) <- c(rep(c("IM SE", "MC Sd"), time = 5 * 3))
# print(sum.table)

write.csv(sum.table, paste0(PATH, "/Result/TableB1.csv"))
