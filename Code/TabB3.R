################################################################################
#
#   Filename    :    TabB3.R
#   Project     :    Biometrics article "Modeling HIV Viral Dynamics Using a
#                    Nonlinear Mixed-Effects Framework for Heavy-Tailed Data
#                    with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#   Date        :    18.08.2026
#   Purpose     :    produce Table B.3 for the simulation study by comparing
#                    AIC and BIC among the MCAR, MAR, and MNAR specifications
#                    of the tNLME model across sample sizes and dropout levels
#
#   Input data files  :  Data_and_Code/Data/simulation/SS-simulation-t25/SIM1-SIM5/
#                         Data_and_Code/Data/simulation/SS-simulation-t50/SIM1-SIM5/
#                         Data_and_Code/Data/simulation/SS-simulation-t75/SIM1-SIM5/
#
#   Output data files :  Data_and_Code/Result/TableB3.csv
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

loglik1 =    as.matrix(read.table(paste(PATH1,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE1 =    as.matrix(read.table(paste(PATH1,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara1 = colMeans(as.matrix(read.table(paste(PATH1,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik2 =    as.matrix(read.table(paste(PATH2,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE2 =    as.matrix(read.table(paste(PATH2,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara2 = colMeans(as.matrix(read.table(paste(PATH2,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik3 =    as.matrix(read.table(paste(PATH3,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE3 =    as.matrix(read.table(paste(PATH3,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara3 = colMeans(as.matrix(read.table(paste(PATH3,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik4 =    as.matrix(read.table(paste(PATH4,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE4 =    as.matrix(read.table(paste(PATH4,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara4 = colMeans(as.matrix(read.table(paste(PATH4,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik5 =    as.matrix(read.table(paste(PATH5,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE5 =    as.matrix(read.table(paste(PATH5,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara5 = colMeans(as.matrix(read.table(paste(PATH5,'realpara.txt',sep=""),na.strings="NA",sep="")))

####################
##########
##########   table1
nn = c(25, 50, 100, 200, 400)
aic1 = aic2 = bic1 = bic2 = NULL
mMNAR = length(realpara1)
mMAR = length(realpara1)-1
mMCAR = length(realpara1)-2

mmm = c(mMNAR, mMAR, mMCAR)

n = nn[1]
AIc = 2*(mmm) -2*t(loglik1[,2:4])
BIc = mmm*log(n) -2*t(loglik1[,2:4])

sum.tabel = rbind(colMeans(t(AIc)),
                  apply(AIc, 1, sd),
                  c(sum(apply(AIc[1:3,], 2, order)[1,]==1), sum(apply(AIc[1:3,], 2, order)[2,]==1), sum(apply(AIc[1:3,], 2, order)[3,]==1)),
                  colMeans(t(BIc)),
                  apply(BIc, 1, sd),
                  c(sum(apply(BIc[1:3,], 2, order)[1,]==1), sum(apply(BIc[1:3,], 2, order)[2,]==1), sum(apply(BIc[1:3,], 2, order)[3,]==1))
                  )


n = nn[2]
AIc = 2*(mmm) -2*t(loglik2[,2:4])
BIc = mmm*log(n) -2*t(loglik2[,2:4])

sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(AIc)),
                        apply(AIc, 1, sd),
                        c(sum(apply(AIc[1:3,], 2, order)[1,]==1), sum(apply(AIc[1:3,], 2, order)[2,]==1), sum(apply(AIc[1:3,], 2, order)[3,]==1)),
                          colMeans(t(BIc)),
                          apply(BIc, 1, sd),
                          c(sum(apply(BIc[1:3,], 2, order)[1,]==1), sum(apply(BIc[1:3,], 2, order)[2,]==1), sum(apply(BIc[1:3,], 2, order)[3,]==1))
                        ))

n = nn[3]
AIc = 2*(mmm) -2*t(loglik3[,2:4])
BIc = mmm*log(n) -2*t(loglik3[,2:4])

sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(AIc)),
                        apply(AIc, 1, sd),
                        c(sum(apply(AIc[1:3,], 2, order)[1,]==1), sum(apply(AIc[1:3,], 2, order)[2,]==1), sum(apply(AIc[1:3,], 2, order)[3,]==1)),
                          colMeans(t(BIc)),
                          apply(BIc, 1, sd),
                          c(sum(apply(BIc[1:3,], 2, order)[1,]==1), sum(apply(BIc[1:3,], 2, order)[2,]==1), sum(apply(BIc[1:3,], 2, order)[3,]==1))
                        ))


n = nn[4]
AIc = 2*(mmm) -2*t(loglik4[,2:4])
BIc = mmm*log(n) -2*t(loglik4[,2:4])

sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(AIc)),
                        apply(AIc, 1, sd),
                        c(sum(apply(AIc[1:3,], 2, order)[1,]==1), sum(apply(AIc[1:3,], 2, order)[2,]==1), sum(apply(AIc[1:3,], 2, order)[3,]==1)),
                          colMeans(t(BIc)),
                          apply(BIc, 1, sd),
                          c(sum(apply(BIc[1:3,], 2, order)[1,]==1), sum(apply(BIc[1:3,], 2, order)[2,]==1), sum(apply(BIc[1:3,], 2, order)[3,]==1))
                        ))

n = nn[2]
AIc = 2*(mmm) -2*t(loglik5[,2:4])
BIc = mmm*log(n) -2*t(loglik5[,2:4])

sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(AIc)),
                        apply(AIc, 1, sd),
                        c(sum(apply(AIc[1:3,], 2, order)[1,]==1), sum(apply(AIc[1:3,], 2, order)[2,]==1), sum(apply(AIc[1:3,], 2, order)[3,]==1)),
                          colMeans(t(BIc)),
                          apply(BIc, 1, sd),
                          c(sum(apply(BIc[1:3,], 2, order)[1,]==1), sum(apply(BIc[1:3,], 2, order)[2,]==1), sum(apply(BIc[1:3,], 2, order)[3,]==1))
                        ))

### MCAR MAR MNAR
sum.tabel.T = round(sum.tabel, 3)[-c(2,5,8,11,14,17,20,23,26,29),][,c(3,2,1)]

################################################################################
#################
#################  50 %; nu=5
#################
PATH0 <- paste(PATH, "/Data/simulation/SS-simulation-t50/", sep = "")
PATH1 <- paste(PATH0, "SIM1/", sep = "")
PATH2 <- paste(PATH0, "SIM2/", sep = "")
PATH3 <- paste(PATH0, "SIM3/", sep = "")
PATH4 <- paste(PATH0, "SIM4/", sep = "")
PATH5 <- paste(PATH0, "SIM5/", sep = "")

loglik1 =    as.matrix(read.table(paste(PATH1,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE1 =    as.matrix(read.table(paste(PATH1,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara1 = colMeans(as.matrix(read.table(paste(PATH1,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik2 =    as.matrix(read.table(paste(PATH2,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE2 =    as.matrix(read.table(paste(PATH2,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara2 = colMeans(as.matrix(read.table(paste(PATH2,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik3 =    as.matrix(read.table(paste(PATH3,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE3 =    as.matrix(read.table(paste(PATH3,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara3 = colMeans(as.matrix(read.table(paste(PATH3,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik4 =    as.matrix(read.table(paste(PATH4,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE4 =    as.matrix(read.table(paste(PATH4,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara4 = colMeans(as.matrix(read.table(paste(PATH4,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik5 =    as.matrix(read.table(paste(PATH5,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE5 =    as.matrix(read.table(paste(PATH5,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara5 = colMeans(as.matrix(read.table(paste(PATH5,'realpara.txt',sep=""),na.strings="NA",sep="")))

####################
##########
##########   table1
nn = c(25, 50, 100, 200, 400)
aic1 = aic2 = bic1 = bic2 = NULL
mMNAR = length(realpara1)
mMAR = length(realpara1)-1
mMCAR = length(realpara1)-2

mmm = c(mMNAR, mMAR, mMCAR)

n = nn[1]
AIc = 2*(mmm) -2*t(loglik1[,2:4])
BIc = mmm*log(n) -2*t(loglik1[,2:4])

sum.tabel = rbind(colMeans(t(AIc)),
                  apply(AIc, 1, sd),
                  c(sum(apply(AIc[1:3,], 2, order)[1,]==1), sum(apply(AIc[1:3,], 2, order)[2,]==1), sum(apply(AIc[1:3,], 2, order)[3,]==1)),
                  colMeans(t(BIc)),
                  apply(BIc, 1, sd),
                  c(sum(apply(BIc[1:3,], 2, order)[1,]==1), sum(apply(BIc[1:3,], 2, order)[2,]==1), sum(apply(BIc[1:3,], 2, order)[3,]==1))
)


n = nn[2]
AIc = 2*(mmm) -2*t(loglik2[,2:4])
BIc = mmm*log(n) -2*t(loglik2[,2:4])

sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(AIc)),
                        apply(AIc, 1, sd),
                        c(sum(apply(AIc[1:3,], 2, order)[1,]==1), sum(apply(AIc[1:3,], 2, order)[2,]==1), sum(apply(AIc[1:3,], 2, order)[3,]==1)),
                        colMeans(t(BIc)),
                        apply(BIc, 1, sd),
                        c(sum(apply(BIc[1:3,], 2, order)[1,]==1), sum(apply(BIc[1:3,], 2, order)[2,]==1), sum(apply(BIc[1:3,], 2, order)[3,]==1))
                  ))

n = nn[3]
AIc = 2*(mmm) -2*t(loglik3[,2:4])
BIc = mmm*log(n) -2*t(loglik3[,2:4])

sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(AIc)),
                        apply(AIc, 1, sd),
                        c(sum(apply(AIc[1:3,], 2, order)[1,]==1), sum(apply(AIc[1:3,], 2, order)[2,]==1), sum(apply(AIc[1:3,], 2, order)[3,]==1)),
                        colMeans(t(BIc)),
                        apply(BIc, 1, sd),
                        c(sum(apply(BIc[1:3,], 2, order)[1,]==1), sum(apply(BIc[1:3,], 2, order)[2,]==1), sum(apply(BIc[1:3,], 2, order)[3,]==1))
                  ))


n = nn[4]
AIc = 2*(mmm) -2*t(loglik4[,2:4])
BIc = mmm*log(n) -2*t(loglik4[,2:4])

sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(AIc)),
                        apply(AIc, 1, sd),
                        c(sum(apply(AIc[1:3,], 2, order)[1,]==1), sum(apply(AIc[1:3,], 2, order)[2,]==1), sum(apply(AIc[1:3,], 2, order)[3,]==1)),
                        colMeans(t(BIc)),
                        apply(BIc, 1, sd),
                        c(sum(apply(BIc[1:3,], 2, order)[1,]==1), sum(apply(BIc[1:3,], 2, order)[2,]==1), sum(apply(BIc[1:3,], 2, order)[3,]==1))
                  ))

n = nn[2]
AIc = 2*(mmm) -2*t(loglik5[,2:4])
BIc = mmm*log(n) -2*t(loglik5[,2:4])

sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(AIc)),
                        apply(AIc, 1, sd),
                        c(sum(apply(AIc[1:3,], 2, order)[1,]==1), sum(apply(AIc[1:3,], 2, order)[2,]==1), sum(apply(AIc[1:3,], 2, order)[3,]==1)),
                        colMeans(t(BIc)),
                        apply(BIc, 1, sd),
                        c(sum(apply(BIc[1:3,], 2, order)[1,]==1), sum(apply(BIc[1:3,], 2, order)[2,]==1), sum(apply(BIc[1:3,], 2, order)[3,]==1))
                  ))

### MCAR MAR MNAR
sum.tabel.T = cbind(sum.tabel.T, round(sum.tabel, 3)[-c(2,5,8,11,14,17,20,23,26,29),][,c(3,2,1)])

################################################################################
#################
#################  75 %; nu=5
#################
PATH0 <- paste(PATH, "/Data/simulation/SS-simulation-t75/", sep = "")
PATH1 <- paste(PATH0, "SIM1/", sep = "")
PATH2 <- paste(PATH0, "SIM2/", sep = "")
PATH3 <- paste(PATH0, "SIM3/", sep = "")
PATH4 <- paste(PATH0, "SIM4/", sep = "")
PATH5 <- paste(PATH0, "SIM5/", sep = "")

loglik1 =    as.matrix(read.table(paste(PATH1,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE1 =    as.matrix(read.table(paste(PATH1,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara1 = colMeans(as.matrix(read.table(paste(PATH1,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik2 =    as.matrix(read.table(paste(PATH2,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE2 =    as.matrix(read.table(paste(PATH2,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara2 = colMeans(as.matrix(read.table(paste(PATH2,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik3 =    as.matrix(read.table(paste(PATH3,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE3 =    as.matrix(read.table(paste(PATH3,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara3 = colMeans(as.matrix(read.table(paste(PATH3,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik4 =    as.matrix(read.table(paste(PATH4,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE4 =    as.matrix(read.table(paste(PATH4,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara4 = colMeans(as.matrix(read.table(paste(PATH4,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik5 =    as.matrix(read.table(paste(PATH5,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE5 =    as.matrix(read.table(paste(PATH5,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara5 = colMeans(as.matrix(read.table(paste(PATH5,'realpara.txt',sep=""),na.strings="NA",sep="")))

####################
##########
##########   table1
nn = c(25, 50, 100, 200, 400)
aic1 = aic2 = bic1 = bic2 = NULL
mMNAR = length(realpara1)
mMAR = length(realpara1)-1
mMCAR = length(realpara1)-2

mmm = c(mMNAR, mMAR, mMCAR)

n = nn[1]
AIc = 2*(mmm) -2*t(loglik1[,2:4])
BIc = mmm*log(n) -2*t(loglik1[,2:4])

sum.tabel = rbind(colMeans(t(AIc)),
                  apply(AIc, 1, sd),
                  c(sum(apply(AIc[1:3,], 2, order)[1,]==1), sum(apply(AIc[1:3,], 2, order)[2,]==1), sum(apply(AIc[1:3,], 2, order)[3,]==1)),
                  colMeans(t(BIc)),
                  apply(BIc, 1, sd),
                  c(sum(apply(BIc[1:3,], 2, order)[1,]==1), sum(apply(BIc[1:3,], 2, order)[2,]==1), sum(apply(BIc[1:3,], 2, order)[3,]==1))
)


n = nn[2]
AIc = 2*(mmm) -2*t(loglik2[,2:4])
BIc = mmm*log(n) -2*t(loglik2[,2:4])

sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(AIc)),
                        apply(AIc, 1, sd),
                        c(sum(apply(AIc[1:3,], 2, order)[1,]==1), sum(apply(AIc[1:3,], 2, order)[2,]==1), sum(apply(AIc[1:3,], 2, order)[3,]==1)),
                        colMeans(t(BIc)),
                        apply(BIc, 1, sd),
                        c(sum(apply(BIc[1:3,], 2, order)[1,]==1), sum(apply(BIc[1:3,], 2, order)[2,]==1), sum(apply(BIc[1:3,], 2, order)[3,]==1))
                  ))

n = nn[3]
AIc = 2*(mmm) -2*t(loglik3[,2:4])
BIc = mmm*log(n) -2*t(loglik3[,2:4])

sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(AIc)),
                        apply(AIc, 1, sd),
                        c(sum(apply(AIc[1:3,], 2, order)[1,]==1), sum(apply(AIc[1:3,], 2, order)[2,]==1), sum(apply(AIc[1:3,], 2, order)[3,]==1)),
                        colMeans(t(BIc)),
                        apply(BIc, 1, sd),
                        c(sum(apply(BIc[1:3,], 2, order)[1,]==1), sum(apply(BIc[1:3,], 2, order)[2,]==1), sum(apply(BIc[1:3,], 2, order)[3,]==1))
                  ))


n = nn[4]
AIc = 2*(mmm) -2*t(loglik4[,2:4])
BIc = mmm*log(n) -2*t(loglik4[,2:4])

sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(AIc)),
                        apply(AIc, 1, sd),
                        c(sum(apply(AIc[1:3,], 2, order)[1,]==1), sum(apply(AIc[1:3,], 2, order)[2,]==1), sum(apply(AIc[1:3,], 2, order)[3,]==1)),
                        colMeans(t(BIc)),
                        apply(BIc, 1, sd),
                        c(sum(apply(BIc[1:3,], 2, order)[1,]==1), sum(apply(BIc[1:3,], 2, order)[2,]==1), sum(apply(BIc[1:3,], 2, order)[3,]==1))
                  ))

n = nn[2]
AIc = 2*(mmm) -2*t(loglik5[,2:4])
BIc = mmm*log(n) -2*t(loglik5[,2:4])

sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(AIc)),
                        apply(AIc, 1, sd),
                        c(sum(apply(AIc[1:3,], 2, order)[1,]==1), sum(apply(AIc[1:3,], 2, order)[2,]==1), sum(apply(AIc[1:3,], 2, order)[3,]==1)),
                        colMeans(t(BIc)),
                        apply(BIc, 1, sd),
                        c(sum(apply(BIc[1:3,], 2, order)[1,]==1), sum(apply(BIc[1:3,], 2, order)[2,]==1), sum(apply(BIc[1:3,], 2, order)[3,]==1))
                  ))

### MCAR MAR MNAR
sum.tabel.T = cbind(sum.tabel.T, round(sum.tabel, 3)[-c(2,5,8,11,14,17,20,23,26,29),][,c(3,2,1)])

colnames(sum.tabel.T) <- rep(c("25% MCAR", "50% MCAR", "75% MCAR"), 3)
row.names(sum.tabel.T) <- c(rep(c("AIC", "(Freq)", "BIC", "(Freq)"), time = 5))

write.csv(sum.tabel.T, paste0(PATH, "/Result/TableB3.csv"))
