################################################################################
#
#   Filename    :    simSAEM75.R
#   Project     :    Biometrics article "Modeling HIV Viral Dynamics Using a
#                    Nonlinear Mixed-Effects Framework for Heavy-Tailed Data
#                    with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#   Date        :    18.08.2026
#   Purpose     :    reproduce the primary simulation results with 75% dropout
#                    across sample sizes N = 25, 50, 100, 200, and 400
#
#   Input code files  :  Data_and_Code/function/simulation/simulate_dropout_study.r
#                         Data_and_Code/function/simulation/tNLMMmissingSAEM.r
#                         Data_and_Code/function/simulation/NLMMmissingSAEM.r
#
#   Output .txt files :  Data_and_Code/Data/simulation/SS-simulation-t75/SIM1-SIM5/
#
#   R Version   :    R-4.6.0
#   Required R packages : mvtnorm; nlme
#
################################################################################
library(mvtnorm)
library(nlme)

PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

source(file.path(PATH, "function", "simulation", "simulate_dropout_study.R"))
source(file.path(PATH, "function", "simulation", "tNLMMmissingSAEM.R"))
source(file.path(PATH, "function", "simulation", "NLMMmissingSAEM.R"))

# SIM1: n=25
# SIM2: n=50
# SIM3: n=100
# SIM4: n=200
# SIM5: 
n=400

# Parameter setting:
p = 3; si = 10; q = 1; g=1
Beta = matrix(c(27, 2, -0.5), ncol=g)
# DD = array(0, dim=c(q,q,g))
# for(i in 1: g) DD[,,i] = matrix(0.3*i, q, q)
DD = matrix(0, q, q)
DD[1,1] = 2
sigma = 5
nu = 5
Phi = 1e-6

### 25 dropout
# alpha = c(-4.91, -0.69, -0.30,  0.30) ## intercept; treat; y_si-1; y_si
### 50 dropout
# alpha = c(-3.82, -0.69, -0.30,  0.30) ## intercept; treat; y_si-1; y_si
### 75 dropout
alpha = c(-2.96, -0.69, -0.30,  0.30) ## intercept; treat; y_si-1; y_si
para = list(alpha=alpha, Beta=Beta, DD=DD, sigma=sigma, Phi=Phi, nu=nu)

cor.type = "UNC"
gen.Data = gen.tnlmm(n, para, cor.type='UNC', si, q)
Data = gen.Data$Data
Ymat = gen.Data$Ymat
b = gen.Data$b

gen.Data.miss = add.MNAR.tNLMM(gen.Data$Data, gen.Data$Ymat, b, si, para$alpha)
Data.miss = gen.Data.miss$Data.na
Data.miss$R = rep(0, nrow(Data.miss))
Data.miss$R[is.na(Data.miss$Var1)] = 1
Data = Data.miss
table(Data$Time)

# library(ggplot2)
# k1 = ggplot()+
#   geom_line(data = Data, aes(x=Time, y=y.c, group=as.factor(Subject))) + 
#   geom_line(data = Data, aes(x=Time, y=Var1,  group=as.factor(Subject))) + xlab("Time") + ylab(NULL) +
#   theme(legend.title=element_blank())
# k1

# ggplot()+
#   geom_boxplot(data = Data, aes(x=as.factor(Time), y=y.c)) + xlab("Time") + ylab(NULL) +
#   theme(legend.title=element_blank())
table(Data$Time)

Data
table(Data$Subject, Data$R)

length(unique(Data$Subject[Data$R==1]))/n
max.iter = 1000

N <- c(25, 50, 100, 200, 400)

Repp <- 100
binrep = 1
seednum = NULL
set.seed(seednum)
for (Rep in binrep:Repp) {
  for (ll in 1:length(N)) {
    n <- N[ll]
    if (n == 25) PATH1 <- file.path(PATH, "Data", "simulation", "SS-simulation-t75", "SIM1")
    if (n == 50) PATH1 <- file.path(PATH, "Data", "simulation", "SS-simulation-t75", "SIM2")
    if (n == 100) PATH1 <- file.path(PATH, "Data", "simulation", "SS-simulation-t75", "SIM3")
    if (n == 200) PATH1 <- file.path(PATH, "Data", "simulation", "SS-simulation-t75", "SIM4")
    if (n == 400) PATH1 <- file.path(PATH, "Data", "simulation", "SS-simulation-t75", "SIM5")
    set.seed(Rep + n*100)
    cat(paste(c(rep('=', 15), rep(' ', 3), 'The ', Rep, ' time simulation: n = ', n, rep(' ', 3), rep('=', 15)), sep = '', collapse = ''), '\n')
    
    ################################
    ################   The new simulation Data
    ################
    cor.type = "UNC"
    gen.Data = gen.tnlmm(n, para, cor.type='UNC', si, q)
    Data = gen.Data$Data
    Ymat = gen.Data$Ymat
    b = gen.Data$b
    
    gen.Data.miss = add.MNAR.tNLMM(gen.Data$Data, gen.Data$Ymat, gen.Data$b, si, para$alpha)
    Data.miss = gen.Data.miss$Data.na
    Data.miss$R = rep(0, nrow(Data.miss))
    Data.miss$R[is.na(Data.miss$Var1)] = 1
    Data = Data.miss
    
    fm1.nlme = try(nlme(mu.fn.nlme, 
                        fixed = beta1 + beta2 + beta3~1,
                        random = beta1~ 1|Subject, 
                        data = Data[-which(is.na(Data$Var1)),], start = c(Beta)), silent = F)
    if(class(fm1.nlme)[1] == "try-error") next;
    
    alpha1 = as.numeric(prediction_ym(Data, fm1.nlme, "UNC")$alpha.hat)
    Data$yc = prediction_ym(Data, fm1.nlme, "UNC")$yc
    
    
    ###################################################
    ###################################################  Run MNAR
    ###################################################
    # mechanism='MNAR'
    cat(rep('=', 15), "Student's t of Linear mixed modles with MNAR missing", cor.type[1], " errors: ", "\n")
    DD1 = matrix(0, q,q)
    DD1[1,1] = as.numeric(VarCorr(fm1.nlme)[1,1])
    init.para = list(Beta=matrix(c(fm1.nlme$coefficients$fixed), ncol=1), DD=DD1, sigma=c(fm1.nlme$sigma^2), Phi=para$Phi, nu=para$nu, alpha = alpha1)
    est.MNAR = try(tNLMM.miss.SAEM(Data, g, init.para, cor.type = c("UNC"), M=10, M.LL = 1000, P=1, tol = 1e-6, max.iter=max.iter, per=500,
                                   mechanism='MNAR'), silent = F)
    if(class(est.MNAR) == "try-error") next;
    
    fit.MNAR = try(NLMM.miss.SAEM(Data, g, init.para, cor.type = c("UNC"), M=10, M.LL = 1000, P=1, tol = 1e-6, max.iter=max.iter, per=500,
                                  mechanism='MNAR'), silent = F)
    if(class(fit.MNAR) == "try-error") next;
    
    ###################################################
    ###################################################  Run MAR
    ###################################################
    cat(rep('=', 15), "Student's t of Linear mixed modles with MAR missing", cor.type[1], " errors: ", "\n")
    init.para$alpha = alpha[-length(alpha)]
    est.MAR = try(tNLMM.miss.SAEM(Data, g, init.para, cor.type = c("UNC"), M=10, M.LL = 1000, P=1, tol = 1e-6, max.iter=max.iter, per=500,
                                  mechanism='MAR'), silent = F)
    if(class(est.MAR) == "try-error") next;
    fit.MAR = try(NLMM.miss.SAEM(Data, g, init.para, cor.type = c("UNC"), M=10, M.LL = 1000, P=1, tol = 1e-6, max.iter=max.iter, per=500,
                                 mechanism='MAR'), silent = F)
    if(class(fit.MAR) == "try-error") next;
    
    ###################################################
    ###################################################  Run MCAR
    ###################################################
    cat(rep('=', 15), "Student's t of Linear mixed modles with MCAR missing", cor.type[1], " errors: ", "\n")
    init.para$alpha = c(-2.51)
    est.MCAR = try(tNLMM.miss.SAEM(Data, g, init.para, cor.type = c("UNC"), M=10, M.LL = 1000, P=1, tol = 1e-6, max.iter=max.iter, per=500,
                                   mechanism='MCAR'), silent = F)
    if(class(est.MCAR) == "try-error") next;
    fit.MCAR = try(NLMM.miss.SAEM(Data, g, init.para, cor.type = c("UNC"), M=10, M.LL = 1000, P=1, tol = 1e-6, max.iter=max.iter, per=500,
                                  mechanism='MCAR'), silent = F)
    if(class(fit.MCAR) == "try-error") next;
    
    
    # ################################################### ###################################################
    # ###################################################  Run MCEM
    # ###################################################  
    # 
    # ###########
    # #############################  Run MNAR
    # ###########
    # # mechanism='MNAR'
    # cat(rep('=', 15), "Student's t of Linear mixed modles with MNAR missing", cor.type[1], " errors: ", "\n")
    # DD1 = matrix(0, q,q)
    # DD1[1,1] = as.numeric(VarCorr(fm1.nlme)[1,1])
    # init.para = list(Beta=matrix(c(fm1.nlme$coefficients$fixed), ncol=1), DD=DD1, sigma=c(fm1.nlme$sigma^2), Phi=para$Phi, nu=para$nu, alpha = alpha1)
    # est.MNAR.MCEM = try(tNLMM.miss.MCEM(Data, g, init.para, cor.type = c("UNC"), M=1000, M.LL = 1000, P=1, tol = 1e-6, max.iter=max.iter, per=500,
    #                                     mechanism='MNAR'), silent = F)
    # if(class(est.MNAR.MCEM) == "try-error") next;
    # 
    # fit.MNAR.MCEM = try(NLMM.miss.MCEM(Data, g, init.para, cor.type = c("UNC"), M=1000, M.LL = 1000, P=1, tol = 1e-6, max.iter=max.iter, per=500,
    #                                    mechanism='MNAR'), silent = F)
    # if(class(fit.MNAR.MCEM) == "try-error") next;
    # 
    # ###########
    # #############################  Run MAR
    # ###########
    # cat(rep('=', 15), "Student's t of Linear mixed modles with MAR missing", cor.type[1], " errors: ", "\n")
    # init.para$alpha = alpha[-length(alpha)]
    # est.MAR.MCEM = try(tNLMM.miss.MCEM(Data, g, init.para, cor.type = c("UNC"), M=1000, M.LL = 1000, P=1, tol = 1e-6, max.iter=max.iter, per=500,
    #                                    mechanism='MAR'), silent = F)
    # if(class(est.MAR.MCEM) == "try-error") next;
    # fit.MAR.MCEM = try(NLMM.miss.MCEM(Data, g, init.para, cor.type = c("UNC"), M=1000, M.LL = 1000, P=1, tol = 1e-6, max.iter=max.iter, per=500,
    #                                   mechanism='MAR'), silent = F)
    # if(class(fit.MAR.MCEM) == "try-error") next;
    # 
    # ###########
    # #############################  Run MCAR
    # ###########
    # cat(rep('=', 15), "Student's t of Linear mixed modles with MCAR missing", cor.type[1], " errors: ", "\n")
    # init.para$alpha = c(-2.51)
    # est.MCAR.MCEM = try(tNLMM.miss.MCEM(Data, g, init.para, cor.type = c("UNC"), M=1000, M.LL = 1000, P=1, tol = 1e-6, max.iter=max.iter, per=500,
    #                                     mechanism='MCAR'), silent = F)
    # if(class(est.MCAR.MCEM) == "try-error") next;
    # fit.MCAR.MCEM = try(NLMM.miss.MCEM(Data, g, init.para, cor.type = c("UNC"), M=1000, M.LL = 1000, P=1, tol = 1e-6, max.iter=max.iter, per=500,
    #                                    mechanism='MCAR'), silent = F)
    # if(class(fit.MCAR.MCEM) == "try-error") next;
    # 
    # ###################################################  
    # ###################################################  end MCEM
    # ################################################### ###################################################
    
    
    if(Rep == 1){
      para.real = c(para$Beta, para$sigma, para$DD[vech.posi(q)], para$nu, para$alpha)
      write(c(para.real), paste(PATH1, 'realpara.txt',sep=""), ncol=1000, append=T)
    }
    if(Rep == 2){
      para.real = c(para$Beta, para$sigma, para$DD[vech.posi(q)], para$nu, para$alpha)
      write(c(para.real), paste(PATH1, 'realpara.txt',sep=""), ncol=1000, append=T)
    }
    
    MNAR.para.est = MAR.para.est = MCAR.para.est = NULL
    MNAR.para.est = c(est.MNAR$para.est$Beta, est.MNAR$para.est$sigma, c(est.MNAR$para.est$D[vech.posi(q)]), est.MNAR$para.est$nu, c(est.MNAR$para.est$alpha))
    MAR.para.est = c(est.MAR$para.est$Beta, est.MAR$para.est$sigma, c(est.MAR$para.est$D[vech.posi(q)]), est.MAR$para.est$nu, c(est.MAR$para.est$alpha))
    MCAR.para.est = c(est.MCAR$para.est$Beta, est.MCAR$para.est$sigma, c(est.MCAR$para.est$D[vech.posi(q)]), est.MCAR$para.est$nu, c(est.MCAR$para.est$alpha))
    
    MNAR.para.fit = MAR.para.fit = MCAR.para.fit = NULL
    MNAR.para.fit = c(fit.MNAR$para.est$Beta, fit.MNAR$para.est$sigma, c(fit.MNAR$para.est$D[vech.posi(q)]), c(fit.MNAR$para.est$alpha))
    MAR.para.fit = c(fit.MAR$para.est$Beta, fit.MAR$para.est$sigma, c(fit.MAR$para.est$D[vech.posi(q)]), c(fit.MAR$para.est$alpha))
    MCAR.para.fit = c(fit.MCAR$para.est$Beta, fit.MCAR$para.est$sigma, c(fit.MCAR$para.est$D[vech.posi(q)]), c(fit.MCAR$para.est$alpha))
    
    loglik = c(est.MNAR$model.inf$loglik, est.MAR$model.inf$loglik, est.MCAR$model.inf$loglik, fit.MNAR$model.inf$loglik, fit.MAR$model.inf$loglik, fit.MCAR$model.inf$loglik)
    aic = c(est.MNAR$model.inf$aic, est.MAR$model.inf$aic, est.MCAR$model.inf$aic, fit.MNAR$model.inf$aic, fit.MAR$model.inf$aic, fit.MCAR$model.inf$aic)
    bic = c(est.MNAR$model.inf$bic, est.MAR$model.inf$bic, est.MCAR$model.inf$bic, fit.MNAR$model.inf$bic, fit.MAR$model.inf$bic, fit.MCAR$model.inf$bic)
    MSE.y = c(est.MNAR$MSE.y, est.MAR$MSE.y, est.MCAR$MSE.y, fit.MNAR$MSE.y, fit.MAR$MSE.y, fit.MCAR$MSE.y)
    MAE.y = c(est.MNAR$MAE.y, est.MAR$MAE.y, est.MCAR$MAE.y, fit.MNAR$MAE.y, fit.MAR$MAE.y, fit.MCAR$MAE.y)
    MAPE.y = c(est.MNAR$MAPE.y, est.MAR$MAPE.y, est.MCAR$MAPE.y, fit.MNAR$MAPE.y, fit.MAR$MAPE.y, fit.MCAR$MAPE.y)
    iter = c(est.MNAR$iter, est.MAR$iter, est.MCAR$iter, fit.MNAR$iter, fit.MAR$iter, fit.MCAR$iter)
    
    MNAR.se.est = est.MNAR$IM$se
    MAR.se.est = est.MAR$IM$se
    MCAR.se.est = est.MCAR$IM$se
    MNAR.se.fit = fit.MNAR$IM$se
    MAR.se.fit = fit.MAR$IM$se
    MCAR.se.fit = fit.MCAR$IM$se
    
    write(c(Rep, loglik), paste(PATH1,'loglik.txt',sep=""), ncol=10000, append=T)
    write(c(Rep, aic), paste(PATH1,'aic.txt',sep=""), ncol=10000, append=T)
    write(c(Rep, bic), paste(PATH1,'bic.txt',sep=""), ncol=10000, append=T)
    write(c(Rep, MSE.y), paste(PATH1,'MSE.txt',sep=""), ncol=10000, append=T)
    write(c(Rep, MAE.y), paste(PATH1,'MAE.txt',sep=""), ncol=10000, append=T)
    write(c(Rep, MAPE.y), paste(PATH1,'MAPE.txt',sep=""), ncol=10000, append=T)
    write(c(Rep, iter), paste(PATH1,'iter.txt',sep=""), ncol=10000, append=T)
    
    write(c(Rep, MNAR.para.est), paste(PATH1,'MNAR.para.est.txt',sep=""), ncol=10000, append=T)
    write(c(Rep, MAR.para.est), paste(PATH1,'MAR.para.est.txt',sep=""), ncol=10000, append=T)
    write(c(Rep, MCAR.para.est), paste(PATH1,'MCAR.para.est.txt',sep=""), ncol=10000, append=T)
    write(c(Rep, MNAR.para.fit), paste(PATH1,'MNAR.para.fit.txt',sep=""), ncol=10000, append=T)
    write(c(Rep, MAR.para.fit), paste(PATH1,'MAR.para.fit.txt',sep=""), ncol=10000, append=T)
    write(c(Rep, MCAR.para.fit), paste(PATH1,'MCAR.para.fit.txt',sep=""), ncol=10000, append=T)
    
    write(c(Rep, MNAR.se.est), paste(PATH1,'MNAR.se.est.txt',sep=""), ncol=10000, append=T)
    write(c(Rep, MAR.se.est), paste(PATH1,'MAR.se.est.txt',sep=""), ncol=10000, append=T)
    write(c(Rep, MCAR.se.est), paste(PATH1,'MCAR.se.est.txt',sep=""), ncol=10000, append=T)
    write(c(Rep, MNAR.se.fit), paste(PATH1,'MNAR.se.fit.txt',sep=""), ncol=10000, append=T)
    write(c(Rep, MAR.se.fit), paste(PATH1,'MAR.se.fit.txt',sep=""), ncol=10000, append=T)
    write(c(Rep, MCAR.se.fit), paste(PATH1,'MCAR.se.fit.txt',sep=""), ncol=10000, append=T)
    
    write(c(Rep, est.MNAR$model.inf$iter.lnL), paste(PATH1,'MNAR.iter.lnL.est.txt',sep=""), ncol=10000, append=T)
    write(c(Rep, est.MAR$model.inf$iter.lnL), paste(PATH1,'MAR.iter.lnL.est.txt',sep=""), ncol=10000, append=T)
    write(c(Rep, est.MCAR$model.inf$iter.lnL), paste(PATH1,'MCAR.iter.lnL.est.txt',sep=""), ncol=10000, append=T)
    write(c(Rep, fit.MNAR$model.inf$iter.lnL), paste(PATH1,'MNAR.iter.lnL.fit.txt',sep=""), ncol=10000, append=T)
    write(c(Rep, fit.MAR$model.inf$iter.lnL), paste(PATH1,'MAR.iter.lnL.fit.txt',sep=""), ncol=10000, append=T)
    write(c(Rep, fit.MCAR$model.inf$iter.lnL), paste(PATH1,'MCAR.iter.lnL.fit.txt',sep=""), ncol=10000, append=T)
    
    # ################################################### ###################################################
    # ###################################################  Run MCEM
    # ###################################################  
    # MNAR.MCEM.para.est = c(est.MNAR.MCEM$para.est$Beta, est.MNAR.MCEM$para.est$sigma, c(est.MNAR.MCEM$para.est$D[vech.posi(q)]), est.MNAR.MCEM$para.est$nu, c(est.MNAR.MCEM$para.est$alpha))
    # MAR.MCEM.para.est = c(est.MAR.MCEM$para.est$Beta, est.MAR.MCEM$para.est$sigma, c(est.MAR.MCEM$para.est$D[vech.posi(q)]), est.MAR.MCEM$para.est$nu, c(est.MAR.MCEM$para.est$alpha))
    # MCAR.MCEM.para.est = c(est.MCAR.MCEM$para.est$Beta, est.MCAR.MCEM$para.est$sigma, c(est.MCAR.MCEM$para.est$D[vech.posi(q)]), est.MCAR.MCEM$para.est$nu, c(est.MCAR.MCEM$para.est$alpha))
    # 
    # MNAR.MCEM.para.fit = c(fit.MNAR.MCEM$para.est$Beta, fit.MNAR.MCEM$para.est$sigma, c(fit.MNAR.MCEM$para.est$D[vech.posi(q)]), c(fit.MNAR.MCEM$para.est$alpha))
    # MAR.MCEM.para.fit = c(fit.MAR.MCEM$para.est$Beta, fit.MAR.MCEM$para.est$sigma, c(fit.MAR.MCEM$para.est$D[vech.posi(q)]), c(fit.MAR.MCEM$para.est$alpha))
    # MCAR.MCEM.para.fit = c(fit.MCAR.MCEM$para.est$Beta, fit.MCAR.MCEM$para.est$sigma, c(fit.MCAR.MCEM$para.est$D[vech.posi(q)]), c(fit.MCAR.MCEM$para.est$alpha))
    # 
    # loglik.MCEM = c(est.MNAR.MCEM$model.inf$loglik, est.MAR.MCEM$model.inf$loglik, est.MCAR.MCEM$model.inf$loglik, fit.MNAR.MCEM$model.inf$loglik, fit.MAR.MCEM$model.inf$loglik, fit.MCAR.MCEM$model.inf$loglik)
    # aic.MCEM = c(est.MNAR.MCEM$model.inf$aic, est.MAR.MCEM$model.inf$aic, est.MCAR.MCEM$model.inf$aic, fit.MNAR.MCEM$model.inf$aic, fit.MAR.MCEM$model.inf$aic, fit.MCAR.MCEM$model.inf$aic)
    # bic.MCEM = c(est.MNAR.MCEM$model.inf$bic, est.MAR.MCEM$model.inf$bic, est.MCAR.MCEM$model.inf$bic, fit.MNAR.MCEM$model.inf$bic, fit.MAR.MCEM$model.inf$bic, fit.MCAR.MCEM$model.inf$bic)
    # MSE.y.MCEM = c(est.MNAR.MCEM$MSE.y, est.MAR.MCEM$MSE.y, est.MCAR.MCEM$MSE.y, fit.MNAR.MCEM$MSE.y, fit.MAR.MCEM$MSE.y, fit.MCAR.MCEM$MSE.y)
    # MAE.y.MCEM = c(est.MNAR.MCEM$MAE.y, est.MAR.MCEM$MAE.y, est.MCAR.MCEM$MAE.y, fit.MNAR.MCEM$MAE.y, fit.MAR.MCEM$MAE.y, fit.MCAR.MCEM$MAE.y)
    # MAPE.y.MCEM = c(est.MNAR.MCEM$MAPE.y, est.MAR.MCEM$MAPE.y, est.MCAR.MCEM$MAPE.y, fit.MNAR.MCEM$MAPE.y, fit.MAR.MCEM$MAPE.y, fit.MCAR.MCEM$MAPE.y)
    # iter.MCEM = c(est.MNAR.MCEM$iter, est.MAR.MCEM$iter, est.MCAR.MCEM$iter, fit.MNAR.MCEM$iter, fit.MAR.MCEM$iter, fit.MCAR.MCEM$iter)
    # 
    # MNAR.MCEM.se.est = est.MNAR.MCEM$IM$se
    # MAR.MCEM.se.est = est.MAR.MCEM$IM$se
    # MCAR.MCEM.se.est = est.MCAR.MCEM$IM$se
    # MNAR.MCEM.se.fit = fit.MNAR.MCEM$IM$se
    # MAR.MCEM.se.fit = fit.MAR.MCEM$IM$se
    # MCAR.MCEM.se.fit = fit.MCAR.MCEM$IM$se
    # 
    # write(c(Rep, loglik.MCEM), paste(PATH1,'loglik.MCEM.txt',sep=""), ncol=10000, append=T)
    # write(c(Rep, aic.MCEM), paste(PATH1,'aic.MCEM.txt',sep=""), ncol=10000, append=T)
    # write(c(Rep, bic.MCEM), paste(PATH1,'bic.MCEM.txt',sep=""), ncol=10000, append=T)
    # write(c(Rep, MSE.y.MCEM), paste(PATH1,'MSE.MCEM.txt',sep=""), ncol=10000, append=T)
    # write(c(Rep, MAE.y.MCEM), paste(PATH1,'MAE.MCEM.txt',sep=""), ncol=10000, append=T)
    # write(c(Rep, MAPE.y.MCEM), paste(PATH1,'MAPE.MCEM.txt',sep=""), ncol=10000, append=T)
    # write(c(Rep, iter.MCEM), paste(PATH1,'iter.MCEM.txt',sep=""), ncol=10000, append=T)
    # 
    # write(c(Rep, MNAR.MCEM.para.est), paste(PATH1,'MNAR.MCEM.para.est.txt',sep=""), ncol=10000, append=T)
    # write(c(Rep, MAR.MCEM.para.est), paste(PATH1,'MAR.MCEM.para.est.txt',sep=""), ncol=10000, append=T)
    # write(c(Rep, MCAR.MCEM.para.est), paste(PATH1,'MCAR.MCEM.para.est.txt',sep=""), ncol=10000, append=T)
    # write(c(Rep, MNAR.MCEM.para.fit), paste(PATH1,'MNAR.MCEM.para.fit.txt',sep=""), ncol=10000, append=T)
    # write(c(Rep, MAR.MCEM.para.fit), paste(PATH1,'MAR.MCEM.para.fit.txt',sep=""), ncol=10000, append=T)
    # write(c(Rep, MCAR.MCEM.para.fit), paste(PATH1,'MCAR.MCEM.para.fit.txt',sep=""), ncol=10000, append=T)
    # 
    # write(c(Rep, MNAR.MCEM.se.est), paste(PATH1,'MNAR.MCEM.se.est.txt',sep=""), ncol=10000, append=T)
    # write(c(Rep, MAR.MCEM.se.est), paste(PATH1,'MAR.MCEM.se.est.txt',sep=""), ncol=10000, append=T)
    # write(c(Rep, MCAR.MCEM.se.est), paste(PATH1,'MCAR.MCEM.se.est.txt',sep=""), ncol=10000, append=T)
    # write(c(Rep, MNAR.MCEM.se.fit), paste(PATH1,'MNAR.MCEM.se.fit.txt',sep=""), ncol=10000, append=T)
    # write(c(Rep, MAR.MCEM.se.fit), paste(PATH1,'MAR.MCEM.se.fit.txt',sep=""), ncol=10000, append=T)
    # write(c(Rep, MCAR.MCEM.se.fit), paste(PATH1,'MCAR.MCEM.se.fit.txt',sep=""), ncol=10000, append=T)
    # 
    # write(c(Rep, est.MNAR.MCEM$model.inf$iter.lnL), paste(PATH1,'MNAR.MCEM.iter.lnL.est.txt',sep=""), ncol=10000, append=T)
    # write(c(Rep, est.MAR.MCEM$model.inf$iter.lnL), paste(PATH1,'MAR.MCEM.iter.lnL.est.txt',sep=""), ncol=10000, append=T)
    # write(c(Rep, est.MCAR.MCEM$model.inf$iter.lnL), paste(PATH1,'MCAR.MCEM.iter.lnL.est.txt',sep=""), ncol=10000, append=T)
    # write(c(Rep, fit.MNAR.MCEM$model.inf$iter.lnL), paste(PATH1,'MNAR.MCEM.iter.lnL.fit.txt',sep=""), ncol=10000, append=T)
    # write(c(Rep, fit.MAR.MCEM$model.inf$iter.lnL), paste(PATH1,'MAR.MCEM.iter.lnL.fit.txt',sep=""), ncol=10000, append=T)
    # write(c(Rep, fit.MCAR.MCEM$model.inf$iter.lnL), paste(PATH1,'MCAR.MCEM.iter.lnL.fit.txt',sep=""), ncol=10000, append=T)
    # 
    # ###################################################  
    # ###################################################  end MCEM
    # ################################################### ###################################################
    
    drop_rate = length(unique(Data$Subject[Data$R==1]))/n
    write(c(Rep, drop_rate), paste(PATH1,'drop_rate.txt',sep=""), ncol=10000, append=T)
  }
}



