####### t NLME model with missing data #######
tNLMM.miss.SAEM = function(Data, g=g, init.para, cor.type = c("UNC", "ARp", "BAND1", "CS"), M = 100, M.LL = 1000, P=1, tol = 1e-6, max.iter=max.iter, per=1, mechanism=c('MNAR',"MAR",'MCAR'))
{
  begin = proc.time()[1]
  #initial values of parameter
  Beta = init.para$Beta
  DD = init.para$DD
  sigma = init.para$sigma
  alpha = init.para$alpha
  Phi = init.para$Phi
  nu = init.para$nu
  p = length(Beta); q = nrow(DD)
  N = length(unique(Data$Subject))
  na.ind = which(is.na(as.vector(t(Data$Var1))))    ## which time point have missing value
  Data.miss = Data
  y = Data.miss$Var1
  R = Data.miss$R
  
  b = mvtnorm::rmvt(N, sigma = DD, df = nu)
  
  A = diag(p)
  B = matrix(c(1,0,0), ncol=q)
  
  ni = numeric(N)
  for(i in 1: N) ni[i] = length(Data$Subject[Data$Subject==i])
  cumsum.ni = cumsum(ni)
  ni.o = numeric(N)
  for(i in 1: N) ni.o[i] = sum(!is.na(Data$Var1[Data$Subject==i]))
  cumsum.ni.o = cumsum(ni.o)
  cumsum.q = cumsum(rep(q, N))
  
  si = max(ni)
  n = sum(ni)
  
  y.na.ind = unique(Data$Subject[which(Data$R==1)])
  Nm = length(y.na.ind)
  mi = numeric(Nm)
  for(i in 1: Nm) mi[i] = sum(Data.miss$R[Data.miss$Subject == y.na.ind[i]])
  cumsum.na = cumsum(mi)
  num.na = length(na.ind)
  na.idx = as.list(N)
  for(i in 1: N) na.idx[[i]] = NA
  na.idx[[y.na.ind[[1]]]] = 1: cumsum.na[1]
  for(i in 2:Nm) na.idx[[y.na.ind[i]]] = (cumsum.na[i-1]+1): cumsum.na[i]
  
  MU = NULL
  TXtilde = 0
  TZtilde = matrix(0, ncol=N*q, nrow=n)
  TLam = TLam.inv = TCor = TCor.inv = matrix(0, n, n)
  for(i in 1: N){
    if(i == 1)idx1 = 1: cumsum.ni[1]
    else idx1 = (cumsum.ni[i-1]+1): cumsum.ni[i]
    eta.ij = A %*% Beta + B %*% b[i,]
    MU = c(MU, mu.fn(eta.ij, Data$Time[Data$Subject == i]))
    dmu.ij = dmu(eta.ij, Data$Time[Data$Subject == i])
    TXtilde = rbind(TXtilde, dmu.ij %*% A)
    Ztilde = dmu.ij %*% B
    TZtilde[idx1, ((i-1)*q+1): (i*q)] = Ztilde
    Lam = Ztilde %*% DD %*% t(Ztilde)+ sigma*cor.fn(Phi, dim=ni[i], type=cor.type[1], Ti=Data$Time[Data$Subject == i], ga=ga)
    TLam[idx1, idx1] = Lam
    TLam.inv[idx1, idx1] = solve(Lam)
    Cor = cor.fn(Phi, dim=ni[i], type=cor.type[1], Ti=Data$Time[Data$Subject == i], ga=ga)
    TCor[idx1,idx1] = Cor
    TCor.inv[idx1,idx1] = solve(Cor)
  }
  TXtilde = TXtilde[-1, ]
  Ytilde = Data$yc - MU + TXtilde %*% Beta + TZtilde %*% as.vector(t(b))
  
  TLam.oo = TLam[-na.ind, -na.ind]
  TLam.mo = TLam[na.ind, -na.ind]
  TLam.mm = TLam[na.ind, na.ind]
  if(num.na == 1) TLam.mo = t(TLam.mo)
  TLam.oo.inv = matrix(0, ncol=sum(ni.o), nrow = sum(ni.o))
  for(i in 1: N){
    if(i == 1) idx1 = 1: cumsum.ni.o[1]
    else idx1 = (cumsum.ni.o[i-1]+1): cumsum.ni.o[i]
    TLam.oo.inv[idx1,idx1] = solve(TLam.oo[idx1,idx1])
  }
  
  ytilde.o = c(Ytilde[-na.ind,])
  ytilde.m = c(Ytilde[na.ind,])
  no = length(ytilde.o)
  Xtilde.o = TXtilde[-na.ind, ]
  Xtilde.m = TXtilde[na.ind, ]
  
  if(mechanism=='MNAR') cat(rep('=', 25), 'tLMM (MNAR) with ', cor.type, ' errors is fitted...; ', 'missing = ', Nm/N*100, '%', rep('=', 25), sep = '', '\n')
  if(mechanism=='MAR') cat(rep('=', 25), 'tLMM (MAR) with ', cor.type, ' errors is fitted...; ', 'missing = ', Nm/N*100, '%', rep('=', 25), sep = '', '\n')
  if(mechanism=='MCAR') cat(rep('=', 25), 'tLMM (MCAR) with ', cor.type, ' errors is fitted...; ', 'missing = ', Nm/N*100, '%', rep('=', 25), sep = '', '\n')
  
  TO = diag(n)[-na.ind, ]
  TM = diag(n)[na.ind, ]
  if(num.na == 1) TM=t(TM)
  vechD = vech.posi(q)
  
  Sig.mm.o.MC = (TLam.mm - TLam.mo %*% TLam.oo.inv %*% t(TLam.mo))
  
  # observed log-likelihood:
  Xbeta = TXtilde %*% Beta
  
  y.samp = matrix(rep(Ytilde, M), nrow=M, ncol=n, byrow=T)
  y.hat = t(TO) %*% ytilde.o + t(TM) %*% ytilde.m
  y.samp[1, ] = y.hat
  
  mc.size = M
  gamma.h = 1/1:max.iter
  gamma.h = c(rep(1, 5), gamma.h)
  b.hat = sum.b2 = sum.Omega2 = E.hat = wij = J.alpha = S.alpha = 0
  IM.H.mean = IM.SS.mean = IM.S.mean = 0
  
  ## A covance matrix for missing propabality
  if(mechanism=='MNAR')
  {
    V.fun = function(y, Data.miss)
    {
      V=NULL
      for(i in 1:N)
      {
        if(i == 1) idx1 = 1: cumsum.ni[1]
        else idx1 = (cumsum.ni[i-1]+1): cumsum.ni[i]
        Data.i = Data.miss[which(Data.miss$Subject==i),]
        Time = Data.i$Time
        k = length(Time)
        fData.i = y[1:(k-1)]
        fData.i = y[idx1][-k]
        fData.i = c(0, fData.i)
        bData.i = y[idx1]
        treat.i = Data.i$treat
        V.i = cbind(1, treat.i, fData.i, bData.i)
        V = rbind(V, V.i)
      }
      return(V)
    }
    V.fun.i = function(y.i, Data.i)
    {
      k = length(y.i)
      fData.i = y.i[1:(k-1)]
      fData.i = y.i[-k]
      fData.i = c(0, fData.i)
      bData.i = y.i
      treat.i = Data.i$treat
      V = cbind(rep(1, k), treat.i, fData.i, bData.i)
      return(V)
    }
    # print(mechanism)
  }
  if(mechanism=='MAR')
  {
    V.fun = function(y, Data.miss)
    {
      V=NULL
      for(i in 1:N)
      {
        if(i == 1) idx1 = 1: cumsum.ni[1]
        else idx1 = (cumsum.ni[i-1]+1): cumsum.ni[i]
        Data.i = Data.miss[which(Data.miss$Subject==i),]
        Time = Data.i$Time
        k = length(Time)
        fData.i = y[1:(k-1)]
        fData.i = y[idx1][-k]
        fData.i = c(0, fData.i)
        bData.i = y[idx1]
        treat.i = Data.i$treat
        V.i = cbind(1, treat.i, fData.i)
        V = rbind(V, V.i)
      }
      return(V)
    }
    V.fun.i = function(y.i, Data.i)
    {
      k = length(y.i)
      fData.i = y.i[1:(k-1)]
      fData.i = y.i[-k]
      fData.i = c(0, fData.i)
      bData.i = y.i
      treat.i = Data.i$treat
      V = cbind(rep(1, k), treat.i, fData.i)
      return(V)
    }
  }
  
  if(mechanism=='MCAR')
  {
    V.fun = function(y, Data.miss)
    {
      V=NULL
      for(i in 1:N)
      {
        if(i == 1) idx1 = 1: cumsum.ni[1]
        else idx1 = (cumsum.ni[i-1]+1): cumsum.ni[i]
        Data.i = Data.miss[which(Data.miss$Subject==i),]
        Time = Data.i$Time
        k = length(Time)
        V.i = cbind(rep(1, length(idx1)))
        V = rbind(V, V.i)
      }
      return(V)
    }
    V.fun.i = function(y.i, Data.i)
    {
      k = length(y.i)
      V = cbind(rep(1, k))
      return(V)
    }
  }
  V = V.fun(y.hat, Data.miss)
  
  yo.cent = ytilde.o - Xtilde.o %*% Beta
  Delta = t(yo.cent[1: cumsum.ni.o[1], ]) %*% TLam.oo.inv[1: cumsum.ni.o[1], 1: cumsum.ni.o[1]] %*% yo.cent[1: cumsum.ni.o[1], ]
  for(i in 2: N) Delta = c(Delta, t(yo.cent[(cumsum.ni.o[i-1]+1): cumsum.ni.o[i], ]) %*% TLam.oo.inv[(cumsum.ni.o[i-1]+1): cumsum.ni.o[i], (cumsum.ni.o[i-1]+1): cumsum.ni.o[i]] %*% yo.cent[(cumsum.ni.o[i-1]+1): cumsum.ni.o[i], ])
  tau = (nu + ni.o)/(nu+Delta)
  
  #### #### #### #### #### #### #### #### #### #### 
  #### log-likelihood
  mu.mo = Xtilde.m%*%Beta + TLam.mo%*% TLam.oo.inv %*% yo.cent
  Sig.mm.o = (TLam.mm - TLam.mo %*% TLam.oo.inv %*% t(TLam.mo))/rep(tau, time=ni)[na.ind]
  wden = numeric(N)
  for(i in 1:N)
  {
    if(i == 1){
      idx0 = 1: cumsum.ni.o[1]
      idx1 = 1: cumsum.ni[1]
    }else{
      idx0 = (cumsum.ni.o[i-1]+1): cumsum.ni.o[i]
      idx1 = (cumsum.ni[i-1]+1): cumsum.ni[i]
    }
    if(sum(i==y.na.ind)==0){
      pvi = c(exp(as.matrix(V[idx1,]) %*% alpha) / (1 + exp(as.matrix(V[idx1,]) %*% alpha)))
      wden[i] = mvtnorm::dmvt(Ytilde[idx1], TXtilde[idx1, ]%*%Beta, TLam[idx1,idx1], df = nu, log=F) * prod(1-pvi)
    }else{
      na.i = which(y.na.ind==i)
      if(na.i == 1){
        idx2 = 1: cumsum.na[1]
      }else{
        idx2 = (cumsum.na[na.i-1]+1): cumsum.na[na.i]
      }
      y.hat.ll = t(mu.mo[idx2,] + t(mvtnorm::rmvt(M.LL, sigma = as.matrix(Sig.mm.o[idx2,idx2]), df = nu+ni.o[y.na.ind][na.i])))
      dt.ym = mvtnorm::dmvt(as.matrix(y.hat.ll), delta = mu.mo[idx2,], sigma = as.matrix(Sig.mm.o[idx2, idx2]), df = nu+ni.o[y.na.ind][na.i], log = F)
      wden.i = NULL
      for(jj in 1:(M.LL/100))
      {
        y.hat.i = c(Ytilde[idx1][R[idx1]==0], y.hat.ll[jj,])
        V.i = V.fun.i(y.i = y.hat.i, Data.miss[idx1,])
        pvi = c(exp(V.i %*% alpha) / (1 + exp(V.i %*% alpha)))
        if(sum(R[idx1])>1) pvi[(sum(R[idx1]==0)+2):length(pvi)] = 1
        if(sum(R[idx1])==1) pvi[length(pvi)] = 1
        fn = mvtnorm::dmvt(y.hat.i, TXtilde[idx1, ]%*%Beta, TLam[idx1,idx1], df = nu, log=F) * prod((1-pvi)^(1-R[idx1]) * pvi^(R[idx1]))
        wden.i = c(wden.i, fn/dt.ym[jj])
      }
      wden[i] = mean(wden.i)
    }
  }
  loglik.old = iter.lnL = sum(log(wden))
  
  theta.old = c(Beta, DD[vechD], sigma, Phi, nu, alpha)
  iter = 0
  cat(paste(rep("=", 50), sep = "", collapse = ""), "\n")
  cat("t nonlinear mixed models with ", cor.type[1], " errors: ", "\n")
  cat("iter = ", iter, ",\t obs.loglik = ", loglik.old, sep = "", "\n")
  Tpara = theta.old; diff.lnL = 10000; diff = 10000
  Taccept.rate = NULL
  repeat
  {
    iter = iter+1
    ##### ##### ##### ##### ##### ##### ##### ##### 
    ######  SA-Step:
    yo.cent = ytilde.o - Xtilde.o %*% Beta
    Delta = t(yo.cent[1: cumsum.ni.o[1], ]) %*% TLam.oo.inv[1: cumsum.ni.o[1], 1: cumsum.ni.o[1]] %*% yo.cent[1: cumsum.ni.o[1], ]
    for(i in 2: N) Delta = c(Delta, t(yo.cent[(cumsum.ni.o[i-1]+1): cumsum.ni.o[i], ]) %*% TLam.oo.inv[(cumsum.ni.o[i-1]+1): cumsum.ni.o[i], (cumsum.ni.o[i-1]+1): cumsum.ni.o[i]] %*% yo.cent[(cumsum.ni.o[i-1]+1): cumsum.ni.o[i], ])
    tau = (nu + ni.o)/(nu+Delta)
    # tau = rep(1, N)
    kappa = digamma((nu + ni.o)/2) - log((nu + Delta)/2)
    
    TD = kronecker(diag(N), DD)
    TSig.b = matrix(0, N*q, N*q)
    for(i in 1:N)
    {
      if(i == 1)
      {
        idx1 = 1: cumsum.ni[1]
        idx2 = 1: cumsum.q[1]
      }
      else
      {
        idx1 = (cumsum.ni[i-1]+1): cumsum.ni[i]
        idx2 = (cumsum.q[i-1]+1): cumsum.q[i]
      }
      TSig.b[idx2,idx2] = solve(t(TZtilde[idx1,idx2]) %*% TCor.inv[idx1,idx1] %*% t(t(TZtilde[idx1,idx2]))/sigma + solve(DD)) /tau[i]
    }
    
    ##### ##### ##### ##### ##### ##### ##### ##### 
    ##### genrate ym by MCMC
    mu.mo = Xtilde.m%*%Beta + TLam.mo%*% TLam.oo.inv %*% yo.cent
    Sig.mm.o = (TLam.mm - TLam.mo %*% TLam.oo.inv %*% t(TLam.mo))/rep(tau, time=ni)[na.ind]
    y.samp = matrix(rep(y.hat, M), nrow=M, ncol=n, byrow=T)
    for(m in 2: M)
    {
      for(i in y.na.ind)
      {
        if(i == 1) idx1 = 1: cumsum.ni[1]
        else idx1 = (cumsum.ni[i-1]+1): cumsum.ni[i]
        y.samp[m, idx1] = MH.y.miss.t(yi=y.samp[(m-1), idx1], Xbeta.i=c(Xbeta)[idx1], Lam.i = TLam[idx1, idx1], tau.i = tau[i], alpha=alpha, 
                                      Vi=as.matrix(V[idx1,]), R.i = R[idx1], n.i=ni[i], mu.mo=mu.mo, Sig.mm.o=Sig.mm.o, 
                                      Sig.mm.o.MC = Sig.mm.o.MC, na.idx.i=na.idx[[i]], Data.miss.i = Data.miss[idx1,],
                                      V.fun.i=V.fun.i)
      }
      # print(m)
    }
    
    accept.rate.iter = NULL
    for(m in 2:M)
    {
      accept.rate.iter = rbind(accept.rate.iter, y.samp[m,] != y.samp[m-1,])
    }
    Taccept.rate = rbind(Taccept.rate, colMeans(accept.rate.iter))
    y.conv = y.samp
    
    #### 
    Xbeta = TXtilde %*% Beta
    b.hat.SA = matrix(0, nrow = N*q)
    sum.b2.SA = matrix(0, N*q, N*q)
    sum.Omega2.SA = matrix(0, n, n)
    wij.SA = rep(0, n)
    
    for(m in 1:mc.size)
    {
      Yhat.cent = c(y.conv[m, ] - Xbeta)
      V = V.fun(y.conv[m,], Data.miss)
      for(i in 1:N)
      {
        if(i == 1){
          idx1 = 1: cumsum.ni[1]
          idx2 = 1: cumsum.q[1]
        }else{
          idx1 = (cumsum.ni[i-1]+1): cumsum.ni[i]
          idx2 = (cumsum.q[i-1]+1): cumsum.q[i]
        }
        b.hat.i = DD %*% t(TZtilde[idx1,idx2]) %*% TLam.inv[idx1,idx1] %*% Yhat.cent[idx1]
        b.hat.SA[idx2] = b.hat.SA[idx2] + b.hat.i
        sum.b2.SA[idx2,idx2] = sum.b2.SA[idx2,idx2] + 
          (DD %*% t(TZtilde[idx1,idx2]) %*% TLam.inv[idx1,idx1] %*% Yhat.cent[idx1] %*% t(DD %*% t(TZtilde[idx1,idx2]) %*% TLam.inv[idx1,idx1] %*% (Yhat.cent[idx1])))
        sum.Omega2.SA[idx1, idx1] = sum.Omega2.SA[idx1, idx1] + (Yhat.cent[idx1] - TZtilde[idx1,idx2] %*% b.hat.i) %*% t(Yhat.cent[idx1] - TZtilde[idx1,idx2] %*% b.hat.i) + 
          t(t(TZtilde[idx1,idx2])) %*% TSig.b[idx2,idx2] %*% t(TZtilde[idx1,idx2])
        wij.SA[idx1] = wij.SA[idx1] + c(exp(as.matrix(V[idx1,]) %*% alpha) / (1 + exp(as.matrix(V[idx1,]) %*% alpha)))
      }
    }
    
    y.hat = y.hat + gamma.h[iter]*(colMeans(y.conv) - y.hat)
    b.hat = b.hat + gamma.h[iter]*(b.hat.SA/mc.size - b.hat)
    sum.b2 = sum.b2 + gamma.h[iter]*(sum.b2.SA/mc.size + TSig.b - sum.b2)
    E.hat = E.hat + gamma.h[iter]*(sum.Omega2.SA/mc.size - E.hat)
    wij = wij + gamma.h[iter]*(wij.SA/mc.size - wij)
    
    
    ##### ##### ##### ##### ##### ##### ##### ##### 
    #####  CM-Step:
    ### Beta
    k1 = k2 = 0
    for(i in 1:N)
    {
      if(i == 1){
        idx1 = 1: cumsum.ni[1]
        idx2 = 1: cumsum.q[1]
      }else{
        idx1 = (cumsum.ni[i-1]+1): cumsum.ni[i]
        idx2 = (cumsum.q[i-1]+1): cumsum.q[i]
      }
      k1 = k1 + t(tau[i]*TXtilde[idx1,]) %*% TCor.inv[idx1,idx1] %*% TXtilde[idx1,]
      k2 = k2 + t(tau[i]*TXtilde[idx1,]) %*% TCor.inv[idx1,idx1] %*% (y.hat[idx1] - t(t(TZtilde[idx1,idx2])) %*% b.hat[idx2])
    }
    Beta = solve(k1) %*% k2
    
    
    ### D
    Nb2 = 0
    for(i in 1: N) Nb2 = Nb2 + tau[i] * sum.b2[((i-1)*q+1):(i*q), ((i-1)*q+1):(i*q)]
    DD = as.matrix(Nb2 / N)
    # DD = init.para$DD
    
    ### sigma
    Ce = 0
    for(i in 1: N)
    {
      if(i == 1) idx1 = 1: cumsum.ni[1]
      else idx1 = (cumsum.ni[i-1]+1): cumsum.ni[i]
      Ce = Ce + sum(tau[i] * TCor.inv[idx1,idx1] * E.hat[idx1, idx1])
    }
    sigma = Ce/sum(ni)
    # sigma = 5
    
    ### Phi
    if(cor.type == 'UNC'){ Phi=1e-6; ga = 1}
    if(cor.type == "CAR1" | cor.type == 'ARp'){
      ga = 1
      Phi = optim(par = Phi, fn = phiga.t.fn, method = "L-BFGS-B", lower = -1+1e-6, upper = 1-1e-6, 
                  sigma=sigma, E.hat1=E.hat, cumsum.ni=cumsum.ni, tau=tau, N=N, ni=ni, Data=Data, cor.type = cor.type[1])$par
    }
    if(cor.type == 'CS'){
      ga = 1
      Phi = optim(par = Phi, fn = phiga.t.fn, method = "L-BFGS-B", lower = 1e-6, upper = 1-1e-6, 
                  sigma=sigma, E.hat1=E.hat, cumsum.ni=cumsum.ni, tau=tau, N=N, ni=ni, Data=Data, cor.type = cor.type[1])$par
    }
    if(cor.type == 'DEC'){
      par.DEC = optim(par = c(Phi, ga), fn = phiga.t.fn, method = "L-BFGS-B", lower = c(1e-6, 1e-6), upper = c(1-1e-6, Inf), 
                      sigma=sigma, E.hat1=E.hat, cumsum.ni=cumsum.ni, tau=tau, N=N, ni=ni, Data=Data, cor.type = cor.type[1])$par
      Phi = par.DEC[1]; ga = par.DEC[2]
    }
    if(cor.type == 'BAND1'){
      ga = 1
      Phi = optim(par = Phi, fn = phiga.t.fn, method = "L-BFGS-B", lower = -1/2, upper = 1/2, 
                  sigma=sigma, E.hat1=E.hat, cumsum.ni=cumsum.ni, tau=tau, N=N, ni=ni, Data=Data, cor.type = cor.type[1])$par
    }
    
    
    ### nu
    nu = optim(par = nu, fn = nu.Q.fn, method = "L-BFGS-B", lower = 2, upper = 200, kappa = kappa, tau = tau, N=N)$par
    
    ### alpha
    J.alpha = S.alpha = 0
    for(i in 1:N)
    {
      if(i == 1) idx1 = 1: cumsum.ni[1]
      else idx1 = (cumsum.ni[i-1]+1): cumsum.ni[i]
      R.i = R[idx1]
      k3 = k4 = 0
      a.i = V.fun.i(y.hat[idx1], Data.miss[idx1,])
      J.alpha = J.alpha + t(a.i)%*%a.i
      S.alpha = S.alpha + t(a.i) %*% (R.i - wij[idx1])
      
    }
    alpha = alpha + 4*solve(J.alpha) %*% S.alpha
    
    
    # evaluate new log-likelihood
    MU = NULL
    TXtilde = 0
    TZtilde = matrix(0, ncol=N*q, nrow=n)
    TLam = TLam.inv = TCor = TCor.inv = matrix(0, n, n)
    for(i in 1: N){
      if(i == 1){
        idx1 = 1: cumsum.ni[1]
        idx2 = 1: cumsum.q[1]
      }else{
        idx1 = (cumsum.ni[i-1]+1): cumsum.ni[i]
        idx2 = (cumsum.q[i-1]+1): cumsum.q[i]
      }
      eta.ij = A %*% Beta + B %*% b.hat[idx2,]
      MU = c(MU, mu.fn(eta.ij, Data$Time[Data$Subject == i]))
      dmu.ij = dmu(eta.ij, Data$Time[Data$Subject == i])
      TXtilde = rbind(TXtilde, dmu.ij %*% A)
      Ztilde = dmu.ij %*% B
      TZtilde[idx1, ((i-1)*q+1): (i*q)] = Ztilde
      Lam = Ztilde %*% DD %*% t(Ztilde)+ sigma*cor.fn(Phi, dim=ni[i], type=cor.type[1], Ti=Data$Time[Data$Subject == i], ga=ga)
      TLam[idx1, idx1] = Lam
      TLam.inv[idx1, idx1] = solve(Lam)
      Cor = cor.fn(Phi, dim=ni[i], type=cor.type[1], Ti=Data$Time[Data$Subject == i], ga=ga)
      TCor[idx1,idx1] = Cor
      TCor.inv[idx1,idx1] = solve(Cor)
    }
    TXtilde = TXtilde[-1, ]
    
    ytilde.o = c(y.hat[-na.ind])
    ytilde.m = c(y.hat[na.ind])
    Xtilde.o = TXtilde[-na.ind, ]
    Xtilde.m = TXtilde[na.ind, ]
    
    TLam.oo = TLam[-na.ind, -na.ind]
    TLam.mo = TLam[na.ind, -na.ind]
    TLam.mm = TLam[na.ind, na.ind]
    if(num.na == 1) TLam.mo = t(TLam.mo)
    TLam.oo.inv = matrix(0, ncol=sum(ni.o), nrow = sum(ni.o))
    for(i in 1: N){
      if(i == 1) idx1 = 1: cumsum.ni.o[1]
      else idx1 = (cumsum.ni.o[i-1]+1): cumsum.ni.o[i]
      TLam.oo.inv[idx1,idx1] = solve(TLam.oo[idx1,idx1])
    }
    
    V = V.fun(y.hat, Data.miss)
    #### #### #### #### #### #### #### #### #### #### 
    #### log-likelihood
    yo.cent = ytilde.o - Xtilde.o %*% Beta
    Delta = t(yo.cent[1: cumsum.ni.o[1], ]) %*% TLam.oo.inv[1: cumsum.ni.o[1], 1: cumsum.ni.o[1]] %*% yo.cent[1: cumsum.ni.o[1], ]
    for(i in 2: N) Delta = c(Delta, t(yo.cent[(cumsum.ni.o[i-1]+1): cumsum.ni.o[i], ]) %*% TLam.oo.inv[(cumsum.ni.o[i-1]+1): cumsum.ni.o[i], (cumsum.ni.o[i-1]+1): cumsum.ni.o[i]] %*% yo.cent[(cumsum.ni.o[i-1]+1): cumsum.ni.o[i], ])
    tau = (nu + ni.o)/(nu+Delta)
    
    mu.mo = Xtilde.m%*%Beta + TLam.mo%*% TLam.oo.inv %*% yo.cent
    Sig.mm.o = (TLam.mm - TLam.mo %*% TLam.oo.inv %*% t(TLam.mo))/rep(tau, time=ni)[na.ind]
    
    if(abs(diff.lnL) < tol || diff < tol || iter >= max.iter){
      M.LL.i = M.LL
    }else{
      M.LL.i = M.LL/100
    }
    wden = numeric(N)
    for(i in 1:N)
    {
      if(i == 1){
        idx0 = 1: cumsum.ni.o[1]
        idx1 = 1: cumsum.ni[1]
      }else{
        idx0 = (cumsum.ni.o[i-1]+1): cumsum.ni.o[i]
        idx1 = (cumsum.ni[i-1]+1): cumsum.ni[i]
      }
      if(sum(i==y.na.ind)==0){
        pvi = c(exp(as.matrix(V[idx1,]) %*% alpha) / (1 + exp(as.matrix(V[idx1,]) %*% alpha)))
        wden[i] = mvtnorm::dmvt(Ytilde[idx1], TXtilde[idx1, ]%*%Beta, TLam[idx1,idx1], df = nu, log=F) * prod(1-pvi)
      }else{
        na.i = which(y.na.ind==i)
        if(na.i == 1){
          idx2 = 1: cumsum.na[1]
        }else{
          idx2 = (cumsum.na[na.i-1]+1): cumsum.na[na.i]
        }
        y.hat.ll = t(mu.mo[idx2,] + t(mvtnorm::rmvt(M.LL, sigma = as.matrix(Sig.mm.o[idx2,idx2]), df = nu+ni.o[y.na.ind][na.i])))
        dt.ym = mvtnorm::dmvt(as.matrix(y.hat.ll), delta = mu.mo[idx2,], sigma = as.matrix(Sig.mm.o[idx2, idx2]), df = nu+ni.o[y.na.ind][na.i], log = F)
        wden.i = NULL
        for(jj in 1:M.LL.i)
        {
          y.hat.i = c(Ytilde[idx1][R[idx1]==0], y.hat.ll[jj,])
          V.i = V.fun.i(y.i = y.hat.i, Data.miss[idx1,])
          pvi = c(exp(V.i %*% alpha) / (1 + exp(V.i %*% alpha)))
          if(sum(R[idx1])>1) pvi[(sum(R[idx1]==0)+2):length(pvi)] = 1
          if(sum(R[idx1])==1) pvi[length(pvi)] = 1
          fn = mvtnorm::dmvt(y.hat.i, TXtilde[idx1, ]%*%Beta, TLam[idx1,idx1], df = nu, log=F) * prod((1-pvi)^(1-R[idx1]) * pvi^(R[idx1]))
          wden.i = c(wden.i, fn/dt.ym[jj])
        }
        wden[i] = mean(wden.i)
      }
    }
    loglik.new = sum(log(wden))
    
    #### Fisher information
    para.est = list(Beta = Beta, sigma = sigma, D = DD, Phi = Phi, nu=nu, alpha=alpha, b=b.hat)
    IM = I.tlmm.missing(para.est, cor.type=cor.type, X=TXtilde, Z = TZtilde, N=N, cumsum.ni=cumsum.ni, cumsum.q=cumsum.q, p=p, ni=ni, mc.size=mc.size, TLam.inv=TLam.inv,
                        y.conv=y.conv, mechanism=mechanism, dropout.idx=dropout.idx, R=R, n=n, q=q, Data.miss=Data.miss)
    IM.H.mean = IM.H.mean + gamma.h[iter]*(IM$H.mean - IM.H.mean)
    IM.SS.mean = IM.SS.mean + gamma.h[iter]*(IM$SS.mean - IM.SS.mean)
    IM.S.mean = IM.S.mean + gamma.h[iter]*(IM$S.mean - IM.S.mean)
    
    iter.lnL = c(iter.lnL, loglik.new)
    diff.lnL = (loglik.new-loglik.old)
    theta.new = c(Beta, DD[vechD], sigma, Phi, nu, alpha)
    Tpara = rbind(Tpara, theta.new)
    diff = mean(abs((theta.new-theta.old)/theta.old)[-6])
    if(iter%%per == 0) cat('iter = ', iter, ',\t obs.loglik = ', loglik.new, ',\t theta.diff = ', diff, ',\t MSE of y = ', sum((y.hat- Data$y.c)^2), sep = ' ', '\n')
    if(iter%%per == 0) cat('sigma = ', sigma, "DD = ", DD[vechD], "phi = ", Phi, "nu = ", nu, sep = ' ', '\n')
    if(abs(diff.lnL) < tol || diff < tol || iter >= max.iter) break
    loglik.old = loglik.new; theta.old = theta.new
  }
  end = proc.time()[1]
  # Parameter estimation
  cat(rep("=", 20), "t nonlinear mixed models with ", cor.type[1], " errors", rep("=", 20), sep = "", "\n")
  cat("It took", end - begin, "seconds.\n")
  cat("iter = ", iter, ",\t obs.loglik = ", loglik.new, sep = "", "\n")
  cat("Beta =", Beta, "\n")
  cat("sigma =", sigma, "\n")
  cat("D =\n")
  print(DD)
  cat("Phi =", Phi, "\n")
  cat("nu =", nu, "\n")
  cat("alpha =", alpha, "\n")
  para.est = list(Beta = Beta, sigma = sigma, D = DD, Phi = Phi, nu=nu, alpha=alpha, b=b.hat)
  if(cor.type == "UNC"){
    EST = c(Beta, DD[vechD], sigma, nu, alpha)
  }else{
    EST = c(Beta, DD[vechD], sigma, Phi, nu, alpha)
  }
  g1 = q*(q+1)/2 ## D
  g2 = r = 1  ## sigma
  g3 = length(Phi) ## phi
  goo = g1+g2+g3
  ma = length(c(alpha))
  I.theta = - IM.H.mean - IM.SS.mean  + IM.S.mean%*%t(IM.S.mean)
  if(cor.type == 'UNC'){
    # V.theta = solve(I.theta[1:(p+g-1),1:(p+g-1)])
    V.theta = solve(I.theta[-(p+goo),-(p+goo)])
    sd.theta = c(sqrt(diag(V.theta)))
  }
  # if(cor.type == 'CS'){
  #   # V.theta = solve(I.theta[1:(p+g-1),1:(p+g-1)])
  #   V.theta = solve(I.theta[-(p+goo),-(p+goo)])
  #   sd.theta = c(sqrt(diag(V.theta)))
  # }
  if(cor.type == "CAR1" | cor.type == 'ARp'| cor.type == "BAND1" | cor.type=="CS"){
    V.theta = solve(I.theta)
    sd.theta = c(sqrt(diag(V.theta)))
  }
  if(cor.type == 'DEC'){
    V.theta = solve(I.theta)
    sd.theta = sqrt(diag(V.theta))
  }
  out = rbind(EST, c(sd.theta))
  if(cor.type=="UNC") colnames(out) = rep(c("beta", "d", "sigma", "nu", "alpha"), c(p, length(DD[vechD]), 1, 1, ma))
  if(cor.type=="ARp" | cor.type=="CS"| cor.type == "BAND1") colnames(out) = rep(c("beta", "d", "sigma", "phi", "nu", "alpha"), c(p, length(DD[vechD]), 1, 1, 1, ma))
  IM = list(out=out, se = c(sd.theta), I.theta=I.theta, V.theta=V.theta)
  IM2 = I.tlmm.missing(para.est, cor.type=cor.type, X=TXtilde, Z = TZtilde, N=N, cumsum.ni=cumsum.ni, cumsum.q=cumsum.q, p=p, ni=ni, mc.size=mc.size, TLam.inv=TLam.inv,
                      y.conv=y.conv, mechanism=mechanism, dropout.idx=dropout.idx, R=R, n=n, q=q, Data.miss=Data.miss)
  
  
  
  
  ma = length(c(alpha))
  if(cor.type == "UNC")   m = g * (p + 1 + q*(q+1)/2 + 1) + ma
  else m = g * (p + 1 + q*(q+1)/2 + 1) + length(as.vector(Phi)) + ma
  aic = 2 * m - 2 * loglik.new
  bic = m * log(N) - 2 * loglik.new
  cat("aic =", aic, "\n")
  cat("bic =", bic, "\n")
  cat("MSE of y =", sum((y.hat- Data$y.c)^2), "\n")
  cat(paste(rep("=", 50), sep = "", collapse = ""), "\n")
  model.inf = list(loglik = loglik.new, iter.lnL = iter.lnL, aic = aic, bic = bic, time = end - begin)
  return(list(model.inf = model.inf, para.est = para.est, iter=iter, y.c=y.hat, IM=IM, IM2=IM2, 
              MSE.y = sum((y.hat- Data$y.c)^2)/sum(Data$R), MAE.y = sum(abs((y.hat- Data$y.c)))/sum(Data$R), MAPE.y = sum(abs((y.hat- Data$y.c)/Data$y.c))/sum(Data$R),
              na.ind = na.ind, Tpara=Tpara, Taccept.rate = Taccept.rate))
}



I.tlmm.missing = function(para.est, cor.type=cor.type, X, Z, N, cumsum.ni, cumsum.q, p, ni, mc.size, TLam.inv, y.conv, mechanism, dropout.idx, R, n, q, Data.miss)
{
  Beta = para.est$Beta
  DD = para.est$D
  sigma = para.est$sigma
  Phi = para.est$Phi
  nu = para.est$nu
  alpha= para.est$alpha
  vechD = vech.posi(q)
  if(cor.type == "UNC"){
    EST = c(Beta, DD[vechD], sigma, nu, alpha)
  }else{
    EST = c(Beta, DD[vechD], sigma, Phi, nu, alpha)
  }
  ## A covance matrix for missing propabality
  if(mechanism=='MNAR')
  {
    V.fun = function(y, Data.miss)
    {
      V=NULL
      for(i in 1:N)
      {
        if(i == 1) idx1 = 1: cumsum.ni[1]
        else idx1 = (cumsum.ni[i-1]+1): cumsum.ni[i]
        Data.i = Data.miss[which(Data.miss$Subject==i),]
        Time = Data.i$Time
        k = length(Time)
        fData.i = y[1:(k-1)]
        fData.i = y[idx1][-k]
        fData.i = c(0, fData.i)
        bData.i = y[idx1]
        treat.i = Data.i$treat
        V.i = cbind(1, treat.i, fData.i, bData.i)
        V = rbind(V, V.i)
      }
      return(V)
    }
    V.fun.i = function(y.i, Data.i)
    {
      k = length(y.i)
      fData.i = y.i[1:(k-1)]
      fData.i = y.i[-k]
      fData.i = c(0, fData.i)
      bData.i = y.i
      treat.i = Data.i$treat
      V = cbind(rep(1, k), treat.i, fData.i, bData.i)
      return(V)
    }
    # print(mechanism)
  }
  if(mechanism=='MAR')
  {
    V.fun = function(y, Data.miss)
    {
      V=NULL
      for(i in 1:N)
      {
        if(i == 1) idx1 = 1: cumsum.ni[1]
        else idx1 = (cumsum.ni[i-1]+1): cumsum.ni[i]
        Data.i = Data.miss[which(Data.miss$Subject==i),]
        Time = Data.i$Time
        k = length(Time)
        fData.i = y[1:(k-1)]
        fData.i = y[idx1][-k]
        fData.i = c(0, fData.i)
        bData.i = y[idx1]
        treat.i = Data.i$treat
        # V.i = cbind(1, Data.i$prevOI, fData.i)
        V.i = cbind(1, treat.i, fData.i)
        V = rbind(V, V.i)
      }
      return(V)
    }
    V.fun.i = function(y.i, Data.i)
    {
      k = length(y.i)
      fData.i = y.i[1:(k-1)]
      fData.i = y.i[-k]
      fData.i = c(0, fData.i)
      bData.i = y.i
      treat.i = Data.i$treat
      V = cbind(rep(1, k), treat.i, fData.i)
      return(V)
    }
  }
  
  if(mechanism=='MCAR')
  {
    V.fun = function(y, Data.miss)
    {
      V=NULL
      for(i in 1:N)
      {
        if(i == 1) idx1 = 1: cumsum.ni[1]
        else idx1 = (cumsum.ni[i-1]+1): cumsum.ni[i]
        Data.i = Data.miss[which(Data.miss$Subject==i),]
        Time = Data.i$Time
        k = length(Time)
        V.i = cbind(rep(1, length(idx1)))
        V = rbind(V, V.i)
      }
      return(V)
    }
    V.fun.i = function(y.i, Data.i)
    {
      k = length(y.i)
      V = cbind(rep(1, k))
      return(V)
    }
  }
  
  
  ## ## ## ## ## ## ## ## ## ## ## ## 
  ## ## ## Information matrix
  g1 = q*(q+1)/2 ## D
  g2 = r = 1  ## sigma
  g3 = length(Phi) ## phi
  g = g1+g2+g3
  ma = length(alpha)
  dot.L = as.list(matrix(0, N, g))
  
  # dot D
  # Tdot.L = array(0, dim=c(n,n,g))
  # for(l in 1:g1)
  # {
  #   TZDZ = matrix(0, ncol=n, nrow=n)
  #   dot.DD = matrix(0, q, q)
  #   dot.DD[matrix(vechD[l, ], 1)] = dot.DD[matrix(rev(vechD[l, ]), 1)] = 1
  #   TZDZ[1: cumsum.ni[1], 1: cumsum.ni[1]] = Z[1: cumsum.ni[1], 1: cumsum.q[1]] %*% dot.DD %*% t(Z[1: cumsum.ni[1], 1: cumsum.q[1]])
  #   for(i in 2: N)
  #   {
  #     rZj = Z[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.q[i-1]+1): cumsum.q[i]]
  #     TZDZ[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]] = rZj %*% dot.DD %*% t(rZj)
  #   }
  #   Tdot.L[,,l] = TZDZ
  # }
  
  # 0 ~ 300 dot D
  for(l in 1: g1)
  {
    dot.DD = matrix(0, q, q)
    dot.DD[matrix(vechD[l, ], 1)] = dot.DD[matrix(rev(vechD[l, ]), 1)] = 1
    for(i in 1: N) 
    {
      if(i == 1){
        idx1 = 1: cumsum.ni[1]
        idx2 = 1: cumsum.q[1]
      }else{
        idx1 = (cumsum.ni[i-1]+1): cumsum.ni[i]
        idx2 = (cumsum.q[i-1]+1): cumsum.q[i]
      }
      dot.L[[(l-1)*N+i]] = Z[idx1, idx2] %*% dot.DD %*% t(Z[idx1, idx2])
    }
  }
  
  # dot sigma
  # if(cor.type=="UNC"){
  #   for(l in 1:g2)
  #   {
  #     Tdot.L[,,g1+l] = diag(n)
  #   }
  # }else{
  #   for(l in 1:g2)
  #   {
  #     Tdot.L[,,g1+l][1: cumsum.ni[1], 1: cumsum.ni[1]] = cor.fn(Phi, dim=ni[1], type=cor.type[1], Ti=Data$Time[Data$Subject == 1], ga=1)
  #     for(i in 2: N) Tdot.L[,,g1+l][(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]] = cor.fn(Phi, dim=ni[i], type=cor.type[1], Ti=Data$Time[Data$Subject == i], ga=ga)
  #   }
  # }
  
  ### 301 ~ 400 dot sigma
  for(l in 1: g2)
  {
    dot.Sig = matrix(0, r, r)
    dot.Sig = 1
    for(i in 1: N) dot.L[[g1*N+(l-1)*N+i]] = kronecker(dot.Sig, cor.fn(Phi, dim=ni[i], type=cor.type[1], Ti=Data$week[Data$Subject == i], ga=ga))
  } 
  
  # dot Phi
  # if(cor.type=="UNC"){
  #   Tdot.L[,,(g1+g2+1)] = diag(0, sum(ni))
  # }
  # if(cor.type == "CAR1"){
  #   for(l in 1: g3){
  #     Tdot.L[,, (g1+g2+l)][1: cumsum.ni[1], 1: cumsum.ni[1]] =  sigma*DEC.dot.phi(phi=Phi, ga=1, Ti=Data$Time[Data$Subject == 1])
  #     for(i in 2: N) Tdot.L[,, (g1+g2+l)][(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]] =  sigma*DEC.dot.phi(phi=Phi, ga=1, Ti=Data$Time[Data$Subject == i])
  #   }
  # }
  # if(cor.type=="ARp"){
  #   for(l in 1: g3){
  #     Tdot.L[,, (g1+g2+l)][1: cumsum.ni[1], 1: cumsum.ni[1]] =  sigma*Arp.Ci.dot(Phi, dim=ni[1])
  #     for(i in 2: N) Tdot.L[,, (g1+g2+l)][(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]] =  sigma*Arp.Ci.dot(Phi, dim=ni[i])
  #   }
  # }
  # if(cor.type == "BAND1"){
  #   for(l in 1: g3){
  #     Tdot.L[,, (g1+g2+l)][1: cumsum.ni[1], 1: cumsum.ni[1]] =  sigma*BAND1.Ci.dot(Phi, dim=ni[1])
  #     for(i in 2: N) Tdot.L[,, (g1+g2+l)][(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]] =  sigma*BAND1.Ci.dot(Phi, dim=ni[i])
  #   }
  # }
  # if(cor.type == "CS"){
  #   for(l in 1: g3){
  #     Tdot.L[,, (g1+g2+l)][1: cumsum.ni[1], 1: cumsum.ni[1]] =  sigma*cs.Ci.dot(Phi, dim=ni[1])
  #     for(i in 2: N) Tdot.L[,, (g1+g2+l)][(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]] =  sigma*cs.Ci.dot(Phi, dim=ni[i])
  #   }
  # }
  # 
  ### 400~ 500 dot phi
  if(cor.type=="UNC"){
    for(i in 1: N) dot.L[[(g1+g2)*N+i]] = diag(0, ni[i])
  }
  if(cor.type == "CAR1" | cor.type=="DEC"){
    for(i in 1: N) dot.L[[(g1+g2)*N+i]] = sigma*DEC.dot.phi(phi=Phi, ga=ga, Ti=Data$week[Data$Subject == i])
  }
  if(cor.type=="ARp"){
    # for(i in 1: N) dot.L[[(g1+g2)*N+i]] = arp.Ci.dot(Phi, dim=ni[i], l)
    for(i in 1: N) dot.L[[(g1+g2)*N+i]] = sigma*Arp.Ci.dot(Phi, dim=ni[i])
  }
  if(cor.type=="CS"){
    for(i in 1: N) dot.L[[(g1+g2)*N+i]] = sigma*cs.Ci.dot(Phi, dim=ni[i])
  }
  if(cor.type == "BAND1"){
    for(i in 1: N) dot.L[[(g1+g2)*N+i]] = sigma*BAND1.Ci.dot(Phi, dim=ni[i])
  }
  
  
  
  
  H = matrix(0, ncol=(p+g+1), nrow=(p+g+1))
  SS = matrix(0, ncol=(p+g+1+ma), nrow=(p+g+1+ma))
  S = numeric((p+g+1 + ma))
  H.alpha = matrix(0, ma, ma)
  H.gamma = 0
  for(m in (1:mc.size))
  {
    # print(m)
    # if(m%%2 == 0) cat(rep("=", 10), "Fisher information ", m/mc.size *100, " % ", rep("=", 10), sep = "", "\n")
    ys.cent = y.conv[m, ] - X%*%Beta
    Delta = t(ys.cent[1: cumsum.ni[1], ]) %*% TLam.inv[1: cumsum.ni[1], 1: cumsum.ni[1]] %*% ys.cent[1: cumsum.ni[1], ]
    for(i in 2: N) Delta = c(Delta, t(ys.cent[(cumsum.ni[i-1]+1): cumsum.ni[i], ]) %*% TLam.inv[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]] %*% ys.cent[(cumsum.ni[i-1]+1): cumsum.ni[i], ])
    tau = (nu + ni)/(nu + Delta)
    ka1 = 1 / (nu + Delta)
    ka2 = (nu + ni)/(nu + Delta)^2
    ### Hesian matrix
    # H beta beta
    # H[1:p, 1:p] = H[1:p, 1:p] - t(rep(tau, ni) * X) %*% TLam.inv %*% X + 
    #   t(rep(2*tau*ka1, ni) * X) %*% TLam.inv %*% ys.cent %*% t(ys.cent) %*% TLam.inv %*% X
    # H beta beta
    for(i in 1:N)
    {
      if(i == 1) idx1 = 1: cumsum.ni[1]
      else idx1 = (cumsum.ni[i-1]+1): cumsum.ni[i]
      H[1:p, 1:p] = H[1:p, 1:p] - t(rep(tau[i], ni[i]) * X[idx1,]) %*% TLam.inv[idx1,idx1] %*% X[idx1,] + 
        t(rep(2*tau[i]*ka1[i], ni[i]) * X[idx1,]) %*% TLam.inv[idx1,idx1] %*% ys.cent[idx1,] %*% t(ys.cent[idx1,]) %*% TLam.inv[idx1,idx1] %*% X[idx1,]
      
    }
    # H beta nu
    H[1:p, p+g+1] = H[1:p, p+g+1] + t(rep(ka1-ka2, ni) * X) %*% TLam.inv %*% ys.cent
    # H nu nu
    H[p+g+1,p+g+1] = H[p+g+1,p+g+1] + (1/4) * 
      sum((trigamma((nu+ni)/2) - trigamma(nu/2) + 2*ni/nu^2 + 4*Delta*ka1/nu - 2*Delta*ka2*(2*nu + Delta)/(nu^2)))
    
    H1 = matrix(0, ncol=(p+g+1), nrow=(p+g+1))
    # TLinv.dotL = as.list(numeric(g1+g2+g3))
    TLinv.dotL = as.list(numeric(N*g))
    for(s in 1: (g1+g2+g3)){
      # TLinv.dotL[[s]] = TLam.inv %*% Tdot.L[,,s]
      for(i in 1:N)
      {
        if(i == 1) idx1 = 1: cumsum.ni[1]
        else idx1 = (cumsum.ni[i-1]+1): cumsum.ni[i]
        TLinv.dotL[[(s-1)*N+i]] = TLam.inv[idx1, idx1] %*% dot.L[[(s-1)*N+i]]
        # H beta gamma
        H1[1:p, (p+s)] = H1[1:p, (p+s)] + t(rep(ka2[i], ni[i]) * X[idx1,]) %*% TLam.inv[idx1,idx1] %*% ys.cent[idx1,] %*% t(ys.cent[idx1,]) %*% TLinv.dotL[[(s-1)*N+i]] %*% TLam.inv[idx1,idx1] %*% ys.cent[idx1,] -
          t(rep(tau[i], ni[i]) * X[idx1,]) %*% TLinv.dotL[[(s-1)*N+i]] %*% TLam.inv[idx1,idx1] %*% ys.cent[idx1,]
        # H gamma gamma
        for(l in 1: s) H1[(p+s),(p+l)] = H1[(p+s),(p+l)] + sum(diag(TLinv.dotL[[(l-1)*N+i]] %*% TLinv.dotL[[(s-1)*N+i]]))/2 + 
          (t(rep(ka2[i], ni[i]) * ys.cent[idx1])%*% (TLinv.dotL[[(l-1)*N+i]] %*% TLam.inv[idx1,idx1] %*% ys.cent[idx1,] %*% t(ys.cent[idx1,]) %*% TLinv.dotL[[(s-1)*N+i]] %*% TLam.inv[idx1,idx1] %*% ys.cent[idx1,])) / 2 -
          (t(rep(tau[i], ni[i]) * ys.cent[idx1])%*% (TLinv.dotL[[(l-1)*N+i]] %*% TLinv.dotL[[(s-1)*N+i]] %*% TLam.inv[idx1,idx1] + TLinv.dotL[[(s-1)*N+i]] %*% TLinv.dotL[[(l-1)*N+i]] %*% TLam.inv[idx1,idx1]) %*% ys.cent[idx1,]) / 2
        # H gamma nu
        H1[(p+s),p+g+1] = H1[(p+s),p+g+1] + (t(rep((ka1[i] - ka2[i]), ni[i]) * ys.cent[idx1,]) %*% TLinv.dotL[[(s-1)*N+i]] %*% TLam.inv[idx1,idx1] %*% ys.cent[idx1,]) / 2
      }
    }
    H.gamma = H.gamma + H1
    
    
    ### Score vector
    S1 = numeric((p+g+1))
    # S beta
    S1[1:p] = S1[1:p] + t(rep(tau, ni) * X)%*%TLam.inv%*%ys.cent
    # # S gamma
    # for(s in 1: g)
    # {
    #   S1[p+s] = S1[p+s] - sum(diag(TLinv.dotL[[s]]))/2 + (t(rep(tau, ni) * ys.cent) %*% TLinv.dotL[[s]] %*% TLam.inv %*%  ys.cent)/2
    # }
    
    for(s in 1: g){
      for(i in 1: N){
        if(i == 1) idx1 = 1: cumsum.ni[1]
        else idx1 = (cumsum.ni[i-1]+1): cumsum.ni[i]
        TLinv.dotL[[(s-1)*N+i]] = TLam.inv[idx1, idx1] %*% dot.L[[(s-1)*N+i]]
        S1[p+s] = S1[p+s] -0.5*sum(diag(TLinv.dotL[[(s-1)*N+i]])) + (t(tau[i] * ys.cent[idx1]) %*% TLinv.dotL[[(s-1)*N+i]] %*% TLam.inv[idx1,idx1] %*%  ys.cent[idx1])/2
      }}
    
    
    # S nu
    S1[p+g+1] = S1[p+g+1] + sum((digamma((nu + ni)/2) - digamma(nu/2) - ni/nu - log(1 + Delta/(nu)) + tau*Delta/nu) / 2)
    # S alpha
    V = V.fun(y.conv[m,], Data.miss)
    pvi = c(exp(V %*% alpha) / (1 + exp(V %*% alpha)))
    ppvi = pvi * (1-pvi)
    S.alpha = t(V) %*% (R-pvi)
    for(i in 1:n)
    {
      H.alpha = H.alpha + t(t(V[i,])) %*% ppvi[i] %*% t(V[i,])
    }
    
    S = S + c(S1, c(S.alpha))
    ### SS
    # SS beta beta
    SS[1:p, 1:p] = SS[1:p, 1:p] + S1[1:p] %*% t(S1[1:p])
    # SS beta gamma
    SS[1:p, (p+1):(p+g)] = SS[1:p, (p+1):(p+g)] + S1[1:p] %*% t(S1[(p+1):(p+g)])
    # SS beta nu
    SS[1:p, (p+g+1)] = SS[1:p, (p+g+1)] + S1[1:p] %*% t(S1[p+g+1])
    ## SS beta alpha
    SS[1:p, (p+g+1+1):(p+g+1+ma)] = SS[1:p, (p+g+1+1):(p+g+1+ma)] + S1[1:p] %*% t(S.alpha)
    # SS gamma gamma
    SS[(p+1):(p+g), (p+1):(p+g)] = SS[(p+1):(p+g), (p+1):(p+g)] + S1[(p+1):(p+g)] %*% t(S1[(p+1):(p+g)])
    # SS gamma nu
    SS[(p+1):(p+g), (p+g+1)] = SS[(p+1):(p+g), (p+g+1)] + S1[(p+1):(p+g)] %*% t(S1[p+g+1])
    # SS gamma alpha
    SS[(p+1):(p+g), (p+g+1+1):(p+g+1+ma)] = SS[(p+1):(p+g), (p+g+1+1):(p+g+1+ma)] + S1[(p+1):(p+g)]  %*% t(S.alpha)
    # SS nu nu
    SS[p+g+1, p+g+1] = SS[p+g+1, p+g+1] + S1[p+g+1]%*% t(S1[p+g+1])
    # SS nu alpha
    SS[p+g+1, (p+g+1+1):(p+g+1+ma)] = SS[p+g+1, (p+g+1+1):(p+g+1+ma)] + S1[p+g+1] %*% t(S.alpha)
    # SS alpha
    SS[(p+g+1+1):(p+g+1+ma), (p+g+1+1):(p+g+1+ma)] = SS[(p+g+1+1):(p+g+1+ma), (p+g+1+1):(p+g+1+ma)] + (S.alpha %*% t(S.alpha))
  }
  
  for(l in (p+1):(p+g-1)) for(s in (l+1):(p+g)) H.gamma[l, s] = H.gamma[s, l]
  H.gamma[(p+g+1), (p+1):(p+g)] = t(H.gamma[(p+1):(p+g), (p+g+1)])
  H.gamma[(p+1):(p+g), 1:p] = t(H.gamma[1:p, (p+1):(p+g)])
  
  H = H + H.gamma
  H[(p+g+1), 1:p] = t(H[1:p, (p+g+1)])
  H = cbind(H, matrix(0, ncol = ma, nrow = (p+g+1)))
  H = rbind(H, matrix(0, ncol = p+g+1+ma, nrow = ma))
  H[(p+g+1+1):(p+g+1+ma), (p+g+1+1):(p+g+1+ma)] = -H.alpha
  H.mean = H/mc.size
  SS[(p+1):(p+g), 1:p] = t(SS[1:p, (p+1):(p+g)])
  SS[(p+g+1), 1:(p+g)] = SS[1:(p+g), (p+g+1)]
  SS[(p+g+1+1):(p+g+1+ma), 1:(p+g+1)] = t(SS[1:(p+g+1), (p+g+1+1):(p+g+1+ma)])
  SS.mean = SS/mc.size
  S.mean = S/mc.size
  
  # I.theta = - H.mean - SS.mean  + S.mean%*%t(S.mean)
  # if(cor.type == 'UNC'){
  #   # V.theta = solve(I.theta[1:(p+g-1),1:(p+g-1)])
  #   V.theta = solve(I.theta[-(p+g),-(p+g)])
  #   sd.theta = c(sqrt(diag(V.theta)))
  # }
  # if(cor.type == 'CS'){
  #   # V.theta = solve(I.theta[1:(p+g-1),1:(p+g-1)])
  #   V.theta = solve(I.theta[-(p+g),-(p+g)])
  #   sd.theta = c(sqrt(diag(V.theta)))
  # }
  # if(cor.type == "CAR1" | cor.type == 'ARp'| cor.type == "BAND1"){
  #   V.theta = solve(I.theta)
  #   sd.theta = c(sqrt(diag(V.theta)))
  # }
  # if(cor.type == 'DEC'){
  #   V.theta = solve(I.theta)
  #   sd.theta = sqrt(diag(V.theta))
  # }
  # out = rbind(EST, c(sd.theta))
  # if(cor.type=="UNC") colnames(out) = rep(c("beta", "d", "sigma", "nu", "alpha"), c(p, length(DD[vechD]), 1, 1, ma))
  # if(cor.type=="ARp" | cor.type=="CS"| cor.type == "BAND1") colnames(out) = rep(c("beta", "d", "sigma", "phi", "nu", "alpha"), c(p, length(DD[vechD]), 1, 1, 1, ma))
  # # SD = list(out=out, se = c(sd.theta), I.theta=I.theta, V.theta=V.theta)
  SD = list(H.mean=H.mean, SS.mean=SS.mean, S.mean=S.mean)
  # }
  return(SD)
}
