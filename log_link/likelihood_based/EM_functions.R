log_lik = expression(-log(2*pi) - log((sigmaU2 + sU2)*(sigmaV2 + sV2)-rho^2*sigmaU2*sigmaV2)/2 - ((sigmaV2 + sV2)*(hUm - muU)^2 - 2*rho*sqrt(sigmaU2*sigmaV2)*(hUm-muU)*(hVm-muV) + (sigmaU2 + sU2)*(hVm-muV)^2)/(2*((sigmaU2 + sU2)*(sigmaV2 + sV2)-rho^2*sigmaU2*sigmaV2)))

#### First Derivatives ####
s_muU = D(log_lik,'muU')
s_muV = D(log_lik,'muV')
s_sigmaU2 = D(log_lik,'sigmaU2')
s_sigmaV2 = D(log_lik,'sigmaV2')
s_rho = D(log_lik,'rho')

#### Second Derivatives ####
s_muU2 = D(s_muU,'muU'); s_muU_muV = D(s_muU,'muV')
s_muU_sigmaU2 = D(s_muU,'sigmaU2'); s_muU_sigmaV2 = D(s_muU,'sigmaV2')
s_muU_rho = D(s_muU,'rho')
s_muV2 = D(s_muV,'muV'); s_muV_sigmaU2 = D(s_muV,'sigmaU2')
s_muV_sigmaV2 = D(s_muV,'sigmaV2'); s_muV_rho = D(s_muV,'rho')
s_sigmaU22 = D(s_sigmaU2,'sigmaU2'); s_sigmaU2_sigmaV2 = D(s_sigmaU2,'sigmaV2')
s_sigmaU2_rho = D(s_sigmaU2,'rho')
s_sigmaV22 = D(s_sigmaV2,'sigmaV2'); s_sigmaV2_rho = D(s_sigmaV2,'rho')
s_rho2 = D(s_rho,'rho')



#### Valid form of log-likelihood ####
#### Generate Data ####
# muU0 = -3.117; muV0 = -2.091; sigmaU20 = 0.7941; sigmaV20 = 1.1237; rho0 = 0.938
# M = 100; s = 0.3; ds = 0.05
# set.seed(1)
# sU2 = rnorm(M,s,ds); sV2 = rnorm(M,s,ds)
# Um = rnorm(M, muU0,sqrt(sigmaU20)); Vm = rnorm(M, muV0,sqrt(sigmaV20))
# hUm = Um + rnorm(M,0,sqrt(sU2)); hVm = Vm + rnorm(M,0,sqrt(sV2))
# library('mvtnorm')
# x = c(hUm[1],hVm[1]); mux = c(muU0,muV0)
# cov = rho0*sqrt(sigmaU20*sigmaV20)
# Sigmax = matrix(c(sigmaU20+ sU2[1], cov, cov, sigmaV20+ sV2[1]),nrow=2)
# dmvnorm(x,mean= mux,sigma = Sigmax,log=T)
# muU = muU0; muV = muV0; sigmaU2 = sigmaU20; sigmaV2 = sigmaV20; rho = rho0
# eval(log_lik)[1]

#### Generate Data ####
gen_data <- function(para0, sU2, sV2){
  muU0 = para0[1]; muV0 = para0[2]; sigmaU20 = para0[3]; sigmaV20 = para0[4]; rho0 = para0[5]
  M = length(sU2)
  Vm = rnorm(M, muV0,sqrt(sigmaV20)); Um = rnorm(M, muU0,sqrt(sigmaU20*(1-rho0^2))) + rho0*sqrt(sigmaU20/sigmaV20)*(Vm-muV0); 
  hUm = Um + rnorm(M,0,sqrt(sU2)); hVm = Vm + rnorm(M,0,sqrt(sV2))
  return(list(hUm = hUm, hVm = hVm, Um = Um, Vm = Vm))
}
gen_data_inc <- function(para0, M){
  muU0 = para0[1]; muV0 = para0[2]; sigmaU20 = para0[3]; sigmaV20 = para0[4]; rho0 = para0[5]
  Vm = rnorm(M, muV0,sqrt(sigmaV20)); Um = rnorm(M, muU0,sqrt(sigmaU20*(1-rho0^2))) + rho0*sqrt(sigmaU20/sigmaV20)*(Vm-muV0)
  n = round(runif(M,2000,5000))
  nU = rbinom(M,n,exp(Um)/(1+exp(Um))); nV = rbinom(M,n,exp(Vm)/(1+exp(Vm)))
  hlambdaU = nU/n; hlambdaV = nV/n
  # hUm = log(hlambdaU/(1-hlambdaU)); hVm = log(hlambdaV/(1-hlambdaV))
  # sU2 = 1/(hlambdaU*(1-hlambdaU)*n); sV2 = 1/(hlambdaV*(1-hlambdaV)*n)
  hUm = log(hlambdaU); hVm = log(hlambdaV)
  sU2 = (1-hlambdaU)/(hlambdaU*n); sV2 =(1-hlambdaV)/(hlambdaV*n)
  
  ## Exclude those with zero incidence estimates ##
  Inc = (hlambdaU>0)
  hUm = hUm[Inc]; hVm = hVm[Inc]; Um = Um[Inc]; Vm = Vm[Inc]; sU2 = sU2[Inc]; sV2 = sV2[Inc]
  return(list(hUm = hUm, hVm = hVm, Um = Um, Vm = Vm, sU2 = sU2, sV2 = sV2))
}

solve_para <- function(data,sU2,sV2){
  hUm = data$hUm; hVm = data$hVm
  muU = mean(hUm); muV = mean(hVm); sigmaU2 = var(hUm)-mean(sU2); sigmaV2 = var(hVm)-mean(sV2)
  rho = cov(hUm,hVm)/sqrt(sigmaU2*sigmaV2)
  d = 1
  for (iter in 1:1000){
    if (d >= 0.001){
      muU_old = muU; muV_old = muV; sigmaU2_old = sigmaU2; sigmaV2_old = sigmaV2; rho_old = rho
      Score = rep(0,5); Hessian = matrix(0,nrow = 5, ncol = 5)
      Score[1] = sum(eval(s_muU)); Score[2] = sum(eval(s_muV))
      Score[3] = sum(eval(s_sigmaU2)); Score[4] = sum(eval(s_sigmaV2))
      Score[5] = sum(eval(s_rho))
      Hessian[1,1] = sum(eval(s_muU2))
      Hessian[1,2] = Hessian[2,1] = sum(eval(s_muU_muV))
      Hessian[1,3] = Hessian[3,1] = sum(eval(s_muU_sigmaU2))
      Hessian[1,4] = Hessian[4,1] = sum(eval(s_muU_sigmaV2))
      Hessian[1,5] = Hessian[5,1] = sum(eval(s_muU_rho))
      Hessian[2,2] = sum(eval(s_muV2))
      Hessian[2,3] = Hessian[3,2] = sum(eval(s_muV_sigmaU2))
      Hessian[2,4] = Hessian[4,2] = sum(eval(s_muV_sigmaV2))
      Hessian[2,5] = Hessian[5,2] = sum(eval(s_muV_rho))
      Hessian[3,3] = sum(eval(s_sigmaU22))
      Hessian[3,4] = Hessian[4,3] = sum(eval(s_sigmaU2_sigmaV2))
      Hessian[3,5] = Hessian[5,3] = sum(eval(s_sigmaU2_rho))
      Hessian[4,4] = sum(eval(s_sigmaV22))
      Hessian[4,5] = Hessian[5,4] = sum(eval(s_sigmaV2_rho))
      Hessian[5,5] = sum(eval(s_rho2))
      m = max(mean(abs(Score))/20,1)
      step = - solve(Hessian)%*%Score/m
      if ((sigmaU2 + step[3]>0)&(sigmaU2 + step[3] <= 10)) sigmaU2 = sigmaU2 + step[3] else {
        if (sigmaU2 + step[3] > 10)  sigmaU2 = 10 else sigmaU2 = sigmaU2/2
      }
      if ((sigmaV2 + step[4]>0)&(sigmaV2 + step[4] <= 10)) sigmaV2 = sigmaV2 + step[4] else {
        if (sigmaV2 + step[4] > 10)  sigmaV2 = 10 else sigmaV2 = sigmaV2/2
      }
      # if (sigmaV2 + step[4]>0) sigmaV2 = sigmaV2 + step[4] else sigmaV2 = sigmaV2/2
      if ((rho + step[5]>-1)&(rho + step[5]<=1)) rho = rho + step[5] else {
        if (rho + step[5]<=-1)  rho = -1 else rho = 1
      }
      
      muU = muU + step[1]; muV = muV + step[2]
      # sigmaU2 = sigmaU2 + step[3]; sigmaV2 = sigmaV2 + step[4]; rho = rho + step[5]
      d = abs(muU-muU_old) + abs(muV - muV_old) + abs(sigmaU2 - sigmaU2_old) + abs(sigmaV2 -sigmaV2_old) + abs(rho - rho_old)
      # print(c(d, mean(eval(log_lik))))
      # print(c(muU,muV,sigmaU2,sigmaV2,rho))
    }
  }
  sd_est = sqrt(diag(-solve(Hessian)))
  para_est = c(muU,muV,sigmaU2,sigmaV2,rho)
  return(list(est = para_est, sd = sd_est,cov = (d<0.001)))
}

solve_para_EM <- function(data,sU2,sV2){
  hUm = data$hUm; hVm = data$hVm
  muU = mean(hUm); muV = mean(hVm); sigmaU2 = var(hUm)-mean(sU2); sigmaV2 = var(hVm)-mean(sV2)
  rho = cov(hUm,hVm)/sqrt(sigmaU2*sigmaV2)
  M = length(hUm)
  d = 1
  for (iter in 1:1000){
    if (d >= 0.001){
      muU_old = muU; muV_old = muV; sigmaU2_old = sigmaU2; sigmaV2_old = sigmaV2; rho_old = rho
      #### E-step ####
      A = matrix(c(sigmaU2,rho*sqrt(sigmaU2*sigmaV2),rho*sqrt(sigmaU2*sigmaV2),sigmaV2), nrow = 2)
      Mean_U = Mean_V = var_U = var_V = cov = 0
      for (m in 1:M){
        Bm = A; Bm[1,1] = Bm[1,1] + sU2[m]; Bm[2,2] = Bm[2,2] + sV2[m]
        mean_cond = c(muU,muV) + A%*%solve(Bm)%*%c(hUm[m]-muU,hVm[m]-muV)
        Mean_U = Mean_U + mean_cond[1]; Mean_V = Mean_V + mean_cond[2]
        var_cond = A - A%*%solve(Bm)%*%A
        var_U = var_U + mean_cond[1]^2 + var_cond[1,1] - 2*mean_cond[1]*muU + muU^2
        var_V = var_V + mean_cond[2]^2 + var_cond[2,2] - 2*mean_cond[2]*muV + muV^2
        cov = cov + mean_cond[1]*mean_cond[2] + var_cond[1,2] - mean_cond[1]*muV - mean_cond[2]*muU + muU*muV
      }
      ### M-step ###
      muU = Mean_U/M; muV = Mean_V/M
      sigmaU2 = var_U/M; sigmaV2 = var_V/M; rho = cov/M/sqrt(sigmaU2* sigmaV2)
      d = abs(muU-muU_old) + abs(muV - muV_old) + abs(sigmaU2 - sigmaU2_old) + abs(sigmaV2 -sigmaV2_old) + abs(rho - rho_old)
      # print(c(d, mean(eval(log_lik))))
      # print(c(muU,muV,sigmaU2,sigmaV2,rho))
    }
  }
  
  Hessian = matrix(0,nrow = 5, ncol = 5)
  Hessian[1,1] = sum(eval(s_muU2))
  Hessian[1,2] = Hessian[2,1] = sum(eval(s_muU_muV))
  Hessian[1,3] = Hessian[3,1] = sum(eval(s_muU_sigmaU2))
  Hessian[1,4] = Hessian[4,1] = sum(eval(s_muU_sigmaV2))
  Hessian[1,5] = Hessian[5,1] = sum(eval(s_muU_rho))
  Hessian[2,2] = sum(eval(s_muV2))
  Hessian[2,3] = Hessian[3,2] = sum(eval(s_muV_sigmaU2))
  Hessian[2,4] = Hessian[4,2] = sum(eval(s_muV_sigmaV2))
  Hessian[2,5] = Hessian[5,2] = sum(eval(s_muV_rho))
  Hessian[3,3] = sum(eval(s_sigmaU22))
  Hessian[3,4] = Hessian[4,3] = sum(eval(s_sigmaU2_sigmaV2))
  Hessian[3,5] = Hessian[5,3] = sum(eval(s_sigmaU2_rho))
  Hessian[4,4] = sum(eval(s_sigmaV22))
  Hessian[4,5] = Hessian[5,4] = sum(eval(s_sigmaV2_rho))
  Hessian[5,5] = sum(eval(s_rho2))
  sd_est = sqrt(diag(-solve(Hessian)))
  para_est = c(muU,muV,sigmaU2,sigmaV2,rho)
  return(list(est = para_est, sd = sd_est,covergence = (d<0.001), cov = -solve(Hessian)))
}


log.lik_incidences<-function(par,est_inc_hiv,est_inc_rgc,py_hiv,py_rgc)
{
  mu_u<-par[1]
  mu_v<-par[2]
  sigma2_u<-par[3]
  sigma2_v<-par[4]
  rho<-par[5]
  
  n<-length(est_inc_hiv)
  
  # est_u<- logit(est_inc_hiv)
  # est_v<- logit(est_inc_rgc)
  # se2_u<- 1 / (est_inc_hiv * (1-est_inc_hiv) * py_hiv)
  # se2_v<- 1 / (est_inc_rgc * (1-est_inc_rgc) * py_rgc)
  est_u<- log(est_inc_hiv)
  est_v<- log(est_inc_rgc)
  se2_u<- (1-est_inc_hiv) / (est_inc_hiv * py_hiv)
  se2_v<- (1-est_inc_rgc) / (est_inc_rgc * py_rgc)
  
  rho_i <- rho * sqrt( sigma2_u * sigma2_v / (sigma2_u + se2_u) / (sigma2_v + se2_v) )
  
  log.lik<-0
  for(i in 1:n)
  {
    log.lik <- log.lik + dmvnorm( c(est_u[i],est_v[i]), mean = c(mu_u,mu_v), 
                                  sigma = matrix(c(sigma2_u + se2_u[i],
                                                   rho_i[i] * sqrt( (sigma2_u + se2_u[i])*(sigma2_v + se2_v[i]) ),
                                                   rho_i[i] * sqrt( (sigma2_u + se2_u[i])*(sigma2_v + se2_v[i]) ),
                                                   sigma2_v + se2_v[i]),2,2),
                                  log = TRUE)
  }
  return(log.lik)
}

sim.incidence<-function(n,log.mu,log.sigma, cohort.sizes = c(200,5000))
{
  # inc.rand<-expit(mvrnorm(n,mu=logit.mu,Sigma=logit.sigma))
  inc.rand<-exp(mvrnorm(n,mu=log.mu,Sigma=log.sigma))
  inc.rand[which(inc.rand>=1)] <- 0.9999
  
  # apply(inc.rand,2,range)
  #apply(inc.rand,2,quantile,c(.1,.9))
  
  #PY_HIV<-floor(runif(n,50,1000))
  PY_HIV<-floor(runif(n,cohort.sizes[1],cohort.sizes[2]))
  PY_RGC <- floor(PY_HIV * rnorm(n,1,.1))
  
  meta.sample<-data.frame(cbind(inc.rand,PY_HIV,PY_RGC))
  colnames(meta.sample)<-c("lambda.HIV","lambda.RGC","PY.HIV","PY.RGC")
  
  case.HIV<-case.RGC<-rep(NA,n)
  est_inc.HIV<-est_inc.RGC<-rep(NA,n)
  se_inc.HIV<-se_inc.RGC<-rep(NA,n)
  for(i in 1:n)
  {
    case.HIV[i] <- rbinom(1,meta.sample$PY.HIV[i],meta.sample$lambda.HIV[i])
    est_inc.HIV[i] <- case.HIV[i]/meta.sample$PY.HIV[i]
    se_inc.HIV[i] <- sqrt(est_inc.HIV[i] * (1-est_inc.HIV[i])/meta.sample$PY.HIV[i])
    
    case.RGC[i] <- rbinom(1,meta.sample$PY.RGC[i],meta.sample$lambda.RGC[i])
    est_inc.RGC[i] <- case.RGC[i]/meta.sample$PY.RGC[i]
    se_inc.RGC[i] <- sqrt(est_inc.RGC[i] * (1-est_inc.RGC[i])/meta.sample$PY.RGC[i])
  }
  
  meta.data <- cbind(meta.sample,case.HIV,est_inc.HIV,se_inc.HIV,case.RGC,est_inc.RGC,se_inc.RGC)
  return(meta.data)
}

sim.incidence.rho_m<-function(n,log.mu,log.sigma, cohort.sizes = c(200,5000), rho_m.range = c(.4,.5))
{
  # inc.rand<-expit(mvrnorm(n,mu=logit.mu,Sigma=logit.sigma))
  inc.rand<-exp(mvrnorm(n,mu=log.mu,Sigma=log.sigma))
  inc.rand[which(inc.rand>=1)] <- 0.9999
  
  #apply(inc.rand,2,range)
  #apply(inc.rand,2,quantile,c(.1,.9))
  
  #PY_HIV<-floor(runif(n,50,1000))
  PY_HIV<-floor(runif(n,cohort.sizes[1],cohort.sizes[2]))
  PY_RGC <- floor(PY_HIV * rnorm(n,1,.1))
  
  meta.sample<-data.frame(cbind(inc.rand,PY_HIV,PY_RGC))
  colnames(meta.sample)<-c("lambda.HIV","lambda.RGC","PY.HIV","PY.RGC")
  
  case.HIV<-case.RGC<-rep(NA,n)
  est_inc.HIV<-est_inc.RGC<-rep(NA,n)
  se_inc.HIV<-se_inc.RGC<-rep(NA,n)
  for(i in 1:n)
  {
    # case.HIV[i] <- rbinom(1,meta.sample$PY.HIV[i],meta.sample$lambda.HIV[i])
    # est_inc.HIV[i] <- case.HIV[i]/meta.sample$PY.HIV[i]
    # se_inc.HIV[i] <- sqrt(est_inc.HIV[i] * (1-est_inc.HIV[i])/meta.sample$PY.HIV[i])
    # 
    # case.RGC[i] <- rbinom(1,meta.sample$PY.RGC[i],meta.sample$lambda.RGC[i])
    # est_inc.RGC[i] <- case.RGC[i]/meta.sample$PY.RGC[i]
    # se_inc.RGC[i] <- sqrt(est_inc.RGC[i] * (1-est_inc.RGC[i])/meta.sample$PY.RGC[i])
    
    # logit.mu.i<-logit(c(meta.sample$lambda.HIV[i],meta.sample$lambda.RGC[i]))
    log.mu.i<-log(c(meta.sample$lambda.HIV[i],meta.sample$lambda.RGC[i]))
    
    # se2_u.i<-1/(meta.sample$lambda.HIV[i] * (1-meta.sample$lambda.HIV[i]) * meta.sample$PY.HIV[i])
    # se2_v.i<-1/(meta.sample$lambda.RGC[i] * (1-meta.sample$lambda.RGC[i]) * meta.sample$PY.RGC[i])
    se2_u.i<-(1-meta.sample$lambda.HIV[i])/(meta.sample$lambda.HIV[i] * meta.sample$PY.HIV[i])
    se2_v.i<-(1-meta.sample$lambda.RGC[i])/(meta.sample$lambda.RGC[i] * meta.sample$PY.RGC[i])
    
    rho.i<-runif(1,rho_m.range[1],rho_m.range[2])
    
    # logit.sigma.i<-matrix(c(se2_u.i,rho.i*sqrt(se2_u.i * se2_v.i),rho.i*sqrt(se2_u.i * se2_v.i),se2_v.i),2,2)
    log.sigma.i<-matrix(c(se2_u.i,rho.i*sqrt(se2_u.i * se2_v.i),rho.i*sqrt(se2_u.i * se2_v.i),se2_v.i),2,2)
    
    # est.inc.i<-expit(mvrnorm(1,mu=logit.mu.i,Sigma=logit.sigma.i))
    est.inc.i<-exp(mvrnorm(1,mu=log.mu.i,Sigma=log.sigma.i))
    est_inc.HIV[i]<-est.inc.i[1]
    se_inc.HIV[i] <- sqrt(est_inc.HIV[i] * (1-est_inc.HIV[i])/meta.sample$PY.HIV[i])
    est_inc.RGC[i]<-est.inc.i[2]
    se_inc.RGC[i] <- sqrt(est_inc.RGC[i] * (1-est_inc.RGC[i])/meta.sample$PY.RGC[i])
  }
  
  meta.data <- cbind(meta.sample,case.HIV,est_inc.HIV,se_inc.HIV,case.RGC,est_inc.RGC,se_inc.RGC)
  return(meta.data)
}

# logit<-function(x){log(x/(1-x))}
# expit<-function(x){ifelse(x>20,1,exp(x - log(1 + exp(x))))}


a.coef = expression(
  muU.est - rho.est * sqrt( sigmaU2.est/sigmaV2.est ) * muV.est
)

#### First Derivatives ####
a.coef.s_muU = D(a.coef,'muU.est')
a.coef.s_muV = D(a.coef,'muV.est')
a.coef.s_sigmaU2 = D(a.coef,'sigmaU2.est')
a.coef.s_sigmaV2 = D(a.coef,'sigmaV2.est')
a.coef.s_rho = D(a.coef,'rho.est')

b.coef = expression(
  rho.est * sqrt( sigmaU2.est / sigmaV2.est )
)

#### First Derivatives ####
b.coef.s_muU = D(b.coef,'muU.est')
b.coef.s_muV = D(b.coef,'muV.est')
b.coef.s_sigmaU2 = D(b.coef,'sigmaU2.est')
b.coef.s_sigmaV2 = D(b.coef,'sigmaV2.est')
b.coef.s_rho = D(b.coef,'rho.est')


a.star = expression(
  muU.est - rho.est * sqrt( sigmaU2.est * sigmaV2.est ) * muV.est / (sigmaV2.est + sV2.est)
)

#### First Derivatives ####
a.star.s_muU = D(a.star,'muU.est')
a.star.s_muV = D(a.star,'muV.est')
a.star.s_sigmaU2 = D(a.star,'sigmaU2.est')
a.star.s_sigmaV2 = D(a.star,'sigmaV2.est')
a.star.s_rho = D(a.star,'rho.est')
a.star.s_sV2 = D(a.star,'sV2.est')

b.star = expression(
  rho.est * sqrt( sigmaU2.est * sigmaV2.est ) / (sigmaV2.est + sV2.est)
)

#### First Derivatives ####
b.star.s_muU = D(b.star,'muU.est')
b.star.s_muV = D(b.star,'muV.est')
b.star.s_sigmaU2 = D(b.star,'sigmaU2.est')
b.star.s_sigmaV2 = D(b.star,'sigmaV2.est')
b.star.s_rho = D(b.star,'rho.est')
b.star.s_sV2 = D(b.star,'sV2.est')

# PE.formula = expression( 1 - exp(u_a - log(1 + exp(u_a)))/exp(u_0 - log(1 + exp(u_0))))
PE.formula = expression( 1 - exp(u_a)/exp(u_0) )
#### First Derivatives ####
PE.s_u_a = D(PE.formula,'u_a')
PE.s_u_0 = D(PE.formula,'u_0')

# U.by.obsV = expression(
#   muU.est + rho.est * sqrt( sigmaU2.est * sigmaV2.est ) / (sigmaV2.est + 1/(trial.size * exp(obs.V.est - log(1 + exp(obs.V.est))) * (1-exp(obs.V.est - log(1 + exp(obs.V.est)))))) * (obs.V.est - muV.est)
# )
U.by.obsV = expression(
  muU.est + rho.est * sqrt( sigmaU2.est * sigmaV2.est ) / (sigmaV2.est + ( 1 - exp(obs.V.est) )/(trial.size * exp(obs.V.est)) ) * (obs.V.est - muV.est)
)
#### First Derivatives ####
U.by.obsV.s_muU = D(U.by.obsV,'muU.est')
U.by.obsV.s_muV = D(U.by.obsV,'muV.est')
U.by.obsV.s_sigmaU2 = D(U.by.obsV,'sigmaU2.est')
U.by.obsV.s_sigmaV2 = D(U.by.obsV,'sigmaV2.est')
U.by.obsV.s_rho = D(U.by.obsV,'rho.est')
U.by.obsV.s_obsV = D(U.by.obsV,'obs.V.est')


fit.counterfactual.inc.EM<-function(ext.cohort, n.active, case.HIV.active, case.RGC.active, alpha=0.05, n.active.RGC = NULL)
{
  meta.sample <- ext.cohort
  cohort.size <- nrow(ext.cohort)
  # trial.size <- n.active
  
  if(is.null(n.active.RGC))
  {
    n.active.HIV <- n.active.RGC <- n.active
    trial.size <- n.active
  }else
  {
    n.active.HIV <- n.active  
    trial.size <- n.active.RGC
  }
  #################################
  #EM estimation
  #################################
  # data.meta.temp<-data.frame(hUm=logit(meta.sample$est_inc.HIV),hVm=logit(meta.sample$est_inc.RGC))
  # se2_u.temp<- 1 / (meta.sample$est_inc.HIV * (1-meta.sample$est_inc.HIV) * meta.sample$PY.HIV)
  # se2_v.temp<- 1 / (meta.sample$est_inc.RGC * (1-meta.sample$est_inc.RGC) * meta.sample$PY.RGC)
  data.meta.temp<-data.frame(hUm=log(meta.sample$est_inc.HIV),hVm=log(meta.sample$est_inc.RGC))
  se2_u.temp<- (1-meta.sample$est_inc.HIV) / (meta.sample$est_inc.HIV * meta.sample$PY.HIV)
  se2_v.temp<- (1-meta.sample$est_inc.RGC) / (meta.sample$est_inc.RGC * meta.sample$PY.RGC)
  
  EM.fit.temp<-solve_para_EM(data.meta.temp,se2_u.temp,se2_v.temp)
  EM_par<-EM.fit.temp$est
  
  muU.est = EM_par[1]; muV.est = EM_par[2]; sigmaU2.est = EM_par[3]; sigmaV2.est = EM_par[4]; rho.est=EM_par[5];
  
  EM.alpha<- EM_par[1] - EM_par[5] * sqrt(EM_par[3]/EM_par[4]) * EM_par[2]
  EM.beta<- EM_par[5] * sqrt(EM_par[3]/EM_par[4])
  
  coef.est_inc_EM <- c(EM.alpha,EM.beta)
  # pred_HIV_inc_EM <- expit(EM.alpha + EM.beta * logit(case.RGC.active/n.active.RGC) )
  pred_HIV_inc_EM <- exp(EM.alpha + EM.beta * log(case.RGC.active/n.active.RGC) )
  
  gradient.a.coef <- c(eval(a.coef.s_muU),eval(a.coef.s_muV),eval(a.coef.s_sigmaU2),eval(a.coef.s_sigmaV2),eval(a.coef.s_rho))
  gradient.b.coef <- c(eval(b.coef.s_muU),eval(b.coef.s_muV),eval(b.coef.s_sigmaU2),eval(b.coef.s_sigmaV2),eval(b.coef.s_rho))
  
  var.a.coef <- t(gradient.a.coef) %*% EM.fit.temp$cov %*% gradient.a.coef
  var.b.coef <- t(gradient.b.coef) %*% EM.fit.temp$cov %*% gradient.b.coef
  cov.a.b.coef <- t(gradient.a.coef) %*% EM.fit.temp$cov %*% gradient.b.coef
  
  # temp.gradident <- c(1, logit(case.RGC.active/n.active.RGC))
  # temp.gradident <- c(1, logit(case.RGC.active/n.active.RGC))
  temp.gradident <- c(1, log(case.RGC.active/n.active.RGC))
  temp.gradident <- c(1, log(case.RGC.active/n.active.RGC))
  temp.se <- sqrt(t(temp.gradident) %*% matrix(c(var.a.coef,cov.a.b.coef,cov.a.b.coef,var.b.coef),2,2) %*% temp.gradident)
  
  # pred_HIV_inc_EM.lwr <- expit(EM.alpha + EM.beta * logit(case.RGC.active/n.active.RGC) - qnorm(1-alpha/2) * temp.se)
  # pred_HIV_inc_EM.upr <- expit(EM.alpha + EM.beta * logit(case.RGC.active/n.active.RGC) + qnorm(1-alpha/2) * temp.se)
  pred_HIV_inc_EM.lwr <- exp(EM.alpha + EM.beta * log(case.RGC.active/n.active.RGC) - qnorm(1-alpha/2) * temp.se)
  pred_HIV_inc_EM.upr <- exp(EM.alpha + EM.beta * log(case.RGC.active/n.active.RGC) + qnorm(1-alpha/2) * temp.se)
  
  
  ########################################################################
  #   prediction of HIV incidence given observed RGC incidence in trt arm
  ########################################################################
  
  # obs.v<-logit(case.RGC.active/n.active.RGC) 
  # se.v<-sqrt(1/(n.active.RGC * expit(obs.v) * (1-expit(obs.v))))
  obs.v<-log(case.RGC.active/n.active.RGC) 
  se.v<-sqrt( (1 - exp(obs.v))/(n.active.RGC * exp(obs.v) ) )
  
  sV2.est = se.v^2
  
  obs.V.est <- obs.v
  pred.u <- eval(U.by.obsV)
  
  
  obs.V.est <- obs.v
  gradient.U.by.obsV <- c(eval(U.by.obsV.s_muU),eval(U.by.obsV.s_muV),eval(U.by.obsV.s_sigmaU2),eval(U.by.obsV.s_sigmaV2),eval(U.by.obsV.s_rho),eval(U.by.obsV.s_obsV)) 
  
  cov.0<-matrix(0,6,6);cov.0[1:5,1:5]<-EM.fit.temp$cov
  cov.0[6,6]<-sV2.est
  
  var.U.by.obsV <- t(gradient.U.by.obsV) %*% cov.0 %*% gradient.U.by.obsV
  
  obs.V.est <- obs.v
  se.pred.u <- sqrt(var.U.by.obsV)
  
  z.score<-qt(1-alpha/2,df = cohort.size - 2)
  pred.u.upr <- pred.u + z.score * se.pred.u
  pred.u.lwr <- pred.u - z.score * se.pred.u
  
  # pred_HIV_inc_trial <- expit(pred.u)
  # pred_HIV_inc_trial.upr <- expit(pred.u.upr)
  # pred_HIV_inc_trial.lwr <- expit(pred.u.lwr)
  pred_HIV_inc_trial <- exp(pred.u)
  pred_HIV_inc_trial.upr <- exp(pred.u.upr)
  pred_HIV_inc_trial.lwr <- exp(pred.u.lwr)
  
  ##########################################
  #  active trial arm HIV inc and CI
  ##########################################
  trial.inc <- case.HIV.active/n.active.HIV
  
  PE.est.EM<- 1 - trial.inc/pred_HIV_inc_trial
  
  # u_a <- logit(trial.inc)
  # u_0 <- logit(pred_HIV_inc_trial)
  u_a <- log(trial.inc)
  u_0 <- log(pred_HIV_inc_trial)
  
  eval(PE.formula)
  
  gradient.VE <- matrix(c(eval(PE.s_u_a),eval(PE.s_u_0)),2,1)
  
  # se.PE.est.EM <- sqrt(t(gradient.VE) %*% 
  #                        diag(c(
  #                          1/ n.active.HIV / trial.inc / (1 - trial.inc),
  #                          se.pred.u^2))  %*% 
  #                        gradient.VE)
  se.PE.est.EM <- sqrt(t(gradient.VE) %*% 
                         diag(c(
                           (1 - trial.inc) / n.active.HIV / trial.inc ,
                           se.pred.u^2))  %*% 
                         gradient.VE)
  
  PE.lwr.EM <- PE.est.EM - qnorm(1-alpha/2) * se.PE.est.EM 
  PE.upr.EM <- PE.est.EM + qnorm(1-alpha/2) * se.PE.est.EM 
  
  return(list(pred.CF_placebo.HIV.inc = pred_HIV_inc_trial,
              CI.CF_placebo.HIV.inc = c(pred_HIV_inc_trial.lwr,pred_HIV_inc_trial.upr),
              est.CF_PE = PE.est.EM,
              CI.CF_PE = c(PE.lwr.EM,PE.upr.EM))
  )
}


##########################################
#simulation wrapper function
##########################################

CF.simulation.EM <- function(model.parameter, inc.HIV.true, cohort.size, trial.size, PE.set, alpha.set = 0.05, n.rep = 5000, PE_0 = 0.3)
{
  # inc.RGC.grid <- expit((logit(inc.HIV.true) - model.parameter$true.alpha)/model.parameter$true.beta )
  inc.RGC.grid <- exp((log(inc.HIV.true) - model.parameter$true.alpha)/model.parameter$true.beta )
  
  est.counterfactual.inc <- rep(NA,n.rep)
  est.ci.counterfactual.inc <- matrix(NA, n.rep,2)
  
  est.PE <- rep(NA,n.rep)
  est.ci.PE <- matrix(NA, n.rep,2)
  
  for(i.rep in 1:n.rep)
  {
    ##simulate external cohorts
    # meta.sample0<-sim.incidence(1000, model.parameter$logit.mu, model.parameter$logit.sigma, c(200,5000))
    meta.sample0<-sim.incidence(cohort.size*2, model.parameter$log.mu, model.parameter$log.sigma, c(200,5000))
    meta.sample1<-subset(meta.sample0,case.HIV > 0 & case.RGC > 0)
    meta.sample<-meta.sample1[sample(1:nrow(meta.sample1),cohort.size),]
    
    ##observed active-controlled trial HIV/RGC events
    case.RGC.obs<-rbinom(1,trial.size, inc.RGC.grid)
    case.HIV.active.obs<-rbinom(1,trial.size, inc.HIV.true*(1-PE.set))
    
    ##estimating counterfactual placebo HIV inc, and PE
    fit.temp<-fit.counterfactual.inc.EM(meta.sample, trial.size, case.HIV.active.obs, case.RGC.obs, alpha.set)
    
    est.counterfactual.inc[i.rep]<- fit.temp$pred.CF_placebo.HIV.inc
    est.ci.counterfactual.inc[i.rep,] <- fit.temp$CI.CF_placebo.HIV.inc
    
    est.PE[i.rep] <- fit.temp$est.CF_PE
    est.ci.PE[i.rep,] <- fit.temp$CI.CF_PE
  }
  
  return(list(CF.HIV.inc.bias.percent = (mean(est.counterfactual.inc, na.rm=T) - inc.HIV.true)/inc.HIV.true*100,
              CF.HIV.inc.cov.percent = mean(est.ci.counterfactual.inc[,1]<inc.HIV.true & est.ci.counterfactual.inc[,2]>inc.HIV.true, na.rm=T)*100,
              CF.HIV.inc.avgCI.percent = apply(est.ci.counterfactual.inc,2,mean,na.rm=T)*100,
              CF.PE.bias.percent = (mean(est.PE, na.rm=T) - PE.set)/PE.set*100,
              CF.PE.cov.percent = mean(est.ci.PE[,1]<PE.set & est.ci.PE[,2]>PE.set, na.rm=T)*100,
              CF.PE.avgCI.percent = apply(est.ci.PE,2,mean, na.rm=T)*100,
              CF.PE.power.PE0.percent = mean(est.ci.PE[,1]>PE_0, na.rm=T)*100))
}

CF.simulation.EM.rho_m <- function(model.parameter, inc.HIV.true, cohort.size, trial.size, PE.set, alpha.set = 0.05, n.rep = 5000, PE_0 = 0.3, rho_m.range = c(0.4,0.5))
{
  # inc.RGC.grid <- expit((logit(inc.HIV.true) - model.parameter$true.alpha)/model.parameter$true.beta )
  inc.RGC.grid <- exp((log(inc.HIV.true) - model.parameter$true.alpha)/model.parameter$true.beta )
  
  est.counterfactual.inc <- rep(NA,n.rep)
  est.ci.counterfactual.inc <- matrix(NA, n.rep,2)
  
  est.PE <- rep(NA,n.rep)
  est.ci.PE <- matrix(NA, n.rep,2)
  
  for(i.rep in 1:n.rep)
  {
    ##simulate external cohorts
    # meta.sample0<-sim.incidence.rho_m(1000,model.parameter$logit.mu, model.parameter$logit.sigma, cohort.sizes = c(200,5000), rho_m.range)
    meta.sample0<-sim.incidence.rho_m(cohort.size*2,model.parameter$log.mu, model.parameter$log.sigma, cohort.sizes = c(200,5000), rho_m.range)
    meta.sample<-meta.sample0[sample(1:nrow(meta.sample0),cohort.size),]
    
    ##observed active-controlled trial HIV/RGC events
    case.RGC.obs<-rbinom(1,trial.size, inc.RGC.grid)
    case.HIV.active.obs<-rbinom(1,trial.size, inc.HIV.true*(1-PE.set))
    
    ##estimating counterfactual placebo HIV inc, and PE
    fit.temp<-fit.counterfactual.inc.EM(meta.sample, trial.size, case.HIV.active.obs, case.RGC.obs, alpha.set)
    
    est.counterfactual.inc[i.rep]<- fit.temp$pred.CF_placebo.HIV.inc
    est.ci.counterfactual.inc[i.rep,] <- fit.temp$CI.CF_placebo.HIV.inc
    
    est.PE[i.rep] <- fit.temp$est.CF_PE
    est.ci.PE[i.rep,] <- fit.temp$CI.CF_PE
  }
  
  return(list(CF.HIV.inc.bias.percent = (mean(est.counterfactual.inc, na.rm=T) - inc.HIV.true)/inc.HIV.true*100,
              CF.HIV.inc.cov.percent = mean(est.ci.counterfactual.inc[,1]<inc.HIV.true & est.ci.counterfactual.inc[,2]>inc.HIV.true, na.rm=T)*100,
              CF.HIV.inc.avgCI.percent = apply(est.ci.counterfactual.inc,2,mean,na.rm=T)*100,
              CF.PE.bias.percent = (mean(est.PE, na.rm=T) - PE.set)/PE.set*100,
              CF.PE.cov.percent = mean(est.ci.PE[,1]<PE.set & est.ci.PE[,2]>PE.set, na.rm=T)*100,
              CF.PE.avgCI.percent = apply(est.ci.PE,2,mean, na.rm=T)*100,
              CF.PE.power.PE0.percent = mean(est.ci.PE[,1]>PE_0, na.rm=T)*100))
}