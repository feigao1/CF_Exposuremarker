log_lik = expression(-log(2*pi) - log((sigmaU2 + sU2)*(sigmaV2 + sV2)-rho^2*sigmaU2*sigmaV2)/2 - ((sigmaV2 + sV2)*(hUm - muU)^2 - 2*rho*sqrt(sigmaU2*sigmaV2)*(hUm-muU)*(hVm-muV) + (sigmaU2 + sU2)*(hVm-muV)^2)/(2*((sigmaU2 + sU2)*(sigmaV2 + sV2)-rho^2*sigmaU2*sigmaV2)))

data = data.meta.temp; sU2= se2_u.temp;sV2=se2_v.temp
solve_para <- function(data,sU2,sV2){
  hUm = data$hUm; hVm = data$hVm
  muU = mean(hUm); muV = mean(hVm); sigmaU2 = var(hUm)-mean(sU2); sigmaV2 = var(hVm)-mean(sV2)
  rho = cov(hUm,hVm)/sqrt(sigmaU2*sigmaV2)
  log_lik <- function(para){
  	muU = para[1]; muV = para[2]; sigmaU2 = para[3]; sigmaV2 = para[4]; rho = para[5]
  	return(sum( -log(2*pi) - log((sigmaU2 + sU2)*(sigmaV2 + sV2)-rho^2*sigmaU2*sigmaV2)/2 - ((sigmaV2 + sV2)*(hUm - muU)^2 - 2*rho*sqrt(sigmaU2*sigmaV2)*(hUm-muU)*(hVm-muV) + (sigmaU2 + sU2)*(hVm-muV)^2)/(2*((sigmaU2 + sU2)*(sigmaV2 + sV2)-rho^2*sigmaU2*sigmaV2))))
  }
  optim(c(0,0,1,1,0),log_lik)
  sd_est = sqrt(diag(-solve(Hessian)))
  para_est = c(muU,muV,sigmaU2,sigmaV2,rho)
  return(list(est = para_est, sd = sd_est,cov = (d<0.001)))
}

fit.counterfactual.inc.directMLE<-function(ext.cohort, n.active, case.HIV.active, case.RGC.active, alpha=0.05, n.active.RGC = NULL)
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
  
  EM.fit.temp<-solve_para(data.meta.temp,se2_u.temp,se2_v.temp)
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