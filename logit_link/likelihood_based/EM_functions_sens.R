fit.counterfactual.inc.EM_sens<-function(ext.cohort, n.active, case.HIV.active, case.RGC.active, alpha=0.05, rho=0.5, n.active.RGC = NULL)
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
  data.meta.temp<-data.frame(hUm=logit(meta.sample$est_inc.HIV),hVm=logit(meta.sample$est_inc.RGC))
  se2_u.temp<- 1 / (meta.sample$est_inc.HIV * (1-meta.sample$est_inc.HIV) * meta.sample$PY.HIV)
  se2_v.temp<- 1 / (meta.sample$est_inc.RGC * (1-meta.sample$est_inc.RGC) * meta.sample$PY.RGC)
  
  EM.fit.temp<-solve_para_EM(data.meta.temp,se2_u.temp,se2_v.temp)
  EM_par<-EM.fit.temp$est
  
  muU.est = EM_par[1]; muV.est = EM_par[2]; sigmaU2.est = EM_par[3]; sigmaV2.est = EM_par[4]; rho.est=rho
  
  EM.alpha<- EM_par[1] - EM_par[5] * sqrt(EM_par[3]/EM_par[4]) * EM_par[2]
  EM.beta<- EM_par[5] * sqrt(EM_par[3]/EM_par[4])
  
  coef.est_inc_EM <- c(EM.alpha,EM.beta)
  pred_HIV_inc_EM <- expit(EM.alpha + EM.beta * logit(case.RGC.active/n.active.RGC) )
  
  gradient.a.coef <- c(eval(a.coef.s_muU),eval(a.coef.s_muV),eval(a.coef.s_sigmaU2),eval(a.coef.s_sigmaV2),eval(a.coef.s_rho))
  gradient.b.coef <- c(eval(b.coef.s_muU),eval(b.coef.s_muV),eval(b.coef.s_sigmaU2),eval(b.coef.s_sigmaV2),eval(b.coef.s_rho))
  
  var.a.coef <- t(gradient.a.coef) %*% EM.fit.temp$cov %*% gradient.a.coef
  var.b.coef <- t(gradient.b.coef) %*% EM.fit.temp$cov %*% gradient.b.coef
  cov.a.b.coef <- t(gradient.a.coef) %*% EM.fit.temp$cov %*% gradient.b.coef
  
  temp.gradident <- c(1, logit(case.RGC.active/n.active.RGC))
  temp.se <- sqrt(t(temp.gradident) %*% matrix(c(var.a.coef,cov.a.b.coef,cov.a.b.coef,var.b.coef),2,2) %*% temp.gradident)
  
  pred_HIV_inc_EM.lwr <- expit(EM.alpha + EM.beta * logit(case.RGC.active/n.active.RGC) - qnorm(1-alpha/2) * temp.se)
  pred_HIV_inc_EM.upr <- expit(EM.alpha + EM.beta * logit(case.RGC.active/n.active.RGC) + qnorm(1-alpha/2) * temp.se)
  
  
  ########################################################################
  #   prediction of HIV incidence given observed RGC incidence in trt arm
  ########################################################################
  
  obs.v<-logit(case.RGC.active/n.active.RGC) 
  se.v<-sqrt(1/(n.active.RGC * expit(obs.v) * (1-expit(obs.v))))
  
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
  
  pred_HIV_inc_trial <- expit(pred.u)
  pred_HIV_inc_trial.upr <- expit(pred.u.upr)
  pred_HIV_inc_trial.lwr <- expit(pred.u.lwr)
  
  ##########################################
  #  active trial arm HIV inc and CI
  ##########################################
  trial.inc <- case.HIV.active/n.active.HIV
  
  PE.est.EM<- 1 - trial.inc/pred_HIV_inc_trial
  
  u_a <- logit(trial.inc)
  u_0 <- logit(pred_HIV_inc_trial)
  
  eval(PE.formula)
  
  gradient.VE <- matrix(c(eval(PE.s_u_a),eval(PE.s_u_0)),2,1)
  
  se.PE.est.EM <- sqrt(t(gradient.VE) %*% 
                         diag(c(
                           1/ n.active.HIV / trial.inc / (1 - trial.inc),
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