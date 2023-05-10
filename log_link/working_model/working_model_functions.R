# logit<-function(x){log(x/(1-x))}
# expit<-function(x){ifelse(x>20,1,exp(x - log(1 + exp(x))))}

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


fit.counterfactual.inc.working_model<-function(ext.cohort, n.active, case.HIV.active, case.RGC.active, alpha=0.05, n.boot=10000, n.active.RGC = NULL)
{
  if(is.null(n.active.RGC))
  {
    n.active.RGC<-n.active.HIV<-n.active
  }else
  {
    n.active.HIV<-n.active;    
  }
  
  meta.sample <- ext.cohort
  cohort.size <- nrow(ext.cohort)
  
  # working model associating observed incidences
  # y1<-logit(meta.sample$est_inc.HIV)
  # x1<-logit(meta.sample$est_inc.RGC)
  y1<-log(meta.sample$est_inc.HIV)
  x1<-log(meta.sample$est_inc.RGC)
  # #fit.logit_est_inc<-lm(y1~x1, weights = meta.sample$PY.HIV)
  # fit.logit_est_inc<-lm(y1~x1, weights = NULL)
  #fit.log_est_inc<-lm(y1~x1, weights = meta.sample$PY.HIV)
  fit.log_est_inc<-lm(y1~x1, weights = NULL)
  
  # coef.est_inc <- coef(fit.logit_est_inc)
  # coef.est_inc.ci <- confint(fit.logit_est_inc)
  # pred_HIV_inc_2 <- expit(coef(fit.logit_est_inc)[1] + coef(fit.logit_est_inc)[2] * logit(case.RGC.active/n.active.RGC)  )
  # pred_HIV_inc_2.upr<-expit(predict(fit.logit_est_inc, newdata = data.frame(x1 = logit(case.RGC.active/n.active.RGC)  ), interval="prediction",level = .95)[,3])
  # pred_HIV_inc_2.lwr<-expit(predict(fit.logit_est_inc, newdata = data.frame(x1 = logit(case.RGC.active/n.active.RGC)  ), interval="prediction",level = .95)[,2])
  # conf_HIV_inc_2.upr<-expit(predict(fit.logit_est_inc, newdata = data.frame(x1 = logit(case.RGC.active/n.active.RGC)  ), interval="confidence",level = .95)[,3])
  # conf_HIV_inc_2.lwr<-expit(predict(fit.logit_est_inc, newdata = data.frame(x1 = logit(case.RGC.active/n.active.RGC)  ), interval="confidence",level = .95)[,2])
  coef.est_inc <- coef(fit.log_est_inc)
  coef.est_inc.ci <- confint(fit.log_est_inc)
  pred_HIV_inc_2 <- exp(coef(fit.log_est_inc)[1] + coef(fit.log_est_inc)[2] * log(case.RGC.active/n.active.RGC)  )
  pred_HIV_inc_2.upr<-exp(predict(fit.log_est_inc, newdata = data.frame(x1 = log(case.RGC.active/n.active.RGC)  ), interval="prediction",level = .95)[,3])
  pred_HIV_inc_2.lwr<-exp(predict(fit.log_est_inc, newdata = data.frame(x1 = log(case.RGC.active/n.active.RGC)  ), interval="prediction",level = .95)[,2])
  conf_HIV_inc_2.upr<-exp(predict(fit.log_est_inc, newdata = data.frame(x1 = log(case.RGC.active/n.active.RGC)  ), interval="confidence",level = .95)[,3])
  conf_HIV_inc_2.lwr<-exp(predict(fit.log_est_inc, newdata = data.frame(x1 = log(case.RGC.active/n.active.RGC)  ), interval="confidence",level = .95)[,2])
  
  #########################################################################
  ##   prediction of HIV incidence given observed RGC incidence in trt arm
  #########################################################################
  
  # obs.v<-logit(case.RGC.active/n.active.RGC) 
  # se.v<-sqrt(1/(n.active.RGC * expit(obs.v) * (1-expit(obs.v))))
  obs.v<-log(case.RGC.active/n.active.RGC) 
  se.v<-sqrt( (1-exp(obs.v))/(n.active.RGC * exp(obs.v) ) )
  
  # pred.u<- coef(fit.logit_est_inc)[1] + coef(fit.logit_est_inc)[2] * obs.v
  pred.u<- coef(fit.log_est_inc)[1] + coef(fit.log_est_inc)[2] * obs.v
  
  v.mat<-obs.v
  
  # se.pred.u <- sqrt(coef(fit.logit_est_inc)[2]^2 * se.v^2 + sigma(fit.logit_est_inc)^2 * (mean(x1^2) - 2 * mean(x1) * v.mat + v.mat^2 + se.v^2)/  (sum(x1^2) - sum(x1)^2/length(x1) ))
  se.pred.u <- sqrt(coef(fit.log_est_inc)[2]^2 * se.v^2 + sigma(fit.log_est_inc)^2 * (mean(x1^2) - 2 * mean(x1) * v.mat + v.mat^2 + se.v^2)/  (sum(x1^2) - sum(x1)^2/length(x1) ))
  
  #z.score<-qnorm(1-alpha/2)
  z.score<-qt(1-alpha/2,df = cohort.size - 2)
  pred.u.upr <- pred.u + z.score * se.pred.u
  pred.u.lwr <- pred.u - z.score * se.pred.u
  
  # pred_HIV_inc_trial <- expit(pred.u)
  # pred_HIV_inc_trial.upr <- expit(pred.u.upr)
  # pred_HIV_inc_trial.lwr <- expit(pred.u.lwr)
  pred_HIV_inc_trial <- exp(pred.u)
  pred_HIV_inc_trial.upr <- exp(pred.u.upr)
  pred_HIV_inc_trial.lwr <- exp(pred.u.lwr)
  
  ###########################################################
  ##  Bootstrap power for PE
  ###########################################################
  
  ###########################################
  ##  lm bootstrap prediction CI
  ###########################################
  # my.reg <- fit.logit_est_inc
  my.reg <- fit.log_est_inc
  
  leverage <- influence(my.reg)$hat
  my.s.resid <- residuals(my.reg)/sqrt(1-leverage)
  my.s.resid <- my.s.resid - mean(my.s.resid)
  
  reg <- my.reg
  s <- my.s.resid
  
  the.replication <- function(reg,s,case.RGC.active,n.active.RGC)
  {
    # Make bootstrap residuals
    ep.star <- sample(s,size=length(reg$residuals),replace=TRUE)
    
    # Make bootstrap Y
    y.star <- fitted(reg)+ep.star
    
    # Do bootstrap regression
    x <- model.frame(reg)[,2]
    bs.reg <- lm(y.star~x)
    
    # Create bootstrapped adjusted residuals
    bs.lev <- influence(bs.reg)$hat
    bs.s   <- residuals(bs.reg)/sqrt(1-bs.lev)
    bs.s   <- bs.s - mean(bs.s)
    
    # Calculate draw on prediction error
    xb.xb <- coef(my.reg)["(Intercept)"] - coef(bs.reg)["(Intercept)"] 
    
    case.RGC.active.boot <- rbinom(1,n.active.RGC,case.RGC.active/n.active.RGC)
    # obs.v.boot<-logit(case.RGC.active.boot/n.active.RGC) 
    obs.v.boot<-log(case.RGC.active.boot/n.active.RGC) 
    
    xb.xb.j <- xb.xb + (coef(my.reg)["x1"] - coef(bs.reg)["x"])*obs.v.boot
    return(unname(xb.xb.j + rep(sample(bs.s,size=1),length=1)))
  }
  
  #############################################################
  obs.RGC.sample <- rbinom(n.active.RGC,1,case.RGC.active/n.active.RGC) 
  
  PE.est.boot.sample<-PE.est.boot.sample2<-PE.est.boot.sample3<-rep(NA,n.boot)
  pred.HIV_placebo_inc.boot.sample<-pred.HIV_placebo_inc.boot.sample2<-pred.HIV_placebo_inc.boot.sample3<-rep(NA,n.boot)
  for(i.boot in 1:n.boot)
  {
    meta.sample.boot <- meta.sample[sample(1:cohort.size,cohort.size,replace = T),]
    rownames(meta.sample.boot) <- NULL
    
    # y1.boot<-logit(meta.sample.boot$est_inc.HIV)
    # x1.boot<-logit(meta.sample.boot$est_inc.RGC)
    # fit.logit_est_inc.boot<-lm(y1.boot~x1.boot, weights = NULL)
    y1.boot<-log(meta.sample.boot$est_inc.HIV)
    x1.boot<-log(meta.sample.boot$est_inc.RGC)
    fit.log_est_inc.boot<-lm(y1.boot~x1.boot, weights = NULL)
    
    # obs.v.boot <- logit(sum(sample(obs.RGC.sample,n.active.RGC,replace=T))/n.active.RGC)
    obs.v.boot <- log(sum(sample(obs.RGC.sample,n.active.RGC,replace=T))/n.active.RGC)
    
    # pred.u.boot<- coef(fit.logit_est_inc.boot)[1] + coef(fit.logit_est_inc.boot)[2] * obs.v.boot
    # pred.HIV_placebo_inc.boot.sample[i.boot] <- expit(pred.u.boot)
    pred.u.boot<- coef(fit.log_est_inc.boot)[1] + coef(fit.log_est_inc.boot)[2] * obs.v.boot
    pred.HIV_placebo_inc.boot.sample[i.boot] <- exp(pred.u.boot)
    
    # pred.u.boot2<- coef(fit.logit_est_inc)[1] + coef(fit.logit_est_inc)[2] * obs.v.boot + the.replication(reg=my.reg,s=my.s.resid,case.RGC.active,n.active.RGC)
    # pred.HIV_placebo_inc.boot.sample2[i.boot] <- expit(pred.u.boot2)
    pred.u.boot2<- coef(fit.log_est_inc.boot)[1] + coef(fit.log_est_inc.boot)[2] * obs.v.boot + the.replication(reg=my.reg,s=my.s.resid,case.RGC.active,n.active.RGC)
    pred.HIV_placebo_inc.boot.sample2[i.boot] <- exp(pred.u.boot2)
    
    # se.pred.u.boot <- sqrt(coef(fit.logit_est_inc.boot)[2]^2 * se.v^2 + sigma(fit.logit_est_inc.boot)^2 * (mean(x1^2) - 2 * mean(x1) * v.mat + v.mat^2 + se.v^2)/  (sum(x1^2) - sum(x1)^2/length(x1) ))
    # pred.u.boot3<- coef(fit.logit_est_inc.boot)[1] + coef(fit.logit_est_inc.boot)[2] * obs.v.boot + rt(1,df=n.active.RGC-2) * se.pred.u.boot
    # pred.HIV_placebo_inc.boot.sample3[i.boot] <- expit(pred.u.boot3)
    se.pred.u.boot <- sqrt(coef(fit.log_est_inc.boot)[2]^2 * se.v^2 + sigma(fit.log_est_inc.boot)^2 * (mean(x1^2) - 2 * mean(x1) * v.mat + v.mat^2 + se.v^2)/  (sum(x1^2) - sum(x1)^2/length(x1) ))
    pred.u.boot3<- coef(fit.log_est_inc.boot)[1] + coef(fit.log_est_inc.boot)[2] * obs.v.boot + rt(1,df=n.active.RGC-2) * se.pred.u.boot
    pred.HIV_placebo_inc.boot.sample3[i.boot] <- exp(pred.u.boot3)
    ###########################################
    ##  active trial arm HIV inc 
    ###########################################
    sim_HIV_inc_active.boot <- sum(sample(c(rep(1,case.HIV.active),rep(0,n.active.HIV-case.HIV.active)),n.active.HIV,replace=T))/n.active.HIV 
    
    # PE.est.boot.sample[i.boot] <- 1 - sim_HIV_inc_active.boot/expit(pred.u.boot)
    # PE.est.boot.sample2[i.boot] <- 1 - sim_HIV_inc_active.boot/expit(pred.u.boot2)
    # PE.est.boot.sample3[i.boot] <- 1 - sim_HIV_inc_active.boot/expit(pred.u.boot3)
    PE.est.boot.sample[i.boot] <- 1 - sim_HIV_inc_active.boot/exp(pred.u.boot)
    PE.est.boot.sample2[i.boot] <- 1 - sim_HIV_inc_active.boot/exp(pred.u.boot2)
    PE.est.boot.sample3[i.boot] <- 1 - sim_HIV_inc_active.boot/exp(pred.u.boot3)
    
  }
  PE.est.boot<-median(PE.est.boot.sample,na.rm=T)
  PE.lwr.boot<-quantile(PE.est.boot.sample,alpha/2,na.rm=T)
  PE.upr.boot<-quantile(PE.est.boot.sample,1-alpha/2,na.rm=T)
  
  PE.est.boot2<-median(PE.est.boot.sample2,na.rm=T)
  PE.lwr.boot2<-quantile(PE.est.boot.sample2,alpha/2,na.rm=T)
  PE.upr.boot2<-quantile(PE.est.boot.sample2,1-alpha/2,na.rm=T)
  
  PE.est.boot3<-median(PE.est.boot.sample3,na.rm=T)
  PE.lwr.boot3<-quantile(PE.est.boot.sample3,alpha/2,na.rm=T)
  PE.upr.boot3<-quantile(PE.est.boot.sample3,1-alpha/2,na.rm=T)
  
  
  pred_HIV_placebo_inc.est.boot<-median(pred.HIV_placebo_inc.boot.sample,na.rm=T)
  pred_HIV_placebo_inc.lwr.boot<-quantile(pred.HIV_placebo_inc.boot.sample,alpha/2,na.rm=T)
  pred_HIV_placebo_inc.upr.boot<-quantile(pred.HIV_placebo_inc.boot.sample,1-alpha/2,na.rm=T)
  
  pred_HIV_placebo_inc.est.boot2<-median(pred.HIV_placebo_inc.boot.sample2,na.rm=T)
  pred_HIV_placebo_inc.lwr.boot2<-quantile(pred.HIV_placebo_inc.boot.sample2,alpha/2,na.rm=T)
  pred_HIV_placebo_inc.upr.boot2<-quantile(pred.HIV_placebo_inc.boot.sample2,1-alpha/2,na.rm=T)
  
  pred_HIV_placebo_inc.est.boot3<-median(pred.HIV_placebo_inc.boot.sample3,na.rm=T)
  pred_HIV_placebo_inc.lwr.boot3<-quantile(pred.HIV_placebo_inc.boot.sample3,alpha/2,na.rm=T)
  pred_HIV_placebo_inc.upr.boot3<-quantile(pred.HIV_placebo_inc.boot.sample3,1-alpha/2,na.rm=T)
  
  
  return(list(pred.CF_placebo.HIV.inc = pred_HIV_inc_trial,
              CI.CF_placebo.HIV.inc = c(pred_HIV_inc_trial.lwr,pred_HIV_inc_trial.upr),
              # pred.CF_placebo.HIV.inc.boot = pred_HIV_placebo_inc.est.boot,
              # CI.CF_placebo.HIV.inc.boot = c(pred_HIV_placebo_inc.lwr.boot,pred_HIV_placebo_inc.upr.boot),
              # est.CF_PE = PE.est.boot2,
              # CI.CF_PE = c(PE.lwr.boot,PE.upr.boot2),
              # est.CF_PE = PE.est.boot3,
              # CI.CF_PE = c(PE.lwr.boot,PE.upr.boot3),
              est.CF_PE = PE.est.boot,
              CI.CF_PE = c(PE.lwr.boot,PE.upr.boot))
  )
}

##########################################
#simulation wrapper function
##########################################

CF.simulation.working_model <- function(model.parameter, inc.HIV.true, cohort.size, trial.size, PE.set, alpha.set = 0.05, n.rep = 5000, n.boot = 10000, PE_0 = 0.3)
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
    fit.temp<-fit.counterfactual.inc.working_model(meta.sample, trial.size, case.HIV.active.obs, case.RGC.obs, alpha.set, n.boot)
    
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

CF.simulation.working_model.rho_m <- function(model.parameter, inc.HIV.true, cohort.size, trial.size, PE.set, alpha.set = 0.05, n.rep = 5000, n.boot = 10000, PE_0 = 0.3, rho_m.range = c(0.4,0.5))
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
    fit.temp<-fit.counterfactual.inc.working_model(meta.sample, trial.size, case.HIV.active.obs, case.RGC.obs, alpha.set, n.boot)
    
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