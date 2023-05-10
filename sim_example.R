##simulation example code
##Author: Yifan Zhu
##Version 1.0 - 10-13-2022
##Version 1.1 - 10-14-2022, added control sequence for link function and approach option, and default options for alpha, PE_0, number of simulation/bootstrap repeats
##Version 1.2 - 10-31-2022, added "windows.directory" option to allow run under windows PC or SCHARP environment

## options: link functions - logit/log
## options: likelihood based/working model approach
## options: CI construction for PE - delta method applied estimate SE for log(RR) or PE - log-link or for PE - logit-link

link.option <- "log"
# link.option <- "logit"

approach.option <- "likelihood"
# approach.option <- "working"

PE.CI.RR.option <- TRUE

alpha.default <- 0.05
PE_0.default <- 0.3
N.rep.default <- 5000
N.boot.default<-10000
windows.directory <- FALSE
####################################################################################
PE.CI.RR.option_flag <- PE.CI.RR.option & (link.option == "log") & (approach.option == "likelihood")

if(windows.directory)
{
	function_source <- paste("T:/vaccine/Holly Janes/HIV counterfactual placebo/R code/", ifelse(link.option=="log", "log_link", "logit_link"), "/", ifelse(approach.option=="likelihood", ifelse(PE.CI.RR.option_flag,"likelihood_based/EM_functions_logRR.R","likelihood_based/EM_functions.R"), "working_model/working_model_functions.R"),sep="")
}else{
	function_source <- paste("/trials/vaccine/Holly Janes/HIV counterfactual placebo/R code/", ifelse(link.option=="log", "log_link", "logit_link"), "/", ifelse(approach.option=="likelihood", ifelse(PE.CI.RR.option_flag,"likelihood_based/EM_functions_logRR.R","likelihood_based/EM_functions.R"), "working_model/working_model_functions.R"),sep="")
}
source(function_source)

## required R packages
require(mvtnorm)
library(MASS)

## Calculating true model parameters from external cohorts, for intial parameter setting
## HIV/RGC incidence summaries from Murray 2019
hiv.inc<-c(2.5,.9,6.6,3.6,3.8,6.4,9,8.3)/100
rgc.inc<-c(3.5,2.3,15.5,10.1,6.2,16.1,33.1,33)/100

hiv.py<-c(943.2,5160,212.1,1000,843.1,50,245,100)
rgc.py<-c(943.2,5160,212.1,1000,726.6,50,596,100)

ini.mu_u <- -3.117
ini.mu_v <- -2.0909
ini.sd2_u <- .7941^2
ini.sd2_v <- 1.1237^2
ini.rho <- .938

ini.par<-c(ini.mu_u,ini.mu_v,ini.sd2_u,ini.sd2_v,ini.rho)

fit.inc_par.MLE<-optim(ini.par,log.lik_incidences,est_inc_hiv = hiv.inc, est_inc_rgc = rgc.inc, py_hiv = hiv.py, py_rgc = rgc.py, method = "L-BFGS-B", lower = c(-Inf,-Inf,0,0,-1), upper = c(Inf,Inf,Inf,Inf,1),control = list(fnscale = -1))

print(fit.inc_par.MLE$par)
# logit-link: [1] -3.1383824 -2.0811034  0.5813719  1.0649200  0.9705414
# log-link: [1] -3.1891158 -2.2454655  0.5365001  0.8143491  0.9800112

# data.ini<-data.frame(hUm=logit(hiv.inc),hVm=logit(rgc.inc))
# se2_u<- 1 / (hiv.inc * (1-hiv.inc) * hiv.py)
# se2_v<- 1 / (rgc.inc * (1-rgc.inc) * rgc.py)
# 
# solve_para_EM(data.ini,se2_u,se2_v)$est
# [1] -3.1382844 -2.0816424  0.5803660  1.0663475  0.9687279

# data.ini<-data.frame(hUm=log(hiv.inc),hVm=log(rgc.inc))
# se2_u<- (1-hiv.inc) / (hiv.inc * hiv.py)
# se2_v<- (1-rgc.inc) / (rgc.inc * rgc.py)
# 
# solve_para_EM(data.ini,se2_u,se2_v)$est
# [1] -3.1886563 -2.2466765  0.5352886  0.8161646  0.9769509


## set true model parameters as the MLE from MM2019
mu_u<-fit.inc_par.MLE$par[1]
mu_v<-fit.inc_par.MLE$par[2]
sd2_u<-fit.inc_par.MLE$par[3]
sd2_v<-fit.inc_par.MLE$par[4]
rho<-fit.inc_par.MLE$par[5]


## scenario selection for simulation parameters
args<-commandArgs(trailingOnly=T)

## (placebo.HIV.inc, N.cohort, rho)
parameter.store<-rbind(
  c(1/100,10, .971),
  c(1/100,20, .971),
  c(1/100,50, .971),
  c(1/100,100, .971),
  c(1/100,1000, .971),
  c(3/100,10, .971),
  c(3/100,20, .971),
  c(3/100,50, .971),
  c(3/100,100, .971),
  c(3/100,1000, .971),
  c(4.5/100,10, .971),
  c(4.5/100,20, .971),
  c(4.5/100,50, .971),
  c(4.5/100,100, .971),
  c(4.5/100,1000, .971),
  c(6/100,10, .971),
  c(6/100,20, .971),
  c(6/100,50, .971),
  c(6/100,100, .971),
  c(6/100,1000, .971),
  c(1/100,10, .5),
  c(1/100,20, .5),
  c(1/100,50, .5),
  c(1/100,100, .5),
  c(1/100,1000, .5),
  c(3/100,10, .5),
  c(3/100,20, .5),
  c(3/100,50, .5),
  c(3/100,100, .5),
  c(3/100,1000, .5),
  c(4.5/100,10, .5),
  c(4.5/100,20, .5),
  c(4.5/100,50, .5),
  c(4.5/100,100, .5),
  c(4.5/100,1000, .5),
  c(6/100,10, .5),
  c(6/100,20, .5),
  c(6/100,50, .5),
  c(6/100,100, .5),
  c(6/100,1000, .5)
)

parameter.set <- parameter.store[as.numeric(args[1]),]

# active.trial.size<-1000
# active.trial.size<-2000
# active.trial.size<-4000

active.trial.size<-as.numeric(args[2])

# PE.set <- .3
# PE.set <- .6
# PE.set <- .75

PE.set <- as.numeric(args[3])

# moderate rho case
if(parameter.set[3]==0.5){rho<-0.5}

######################################################################
# run simulation
######################################################################
ifelse(link.option == "log",
		par.set<-list(log.mu = c(mu_u,mu_v),
					  log.sigma = matrix(c(sd2_u,rho*sqrt(sd2_u * sd2_v),rho*sqrt(sd2_u * sd2_v),sd2_v),2,2),
					  true.alpha = mu_u - rho * sqrt(sd2_u/sd2_v) * mu_v,
					  true.beta = rho * sqrt(sd2_u/sd2_v),
					  true.sigma = sqrt(sd2_u * (1 - rho^2))),
		par.set<-list(logit.mu = c(mu_u,mu_v),
					  logit.sigma = matrix(c(sd2_u,rho*sqrt(sd2_u * sd2_v),rho*sqrt(sd2_u * sd2_v),sd2_v),2,2),
					  true.alpha = mu_u - rho * sqrt(sd2_u/sd2_v) * mu_v,
					  true.beta = rho * sqrt(sd2_u/sd2_v),
					  true.sigma = sqrt(sd2_u * (1 - rho^2)))
)

##run simulations
set.seed(1)

ifelse(approach.option == "likelihood",
test.sim<- CF.simulation.EM(model.parameter = par.set, inc.HIV.true = parameter.set[1], cohort.size = parameter.set[2], trial.size = active.trial.size, PE.set = PE.set, alpha.set = alpha.default, n.rep = N.rep.default, PE_0 = PE_0.default),
test.sim<- CF.simulation.working_model(model.parameter = par.set, inc.HIV.true = parameter.set[1], cohort.size = parameter.set[2], trial.size = active.trial.size, PE.set = PE.set, alpha.set = alpha.default, n.rep = N.rep.default, n.boot = N.boot.default, PE_0 = PE_0.default)
)

# print results
test.sim

if(PE.set==0){PE.write<-0}
if(PE.set==0.3){PE.write<-30}
if(PE.set==0.6){PE.write<-60}
if(PE.set==0.75){PE.write<-75}
if(PE.set==0.9){PE.write<-90}

if(active.trial.size==500){active.trial.size.write<-"500"}
if(active.trial.size==1000){active.trial.size.write<-"1k"}
if(active.trial.size==2000){active.trial.size.write<-"2k"}
if(active.trial.size==4000){active.trial.size.write<-"4k"}
if(active.trial.size==6000){active.trial.size.write<-"6k"}
if(active.trial.size==10000){active.trial.size.write<-"10k"}
if(active.trial.size==20000){active.trial.size.write<-"20k"}

# sample output directory
if(windows.directory)
{
	output.directory <- paste("T:/vaccine/Holly Janes/HIV counterfactual placebo/Output/",ifelse(link.option=="log", "log_link", "logit_link"), "/", ifelse(approach.option=="likelihood", ifelse(PE.CI.RR.option,"likelihood_based_logRR/","likelihood_based/"), "working_model/"),sep="")
}else{
	output.directory <- paste("/trials/vaccine/Holly Janes/HIV counterfactual placebo/Output/",ifelse(link.option=="log", "log_link", "logit_link"), "/", ifelse(approach.option=="likelihood", ifelse(PE.CI.RR.option,"likelihood_based_logRR/","likelihood_based/"), "working_model/"),sep="")
}

rdata.path <- paste(output.directory,"Rdata_HIV_inc_",parameter.set[1],"_rho_",ifelse(rho==.5,5,97),"_M_",parameter.set[2],"_VE_",PE.write,"_",active.trial.size.write,".Rdata",sep="")

save.image(rdata.path)
