## required R library and custom functions
source("working_model_functions.R")
library(MASS)

##test with Murray2019 external cohorts and a active-controlled trial

## HIV/RGC incidence summaries from Murray 2019
hiv.inc<-c(2.5,.9,6.6,3.6,3.8,6.4,9,8.3)/100
rgc.inc<-c(3.5,2.3,15.5,10.1,6.2,16.1,33.1,33)/100
hiv.py<-c(943.2,5160,212.1,1000,843.1,50,245,100)
rgc.py<-c(943.2,5160,212.1,1000,726.6,50,596,100)

MM2019.cohort<-data.frame(est_inc.HIV = hiv.inc,
                          est_inc.RGC = rgc.inc,
                          PY.HIV = hiv.py,
                          PY.RGC = rgc.py)

fit.counterfactual.inc.working_model(MM2019.cohort, 2000, 15, 120, .05,1000)
fit.counterfactual.inc.working_model(MM2019.cohort, 2000, 15, 120, .05, 1000, 2000)
fit.counterfactual.inc.working_model(MM2019.cohort, 2000, 15, 120, .05, 1000, 2500)
fit.counterfactual.inc.working_model(MM2019.cohort, 2000, 30, 205, .05)
fit.counterfactual.inc.working_model(MM2019.cohort, 2000, 70, 205, .05)

#Discover trial 
#F/TAF 651/3014 PYRs, inc=21.6/100 PYRs
#F/TDF 662/3229 PYRs, inc=20.5/100 PYRs
#total 1313/6243 PYRs, inc=21.0/100 PYRs

fit.counterfactual.inc.working_model(MM2019.cohort, 6243, 70, 1313, .05)
# $pred.CF_placebo.HIV.inc
# (Intercept) 
# 0.0686784 
# 
# $CI.CF_placebo.HIV.inc
# (Intercept) (Intercept) 
# 0.05078126  0.09227000 
# 
# $est.CF_PE
# [1] 0.8386175
# 
# $CI.CF_PE
# 2.5%     97.5% 
#   0.7904965 0.8868269 