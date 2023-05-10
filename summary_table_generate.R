#result summary, log-link/log_link_RR/logit-link, likelihood/working model based estimation
##Author: Yifan Zhu
##Version 1.0 - 10-31-2022

link.option <- "log"
# link.option <- "logit"

approach.option <- "likelihood"
# approach.option <- "working"

PE.CI.RR.option <- TRUE
windows.directory.set <- TRUE
####################################################################################

result.table.inc<-result.table.PE<-matrix(NA,15,16)
result.table.power<-matrix(NA,5,8)

parameter.store0<-rbind(
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

for(i in 1:40)
{
  parameter.set <- parameter.store0[i,]

  # active.trial.size<-500
  # active.trial.size<-1000
  # active.trial.size<-2000
  active.trial.size<-4000
  # active.trial.size<-10000
  # active.trial.size<-20000

  # PE.set <- .3
  PE.set <- .6
  # PE.set <- .75

  # moderate rho case
  if(parameter.set[3]==0.5){rho<-0.5}else{rho=0.971}

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


	if(windows.directory.set)
	{
		output.directory <- paste("T:/vaccine/Holly Janes/HIV counterfactual placebo/Output/",ifelse(link.option=="log", "log_link", "logit_link"), "/", ifelse(approach.option=="likelihood", ifelse(PE.CI.RR.option,"likelihood_based_logRR/","likelihood_based/"), "working_model/"),sep="")
	}else{
		output.directory <- paste("/trials/vaccine/Holly Janes/HIV counterfactual placebo/Output/",ifelse(link.option=="log", "log_link", "logit_link"), "/", ifelse(approach.option=="likelihood", ifelse(PE.CI.RR.option,"likelihood_based_logRR/","likelihood_based/"), "working_model/"),sep="")
	}

	rdata.path <- paste(output.directory,"Rdata_HIV_inc_",parameter.set[1],"_rho_",ifelse(rho==.5,5,97),"_M_",parameter.set[2],"_VE_",PE.write,"_",active.trial.size.write,".Rdata",sep="")

	print(parameter.set)
	print(rdata.path)

  test.sim<-NULL
  try(load(rdata.path), silent=T)
  if(!is.null(test.sim))
  {
    col.id<-floor((i-1)/5) + 1
    row.id<-3*(i - 5*floor((i-1)/5)-1) + 1

    result.table.inc[row.id, 2*(col.id-1) + 1] <- round(test.sim$CF.HIV.inc.bias.percent,digits=2)
    result.table.inc[row.id+1, 2*(col.id-1) + 1] <- round(test.sim$CF.HIV.inc.cov.percent,digits=1)
    result.table.inc[row.id+2, 2*(col.id-1) + c(1,2)] <- round(test.sim$CF.HIV.inc.avgCI.percent,digits=2)

    result.table.PE[row.id, 2*(col.id-1) + 1] <- round(test.sim$CF.PE.bias.percent,digits=2)
    result.table.PE[row.id+1, 2*(col.id-1) + 1] <- round(test.sim$CF.PE.cov.percent,digits=1)
    result.table.PE[row.id+2, 2*(col.id-1) + c(1,2)] <- round(test.sim$CF.PE.avgCI.percent,digits=2)

    row.id2<-i - 5*floor((i-1)/5)
    col.id2<-floor((i-1)/5) + 1
    # print(c(row.id2,col.id2))
    result.table.power[row.id2, col.id2] <- round(test.sim$CF.PE.power.PE0.percent,digits=2)
  }
}

View(result.table.inc)

View(result.table.PE)

View(result.table.power)



