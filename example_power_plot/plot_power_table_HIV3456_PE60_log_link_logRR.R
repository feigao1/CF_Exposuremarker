
#################################################
#power table plot
#################################################
#load data from simulation calculated power
load("T:/vaccine/Holly Janes/HIV counterfactual placebo/R code/example_power_plot/power_table_HIV3456_PE60_log_link_logRR.Rdata")

####################################################################
#PE=60%, rho=0.980, n.trial.v=1k & 2k
####################################################################
#power table plot
x11(16,9)
par(mar=c(4,7,4,0))

par(fig=c(0.5,.95,0.2,1))
plot(0,0,cex=0,main=bquote(rho == 0.5~","~n[x] == 2000),ylab="power %",cex.lab=1.8,cex.main=1.7,xaxt="n",cex.axis=1.8,axes=F,xlim=c(1,5.5),ylim=c(40,100),xlab="Number of external cohorts")
axis(1,at=1:5,labels=c("10","20","50","100","1000"),cex.axis=1.8)
axis(2,cex.axis=1.8,at=seq(40,100,by=10),labels=seq(40,100,by=10))

#simulation based, from power_Trial_with_placebo.R
n.power.3<-c(804,1340,2528,7467,"++")
n.power.45<-c(395,612,881,1494,3311)
n.power.6<-c(267,376,535,806,1433)

axis(4,cex.axis=1.4,at=seq(40,100,by=10),labels=NA)
text(5+0.5,y=c(50,60,70,80,90)+2,labels=n.power.3,col=1,cex=1.2)
text(5+0.5,y=c(50,60,70,80,90),labels=n.power.45,col=2,cex=1.2)
text(5+0.5,y=c(50,60,70,80,90)-2,labels=n.power.6,col=3,cex=1.2)
mtext("Required PYRs for placebo-controlled RCT",side=4,line=1,cex=1.4)

polygon(c(0,10,10,0), c(90,90,80,80),col=yarrr::transparent("green", trans.val = .8),border=NA)  
polygon(c(0,10,10,0), c(100,100,90,90),col=yarrr::transparent("green", trans.val = .6),border=NA) 
abline(h=c(50,60,70,80,90),lty=3)

for(j in 1:3+3)
{
  # points(1:5, pow_tab_CF_2k[1:5,j],pch=16,cex=2,col=j-3)
  # lines(1:5, pow_tab_CF_2k[1:5,j],lwd=2,col=j-3,lty=2)
  points(2:5, pow_tab_CF_2k[2:5,j],pch=16,cex=2,col=j-3)
  lines(2:5, pow_tab_CF_2k[2:5,j],lwd=2,col=j-3,lty=2)
  
  # points(1:5, pow_tab_CF_2k[1:5+5,j],pch=16,cex=2,col=j+3)
  # lines(1:5, pow_tab_CF_2k[1:5+5,j],lwd=2,col=j+3,lty=2)
  
  # lines(2:3, c(pow_tab_working_2k[2,j],pow_tab_CF_2k[3,j]),lwd=2,col=j-3,lty=2)
  
  # points(1:5, pow_tab_working_2k[1:5,j],pch=4,cex=2,col=j-3)
  points(1:2, pow_tab_working_2k[1:2,j],pch=4,cex=2,col=j-3)
  # lines(1:2, pow_tab_working_2k[1:2,j],lwd=2,col=j-3,lty=2)
}
legend("bottomright",legend=c("Working reg.","MLE"),pch=c(4,16), bty="n",cex=1.5)

##################################################################

par(fig=c(0,0.45,0.2,1),new=T)
plot(0,0,cex=0,main=bquote(rho == 0.980~","~n[x] == 2000),ylab="power %",cex.lab=1.8,cex.main=1.7,xaxt="n",cex.axis=1.8,axes=F,xlim=c(1,5.5),ylim=c(40,100),xlab="Number of external cohorts")
axis(1,at=1:5,labels=c("10","20","50","100","1000"),cex.axis=1.8)
axis(2,cex.axis=1.8,at=seq(40,100,by=10),labels=seq(40,100,by=10))

#simulation based, from power_Trial_with_placebo.R
n.power.3<-c(804,1340,2528,7467,"++")
n.power.45<-c(395,612,881,1494,3311)
n.power.6<-c(267,376,535,806,1433)

axis(4,cex.axis=1.4,at=seq(40,100,by=10),labels=NA)
text(5+0.5,y=c(50,60,70,80,90)+2,labels=n.power.3,col=1,cex=1.2)
text(5+0.5,y=c(50,60,70,80,90),labels=n.power.45,col=2,cex=1.2)
text(5+0.5,y=c(50,60,70,80,90)-2,labels=n.power.6,col=3,cex=1.2)
mtext("Required PYRs for placebo-controlled RCT",side=4,line=1,cex=1.4)

polygon(c(0,10,10,0), c(90,90,80,80),col=yarrr::transparent("green", trans.val = .8),border=NA)  
polygon(c(0,10,10,0), c(100,100,90,90),col=yarrr::transparent("green", trans.val = .6),border=NA) 
abline(h=c(50,60,70,80,90),lty=3)

for(j in 1:3)
{
  # points(1:5, pow_tab_CF_2k[1:5,j],pch=16,cex=2,col=j)
  # lines(1:5, pow_tab_CF_2k[1:5,j],lwd=2,col=j,lty=2)
  points(2:5, pow_tab_CF_2k[2:5,j],pch=16,cex=2,col=j)
  lines(2:5, pow_tab_CF_2k[2:5,j],lwd=2,col=j,lty=2)
  
  
  # points(1:5, pow_tab_CF_2k[1:5+5,j],pch=16,cex=2,col=j+3)
  # lines(1:5, pow_tab_CF_2k[1:5+5,j],lwd=2,col=j+3,lty=2)
  
  # lines(2:3, c(pow_tab_working_2k[2,j],pow_tab_CF_2k[3,j]),lwd=2,col=j,lty=2)
  
  # points(1:5, pow_tab_working_2k[1:5,j],pch=4,cex=2,col=j)
  points(1:2, pow_tab_working_2k[1:2,j],pch=4,cex=2,col=j)
  # lines(1:2, pow_tab_working_2k[1:2,j],lwd=2,col=j,lty=2)
}
legend("bottomright",legend=c("Working reg.","MLE"),pch=c(4,16), bty="n",cex=1.5)

par(fig=c(0.22,0.42,0,0.3),new=T)
plot(0,0,cex=0,main="",ylab="",xaxt="n",axes=F,xlim=c(1,5),ylim=c(0,100),xlab="")
legend("bottom",legend=c(expression(paste(lambda[HIV]==3,"%",sep=""))),pch=18,lty=2,col=1,bty="n",cex=1.8)

par(fig=c(.38,0.58,0,0.3),new=T)
plot(0,0,cex=0,main="",ylab="",xaxt="n",axes=F,xlim=c(1,5),ylim=c(0,100),xlab="")
legend("bottom",legend=c(expression(paste(lambda[HIV]==4.5,"%",sep=""))),pch=18,lty=2,col=2,bty="n",cex=1.8)

par(fig=c(0.54,0.74,0,0.3),new=T)
plot(0,0,cex=0,main="",ylab="",xaxt="n",axes=F,xlim=c(1,5),ylim=c(0,100),xlab="")
legend("bottom",legend=c(expression(paste(lambda[HIV]==6,"%",sep=""))),pch=18,lty=2,col=3,bty="n",cex=1.8)
##################################################################
