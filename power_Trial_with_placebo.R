##placebo-controlled sample size calculation code
##Author: Yifan Zhu
##Version 1.0 - 10-14-2022

logit<-function(x){log(x/(1-x))}
expit<-function(x){ifelse(x>20,1,exp(x - log(1 + exp(x))))}

#logit-link
get.power<-function(n.p, n.v,PE_a,PE_0, HIV.inc.0,alpha=0.05)
{
  HIV.inc.v <- HIV.inc.0 * (1-PE_a)
  set.seed(1)
  
  n.rep<-20000
  PE.hat.p1.sum <- rep(NA,n.rep)
  CI.PE.p1.sum <- matrix(NA,2,n.rep)
  
  for(i.rep in 1:n.rep)
  {
    p.hat.v <- rbinom(length(HIV.inc.0),n.v,HIV.inc.v)/n.v
    se.logit.p.hat.v <- sqrt(p.hat.v*(1-p.hat.v)/n.v) / (p.hat.v*(1-p.hat.v))
    
    p.hat.p1 <- rbinom(length(HIV.inc.0),n.p,HIV.inc.0)/n.p
    se.logit.p.hat.p1 <- sqrt(p.hat.p1*(1-p.hat.p1)/n.p) / (p.hat.p1*(1-p.hat.p1))
    
    PE.hat.p1 <- 1 - p.hat.v/p.hat.p1
    se.PE.hat.p1 <- sqrt( se.logit.p.hat.v^2 * (1-p.hat.v)^2 * p.hat.v^2 / p.hat.p1^2 
                          + se.logit.p.hat.p1^2 * p.hat.v^2 * (1-p.hat.p1)^2 * p.hat.p1^2 / p.hat.p1^4  )
    
    CI.PE.hat.p1 <- rbind(PE.hat.p1 - qnorm(1-alpha/2) * se.PE.hat.p1, PE.hat.p1 + qnorm(1-alpha/2) * se.PE.hat.p1)
    
    PE.hat.p1.sum[i.rep] <- PE.hat.p1
    CI.PE.p1.sum[,i.rep] <- CI.PE.hat.p1
  }
  return(round(mean(CI.PE.p1.sum[1,] > PE_0,na.rm=T)*100, digits=1))
}

# get.power(1000,1000,.6,.3,.06)

get.samplesize<-function(beta, n.v,PE_a,PE_0, HIV.inc.0,alpha=0.05)
{
  ini.n.p<-10^(1:6)
  ini.power <- rep(NA,length(ini.n.p))
  for(j in 1:length(ini.power))
  {
    ini.power[j]<-get.power(ini.n.p[j], n.v,PE_a,PE_0, HIV.inc.0,alpha)
  }
 
  if(ini.power[length(ini.n.p)]<beta){return(Inf)}else
  {
    n1<-ini.n.p[max(which(ini.power<=beta))]
    n2<-ini.n.p[max(which(ini.power<=beta))+1]
  
    while(n2>n1+1)
    {
      n3<-max(floor((n1+n2)/2),n1)
      power.new<-get.power(n3, n.v,PE_a,PE_0, HIV.inc.0,alpha)
      if(power.new>=beta)
      {
        n2 <- n3
      }else
      {
        n1<-n3
      }
    }
    return(n2)
  }
}

######log-link, PE based CI
get.power.loglink<-function(n.p, n.v,PE_a,PE_0, HIV.inc.0,alpha=0.05)
{
  HIV.inc.v <- HIV.inc.0 * (1-PE_a)
  set.seed(1)
  
  n.rep<-20000
  PE.hat.p1.sum <- rep(NA,n.rep)
  CI.PE.p1.sum <- matrix(NA,2,n.rep)
  
  for(i.rep in 1:n.rep)
  {
    p.hat.v <- rbinom(length(HIV.inc.0),n.v,HIV.inc.v)/n.v
    # se.logit.p.hat.v <- sqrt(p.hat.v*(1-p.hat.v)/n.v) / (p.hat.v*(1-p.hat.v))
    se.log.p.hat.v <- sqrt((1-p.hat.v)/n.v/p.hat.v)
    
    p.hat.p1 <- rbinom(length(HIV.inc.0),n.p,HIV.inc.0)/n.p
    # se.logit.p.hat.p1 <- sqrt(p.hat.p1*(1-p.hat.p1)/n.p) / (p.hat.p1*(1-p.hat.p1))
    se.log.p.hat.p1 <- sqrt((1-p.hat.p1)/n.p/p.hat.p1)
    
    PE.hat.p1 <- 1 - p.hat.v/p.hat.p1
    # se.PE.hat.p1 <- sqrt( se.logit.p.hat.v^2 * (1-p.hat.v)^2 * p.hat.v^2 / p.hat.p1^2 
    #                       + se.logit.p.hat.p1^2 * p.hat.v^2 * (1-p.hat.p1)^2 * p.hat.p1^2 / p.hat.p1^4  )
    se.PE.hat.p1 <- sqrt( se.log.p.hat.v^2 * p.hat.v^2 / p.hat.p1^2 
                          + se.log.p.hat.p1^2 * p.hat.v^2 / p.hat.p1^2 )
    
    CI.PE.hat.p1 <- rbind(PE.hat.p1 - qnorm(1-alpha/2) * se.PE.hat.p1, PE.hat.p1 + qnorm(1-alpha/2) * se.PE.hat.p1)
    
    PE.hat.p1.sum[i.rep] <- PE.hat.p1
    CI.PE.p1.sum[,i.rep] <- CI.PE.hat.p1
  }
  return(round(mean(CI.PE.p1.sum[1,] > PE_0,na.rm=T)*100, digits=1))
}

##log-RR based CI
get.power.loglink2<-function(n.p, n.v,PE_a,PE_0, HIV.inc.0,alpha=0.05)
{
  HIV.inc.v <- HIV.inc.0 * (1-PE_a)
  set.seed(1)
  
  n.rep<-20000
  PE.hat.p1.sum <- rep(NA,n.rep)
  CI.PE.p1.sum <- matrix(NA,2,n.rep)
  CI.log.RR.p1.sum <- matrix(NA,2,n.rep)
  
  for(i.rep in 1:n.rep)
  {
    p.hat.v <- rbinom(length(HIV.inc.0),n.v,HIV.inc.v)/n.v
    # se.logit.p.hat.v <- sqrt(p.hat.v*(1-p.hat.v)/n.v) / (p.hat.v*(1-p.hat.v))
    se.log.p.hat.v <- sqrt((1-p.hat.v)/n.v/p.hat.v)
    
    p.hat.p1 <- rbinom(length(HIV.inc.0),n.p,HIV.inc.0)/n.p
    # se.logit.p.hat.p1 <- sqrt(p.hat.p1*(1-p.hat.p1)/n.p) / (p.hat.p1*(1-p.hat.p1))
    se.log.p.hat.p1 <- sqrt((1-p.hat.p1)/n.p/p.hat.p1)
    
    log.RR.hat <- log(p.hat.v) - log(p.hat.p1)
    se.log.RR.hat <- sqrt(se.log.p.hat.v^2 + se.log.p.hat.p1^2)
    CI.log.RR.p1 <- c(log.RR.hat - qnorm(1-alpha/2) * se.log.RR.hat, log.RR.hat + qnorm(1-alpha/2) * se.log.RR.hat)
    
    
    PE.hat.p1.sum[i.rep] <- 1 - exp(log.RR.hat)
    CI.log.RR.p1.sum[,i.rep] <- CI.log.RR.p1
    CI.PE.p1.sum[,i.rep] <- rev(1 - exp(CI.log.RR.p1))
  }
  return(round(mean(CI.PE.p1.sum[1,] > PE_0,na.rm=T)*100, digits=1))
}

# get.power.loglink(1000,1000,.6,.3,.06)

get.samplesize.loglink<-function(beta, n.v,PE_a,PE_0, HIV.inc.0,alpha=0.05)
{
  ini.n.p<-10^(1:6)
  ini.power <- rep(NA,length(ini.n.p))
  for(j in 1:length(ini.power))
  {
    ini.power[j]<-get.power.loglink(ini.n.p[j], n.v,PE_a,PE_0, HIV.inc.0,alpha)
  }
  
  if(ini.power[length(ini.n.p)]<beta){return(Inf)}else
  {
    n1<-ini.n.p[max(which(ini.power<=beta))]
    n2<-ini.n.p[max(which(ini.power<=beta))+1]
    
    while(n2>n1+1)
    {
      n3<-max(floor((n1+n2)/2),n1)
      power.new<-get.power.loglink(n3, n.v,PE_a,PE_0, HIV.inc.0,alpha)
      if(power.new>=beta)
      {
        n2 <- n3
      }else
      {
        n1<-n3
      }
    }
    return(n2)
  }
}

get.samplesize.loglink2<-function(beta, n.v,PE_a,PE_0, HIV.inc.0,alpha=0.05)
{
  ini.n.p<-10^(1:6)
  ini.power <- rep(NA,length(ini.n.p))
  for(j in 1:length(ini.power))
  {
    ini.power[j]<-get.power.loglink2(ini.n.p[j], n.v,PE_a,PE_0, HIV.inc.0,alpha)
  }
  
  if(ini.power[length(ini.n.p)]<beta){return(Inf)}else
  {
    n1<-ini.n.p[max(which(ini.power<=beta))]
    n2<-ini.n.p[max(which(ini.power<=beta))+1]
    
    while(n2>n1+1)
    {
      n3<-max(floor((n1+n2)/2),n1)
      power.new<-get.power.loglink2(n3, n.v,PE_a,PE_0, HIV.inc.0,alpha)
      if(power.new>=beta)
      {
        n2 <- n3
      }else
      {
        n1<-n3
      }
    }
    return(n2)
  }
}


##analytical, log-link, log-RR based CI
get.samplesize.log<-function(beta, n.v,PE_a,PE_0, HIV.inc.0,alpha=0.05)
{
  r0<-1-PE_0
  r1<-1-PE_a
  temp <- ( (log(r0)-log(r1))/(qnorm(1-alpha/2) + qnorm(beta/100)) )^2 - (1 - r1*HIV.inc.0)/( n.v * r1 * HIV.inc.0 )
  
  n.a <- (1-HIV.inc.0) / HIV.inc.0 / temp
  if(n.a<0){n.a<-Inf}
  return(ceiling(n.a))
}

##analytical, log-link, PE based CI
get.samplesize.logPE<-function(beta, n.v,PE_a,PE_0, HIV.inc.0,alpha=0.05)
{
  HIV.inc.v = (1 - PE_0)*HIV.inc.0
  
  temp<- ((PE_a-PE_0)/(qnorm(1-alpha/2) + qnorm(beta/100)))^2 - HIV.inc.v*(1-HIV.inc.v)/n.v/HIV.inc.0^2
  
  n.a <- HIV.inc.v^2 * (1-HIV.inc.0)/ temp / HIV.inc.0^3
  if(n.a<0){n.a<-Inf}
  return(ceiling(n.a))
}

n.active<-2000
PE_a.in<-.6
PE_0.in<-.3
# HIV.inc.p<-3/100
# HIV.inc.p<-4.5/100
HIV.inc.p<-6/100
alpha.in<-.05

# get.samplesize(60,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)
get.samplesize(70,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)
get.samplesize(80,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)
get.samplesize(90,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)

##PE CI
get.samplesize.loglink(60,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)
get.samplesize.loglink(70,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)
get.samplesize.loglink(80,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)
get.samplesize.loglink(90,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)

get.samplesize.logPE(60,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)
get.samplesize.logPE(70,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)
get.samplesize.logPE(80,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)
get.samplesize.logPE(90,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)


##PE CI = 1- exp(logRR CI)
# get.samplesize.loglink2(50,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)
# get.samplesize.loglink2(60,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)
get.samplesize.loglink2(70,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)
get.samplesize.loglink2(80,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)
get.samplesize.loglink2(90,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)

# get.samplesize.log(60,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)
get.samplesize.log(70,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)
get.samplesize.log(80,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)
get.samplesize.log(90,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)

# get.samplesize(99,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)
# get.samplesize.log(99,n.active,PE_a.in,PE_0.in,HIV.inc.p,alpha.in)
