critic1.cauchy=function(n,mu_A,mu_B,s_A,s_B,alpha)
{
  T1=NULL
  T2=NULL
  t.stat=NULL
  for(j in 1:5000)
  {
    u=sample(c(1:n),n,replace=T)
    c=45
    del=c(rep(1,5),rep(0,5))
    X=c(rcauchy(5,location=mu_A,scale=s_A),rcauchy(5,location=mu_B,scale=s_B))
    
    
    library(univariateML)
    mu.A.hat=mlcauchy(X[1:5])[1]
    mu.B.hat=mlcauchy(X[6:10])[1]
    
    s.A.hat=mlcauchy(X[1:5])[2]
    s.B.hat=mlcauchy(X[6:10])[2]
    
    m.hat.A=c()
    m.hat.B=c()
    s.A=c()
    s.B=c()
    m.hat.A[10]=mu.A.hat
    m.hat.B[10]=mu.B.hat
    s.A[10]=s.A.hat
    s.B[10]=s.B.hat
    prob.A=c()
    for(i in c(11:n))
    {
      if(u[i]>=c)
      {
        prob.A[i]=max(pcauchy((m.hat.A[i-1]-m.hat.B[i-1])/(s.A[i-1]+s.B[i-1])),integrate(function (x) 2*log(x^2)/(pi^2*(x^2-1)), lower=0, upper=s.A[i-1]/s.B[i-1])$value)
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rcauchy(1,location=mu_A,scale=s_A)
          del[i]=1
        }
        else
        {
          X[i]=rcauchy(1,location=mu_B,scale=s_B)
          del[i]=0
        }
      }
      else
      {
        prob.A[i]=min(pcauchy((m.hat.A[i-1]-m.hat.B[i-1])/(s.A[i-1]+s.B[i-1])),integrate(function (x) 2*log(x^2)/(pi^2*(x^2-1)), lower=0, upper=s.A[i-1]/s.B[i-1])$value)
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rcauchy(1,location=mu_A,scale=s_A)
          del[i]=1
        }
        else
        {
          X[i]=rcauchy(1,location=mu_B,scale=s_B)
          del[i]=0
        }
      }
      m.hat.A[i]=mlcauchy(X[del==1])[1]
      s.A[i]=mlcauchy(X[del==1])[2]
      m.hat.B[i]=mlcauchy(X[del==0])[1]
      s.B[i]=mlcauchy(X[del==0])[2]
    }
    T1[j]=(m.hat.A[n]-m.hat.B[n])
    T2[j]=(s.A[n]-s.B[n])
    t.stat[j]=max(T1[j],T2[j])
  }
  sorted=sort(t.stat)
  critic.val=sorted[5000-(alpha*5000)]
  return(critic.val)
}

power1.cauchy=function(n,mu_A,mu_B,s_A,s_B,critic.val)
{
  T1=NULL
  T2=NULL
  t.stat=NULL
  for(j in 1:5000)
  {
    u=sample(c(1:n),n,replace=T)
    c=45
    del=c(rep(1,5),rep(0,5))
    X=c(rcauchy(5,location=mu_A,scale=s_A),rcauchy(5,location=mu_B,scale=s_B))
    
    
    library(univariateML)
    mu.A.hat=mlcauchy(X[1:5])[1]
    mu.B.hat=mlcauchy(X[6:10])[1]
    
    s.A.hat=mlcauchy(X[1:5])[2]
    s.B.hat=mlcauchy(X[6:10])[2]
    
    m.hat.A=c()
    m.hat.B=c()
    s.A=c()
    s.B=c()
    m.hat.A[10]=mu.A.hat
    m.hat.B[10]=mu.B.hat
    s.A[10]=s.A.hat
    s.B[10]=s.B.hat
    prob.A=c()
    for(i in c(11:n))
    {
      if(u[i]>=c)
      {
        prob.A[i]=max(pcauchy((m.hat.A[i-1]-m.hat.B[i-1])/(s.A[i-1]+s.B[i-1])),integrate(function (x) 2*log(x^2)/(pi^2*(x^2-1)), lower=0, upper=s.A[i-1]/s.B[i-1])$value)
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rcauchy(1,location=mu_A,scale=s_A)
          del[i]=1
        }
        else
        {
          X[i]=rcauchy(1,location=mu_B,scale=s_B)
          del[i]=0
        }
      }
      else
      {
        prob.A[i]=min(pcauchy((m.hat.A[i-1]-m.hat.B[i-1])/(s.A[i-1]+s.B[i-1])),integrate(function (x) 2*log(x^2)/(pi^2*(x^2-1)), lower=0, upper=s.A[i-1]/s.B[i-1])$value)
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rcauchy(1,location=mu_A,scale=s_A)
          del[i]=1
        }
        else
        {
          X[i]=rcauchy(1,location=mu_B,scale=s_B)
          del[i]=0
        }
      }
      m.hat.A[i]=mlcauchy(X[del==1])[1]
      s.A[i]=mlcauchy(X[del==1])[2]
      m.hat.B[i]=mlcauchy(X[del==0])[1]
      s.B[i]=mlcauchy(X[del==0])[2]
    }
    T1[j]=(m.hat.A[n]-m.hat.B[n])
    T2[j]=(s.A[n]-s.B[n])
    t.stat[j]=max(T1[j],T2[j])
  }
  sorted=sort(t.stat)
  power=mean(sorted>critic.val)
  return(power)
}

critic2.cauchy=function(n,mu_A,mu_B,s_A,s_B,alpha)
{
  T1=NULL
  T2=NULL
  t.stat=NULL
  for(j in 1:5000)
  {
    u=sample(c(1:n),n,replace=T)
    c=45
    del=c(rep(1,5),rep(0,5))
    X=c(rcauchy(5,location=mu_A,scale=s_A),rcauchy(5,location=mu_B,scale=s_B))
    
    
    library(univariateML)
    mu.A.hat=mlcauchy(X[1:5])[1]
    mu.B.hat=mlcauchy(X[6:10])[1]
    
    s.A.hat=mlcauchy(X[1:5])[2]
    s.B.hat=mlcauchy(X[6:10])[2]
    
    m.hat.A=c()
    m.hat.B=c()
    s.A=c()
    s.B=c()
    m.hat.A[10]=mu.A.hat
    m.hat.B[10]=mu.B.hat
    s.A[10]=s.A.hat
    s.B[10]=s.B.hat
    prob.A=c()
    for(i in c(11:n))
    {
      if(u[i]>=c)
      {
        prob.A[i]=max(pcauchy((m.hat.A[i-1]-m.hat.B[i-1])/(s.A[i-1]+s.B[i-1])),integrate(function (x) 2*log(x^2)/(pi^2*(x^2-1)), lower=0, upper=s.A[i-1]/s.B[i-1])$value)
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rcauchy(1,location=mu_A,scale=s_A)
          del[i]=1
        }
        else
        {
          X[i]=rcauchy(1,location=mu_B,scale=s_B)
          del[i]=0
        }
      }
      else
      {
        prob.A[i]=min(pcauchy((m.hat.A[i-1]-m.hat.B[i-1])/(s.A[i-1]+s.B[i-1])),integrate(function (x) 2*log(x^2)/(pi^2*(x^2-1)), lower=0, upper=s.A[i-1]/s.B[i-1])$value)
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rcauchy(1,location=mu_A,scale=s_A)
          del[i]=1
        }
        else
        {
          X[i]=rcauchy(1,location=mu_B,scale=s_B)
          del[i]=0
        }
      }
      m.hat.A[i]=mlcauchy(X[del==1])[1]
      s.A[i]=mlcauchy(X[del==1])[2]
      m.hat.B[i]=mlcauchy(X[del==0])[1]
      s.B[i]=mlcauchy(X[del==0])[2]
    }
    T1[j]=(m.hat.A[n]-m.hat.B[n])
    T2[j]=(s.A[n]-s.B[n])
    t.stat[j]=(T1[j]+T2[j])/2
  }
  sorted=sort(t.stat)
  critic.val=sorted[5000-(alpha*5000)]
  return(critic.val)
}

power2.cauchy=function(n,mu_A,mu_B,s_A,s_B,critic.val)
{
  T1=NULL
  T2=NULL
  t.stat=NULL
  for(j in 1:5000)
  {
    u=sample(c(1:n),n,replace=T)
    c=45
    del=c(rep(1,5),rep(0,5))
    X=c(rcauchy(5,location=mu_A,scale=s_A),rcauchy(5,location=mu_B,scale=s_B))
    
    
    library(univariateML)
    mu.A.hat=mlcauchy(X[1:5])[1]
    mu.B.hat=mlcauchy(X[6:10])[1]
    
    s.A.hat=mlcauchy(X[1:5])[2]
    s.B.hat=mlcauchy(X[6:10])[2]
    
    m.hat.A=c()
    m.hat.B=c()
    s.A=c()
    s.B=c()
    m.hat.A[10]=mu.A.hat
    m.hat.B[10]=mu.B.hat
    s.A[10]=s.A.hat
    s.B[10]=s.B.hat
    prob.A=c()
    for(i in c(11:n))
    {
      if(u[i]>=c)
      {
        prob.A[i]=max(pcauchy((m.hat.A[i-1]-m.hat.B[i-1])/(s.A[i-1]+s.B[i-1])),integrate(function (x) 2*log(x^2)/(pi^2*(x^2-1)), lower=0, upper=s.A[i-1]/s.B[i-1])$value)
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rcauchy(1,location=mu_A,scale=s_A)
          del[i]=1
        }
        else
        {
          X[i]=rcauchy(1,location=mu_B,scale=s_B)
          del[i]=0
        }
      }
      else
      {
        prob.A[i]=min(pcauchy((m.hat.A[i-1]-m.hat.B[i-1])/(s.A[i-1]+s.B[i-1])),integrate(function (x) 2*log(x^2)/(pi^2*(x^2-1)), lower=0, upper=s.A[i-1]/s.B[i-1])$value)
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rcauchy(1,location=mu_A,scale=s_A)
          del[i]=1
        }
        else
        {
          X[i]=rcauchy(1,location=mu_B,scale=s_B)
          del[i]=0
        }
      }
      m.hat.A[i]=mlcauchy(X[del==1])[1]
      s.A[i]=mlcauchy(X[del==1])[2]
      m.hat.B[i]=mlcauchy(X[del==0])[1]
      s.B[i]=mlcauchy(X[del==0])[2]
    }
    T1[j]=(m.hat.A[n]-m.hat.B[n])
    T2[j]=(s.A[n]-s.B[n])
    t.stat[j]=(T1[j]+T2[j])/2
  }
  sorted=sort(t.stat)
  power=mean(sorted>critic.val)
  return(power)
}

critic3.cauchy=function(n,mu_A,mu_B,s_A,s_B,alpha)
{
  Q=NULL
  for(j in 1:5000)
  {
    u=sample(c(1:n),n,replace=T)
    c=45
    del=c(rep(1,5),rep(0,5))
    X=c(rcauchy(5,location=mu_A,scale=s_A),rcauchy(5,location=mu_B,scale=s_B))
    
    
    library(univariateML)
    mu.A.hat=mlcauchy(X[1:5])[1]
    mu.B.hat=mlcauchy(X[6:10])[1]
    
    s.A.hat=mlcauchy(X[1:5])[2]
    s.B.hat=mlcauchy(X[6:10])[2]
    
    m.hat.A=c()
    m.hat.B=c()
    s.A=c()
    s.B=c()
    m.hat.A[10]=mu.A.hat
    m.hat.B[10]=mu.B.hat
    s.A[10]=s.A.hat
    s.B[10]=s.B.hat
    prob.A=c()
    for(i in c(11:n))
    {
      if(u[i]>=c)
      {
        prob.A[i]=max(pcauchy((m.hat.A[i-1]-m.hat.B[i-1])/(s.A[i-1]+s.B[i-1])),integrate(function (x) 2*log(x^2)/(pi^2*(x^2-1)), lower=0, upper=s.A[i-1]/s.B[i-1])$value)
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rcauchy(1,location=mu_A,scale=s_A)
          del[i]=1
        }
        else
        {
          X[i]=rcauchy(1,location=mu_B,scale=s_B)
          del[i]=0
        }
      }
      else
      {
        prob.A[i]=min(pcauchy((m.hat.A[i-1]-m.hat.B[i-1])/(s.A[i-1]+s.B[i-1])),integrate(function (x) 2*log(x^2)/(pi^2*(x^2-1)), lower=0, upper=s.A[i-1]/s.B[i-1])$value)
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rcauchy(1,location=mu_A,scale=s_A)
          del[i]=1
        }
        else
        {
          X[i]=rcauchy(1,location=mu_B,scale=s_B)
          del[i]=0
        }
      }
      m.hat.A[i]=mlcauchy(X[del==1])[1]
      s.A[i]=mlcauchy(X[del==1])[2]
      m.hat.B[i]=mlcauchy(X[del==0])[1]
      s.B[i]=mlcauchy(X[del==0])[2]
    }
    D1=(m.hat.A[n]-m.hat.B[n])/sd(m.hat.A[11:n]-m.hat.B[11:n])
    D2=(s.A[n]-s.B[n])/sd(s.A[11:n]-s.B[11:n])
    r=cor(m.hat.A[11:n]-m.hat.B[11:n],s.A[11:n]-s.B[11:n])
    
    if(D1>0 & D2>0){
      Q[j]=(D1^2+D2^2-2*r*D1*D2)/sqrt(1-r^2)
    }else if(D2>D1 & D1<0){
      Q[j]=(D2-r*D1)/sqrt(1-r^2)
    }else if(D1>D2 & D2<0){
      Q[j]=(D1-r*D2)/sqrt(1-r^2)
    }else{
      Q[j]=0}
    
  }
  sorted=sort(Q)
  critic.val=sorted[5000-(alpha*5000)]
  return(critic.val)
}

power3.cauchy=function(n,mu_A,mu_B,s_A,s_B,critic.val)
{
  Q=NULL
  for(j in 1:5000)
  {
    u=sample(c(1:n),n,replace=T)
    c=45
    del=c(rep(1,5),rep(0,5))
    X=c(rcauchy(5,location=mu_A,scale=s_A),rcauchy(5,location=mu_B,scale=s_B))
    
    
    library(univariateML)
    mu.A.hat=mlcauchy(X[1:5])[1]
    mu.B.hat=mlcauchy(X[6:10])[1]
    
    s.A.hat=mlcauchy(X[1:5])[2]
    s.B.hat=mlcauchy(X[6:10])[2]
    
    m.hat.A=c()
    m.hat.B=c()
    s.A=c()
    s.B=c()
    m.hat.A[10]=mu.A.hat
    m.hat.B[10]=mu.B.hat
    s.A[10]=s.A.hat
    s.B[10]=s.B.hat
    prob.A=c()
    for(i in c(11:n))
    {
      if(u[i]>=c)
      {
        prob.A[i]=max(pcauchy((m.hat.A[i-1]-m.hat.B[i-1])/(s.A[i-1]+s.B[i-1])),integrate(function (x) 2*log(x^2)/(pi^2*(x^2-1)), lower=0, upper=s.A[i-1]/s.B[i-1])$value)
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rcauchy(1,location=mu_A,scale=s_A)
          del[i]=1
        }
        else
        {
          X[i]=rcauchy(1,location=mu_B,scale=s_B)
          del[i]=0
        }
      }
      else
      {
        prob.A[i]=min(pcauchy((m.hat.A[i-1]-m.hat.B[i-1])/(s.A[i-1]+s.B[i-1])),integrate(function (x) 2*log(x^2)/(pi^2*(x^2-1)), lower=0, upper=s.A[i-1]/s.B[i-1])$value)
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rcauchy(1,location=mu_A,scale=s_A)
          del[i]=1
        }
        else
        {
          X[i]=rcauchy(1,location=mu_B,scale=s_B)
          del[i]=0
        }
      }
      m.hat.A[i]=mlcauchy(X[del==1])[1]
      s.A[i]=mlcauchy(X[del==1])[2]
      m.hat.B[i]=mlcauchy(X[del==0])[1]
      s.B[i]=mlcauchy(X[del==0])[2]
    }
    D1=(m.hat.A[n]-m.hat.B[n])/sd(m.hat.A[11:n]-m.hat.B[11:n])
    D2=(s.A[n]-s.B[n])/sd(s.A[11:n]-s.B[11:n])
    r=cor(m.hat.A[11:n]-m.hat.B[11:n],s.A[11:n]-s.B[11:n])
    
    if(D1>0 & D2>0){
      Q[j]=(D1^2+D2^2-2*r*D1*D2)/sqrt(1-r^2)
    }else if(D2>D1 & D1<0){
      Q[j]=(D2-r*D1)/sqrt(1-r^2)
    }else if(D1>D2 & D2<0){
      Q[j]=(D1-r*D2)/sqrt(1-r^2)
    }else{
      Q[j]=0}
    
  }
  sorted=sort(Q)
  power=mean(sorted>critic.val)
  return(power)
}
#max statistic
set.seed(100)
cr.val=critic1.cauchy(100,5,5,2,2,0.05)
cr.val=1.105023
mA=c(5,6,7,8)
sA=seq(2,4,by=0.1)
po=NULL
for(i in c(1:length(mA)))
{
  set.seed(100)
  po[i]=power1.cauchy(100,mA[i],5,2,2,cr.val)
}
po_val=c(0.0500,0.4554,0.9328,0.9992)
#Average_stat
cr.val2=critic2.cauchy(100,5,5,2,2,0.05)
cr.val2=0.6715186
po2=NULL
for(i in c(1:length(mA)))
{
  set.seed(100)
  po2[i]=power2.cauchy(100,mA[i],5,2,2,cr.val2)
}
po2_val=c(0.0492,0.3480,0.7988,0.9716)
#Q_Stat
cr.val3=critic3.cauchy(100,5,5,2,2,0.05)
cr.val3=6.605703
po3=NULL
for(i in c(1:length(mA)))
{
  set.seed(100)
  po3[i]=power3.cauchy(100,mA[i],5,2,2,cr.val3)
}
po3_val=c(0.0406,0.2434,0.5306,0.7260)

plot(mA-5,po_val,type="l",lty=1,main="Power Curve",xlab="Difference in mu",ylab="Power",lwd=2)
lines(mA-5,po2_val,lty=2,lwd=2)
lines(mA-5,po3_val,lty=15,lwd=2)
legend("bottomright", c("Max","Average","Q-Stat"),lty=c(1,2,15),lwd=3)