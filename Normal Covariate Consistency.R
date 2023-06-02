factor.age=function(age)
{
  n=length(age)
  f=NULL
  for(i in 1:n)  
  {
    if(1<=age[i] & age[i]<=18)
      f[i]=0
    else if (19<=age[i] & age[i]<=45)
      f[i]=1
    else if (46<=age[i] & age[i]<=75)
      f[i]=2
    else
      f[i]=3
  }
  return(f)
}


critic.cov.norm1=function(n=100,a1,a2,b1,b2,s_A,s_B,alpha)
{
  u=sample(c(1:100),size=n,replace=T)
  T1=NULL
  T2=NULL
  T3=NULL
  t.stat=NULL
  for(j in 1:5000)
  {
    count.A=5
    count.B=5
    fac=factor.age(u)
    muA=NULL
    muB=NULL
    
    for(i in 1:n)
    {
      if(fac[i]==0)
      {muA[i]=a1
      muB[i]=a2
      }
      else if(fac[i]==1)
      {
        muA[i]=a1+b1
        muB[i]=a2+b2
      }
      else if(fac[i]==2)
      {
        muA[i]=a1+2*b1
        muB[i]=a2+2*b2
      }
      else
      {
        muA[i]=a1+3*b1
        muB[i]=a2+3*b2
        
      }
    }
    c=45
    del=NULL
    del[1:10]=c(rep(1,5),rep(0,5))
    
    
    X=NULL
    for(i in 1:5)
      X[i]=rnorm(1,muA[i],sqrt(s_A))
    for(i in 6:10)
      X[i]=rnorm(1,muB[i],sqrt(s_B))
    
    #estimates of a alphas and betas
    a1.hat=NULL
    a2.hat=NULL
    b1.hat=NULL
    b2.hat=NULL
    s.A=NULL
    s.B=NULL
    
    a1.hat[10]=lm(X[which(del==1)]~fac[which(del==1)])$coefficients[1]
    a2.hat[10]=lm(X[which(del==0)]~fac[which(del==0)])$coefficients[1]
    b1.hat[10]=lm(X[which(del==1)]~fac[which(del==1)])$coefficients[2]
    b2.hat[10]=lm(X[which(del==0)]~fac[which(del==0)])$coefficients[2]
    
    s.A[10]=(1/(count.A*(count.A-1)))*sum((as.vector(dist(X[del==1],method="manhattan")))^2)
    s.B[10]=(1/(count.B*(count.B-1)))*sum((as.vector(dist(X[del==0],method="manhattan")))^2)
    
    prob.A=c()
    for(i in c(11:n))
    {
      if(u[i]>=c)
      {
        prob.A[i]=max(pnorm((a1.hat[i-1]-a2.hat[i-1]+(b1.hat[i-1]-b2.hat[i-1])*fac[i])/sqrt(s.A[i-1]+s.B[i-1])),(2/pi)*asin(sqrt(s.A[i-1]/(s.A[i-1]+s.B[i-1]))))
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rnorm(1,mean=muA[i],sd=sqrt(s_A))
          count.A=count.A+1
          del[i]=1
        }
        else
        {
          X[i]=rnorm(1,mean=muB[i],sd=sqrt(s_B))
          count.B=count.B+1
          del[i]=0
        }
      }
      else
      {
        prob.A[i]=min(pnorm((a1.hat[i-1]-a2.hat[i-1]+(b1.hat[i-1]-b2.hat[i-1])*fac[i])/sqrt(s.A[i-1]+s.B[i-1])),(2/pi)*asin(sqrt(s.A[i-1]/(s.A[i-1]+s.B[i-1]))))
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rnorm(1,mean=muA[i],sd=sqrt(s_A))
          count.A=count.A+1
          del[i]=1
        }
        else
        {
          X[i]=rnorm(1,mean=muB[i],sd=sqrt(s_B))
          count.B=count.B+1
          del[i]=0
        }
      }
      
      a1.hat[i]=lm(X[which(del==1)]~fac[which(del==1)])$coefficients[1]
      a2.hat[i]=lm(X[which(del==0)]~fac[which(del==0)])$coefficients[1]
      b1.hat[i]=lm(X[which(del==1)]~fac[which(del==1)])$coefficients[2]
      b2.hat[i]=lm(X[which(del==0)]~fac[which(del==0)])$coefficients[2]
      
      s.A[i]=(1/(count.A*(count.A-1)))*sum((as.vector(dist(X[del==1],method="manhattan")))^2)
      s.B[i]=(1/(count.B*(count.B-1)))*sum((as.vector(dist(X[del==0],method="manhattan")))^2)
    }
    T1[j]=(a1.hat[n]-a2.hat[n])
    T2[j]=(b1.hat[n]-b2.hat[n])
    T3[j]=(s.A[n]-s.B[n])
    t.stat[j]=max(T1[j],T2[j],T3[j])
  }
  sorted=sort(t.stat)
  critic.val=sorted[5000-(alpha*5000)]
  return(critic.val)
}

power.cov.norm1=function(n=100,a1,a2,b1,b2,s_A,s_B,critic.val)
{
  u=sample(c(1:100),size=n,replace=T)
  T1=NULL
  T2=NULL
  T3=NULL
  t.stat=NULL
  for(j in 1:5000)
  {
    count.A=5
    count.B=5
    fac=factor.age(u)
    muA=NULL
    muB=NULL
    
    for(i in 1:n)
    {
      if(fac[i]==0)
      {muA[i]=a1
      muB[i]=a2
      }
      else if(fac[i]==1)
      {
        muA[i]=a1+b1
        muB[i]=a2+b2
      }
      else if(fac[i]==2)
      {
        muA[i]=a1+2*b1
        muB[i]=a2+2*b2
      }
      else
      {
        muA[i]=a1+3*b1
        muB[i]=a2+3*b2
        
      }
    }
    c=45
    del=NULL
    del[1:10]=c(rep(1,5),rep(0,5))
    
    
    X=NULL
    for(i in 1:5)
      X[i]=rnorm(1,muA[i],sqrt(s_A))
    for(i in 6:10)
      X[i]=rnorm(1,muB[i],sqrt(s_B))
    
    #estimates of a alphas and betas
    a1.hat=NULL
    a2.hat=NULL
    b1.hat=NULL
    b2.hat=NULL
    s.A=NULL
    s.B=NULL
    
    a1.hat[10]=lm(X[which(del==1)]~fac[which(del==1)])$coefficients[1]
    a2.hat[10]=lm(X[which(del==0)]~fac[which(del==0)])$coefficients[1]
    b1.hat[10]=lm(X[which(del==1)]~fac[which(del==1)])$coefficients[2]
    b2.hat[10]=lm(X[which(del==0)]~fac[which(del==0)])$coefficients[2]
    
    s.A[10]=(1/(count.A*(count.A-1)))*sum((as.vector(dist(X[del==1],method="manhattan")))^2)
    s.B[10]=(1/(count.B*(count.B-1)))*sum((as.vector(dist(X[del==0],method="manhattan")))^2)
    
    prob.A=c()
    for(i in c(11:n))
    {
      if(u[i]>=c)
      {
        prob.A[i]=max(pnorm((a1.hat[i-1]-a2.hat[i-1]+(b1.hat[i-1]-b2.hat[i-1])*fac[i])/sqrt(s.A[i-1]+s.B[i-1])),(2/pi)*asin(sqrt(s.A[i-1]/(s.A[i-1]+s.B[i-1]))))
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rnorm(1,mean=muA[i],sd=sqrt(s_A))
          count.A=count.A+1
          del[i]=1
        }
        else
        {
          X[i]=rnorm(1,mean=muB[i],sd=sqrt(s_B))
          count.B=count.B+1
          del[i]=0
        }
      }
      else
      {
        prob.A[i]=min(pnorm((a1.hat[i-1]-a2.hat[i-1]+(b1.hat[i-1]-b2.hat[i-1])*fac[i])/sqrt(s.A[i-1]+s.B[i-1])),(2/pi)*asin(sqrt(s.A[i-1]/(s.A[i-1]+s.B[i-1]))))
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rnorm(1,mean=muA[i],sd=sqrt(s_A))
          count.A=count.A+1
          del[i]=1
        }
        else
        {
          X[i]=rnorm(1,mean=muB[i],sd=sqrt(s_B))
          count.B=count.B+1
          del[i]=0
        }
      }
      
      a1.hat[i]=lm(X[which(del==1)]~fac[which(del==1)])$coefficients[1]
      a2.hat[i]=lm(X[which(del==0)]~fac[which(del==0)])$coefficients[1]
      b1.hat[i]=lm(X[which(del==1)]~fac[which(del==1)])$coefficients[2]
      b2.hat[i]=lm(X[which(del==0)]~fac[which(del==0)])$coefficients[2]
      
      s.A[i]=(1/(count.A*(count.A-1)))*sum((as.vector(dist(X[del==1],method="manhattan")))^2)
      s.B[i]=(1/(count.B*(count.B-1)))*sum((as.vector(dist(X[del==0],method="manhattan")))^2)
    }
    T1[j]=(a1.hat[n]-a2.hat[n])
    T2[j]=(b1.hat[n]-b2.hat[n])
    T3[j]=(s.A[n]-s.B[n])
    t.stat[j]=max(T1[j],T2[j],T3[j])
  }
  sorted=sort(t.stat)
  power=mean(sorted>critic.val)
  return(power)
}

critic.cov.norm2=function(n=100,a1,a2,b1,b2,s_A,s_B,alpha)
{
  u=sample(c(1:100),n,replace=T)
  T1=NULL
  T2=NULL
  T3=NULL
  t.stat=NULL
  for(j in 1:5000)
  {
    count.A=5
    count.B=5
    fac=factor.age(u)
    muA=NULL
    muB=NULL
    
    for(i in 1:n)
    {
      if(fac[i]==0)
      {muA[i]=a1
      muB[i]=a2
      }
      else if(fac[i]==1)
      {
        muA[i]=a1+b1
        muB[i]=a2+b2
      }
      else if(fac[i]==2)
      {
        muA[i]=a1+2*b1
        muB[i]=a2+2*b2
      }
      else
      {
        muA[i]=a1+3*b1
        muB[i]=a2+3*b2
        
      }
    }
    c=45
    del=NULL
    del[1:10]=c(rep(1,5),rep(0,5))
    
    
    X=NULL
    for(i in 1:5)
      X[i]=rnorm(1,muA[i],sqrt(s_A))
    for(i in 6:10)
      X[i]=rnorm(1,muB[i],sqrt(s_B))
    
    #estimates of a alphas and betas
    a1.hat=NULL
    a2.hat=NULL
    b1.hat=NULL
    b2.hat=NULL
    s.A=NULL
    s.B=NULL
    
    a1.hat[10]=lm(X[which(del==1)]~fac[which(del==1)])$coefficients[1]
    a2.hat[10]=lm(X[which(del==0)]~fac[which(del==0)])$coefficients[1]
    b1.hat[10]=lm(X[which(del==1)]~fac[which(del==1)])$coefficients[2]
    b2.hat[10]=lm(X[which(del==0)]~fac[which(del==0)])$coefficients[2]
    
    s.A[10]=(1/(count.A*(count.A-1)))*sum((as.vector(dist(X[del==1],method="manhattan")))^2)
    s.B[10]=(1/(count.B*(count.B-1)))*sum((as.vector(dist(X[del==0],method="manhattan")))^2)
    
    prob.A=c()
    for(i in c(11:n))
    {
      if(u[i]>=c)
      {
        prob.A[i]=max(pnorm((a1.hat[i-1]-a2.hat[i-1]+(b1.hat[i-1]-b2.hat[i-1])*fac[i])/sqrt(s.A[i-1]+s.B[i-1])),(2/pi)*asin(sqrt(s.A[i-1]/(s.A[i-1]+s.B[i-1]))))
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rnorm(1,mean=muA[i],sd=sqrt(s_A))
          count.A=count.A+1
          del[i]=1
        }
        else
        {
          X[i]=rnorm(1,mean=muB[i],sd=sqrt(s_B))
          count.B=count.B+1
          del[i]=0
        }
      }
      else
      {
        prob.A[i]=min(pnorm((a1.hat[i-1]-a2.hat[i-1]+(b1.hat[i-1]-b2.hat[i-1])*fac[i])/sqrt(s.A[i-1]+s.B[i-1])),(2/pi)*asin(sqrt(s.A[i-1]/(s.A[i-1]+s.B[i-1]))))
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rnorm(1,mean=muA[i],sd=sqrt(s_A))
          count.A=count.A+1
          del[i]=1
        }
        else
        {
          X[i]=rnorm(1,mean=muB[i],sd=sqrt(s_B))
          count.B=count.B+1
          del[i]=0
        }
      }
      
      a1.hat[i]=lm(X[which(del==1)]~fac[which(del==1)])$coefficients[1]
      a2.hat[i]=lm(X[which(del==0)]~fac[which(del==0)])$coefficients[1]
      b1.hat[i]=lm(X[which(del==1)]~fac[which(del==1)])$coefficients[2]
      b2.hat[i]=lm(X[which(del==0)]~fac[which(del==0)])$coefficients[2]
      
      s.A[i]=(1/(count.A*(count.A-1)))*sum((as.vector(dist(X[del==1],method="manhattan")))^2)
      s.B[i]=(1/(count.B*(count.B-1)))*sum((as.vector(dist(X[del==0],method="manhattan")))^2)
    }
    T1[j]=(a1.hat[n]-a2.hat[n])
    T2[j]=(b1.hat[n]-b2.hat[n])
    T3[j]=(s.A[n]-s.B[n])
    t.stat[j]=(T1[j]+T2[j]+T3[j])/3
  }
  sorted=sort(t.stat)
  critic.val=sorted[5000-(alpha*5000)]
  return(critic.val)
}

power.cov.norm2=function(n=100,a1,a2,b1,b2,s_A,s_B,critic.val)
{
  u=sample(c(1:100),n,replace=T)
  T1=NULL
  T2=NULL
  T3=NULL
  t.stat=NULL
  for(j in 1:5000)
  {
    count.A=5
    count.B=5
    fac=factor.age(u)
    muA=NULL
    muB=NULL
    
    for(i in 1:n)
    {
      if(fac[i]==0)
      {muA[i]=a1
      muB[i]=a2
      }
      else if(fac[i]==1)
      {
        muA[i]=a1+b1
        muB[i]=a2+b2
      }
      else if(fac[i]==2)
      {
        muA[i]=a1+2*b1
        muB[i]=a2+2*b2
      }
      else
      {
        muA[i]=a1+3*b1
        muB[i]=a2+3*b2
        
      }
    }
    c=45
    del=NULL
    del[1:10]=c(rep(1,5),rep(0,5))
    
    
    X=NULL
    for(i in 1:5)
      X[i]=rnorm(1,muA[i],sqrt(s_A))
    for(i in 6:10)
      X[i]=rnorm(1,muB[i],sqrt(s_B))
    
    #estimates of a alphas and betas
    a1.hat=NULL
    a2.hat=NULL
    b1.hat=NULL
    b2.hat=NULL
    s.A=NULL
    s.B=NULL
    
    a1.hat[10]=lm(X[which(del==1)]~fac[which(del==1)])$coefficients[1]
    a2.hat[10]=lm(X[which(del==0)]~fac[which(del==0)])$coefficients[1]
    b1.hat[10]=lm(X[which(del==1)]~fac[which(del==1)])$coefficients[2]
    b2.hat[10]=lm(X[which(del==0)]~fac[which(del==0)])$coefficients[2]
    
    s.A[10]=(1/(count.A*(count.A-1)))*sum((as.vector(dist(X[del==1],method="manhattan")))^2)
    s.B[10]=(1/(count.B*(count.B-1)))*sum((as.vector(dist(X[del==0],method="manhattan")))^2)
    
    prob.A=c()
    for(i in c(11:n))
    {
      if(u[i]>=c)
      {
        prob.A[i]=max(pnorm((a1.hat[i-1]-a2.hat[i-1]+(b1.hat[i-1]-b2.hat[i-1])*fac[i])/sqrt(s.A[i-1]+s.B[i-1])),(2/pi)*asin(sqrt(s.A[i-1]/(s.A[i-1]+s.B[i-1]))))
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rnorm(1,mean=muA[i],sd=sqrt(s_A))
          count.A=count.A+1
          del[i]=1
        }
        else
        {
          X[i]=rnorm(1,mean=muB[i],sd=sqrt(s_B))
          count.B=count.B+1
          del[i]=0
        }
      }
      else
      {
        prob.A[i]=min(pnorm((a1.hat[i-1]-a2.hat[i-1]+(b1.hat[i-1]-b2.hat[i-1])*fac[i])/sqrt(s.A[i-1]+s.B[i-1])),(2/pi)*asin(sqrt(s.A[i-1]/(s.A[i-1]+s.B[i-1]))))
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rnorm(1,mean=muA[i],sd=sqrt(s_A))
          count.A=count.A+1
          del[i]=1
        }
        else
        {
          X[i]=rnorm(1,mean=muB[i],sd=sqrt(s_B))
          count.B=count.B+1
          del[i]=0
        }
      }
      
      a1.hat[i]=lm(X[which(del==1)]~fac[which(del==1)])$coefficients[1]
      a2.hat[i]=lm(X[which(del==0)]~fac[which(del==0)])$coefficients[1]
      b1.hat[i]=lm(X[which(del==1)]~fac[which(del==1)])$coefficients[2]
      b2.hat[i]=lm(X[which(del==0)]~fac[which(del==0)])$coefficients[2]
      
      s.A[i]=(1/(count.A*(count.A-1)))*sum((as.vector(dist(X[del==1],method="manhattan")))^2)
      s.B[i]=(1/(count.B*(count.B-1)))*sum((as.vector(dist(X[del==0],method="manhattan")))^2)
    }
    T1[j]=(a1.hat[n]-a2.hat[n])
    T2[j]=(b1.hat[n]-b2.hat[n])
    T3[j]=(s.A[n]-s.B[n])
    t.stat[j]=(T1[j]+T2[j]+T3[j])/3
  }
  sorted=sort(t.stat)
  power=mean(sorted>critic.val)
  return(power)
}
s=c(200,300)
cr1=NULL
p1=NULL
cr2=NULL
p2=NULL
cr1=c(5.258819,3.569499,2.655681)
p1=c(0.0396,0.2030,0.3188)
cr2=c(1.553204,0.9945917,0.7372290)
p2=c(0.08,0.3154,0.4252)
cr2[3]=critic.cov.norm2(n=300,a1=2,a2=2,b1=3,b2=3,s_A=4,s_B=4,alpha=0.05)
p2[3]=power.cov.norm2(n=300,a1=2,a2=2,b1=3,b2=3,s_A=6,s_B=4,cr2[3])
plot(c(100,200,300),p1,type="l",lty=1,lwd=2,main="Consistency of Test Statistics",xlab="Sample Size",ylab="Power",ylim=c(0,0.5))
lines(c(100,200,300),p2,lty=2,lwd=2)
legend("bottomright", c("Max","Average"),lty=c(1,2),lwd=3)
