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


critic.cov.norm1=function(a1,a2,b1,b2,s_A,s_B,u=sample(c(1:100),100,replace=T),alpha)
{
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
    
    for(i in 1:100)
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
    for(i in c(11:100))
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
    T1[j]=(a1.hat[100]-a2.hat[100])
    T2[j]=(b1.hat[100]-b2.hat[100])
    T3[j]=(s.A[100]-s.B[100])
    t.stat[j]=max(T1[j],T2[j],T3[j])
  }
  sorted=sort(t.stat)
  critic.val=sorted[5000-(alpha*5000)]
  return(critic.val)
}

power.cov.norm1=function(a1,a2,b1,b2,s_A,s_B,u=sample(c(1:100),100,replace=T),critic.val)
{
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
    
    for(i in 1:100)
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
    for(i in c(11:100))
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
    T1[j]=(a1.hat[100]-a2.hat[100])
    T2[j]=(b1.hat[100]-b2.hat[100])
    T3[j]=(s.A[100]-s.B[100])
    t.stat[j]=max(T1[j],T2[j],T3[j])
  }
  sorted=sort(t.stat)
  power=mean(sorted>critic.val)
  return(power)
}

critic.cov.norm2=function(a1,a2,b1,b2,s_A,s_B,u=sample(c(1:100),100,replace=T),alpha)
{
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
    
    for(i in 1:100)
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
    for(i in c(11:100))
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
    T1[j]=(a1.hat[100]-a2.hat[100])
    T2[j]=(b1.hat[100]-b2.hat[100])
    T3[j]=(s.A[100]-s.B[100])
    t.stat[j]=(T1[j]+T2[j]+T3[j])/3
  }
  sorted=sort(t.stat)
  critic.val=sorted[5000-(alpha*5000)]
  return(critic.val)
}

power.cov.norm2=function(a1,a2,b1,b2,s_A,s_B,u=sample(c(1:100),100,replace=T),critic.val)
{
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
    
    for(i in 1:100)
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
    for(i in c(11:100))
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
    T1[j]=(a1.hat[100]-a2.hat[100])
    T2[j]=(b1.hat[100]-b2.hat[100])
    T3[j]=(s.A[100]-s.B[100])
    t.stat[j]=(T1[j]+T2[j]+T3[j])/3
  }
  sorted=sort(t.stat)
  power=mean(sorted>critic.val)
  return(power)
}

#Power_Curve
#max_statistic
set.seed(100)
cr.val1=critic.cov.norm1(a1=2,a2=2,b1=3,b2=3,s_A=4,s_B=4,alpha=0.05)
cr.val1=5.258819
po10=power.cov.norm1(2,2,3,3,s_A=4,s_B=4,critic.val=cr.val1)
po11=power.cov.norm1(2,2,3,3,s_A=6,s_B=4,critic.val=cr.val1)
po12=power.cov.norm1(2,2,3,3,s_A=8,s_B=4,critic.val=cr.val1)
po13=power.cov.norm1(2,2,3,3,s_A=10,s_B=4,critic.val=cr.val1)
p10=0.0396
p11=0.1554
p12=0.3849
p13=0.53
plot(c(0,2,4,6),c(p10,p11,p12,p13),type="l",lty=1,lwd=2,xlab="Difference in sigma",ylab="Power",main="Power Function")

#Average_statistic
set.seed(100)
cr.val2=critic.cov.norm2(a1=2,a2=2,b1=3,b2=3,s_A=4,s_B=4,alpha=0.05)
cr.val2=1.553204
po20=power.cov.norm2(2,2,3,3,s_A=4,s_B=4,critic.val=cr.val2)
po21=power.cov.norm2(2,2,3,3,s_A=6,s_B=4,critic.val=cr.val2)
po22=power.cov.norm2(2,2,3,3,s_A=8,s_B=4,critic.val=cr.val2)
po23=power.cov.norm2(2,2,3,3,s_A=10,s_B=4,critic.val=cr.val2)
p20=0.08
p21=0.186
p22=0.4514
p23=0.583
lines(c(0,2,4,6),c(p20,p21,p22,p23),lty=2,lwd=2)
legend("bottomright",legend=c("Max stat","Average stat"),lty=c(1,2),lwd=2)
