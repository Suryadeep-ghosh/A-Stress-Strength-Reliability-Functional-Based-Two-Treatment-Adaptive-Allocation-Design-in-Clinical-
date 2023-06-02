factor.age=function(age)
{
  n=length(age)
  f=NULL
  for(i in 1:n)  
  {
    if(1<=age[i] & age[i]<=35)
      f[i]=0
    else if (36<=age[i] & age[i]<=70)
      f[i]=1
    else
      f[i]=2
  }
  return(f)
}
library(MASS)
func1=function(a1,a2,b1,b2,s_A,s_B,n)
{
  u=sample(1:n,n,replace=T)
  count.A=10
  count.B=10
  fac=factor.age(u)
  muA=NULL
  muB=NULL
  
  for(i in 1:n)
  {
    muA[i]=a1+(b1*fac[i])
    muB[i]=a2+(b2*fac[i])
  }
  c=45
  del=NULL
  del=c(rep(1,10),rep(0,10))
  
  X=NULL
  for(i in 1:10)
    X[i]=rcauchy(1,muA[i],s_A)
  for(i in 11:20)
    X[i]=rcauchy(1,muB[i],s_B)
  
  
  scale.A=NULL
  scale.A[20]=MASS::rlm(X[del==1]~fac[1:20][which(del==1)])$s
  scale.B=NULL
  scale.B[20]=MASS::rlm(X[del==0]~fac[1:20][which(del==0)])$s
  
  
  
  a1.hat=NULL
  a2.hat=NULL
  b1.hat=NULL
  b2.hat=NULL
  
  a1.hat[20]=MASS::rlm(X[del==1]~fac[1:20][which(del==1)])$coefficients[1]
  b1.hat[20]=MASS::rlm(X[del==1]~fac[1:20][which(del==1)])$coefficients[2]
  a2.hat[20]=MASS::rlm(X[del==0]~fac[1:20][which(del==0)])$coefficients[1]
  b2.hat[20]=MASS::rlm(X[del==0]~fac[1:20][which(del==0)])$coefficients[2]
  
  prob.A=NULL
  for(i in c(21:n))
  {
    if(u[i]>=c)
    {
      prob.A[i]=max(pcauchy((a1.hat[i-1]-a2.hat[i-1]+(b1.hat[i-1]-b2.hat[i-1])*fac[i])/sqrt(scale.A[i-1]+scale.B[i-1])),integrate(function (x) 2*log(x^2)/(pi^2*(x^2-1)), lower=0, upper=scale.A[i-1]/scale.B[i-1])$value)
      if(isTRUE(runif(1)<prob.A[i]))
      {
        X[i]=rcauchy(1,location=muA[i],scale=s_A)
        count.A=count.A+1
        del[i]=1
      }
      else
      {
        X[i]=rcauchy(1,location=muB[i],scale=s_B)
        count.B=count.B+1
        del[i]=0
      }
    }
    else
    {
      prob.A[i]=min(pcauchy((a1.hat[i-1]-a2.hat[i-1]+(b1.hat[i-1]-b2.hat[i-1])*fac[i])/sqrt(scale.A[i-1]+scale.B[i-1])),integrate(function (x) 2*log(x^2)/(pi^2*(x^2-1)), lower=0, upper=scale.A[i-1]/scale.B[i-1])$value)
      if(isTRUE(runif(1)<prob.A[i]))
      {
        X[i]=rcauchy(1,location=muA[i],scale=s_A)
        count.A=count.A+1
        del[i]=1
      }
      else
      {
        X[i]=rcauchy(1,location=muB[i],scale=s_B)
        count.B=count.B+1
        del[i]=0
      }
    }
    scale.A[i]=MASS::rlm(X[del==1]~fac[1:i][which(del==1)])$s
    scale.B[i]=MASS::rlm(X[del==0]~fac[1:i][which(del==0)])$s
    
    a1.hat[i]=MASS::rlm(X[del==1]~fac[1:i][which(del==1)])$coefficients[1]
    b1.hat[i]=MASS::rlm(X[del==1]~fac[1:i][which(del==1)])$coefficients[2]
    a2.hat[i]=MASS::rlm(X[del==0]~fac[1:i][which(del==0)])$coefficients[1]
    b2.hat[i]=MASS::rlm(X[del==0]~fac[1:i][which(del==0)])$coefficients[2]
    
  }
  l=list(a1.hat[100],a2.hat[100],b1.hat[100],b2.hat[100],scale.A[100],scale.B[100])
  names(l)=c("a1.h","a2.h","b1.h","b2.h","sA.h","sB.h")
  return(l)
}
t.stat.av=NULL
set.seed(100)
for(i in 2985:3000)
{
  f=func1(2,2,4,4,0.5,0.5,200)
  T1=f$a1.h-f$a2.h
  T2=f$b1.h-f$b2.h
  T3=f$sA.h-f$sB.h
  t.stat.av[i]=(T1+T2+T3)/3
}
alpha=0.05
sorted=sort(t.stat.av)
critic.val.av=sorted[3000-(alpha*3000)]


length(t.stat.av)
critic.val.av.150=0.1810199
critic.val.av.200=0.1935402