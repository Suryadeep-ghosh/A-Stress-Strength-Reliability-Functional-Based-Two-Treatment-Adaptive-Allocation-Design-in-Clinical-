clinical.int.norm=function(a1=6,a2=6,b1=2,b2=2,s_A=4,s_B=6,u=sample(c(1:100),100,replace=T))
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
  return(c(sum(del)/100,sum(1-del)/100))
}

result=matrix(0,nrow=5000,ncol=2)
for(i in 1:5000)
  result[i,]=clinical.int.norm()

propA=mean(result[,1])
sdA=sd(result[,1])
l=list(propA,sdA)
names(l)=c("Mean Ratio","Standard Error")
l
