critic1.norm=function(n,mu_A,mu_B,s_A,s_B,alpha)
{
  T1=NULL
  T2=NULL
  t.stat=NULL
  for(j in 1:5000)
  {
    count.A=5
    count.B=5
    u=sample(c(1:n),n,replace=T)
    c=45
    del=c(rep(1,5),rep(0,5))
    X=c(rnorm(5,mean=mu_A,sd=sqrt(s_A)),rnorm(5,mean = mu_B,sd= sqrt(s_B)))
    
    #estimate of the mean response and the variation of treatment A
    mu.hat.A=sum(del*X)/count.A
    s.hat.A=(sum(del*X^2)/count.A)-(mu.hat.A)^2
    mu.hat.B=sum((1-del)*X)/count.B
    s.hat.B=(sum((1-del)*X^2)/count.B)-(mu.hat.B)^2
    
    m.hat.A=c()
    m.hat.B=c()
    s.A=c()
    s.B=c()
    m.hat.A[10]=mu.hat.A
    m.hat.B[10]=mu.hat.B
    s.A[10]=s.hat.A
    s.B[10]=s.hat.B
    prob.A=c()
    for(i in c(11:n))
    {
      if(u[i]>=c)
      {
        prob.A[i]=max(pnorm((m.hat.A[i-1]-m.hat.B[i-1])/sqrt(s.A[i-1]+s.B[i-1])),(2/pi)*asin(sqrt(s.A[i-1]/(s.A[i-1]+s.B[i-1]))))
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rnorm(1,mean=mu_A,sd=sqrt(s_A))
          count.A=count.A+1
          del[i]=1
        }
        else
        {
          X[i]=rnorm(1,mean=mu_B,sd=sqrt(s_B))
          count.B=count.B+1
          del[i]=0
        }
      }
      else
      {
        prob.A[i]=min(pnorm((m.hat.A[i-1]-m.hat.B[i-1])/sqrt(s.A[i-1]+s.B[i-1])),(2/pi)*asin(sqrt(s.A[i-1]/(s.A[i-1]+s.B[i-1]))))
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rnorm(1,mean=mu_A,sd=sqrt(s_A))
          count.A=count.A+1
          del[i]=1
        }
        else
        {
          X[i]=rnorm(1,mean=mu_B,sd=sqrt(s_B))
          count.B=count.B+1
          del[i]=0
        }
      }
      
      m.hat.A[i]=sum(del*X)/count.A
      s.A[i]=(sum(del*X^2)/count.A)-(m.hat.A[i])^2
      m.hat.B[i]=sum((1-del)*X)/count.B
      s.B[i]=(sum((1-del)*X^2)/count.B)-(m.hat.B[i])^2
    }
    T1[j]=(m.hat.A[n]-m.hat.B[n])
    T2[j]=(s.A[n]-s.B[n])
    t.stat[j]=max(T1[j],T2[j])
  }
  sorted=sort(t.stat)
  critic.val=sorted[5000-(alpha*5000)]
  return(critic.val)
}


power1.norm=function(n,mu_A,mu_B,s_A,s_B,critic.val)
{
  T1=NULL
  T2=NULL
  t.stat=NULL
  for(j in 1:5000)
  {
    count.A=5
    count.B=5
    u=sample(c(1:n),n,replace=T)
    c=45
    del=c(rep(1,5),rep(0,5))
    X=c(rnorm(5,mean=mu_A,sd=sqrt(s_A)),rnorm(5,mean = mu_B,sd= sqrt(s_B)))
    
    #estimate of the mean response and the variation of treatment A
    mu.hat.A=sum(del*X)/count.A
    s.hat.A=(sum(del*X^2)/count.A)-(mu.hat.A)^2
    mu.hat.B=sum((1-del)*X)/count.B
    s.hat.B=(sum((1-del)*X^2)/count.B)-(mu.hat.B)^2
    
    m.hat.A=c()
    m.hat.B=c()
    s.A=c()
    s.B=c()
    m.hat.A[10]=mu.hat.A
    m.hat.B[10]=mu.hat.B
    s.A[10]=s.hat.A
    s.B[10]=s.hat.B
    prob.A=c()
    for(i in c(11:n))
    {
      if(isTRUE(u[i]>=c))
      {
        prob.A[i]=max(pnorm((m.hat.A[i-1]-m.hat.B[i-1])/sqrt(s.A[i-1]+s.B[i-1])),(2/pi)*asin(sqrt(s.A[i-1]/(s.A[i-1]+s.B[i-1]))))
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rnorm(1,mean=mu_A,sd=sqrt(s_A))
          count.A=count.A+1
          del[i]=1
        }
        else
        {
          X[i]=rnorm(1,mean=mu_B,sd=sqrt(s_B))
          count.B=count.B+1
          del[i]=0
        }
      }
      else
      {
        prob.A[i]=min(pnorm((m.hat.A[i-1]-m.hat.B[i-1])/sqrt(s.A[i-1]+s.B[i-1])),(2/pi)*asin(sqrt(s.A[i-1]/(s.A[i-1]+s.B[i-1]))))
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rnorm(1,mean=mu_A,sd=sqrt(s_A))
          count.A=count.A+1
          del[i]=1
        }
        else
        {
          X[i]=rnorm(1,mean=mu_B,sd=sqrt(s_B))
          count.B=count.B+1
          del[i]=0
        }
      }
      
      m.hat.A[i]=sum(del*X)/count.A
      s.A[i]=(sum(del*X^2)/count.A)-(m.hat.A[i])^2
      m.hat.B[i]=sum((1-del)*X)/count.B
      s.B[i]=(sum((1-del)*X^2)/count.B)-(m.hat.B[i])^2
    }
    T1[j]=(m.hat.A[n]-m.hat.B[n])
    T2[j]=(s.A[n]-s.B[n])
    t.stat[j]=max(T1[j],T2[j])
  }
  sorted=sort(t.stat)
  power=mean(sorted>critic.val)
  return(power)
}


critic4.norm=function(n,mu_A,mu_B,s_A,s_B,alpha)
{
  Q=NULL
  for(j in 1:5000)
  {
    count.A=5
    count.B=5
    u=sample(c(1:n),n,replace=T)
    c=45
    del=c(rep(1,5),rep(0,5))
    X=c(rnorm(5,mean=mu_A,sd=sqrt(s_A)),rnorm(5,mean = mu_B,sd= sqrt(s_B)))
    
    #estimate of the mean response and the variation of treatment A
    mu.hat.A=sum(del*X)/count.A
    s.hat.A=(sum(del*X^2)/count.A)-(mu.hat.A)^2
    mu.hat.B=sum((1-del)*X)/count.B
    s.hat.B=(sum((1-del)*X^2)/count.B)-(mu.hat.B)^2
    
    m.hat.A=c()
    m.hat.B=c()
    s.A=c()
    s.B=c()
    m.hat.A[10]=mu.hat.A
    m.hat.B[10]=mu.hat.B
    s.A[10]=s.hat.A
    s.B[10]=s.hat.B
    prob.A=c()
    for(i in c(11:n))
    {
      if(u[i]>=c)
      {
        prob.A[i]=max(pnorm((m.hat.A[i-1]-m.hat.B[i-1])/sqrt(s.A[i-1]+s.B[i-1])),(2/pi)*asin(sqrt(s.A[i-1]/(s.A[i-1]+s.B[i-1]))))
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rnorm(1,mean=mu_A,sd=sqrt(s_A))
          count.A=count.A+1
          del[i]=1
        }
        else
        {
          X[i]=rnorm(1,mean=mu_B,sd=sqrt(s_B))
          count.B=count.B+1
          del[i]=0
        }
      }
      else
      {
        prob.A[i]=min(pnorm((m.hat.A[i-1]-m.hat.B[i-1])/sqrt(s.A[i-1]+s.B[i-1])),(2/pi)*asin(sqrt(s.A[i-1]/(s.A[i-1]+s.B[i-1]))))
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rnorm(1,mean=mu_A,sd=sqrt(s_A))
          count.A=count.A+1
          del[i]=1
        }
        else
        {
          X[i]=rnorm(1,mean=mu_B,sd=sqrt(s_B))
          count.B=count.B+1
          del[i]=0
        }
      }
      
      m.hat.A[i]=sum(del*X)/count.A
      s.A[i]=(sum(del*X^2)/count.A)-(m.hat.A[i])^2
      m.hat.B[i]=sum((1-del)*X)/count.B
      s.B[i]=(sum((1-del)*X^2)/count.B)-(m.hat.B[i])^2
    }
    D1=(m.hat.A[n]-m.hat.B[n])/sd(m.hat.A[11:n]-m.hat.B[11:n])
    D2=(s.A[n]-s.B[n])/sd(s.A[11:n]-s.B[11:n])
    r=cor(m.hat.A[11:n]-m.hat.B[11:n],s.A[11:n]-s.B[11:n])
    
    if(D1>0 & D2>0){
      Q[j]=sqrt(D1^2+D2^2-2*r*D1*D2)/sqrt(1-r^2)
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


power4.norm=function(n,mu_A,mu_B,s_A,s_B,critic.val)
{
  Q=NULL
  for(j in 1:5000)
  {
    count.A=5
    count.B=5
    u=sample(c(1:n),n,replace=T)
    c=45
    del=c(rep(1,5),rep(0,5))
    X=c(rnorm(5,mean=mu_A,sd=sqrt(s_A)),rnorm(5,mean = mu_B,sd= sqrt(s_B)))
    
    #estimate of the mean response and the variation of treatment A
    mu.hat.A=sum(del*X)/count.A
    s.hat.A=(sum(del*X^2)/count.A)-(mu.hat.A)^2
    mu.hat.B=sum((1-del)*X)/count.B
    s.hat.B=(sum((1-del)*X^2)/count.B)-(mu.hat.B)^2
    
    m.hat.A=c()
    m.hat.B=c()
    s.A=c()
    s.B=c()
    m.hat.A[10]=mu.hat.A
    m.hat.B[10]=mu.hat.B
    s.A[10]=s.hat.A
    s.B[10]=s.hat.B
    prob.A=c()
    for(i in c(11:n))
    {
      if(u[i]>=c)
      {
        prob.A[i]=max(pnorm((m.hat.A[i-1]-m.hat.B[i-1])/sqrt(s.A[i-1]+s.B[i-1])),(2/pi)*asin(sqrt(s.A[i-1]/(s.A[i-1]+s.B[i-1]))))
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rnorm(1,mean=mu_A,sd=sqrt(s_A))
          count.A=count.A+1
          del[i]=1
        }
        else
        {
          X[i]=rnorm(1,mean=mu_B,sd=sqrt(s_B))
          count.B=count.B+1
          del[i]=0
        }
      }
      else
      {
        prob.A[i]=min(pnorm((m.hat.A[i-1]-m.hat.B[i-1])/sqrt(s.A[i-1]+s.B[i-1])),(2/pi)*asin(sqrt(s.A[i-1]/(s.A[i-1]+s.B[i-1]))))
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rnorm(1,mean=mu_A,sd=sqrt(s_A))
          count.A=count.A+1
          del[i]=1
        }
        else
        {
          X[i]=rnorm(1,mean=mu_B,sd=sqrt(s_B))
          count.B=count.B+1
          del[i]=0
        }
      }
      
      m.hat.A[i]=sum(del*X)/count.A
      s.A[i]=(sum(del*X^2)/count.A)-(m.hat.A[i])^2
      m.hat.B[i]=sum((1-del)*X)/count.B
      s.B[i]=(sum((1-del)*X^2)/count.B)-(m.hat.B[i])^2
    }
    D1=(m.hat.A[n]-m.hat.B[n])/sd(m.hat.A[11:n]-m.hat.B[11:n])
    D2=(s.A[n]-s.B[n])/sd(s.A[11:n]-s.B[11:n])
    r=cor(m.hat.A[11:n]-m.hat.B[11:n],s.A[11:n]-s.B[11:n])
    
    if(D1>0 & D2>0){
      Q[j]=sqrt(D1^2+D2^2-2*r*D1*D2)/sqrt(1-r^2)
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

critic3.norm=function(n,mu_A,mu_B,s_A,s_B,alpha)
{
  T1=NULL
  T2=NULL
  t.stat=NULL
  for(j in 1:5000)
  {
    count.A=5
    count.B=5
    u=sample(c(1:n),n,replace=T)
    c=45
    del=c(rep(1,5),rep(0,5))
    X=c(rnorm(5,mean=mu_A,sd=sqrt(s_A)),rnorm(5,mean = mu_B,sd= sqrt(s_B)))
    
    #estimate of the mean response and the variation of treatment A
    mu.hat.A=sum(del*X)/count.A
    s.hat.A=(sum(del*X^2)/count.A)-(mu.hat.A)^2
    mu.hat.B=sum((1-del)*X)/count.B
    s.hat.B=(sum((1-del)*X^2)/count.B)-(mu.hat.B)^2
    
    m.hat.A=c()
    m.hat.B=c()
    s.A=c()
    s.B=c()
    m.hat.A[10]=mu.hat.A
    m.hat.B[10]=mu.hat.B
    s.A[10]=s.hat.A
    s.B[10]=s.hat.B
    prob.A=c()
    for(i in c(11:n))
    {
      if(u[i]>=c)
      {
        prob.A[i]=max(pnorm((m.hat.A[i-1]-m.hat.B[i-1])/sqrt(s.A[i-1]+s.B[i-1])),(2/pi)*asin(sqrt(s.A[i-1]/(s.A[i-1]+s.B[i-1]))))
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rnorm(1,mean=mu_A,sd=sqrt(s_A))
          count.A=count.A+1
          del[i]=1
        }
        else
        {
          X[i]=rnorm(1,mean=mu_B,sd=sqrt(s_B))
          count.B=count.B+1
          del[i]=0
        }
      }
      else
      {
        prob.A[i]=min(pnorm((m.hat.A[i-1]-m.hat.B[i-1])/sqrt(s.A[i-1]+s.B[i-1])),(2/pi)*asin(sqrt(s.A[i-1]/(s.A[i-1]+s.B[i-1]))))
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rnorm(1,mean=mu_A,sd=sqrt(s_A))
          count.A=count.A+1
          del[i]=1
        }
        else
        {
          X[i]=rnorm(1,mean=mu_B,sd=sqrt(s_B))
          count.B=count.B+1
          del[i]=0
        }
      }
      
      m.hat.A[i]=sum(del*X)/count.A
      s.A[i]=(sum(del*X^2)/count.A)-(m.hat.A[i])^2
      m.hat.B[i]=sum((1-del)*X)/count.B
      s.B[i]=(sum((1-del)*X^2)/count.B)-(m.hat.B[i])^2
    }
    T1[j]=m.hat.A[n]-m.hat.B[n]
    T2[j]=s.A[n]-s.B[n]
    t.stat[j]=(T1[j]+T2[j])/2
  }
  sorted=sort(t.stat)
  critic.val=sorted[5000-(alpha*5000)]
  return(critic.val)
}

power3.norm=function(n,mu_A,mu_B,s_A,s_B,critic.val)
{
  T1=NULL
  T2=NULL
  t.stat=NULL
  for(j in 1:5000)
  {
    count.A=5
    count.B=5
    u=sample(c(1:n),n,replace=T)
    c=45
    del=c(rep(1,5),rep(0,5))
    X=c(rnorm(5,mean=mu_A,sd=sqrt(s_A)),rnorm(5,mean = mu_B,sd= sqrt(s_B)))
    
    #estimate of the mean response and the variation of treatment A
    mu.hat.A=sum(del*X)/count.A
    s.hat.A=(sum(del*X^2)/count.A)-(mu.hat.A)^2
    mu.hat.B=sum((1-del)*X)/count.B
    s.hat.B=(sum((1-del)*X^2)/count.B)-(mu.hat.B)^2
    
    m.hat.A=c()
    m.hat.B=c()
    s.A=c()
    s.B=c()
    m.hat.A[10]=mu.hat.A
    m.hat.B[10]=mu.hat.B
    s.A[10]=s.hat.A
    s.B[10]=s.hat.B
    prob.A=c()
    for(i in c(11:n))
    {
      if(u[i]>=c)
      {
        prob.A[i]=max(pnorm((m.hat.A[i-1]-m.hat.B[i-1])/sqrt(s.A[i-1]+s.B[i-1])),(2/pi)*asin(sqrt(s.A[i-1]/(s.A[i-1]+s.B[i-1]))))
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rnorm(1,mean=mu_A,sd=sqrt(s_A))
          count.A=count.A+1
          del[i]=1
        }
        else
        {
          X[i]=rnorm(1,mean=mu_B,sd=sqrt(s_B))
          count.B=count.B+1
          del[i]=0
        }
      }
      else
      {
        prob.A[i]=min(pnorm((m.hat.A[i-1]-m.hat.B[i-1])/sqrt(s.A[i-1]+s.B[i-1])),(2/pi)*asin(sqrt(s.A[i-1]/(s.A[i-1]+s.B[i-1]))))
        if(isTRUE(runif(1)<prob.A[i]))
        {
          X[i]=rnorm(1,mean=mu_A,sd=sqrt(s_A))
          count.A=count.A+1
          del[i]=1
        }
        else
        {
          X[i]=rnorm(1,mean=mu_B,sd=sqrt(s_B))
          count.B=count.B+1
          del[i]=0
        }
      }
      
      m.hat.A[i]=sum(del*X)/count.A
      s.A[i]=(sum(del*X^2)/count.A)-(m.hat.A[i])^2
      m.hat.B[i]=sum((1-del)*X)/count.B
      s.B[i]=(sum((1-del)*X^2)/count.B)-(m.hat.B[i])^2
    }
    T1[j]=m.hat.A[n]-m.hat.B[n]
    T2[j]=s.A[n]-s.B[n]
    t.stat[j]=(T1[j]+T2[j])/2
  }
  sorted=sort(t.stat)
  power=mean(sorted>critic.val)
  return(power)
}
#power-curve

#max statistic
set.seed(100)
cr.val=critic1.norm(100,5,5,2,2,0.05)

mA=seq(5,8,by=0.1)
sA=seq(2,9,by=0.1)
po=NULL
for(i in c(1:length(sA)))
{
  set.seed(100)
  po[i]=power1.norm(100,5,5,sA[i],2,cr.val)
}
plot(sA-2,po,type="l",lwd=2,col="red",xlab="difference in sigma",ylab="Power",main="Power Curve")

#Qn statistic
set.seed(100)
cr.val.Q=critic4.norm(100,5,5,2,2,0.05)
po2=NULL
for(i in c(1:length(sA)))
{
  set.seed(100)
  po2[i]=power4.norm(100,5,5,sA[i],2,cr.val.Q)
}
lines(sA-2,po2,col="blue",lwd=2)


#Average Statistic
set.seed(100)
cr.val.A=critic3.norm(100,5,5,2,2,0.05)
po3=NULL
for(i in c(1:length(sA)))
{
  set.seed(100)
  po3[i]=power3.norm(100,5,5,sA[i],2,cr.val.A)
}
lines(sA-2,po3,col="orange",lwd=2)
legend("bottomright",legend=c("Max stat","Q stat","Average stat"),col=c("red","blue","orange"),lty=1:1,lwd=2)


plot(sA-2,po,type="l",lwd=2,lty=1,xlab="difference in sigma",ylab="Power",main="Power Curve")
lines(sA-2,po2,lty=2,lwd=2)
lines(sA-2,po3,lty=15,lwd=2)
legend("bottomright",legend=c("Max stat","Q stat","Average stat"),lty=c(1,2,15),lwd=2)




