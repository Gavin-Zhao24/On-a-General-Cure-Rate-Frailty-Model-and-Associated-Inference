library(stats4)
library(SuppDists)
cut_1=cut_data[cut_data$x==1,]
cut_2=cut_data[cut_data$x==2,]
cut_3=cut_data[cut_data$x==3,]
cut_4=cut_data[cut_data$x==4,]

h0<-function(t,k,lambda) k/lambda*(t/lambda)^(k-1)*exp(-(t/lambda)^k)

E<-function(k,lambda,beta){
  sum0<-0
  for (i in 1:nrow(cut_1)){
    sum0=sum0+(1-exp(-(cut_1[i,1]/lambda)^k))*exp(beta*cut_1[i,3])
  }
  sum0
}

###when d=1
llh1<-function(beta1,k,lambda){
  sum1<-0
  for (i in 1:nrow(cut_1)){
    sum1<-sum1+beta1*cut_1[i,3]+log(h0(cut_1[i,1],k,lambda))
  }
  sum1
}

###when d=0
llh2<-function(ind,p=0.56){
  sum2<-0
  for (i in 1:nrow(cut_0)){
    sum2<-sum2+(1-ind)*log(p)+ind*log(1-p)
  }
  sum2
}


llh3<-function(ind,m=nrow(cut_1),p=0.56){
  m*(1-ind)*log(p)+m*ind*log(1-p)
}


llh4<-function(gam,lambda,k,beta1,m=nrow(cut_1)){
  gam*(1-sqrt(1+2*E(k,lambda,beta1)/gam))+m*log(pi*gam/2)
}


llh<-function(beta1,k,lambda,gam){
  lik<-llh1(beta1,k,lambda)+llh2(1)+llh2(0)+llh3(1)+llh3(0)+llh4(gam,lambda,k,beta1)
  lik
}
llh(0.13,1.75,10,1)


parameter<-function(n){
  k<-0
  ma<--Inf
  count<-0
  while (count<=n) {
    a<-runif(1,min=0,max=0.5)
    b<-runif(1,min=0,max=3)
    c<-runif(1,min=0,max=15)
    d<-runif(1,min=0,max=50)
    k<-llh(a,b,c,d)
    if (is.na(k>=ma)){
      next
    } else if (k>=ma){
      ma<-k
      coe<-c(a,b,c,d)
    }
    count<-count+1
  }
  return(c(coe,ma))
}  
a<-parameter(200)

iteration<-function(n){
  beta<-c(1:n)
  k<-c(1:n)
  lambda<-c(1:n)
  gam<-c(1:n)
  maxi<-c(1:n)
  for (i in 1:n){
    onetime<-parameter(100)
    beta[i]<-onetime[1]
    k[i]<-onetime[2]
    lambda[i]<-onetime[3]
    gam[i]<-onetime[4]
    maxi[i]<-onetime[5]
  }
  result<-data.frame(beta,k,lambda,gam,maxi)
  return(result)
}

estimation<-iteration(1000)
estimation1<-iteration(1000)
estimation2<-iteration(1000)
estimation3<-iteration(1000)

mean(estimation2$beta)
median(estimation2$beta)
sd(estimation2$beta)

mean(estimation2$k)
median(estimation2$k)
sd(estimation2$k)

mean(estimation2$lambda)
median(estimation2$lambda)
sd(estimation2$lambda)

mean(estimation2$gam)
median(estimation2$gam)
sd(estimation2$gam)

plot(estimation1$beta,ylab="Estimation of Beta")
plot(estimation1$k,ylab="Estimation of k")
plot(estimation1$lambda,ylab="Estimation of lambda")
plot(estimation1$gam,ylab="Estimation of gam")

plot(estimation2$beta)
plot(estimation2$k)
plot(estimation2$lambda)
plot(estimation2$gam)

y<-rinvGauss(1,nu = 1,lambda = 0.1)
y<-rinvGauss(1,nu = 1,lambda = 1)
surv<-function(p,t,d,x,y,lambda,k,beta){
  prob<-p+(1-p)*exp(-y*(1-exp(-(t/lambda)^k))*exp(beta*x))
  return(prob)
}
y<-rinvGauss(1,1,50)
plot(cut_1$t,surv(0.676,cut_1$t,cut_1$d,cut_1$x,y,2.9,1.76,0.5),xlab="Survival Time",
     ylab="Proportion Surviving",type = "p")
plot(cut_2$t,surv(0.613,cut_2$t,cut_2$d,cut_2$x,y,2.9,1.76,0.5),xlab="Survival Time",
     ylab="Proportion Surviving",type = "p")
plot(cut_3$t,surv(0.529,cut_3$t,cut_3$d,cut_3$x,y,2.9,1.76,0.5),xlab="Survival Time",
     ylab="Proportion Surviving",type = "p")
plot(cut_4$t,surv(0.329,cut_4$t,cut_4$d,cut_4$x,y,2.9,1.76,0.5),xlab="Survival Time",
     ylab="Proportion Surviving",type = "p")






###it is link to the covariate between 0 and 1,cure rate=1/(1+exp(beta1*x)) 0.65 0.25


beta_esti<-function(k,lambda){
  beta<-0.01
  gam<-0.1
  be_vec<-c(0:50,0.01)
  ga_vec<-c(1:500)
  ma<-c(1:500)
  for (i in 1:50){
    for (j in 1:50){
      x<-llh(beta,k,lambda,gam)
      ma[i]<-x
      be_vec[i]<-beta
      beta<-beta+0.01
    }
    
    ma[i]<-x
    be_vec[i]<-beta
    beta<-beta+0.01
  }
  return(data.frame(be_vec,ma))
}  
s=beta_esti(1.76,2.90,43.255)
llh(0.4,1.76,2.9,43.255)

be_vec<-c(1:50)*0.01
ga_vec<-c(1:500)*0.1
ma<-c(1:25000)
m<-c(1:25000)
n<-c(1:25000)
coun<-1
for (i in be_vec){
  for (j in ga_vec){
    s<-llh(i,1.76,2.9,j)
    ma[coun]<-s
    m[coun]<-i
    n[coun]<-j
    coun<-coun+1
  }
}
kkk<-data.frame(ma,m,n)