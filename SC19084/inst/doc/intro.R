## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval=TRUE----------------------------------------------------------------
library(SC19084)
set.seed(123)
n<-10;p<-15
b0<-matrix(0,p,1)
b0[1:7]<-c(1,-0.5,0.7,-1.2,-0.9,0.3,0.55)
epsi <-0.01
cormat=matrix(1,p,p)
corval = 0.5
for (i in 1:(p-1)){
  for (j in (i+1):p){
    cormat[i,j]=corval^abs(i-j)
    cormat[j,i]=cormat[i,j]
 }
}
sigmat=cormat^0.5
X<-matrix(rnorm(n*p),n,p)%*%sigmat
y<-X%*%b0+0.4*matrix(rnorm(n),n,1)
Xsca <- sqrt(colSums(X^2))/sqrt(n)
X <- X/(matrix(1,n,1)%*%Xsca)

## -----------------------------------------------------------------------------
XXmat <- n^(-1)*t(X)%*%X
cvec <- n^(-1)*t(X)%*%y
lammax = max(abs(cvec))
lammin = epsi*lammax
lamvec = lammax - (lammax - lammin)/(50 - 1)*c(0:(50 - 1))
inival<-rep(0,p)
lambda<-lamvec[50]
varset<-vector(mode="numeric")

## -----------------------------------------------------------------------------
ol1(XXmat, cvec, inival, lambda, varset, 50, 1e-4)

## -----------------------------------------------------------------------------
set.seed(321)
n<-10;p<-15
b0<-matrix(0,p,1)
b0[1:7]<-c(1,-0.5,0.7,-1.2,-0.9,0.3,0.55)
epsi <-0.01
cormat=matrix(1,p,p)
corval = 0.5
for (i in 1:(p-1)){
  for (j in (i+1):p){
    cormat[i,j]=corval^abs(i-j)
    cormat[j,i]=cormat[i,j]
 }
}
sigmat=cormat^0.5
X<-matrix(rnorm(n*p),n,p)%*%sigmat
y<-X%*%b0+0.4*matrix(rnorm(n),n,1)
Xsca <- sqrt(colSums(X^2))/sqrt(n)
X <- X/(matrix(1,n,1)%*%Xsca)

## -----------------------------------------------------------------------------
XXmat <- n^(-1)*t(X)%*%X
z1 = XXmat[1, 1]
Lam<-1/z1
cvec <- n^(-1)*t(X)%*%y
z<-cvec[1]/z1
lammax = max(abs(cvec))
lammin = epsi*lammax
lamvec = lammax - (lammax - lammin)/(50 - 1)*c(0:(50 - 1))
lambda0<-lamvec[49]
lambda<-lamvec[50]

## -----------------------------------------------------------------------------
uhard(z, Lam, lambda0, lambda)

## -----------------------------------------------------------------------------
#par(bg="ivory3")
set.seed(123)
x=rnorm(100,mean=100,sd=10)
hist(x,breaks=20,col=2,main="Simple Histogram")

## -----------------------------------------------------------------------------
A=head(freeny)
M=A[-2]
knitr::kable(M)

## -----------------------------------------------------------------------------
x=rbinom(50,3000,0.001)
y=rpois(50,3)
lm.cf=lm(x~y)
summary(lm.cf)$coef

## -----------------------------------------------------------------------------
n<-1e4
myrn<-function(sigma){
u<-runif(n)
x<-(-2*(sigma^2)*log(1-u))^(1/2)
hist(x,breaks=15, prob = TRUE,#main=expression(f(x)==x/sigma^2*exp(-(x^2)/(2*sigma^2))))
      main=paste("sigma=",sigma))
a <- seq(.001, 10, .001);y<-a/(sigma^2)*exp(-(a^2)/(2*sigma^2))
lines(a, y)
}

## -----------------------------------------------------------------------------
#par(mfrow=c(2,2))
myrn(.5);myrn(1);myrn(2);myrn(4)

## -----------------------------------------------------------------------------
n <- 1e3
x1 <- rnorm(n,0,1)
x2 <- rnorm(n,3,1)
myz<-function(r){  # function to generate different z as value of p1 changes
z <- r*x1+(1-r)*x2
hist(z,breaks=30,main=paste("z=",r,"X1+",(1-r),"X2"))
}
#par(mfrow=c(2,2))
myz(.75);myz(.5);myz(.25);myz(.1)#show histograms where p1=0,75,0.5,0.25,0.1

## -----------------------------------------------------------------------------
myW=function(n){
  sigma<-matrix(c(9,18,3,18,40,16,2,16,27),3,3)#set a three dimentional scale matrix sigma
  a<-rep(0,9);T<-matrix(a,3,3)
  T[1,1]<-sqrt(rchisq(1,n));T[2,1]<-rnorm(1);T[2,2]<-sqrt(rchisq(1,n-1))
  T[3,1]<-rnorm(1);T[3,2]<-rnorm(1);T[3,3]<-sqrt(rchisq(1,n-2))#set matrix T
  A<-T%*%t(T) #obtain matrix A
  L=matrix(c(3,6,1,0,2,5,0,0,1),3,3)#show the Choleski factorization of sigma
  X<-L%*%A%*%t(L) # obtain random sample X
  X
}

## -----------------------------------------------------------------------------
myW(30)

## -----------------------------------------------------------------------------
set.seed(12345)   
# for reproducible research 
m <-1e5
x <-runif(m,0,pi/3)
# X with uniform distribution
g <-(pi/3)*sin(x)
est1 <-mean(g)
# Monte Carlo estimate of the integral
print(c(est1,cos(0)-cos(pi/3)))
# compare with the exact value

## -----------------------------------------------------------------------------
d <-sd(g)/sqrt(m)
# standard error of Monte Carlo estimate
c(est1 - 1.96 * d, est1 + 1.96 * d)
# 95% confidence intervals

## -----------------------------------------------------------------------------
MC <-function(n){
  x <-runif(n)
  g <-exp(-x)/(1+x^2)
  est101 <-mean(g)
  est101
}
# function of Monte Carlo estimate without antithetic variables
MC.anti<-function(n){
  u <-runif(n/2)
  v <-1-u
  g1 <-exp(-u)/(1+u^2)
  g2 <-exp(-v)/(1+v^2)
# antithetic variables
  est102 <-mean(g1+g2)/2
  est102
}
# function of Monte Carlo estimate with antithetic variables
m <-1000
set.seed(321)
MC1 <-MC2 <-numeric(m)
for (i in 1:m){
  MC1[i] <-MC(m)
  MC2[i] <-MC.anti(m/2)
}

## -----------------------------------------------------------------------------
print((var(MC1) - var(MC2))/var(MC1))

## -----------------------------------------------------------------------------
cdf <-function(p){
  uniroot(function(x) exp(1)*(1-exp(-x))-p*(exp(1)-1),c(0,1),tol=1e-9)$root
}
# funtion to compute quantiles
x <-numeric(6)
p <-rep(.2,5)
cp <-cumsum(p)
for(i in 1:5){
  x[i+1]<-cdf(cp[i])
}
# five subintervals
m <-1e4
# number of replicates
k <-5
# number of strata
r <- m/k
# replicates per stratum
n <- 50 
# number of times to repeat the estimation
T <- numeric(k)
estimates <- matrix(0,n,2)
g <- function(x){
  exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
}
F <- function(x){
  exp(-x)*5/(exp(-1)-1) * (x > 0) * (x < 1)
}
for(i in 1:n) {
  v<- runif(m)
  # inverse transform method
  s <- - log(1 - v * (1 - exp(-1)))
  # random number generated from f3(x) in example 5.10
  estimates[i, 1] <- mean(g(s) / (exp(-s) / (1 - exp(-1))))
  # MC estimator with importance sampling method
for (j in 1:k){
  u<-runif(r)
  # inverse transform method
  #t<-uniroot(function(y) F(y)-F(x[j])-u,c(0,1),tol=1e-9)$root
  t<- -log(exp(-x[j])-u*(1-exp(-1))/5)
  # random number generated from f(x)
  gi<-g(t)
  fi<-exp(-t)*5/(1-exp(-1))
  Y <- gi/fi
  T[j]<-mean(Y)
  # MC estimator with importance sampling method in jth interval
}
  estimates[i, 2] <- sum(T)
  # stratified importance sampling estimate
}

## -----------------------------------------------------------------------------
 apply(estimates, 2, mean)

## -----------------------------------------------------------------------------
 apply(estimates, 2, var)

## -----------------------------------------------------------------------------
xt <-numeric(11)
pt <-rep(.1,10)
cpt <-cumsum(pt)
for(i in 1:10){
  xt[i+1]<-cdf(cpt[i])
}
kt <-10
# number of strata
rt <- m/kt
# replicates per stratum
F <- numeric(kt)
test <- matrix(0,n,1)
for(i in 1:n) {
  for (j in 1:kt){
    ut<-runif(rt)
    tt<- -log(exp(-xt[j])-ut*(1-exp(-1))/10)
    git<-g(tt)
    fit<-exp(-tt)*10/(1-exp(-1))
    Yt <- git/fit
    F[j]<-mean(Yt)
  }
  # MC estimator with importance sampling method in jth interval
  test[i] <- sum(F)
  # stratified importance sampling estimate
}

## -----------------------------------------------------------------------------
 mean(test)

## -----------------------------------------------------------------------------
 var(test)

## -----------------------------------------------------------------------------
rm(list=ls())
# remove variables
set.seed(123)
# for repruducible research
n<-20;alpha<-.05
CI_1<-function(a,b){
  x<-rchisq(a,2)
  UCL<-mean(x)+sqrt(var(x)/a)*qt(1-b,df=a-1)
  return(UCL)
}
# function to obtain confidence interval appling t-interval
CI_2<-function(c,d){
  y <- rnorm(c, mean=0, sd=2)
  UCL<-(c-1) * var(y) / qchisq(d, df=c-1)
  return(UCL)
}
# function to obtain confidence interval in example 6.4

## ----results='asis'-----------------------------------------------------------
m<-1e4
# number of replication
ci_1<-replicate(m,expr=CI_1(n,alpha))
# upper confidence limit appling t-interval
ci_2<-replicate(m,expr=CI_2(n,alpha))
# upper confidence limit in example 6.4
c<-matrix(c(mean(ci_1>2),mean(ci_2>4)),1,2)
#Monte Carlo estimates of confidence level
colnames(c)<-c("ECL_1","ECL_2")
pander::pandoc.table(c)

## ----results='asis'-----------------------------------------------------------
set.seed(002)
n_2<-40;n_3<-80;n_4<-160
# a series of sample size
ci_11<-replicate(m,expr=CI_1(n,alpha))
ci_12<-replicate(m,expr=CI_1(n_2,alpha))
ci_13<-replicate(m,expr=CI_1(n_3,alpha))
ci_14<-replicate(m,expr=CI_1(n_4,alpha))
d<-matrix(c(mean(ci_11>2),mean(ci_12>2),mean(ci_13>2),mean(ci_14>2)),1,4)
#Monte Carlo estimates of confidence level
colnames(d)<-c("n=20","n=40","n=80","n=160")
pander::pandoc.table(d)

## ----results='asis'-----------------------------------------------------------
rm(list=ls())
# remove variables
set.seed(0321)
n<-1e3
# sample size
m<-1e4
# number of replication
sk <- function(x) {
#computes the sample skewness.
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
q1<-quantile(replicate(m,expr=sk(rnorm(n))),c(.025,.05,.95,.975))
# quantiles of the skewness under normality by a Monte Carlo experiment
q2<-qnorm(c(.025,.05,.95,.975),mean=0,sd=sqrt(6/n))
# quantiles of the large sample approximation
x<-matrix(c(q1,q2),2,4,byrow=T)
rownames(x)<-c("eatimate","approximation")
colnames(x)<-c("0.025","0.05","0.95","0.975")
pander::pandoc.table(x)

## ----results='asis'-----------------------------------------------------------
set.seed(007)
sd<-numeric(4)
c<-c(.025,.05,.95,.975)
for (i in 1:4){
  sd[i]<-sqrt(c[i]*(1-c[i])/(n*(dnorm(x[1,i],mean=0,sd=sqrt(6*(n-2)/((n+1)*(n+3))))^2)))
}
# standard error of the estimates
y<-matrix(sd,1,4)
rownames(y)<-"standard error"
colnames(y)<-c("0.025","0.05","0.95","0.975")
pander::pandoc.table(y)

## ----results='asis'-----------------------------------------------------------
n_2<-1e2
# decreased sample size
sd_2<-numeric(4)
for (i in 1:4){
  sd_2[i]<-sqrt(c[i]*(1-c[i])/(n_2*(dnorm(x[1,i],mean=0,sd=sqrt(6*(n_2-2)/((n_2+1)*(n_2+3))))^2)))
}
z<-matrix(sd_2,1,4)
rownames(z)<-"standard error 2"
colnames(z)<-c("0.025","0.05","0.95","0.975")
pander::pandoc.table(z)

## -----------------------------------------------------------------------------
rm(list=ls())
# remove variables
#library(ggplot2)
# to display scatter plots
library(bootstrap)
# to obtain data
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * (1 + r) / 2)
}
# put (absolute) correlations on the upper panels with size proportional to the correlations
pairs(scor,upper.panel = panel.cor,gap=0, row1attop=FALSE)

## ----results='asis'-----------------------------------------------------------
B <- 1e3
set.seed(1234)
r12 <- function(x, i) {
#want correlation of columns 1 and 2
cor(x[i,1], x[i,2])
}
r34 <- function(x, i) {
#want correlation of columns 3 and 4
cor(x[i,3], x[i,4])
}
r35 <- function(x, i) {
#want correlation of columns 3 and 5
cor(x[i,3], x[i,5])
}
r45 <- function(x, i) {
#want correlation of columns 4 and 5
cor(x[i,4], x[i,5])
}
library(boot) 
# for boot function
sd_boot<-numeric(4)
#for (i in 1:4){
sd_boot[1]<-sd(boot(data=scor,statistic=r12,R=2000)$t)
sd_boot[2]<-sd(boot(data=scor,statistic=r34,R=2000)$t)
sd_boot[3]<-sd(boot(data=scor,statistic=r35,R=2000)$t)
sd_boot[4]<-sd(boot(data=scor,statistic=r45,R=2000)$t)
#}
a<-matrix(sd_boot,1,4)
colnames(a)<-c("Cor(1,2)","Cor(3,4)","Cor(3,5)","Cor(4,5)")
rownames(a)<-"sd.boot"
pander::pandoc.table(a)

## ----results='asis'-----------------------------------------------------------
rm(list=ls())
# remove variables
library(boot)
set.seed(0321)
sk_nor<-0;sk_chi<-2/sqrt(5/2)
# skewnesses with normal distribution and chi-square distribution
n <- 20
# sample size for statistic
m<-10
# sample size for confidence interval
M<-1e2
# number of replicate for empirical coverage rates
sk.boot <- function(x,i) {
  # function to compute the sample skewness statistic
  xbar <- mean(x[i,])
  m3 <- mean((x[i,] - xbar)^3)
  m2 <- mean((x[i,] - xbar)^2)
  return( m3 / m2^1.5 )
}
ci.norm_nor<-ci.norm_chi<-ci.basic_nor<-ci.basic_chi<-ci.perc_nor<-ci.perc_chi<-matrix(NA,M,2)
for(i in 1:M){
  x<-matrix(rnorm(m*n),m,n)
  # sample from normal population
  y<-matrix(rchisq(m*n,5),m,n)
  # sample from chi-square population
  de_nor <- boot(x, statistic = sk.boot, R = 1000)
  ci_nor <- boot.ci(de_nor,type=c("norm","basic","perc"))
  # bootstrap CI for normality
  de_chi <- boot(y, statistic = sk.boot, R = 1000)
  ci_chi <- boot.ci(de_chi,type=c("norm","basic","perc"))
  # bootstrap CI for chi-square
  ci.norm_nor[i,]<-ci_nor$norm[2:3]
  ci.basic_nor[i,]<-ci_nor$basic[4:5]
  ci.perc_nor[i,]<-ci_nor$percent[4:5]
  # extract three sorts of CI for normality
  ci.norm_chi[i,]<-ci_chi$norm[2:3]
  ci.basic_chi[i,]<-ci_chi$basic[4:5]
  ci.perc_chi[i,]<-ci_chi$percent[4:5]
  # extract three sorts of CI for chi-square
}
CR_nor <- c(mean(ci.norm_nor[,1]<=sk_nor & ci.norm_nor[,2]>=sk_nor),mean(ci.basic_nor[,1]<=sk_nor & ci.basic_nor[,2]>=sk_nor),mean(ci.perc_nor[,1]<=sk_nor & ci.perc_nor[,2]>=sk_nor),mean(ci.norm_nor[,1]>=sk_nor ),mean(ci.basic_nor[,1]>=sk_nor),mean(ci.perc_nor[,1]>=sk_nor),mean(ci.norm_nor[,2]<=sk_nor ),mean(ci.basic_nor[,2]<=sk_nor),mean(ci.perc_nor[,2]<=sk_nor))
# Monte Carlo estimates of CI and miss rates of both sides for normality
b<-matrix(CR_nor,3,3,byrow=T)
colnames(b)<-c("norm","basic","perc")
rownames(b)<-c("CR for normality","miss rates on the left","miss rates on the right")

CR_chi <- c(mean(ci.norm_chi[,1]<=sk_chi & ci.norm_chi[,2]>=sk_chi),mean(ci.basic_chi[,1]<=sk_chi & ci.basic_chi[,2]>=sk_chi),mean(ci.perc_chi[,1]<=sk_chi & ci.perc_chi[,2]>=sk_chi),mean(ci.norm_chi[,1]>=sk_chi ),mean(ci.basic_chi[,1]>=sk_chi),mean(ci.perc_chi[,1]>=sk_chi),mean(ci.norm_chi[,2]<=sk_chi ),mean(ci.basic_chi[,2]<=sk_chi),mean(ci.perc_chi[,2]<=sk_chi))
#Monte Carlo estimates of CI and miss rates of both sides for chi-square
c<-matrix(CR_chi,3,3,byrow=T)
colnames(c)<-c("norm","basic","perc")
rownames(c)<-c("CR for chi-square","miss rates on the left","miss rates on the right")

pander::pandoc.table(b)
pander::pandoc.table(c)

## ----results='asis'-----------------------------------------------------------
set.seed(123)
sk_exp<-2
# skewnesses with exponential distribution
ci.norm_exp<-ci.perc_exp<-ci.basic_exp<-matrix(NA,M,2)
for(i in 1:M){
  x<-matrix(rexp(m*n,2),m,n)
  # sample from exp(2) population
  de_exp <- boot(x, statistic = sk.boot, R = 1000)
  ci_exp <- boot.ci(de_exp,type=c("norm","basic","perc"))
  # bootstrap CI for exponential distribution
  ci.norm_exp[i,]<-ci_exp$norm[2:3]
  ci.basic_exp[i,]<-ci_exp$basic[4:5]
  ci.perc_exp[i,]<-ci_exp$percent[4:5]
  # extract three sorts of CI for exponential distribution
}
  CR_exp <- c(mean(ci.norm_exp[,1]<=sk_exp & ci.norm_exp[,2]>=sk_exp),mean(ci.basic_exp[,1]<=sk_exp & ci.basic_exp[,2]>=sk_exp),mean(ci.perc_exp[,1]<=sk_exp & ci.perc_exp[,2]>=sk_exp),mean(ci.norm_exp[,1]>=sk_exp ),mean(ci.basic_exp[,1]>=sk_exp),mean(ci.perc_exp[,1]>=sk_exp),mean(ci.norm_exp[,2]<=sk_exp ),mean(ci.basic_exp[,2]<=sk_exp),mean(ci.perc_exp[,2]<=sk_exp))
# Monte Carlo estimates of CI and miss rates of both sides for exponential
d<-matrix(CR_exp,3,3,byrow=T)
colnames(d)<-c("norm","basic","perc")
rownames(d)<-c("CR for exp(2)","miss rates on the left","miss rates on the right")
pander::pandoc.table(d)


## -----------------------------------------------------------------------------
rm(list=ls())
# remove variables
set.seed(123)
alpha <- .1
n <- 30
m <- 2500
a<-seq(1,20) # parameter of Beta distribution
#epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
#N <- length(epsilon)
N <- length(a)
power <- numeric(N)
sk <- function(x) {
  # computes the sample skewness coeff.
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
# critical value for the skewness test
for (j in 1:N) { #for each parameter of Beta distribution
e <- a[j]
stests <- numeric(m)
for (i in 1:m) { #for each replicate
  x<- rbeta(n,e,e)
  stests[i] <- as.integer(abs(sk(x)) >= cv)
}
power[j] <- mean(stests)
}
#plot power vs parameter of Beta distribution
plot(a, power, type = "b", xlab = bquote(alpha), ylim = c(0,1))
abline(h = .1, lty = 3)
se <- sqrt(power * (1-power) / m) #add standard errors
lines(a, power+se, lty = 3)
lines(a, power-se, lty = 3)

## ----results='asis'-----------------------------------------------------------
power_2 <- power_3 <- numeric(N)
for (j in 1:N) { #for each parameter of t distribution
e <- a[j]
stests_2 <- numeric(m)
for (i in 1:m) { #for each replicate
  x <- rt(n,e)
  stests_2[i] <- as.integer(abs(sk(x)) >= cv)
}
power_2[j] <- mean(stests_2)
}
for (j in 1:N) { #for each parameter of uniform distribution
e <- a[j]
stests_3 <- numeric(m)
for (i in 1:m) { #for each replicate
  x <- runif(n,0,e)
  stests_3[i] <- as.integer(abs(sk(x)) >= cv)
}
power_3[j] <- mean(stests_3)
}
plot(a,power, type = "l", xlab = "parameter", ylim = c(0,1),ylab="power")
lines(a,power_2,lty=2)
lines(a,power_3,lty=3)
abline(h = alpha, lty = 4)
legend("topright",1,c("normal","t","unif"),lty=c(1,2,3),inset=.02)

## ----results='asis'-----------------------------------------------------------
rm(list=ls())
# remove variables
set.seed(0321)
n <- 20
# sample size
m<-1e4
# number of replication
p <- matrix(0,nrow=m,ncol=3)
p.hat <- numeric(3)
alpha<-.05
mu0 <- 1

for(j in 1:m){
  x <- matrix(0,n,3)
  x[,1] <- rchisq(n,1) # random numbers from chi-square distribution
  x[,2] <- runif(n,0,2) # random numbers from uniform distribution
  x[,3] <- rexp(n,1) # random numbers from exponential distribution
  for (k in 1:3){
    ttest <- t.test(x[,k],alternative="greater",mu=mu0)
    p[j,k] <- ttest$p.value
  }
}
p.hat <- colMeans(p<alpha)
c<-matrix(p.hat,1,3)
#Monte Carlo estimates of confidence level
colnames(c)<-c("Chi-square","Uniform","Exponential")
rownames(c)<-"p.hat"
pander::pandoc.table(c)

## ----results='asis'-----------------------------------------------------------
rm(list=ls())
# remove variables
library(bootstrap)
# to obtain data
n<-nrow(scor)
# sample size for jackknife
theta.j<-numeric(n)
lamda.h<-eigen((n-1)*cov(scor)/n,only.values=T)$values
theta.h<-max(lamda.h)/sum(lamda.h)
# value of theta.hat
for(i in 1:n){
  y<-scor[-i,]
  x<-(n-2)*cov(y)/(n-1)
  # MLE of sigma
  lamda<-eigen(x,only.values=T)$values
  theta.j[i]<-max(lamda)/sum(lamda)
}
# theta.hat with jackknife
bias.j<-(n-1)*(mean(theta.j)-theta.h)
# bias of theta.hat with jackknife
se.j <- sqrt((n-1)*mean((theta.j-mean(theta.j))^2))
# standard error of theta.hat with jackknife
c<-matrix(round(c(original=theta.h,bias.jack=bias.j,se.jack=se.j),3),1,3)
colnames(c)<-c("original","basic.jack","se.jack")
pander::pandoc.table(c)

## ----results='asis'-----------------------------------------------------------
rm(list=ls())
# remove variables
library(DAAG);attach(ironslag)
n <- length(magnetic) #in DAAG ironslag
er1 <- er2 <- er3 <- er4 <- numeric(n)
# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:n) {
  y <- magnetic[-k]
  x <- chemical[-k]
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
  er1[k] <- magnetic[k] - yhat1
  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +J2$coef[3] * chemical[k]^2
  er2[k] <- magnetic[k] - yhat2
  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
  yhat3 <- exp(logyhat3)
  er3[k] <- magnetic[k] - yhat3
  J4 <- lm(y ~ x+I(x^2)+I(x^3))
  yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k]+J4$coef[3] * chemical[k]^2+J4$coef[4] * chemical[k]^3
  er4[k] <- magnetic[k] - yhat4
}
b<-matrix(c(mean(er1^2), mean(er2^2), mean(er3^2), mean(er4^2)),1,4)
colnames(b)<-c("Linear","Quadratic","Exponential","Cubic polynomial model")
rownames(b)<-"Error"
pander::pandoc.table(b)

## -----------------------------------------------------------------------------
#par(mfrow=c(2,2))
data(ironslag,package="DAAG")
a <- seq(10, 40, .1) #sequence for plotting fits
L1 <- lm(magnetic ~ chemical)
plot(chemical, magnetic, main="Linear", pch=16)
yhat1 <- L1$coef[1] + L1$coef[2] * a
lines(a, yhat1, lwd=2)
L2 <- lm(magnetic ~ chemical + I(chemical^2))
plot(chemical, magnetic, main="Quadratic", pch=16)
yhat2 <- L2$coef[1] + L2$coef[2] * a + L2$coef[3] * a^2
lines(a, yhat2, lwd=2)
L3 <- lm(log(magnetic) ~ chemical)
plot(chemical, magnetic, main="Exponential", pch=16)
logyhat3 <- L3$coef[1] + L3$coef[2] * a
yhat3 <- exp(logyhat3)
lines(a, yhat3, lwd=2)
L4 <- lm(magnetic ~ chemical + I(chemical^2)+I(chemical^3))
plot(chemical, magnetic, main="Cubic polynomial model", pch=16)
yhat4 <- L4$coef[1] + L4$coef[2] * a + L4$coef[3] * a^2+L4$coef[4] * a^3
lines(a, yhat4, lwd=2)

## ----results='asis'-----------------------------------------------------------
AR<-numeric(4)
AR[1]<-summary(L1)[9]
AR[2]<-summary(L2)[9]
AR[3]<-summary(L3)[9]
AR[4]<-summary(L4)[9]
ajust_R2<-round(as.numeric(AR),3)
c<-matrix(ajust_R2,1,4)
colnames(c)<-c("Linear","Quadratic","Exponential","Cubic polynomial model")
rownames(c)<-"adj.r.squared"
pander::pandoc.table(c)


## ----results='asis'-----------------------------------------------------------
rm(list=ls())
# remove variables
library(boot);library(Ball)
# to obtain function
set.seed(234)
m<-1e3;R<-999
n1<-20;n2<-30
# different sample sizes
N<-c(n1,n2)
mu1<-mu2<-0
sigma1<-sigma2<-1
# equal variance
sigma3<-2
# unequal variance
Tn <- function(d,ix,sizes){
  n1<-sizes[1];n2<-sizes[2]
  n<-n1+n2
  ii<-ix[1:n1];jj<-ix[(n1+1):n]
  outx<-sum(d[ii]>max(d[jj]))+sum(d[ii]<min(d[jj]))
  outy<-sum(d[jj]>max(d[ii]))+sum(d[jj]<min(d[ii]))
  return(max(outx,outy))
}
# function to obtain maximum number of extreme points
myeq.test<-function(d){
  boot.obj<-boot(data=d,R=R,statistic=Tn,sim="permutation",sizes=N)
  ts<-c(boot.obj$t0,boot.obj$t)
  p.value<-mean(ts >= ts[1])
  return(p.value)
}
# function to obtain p.value of permutation test
p.values<-matrix(NA,m,4)
for(i in 1:m){
  x<-rnorm(n1,mu1,sigma1)
  y<-rnorm(n2,mu2,sigma2)
  z<-rnorm(n2,mu2,sigma3)
  d1<-c(x,y);d2<-c(x,z)
  p.values[i,1]<-bd.test(x=x,y=y,R=999,seed=i*123)$p.value
  # p.value of Ball test with H0 being TURE 
  p.values[i,2]<-myeq.test(d1)
  # p.value of permutation test with H0 being TURE 
  p.values[i,3]<-bd.test(x=x,y=z,R=999,seed=i*123)$p.value
  # p.value of Ball test with H0 being FALSE
  p.values[i,4]<-myeq.test(d2)
  # p.value of permutation test with H0 being FALSE
}
c<-matrix(colMeans(p.values),1,4)
colnames(c)<-c("bd.test(H0=T)","myeq.test(H0=T)","bd.test(H0=F)","myeq.test(H0=F)")
rownames(c)<-"p.value"
pander::pandoc.table(c)

## -----------------------------------------------------------------------------
set.seed(111)
ox<-rnorm(n1,mu1,sigma1)
oy<-rnorm(n2,mu2,sigma2)
od<-c(ox,oy)
obj.o<-boot(data=od,R=999,statistic=Tn,sim="permutation",sizes=N)
tb <- c(obj.o$t,obj.o$t0)
hist(tb, freq=FALSE, main="",xlab="replicates of maximum number of extreme points")
abline(v=obj.o$t0,col='red',lwd=2)

## ----results='asis'-----------------------------------------------------------
#rm(list=ls())
## remove variables
#library(Ball);library(MASS);library(boot)
#par(mfrow=c(1,2))
##set.seed(123)
#M<-100
##sample size for power
#alpha <- 0.1
#sigma<-matrix(c(1,0,0,1),2,2)
#dCov <- function(x, y){
#  x <- as.matrix(x)
#  y <- as.matrix(y)
#  n <- nrow(x)
#  m <- nrow(y)
#  if (n != m || n < 2) stop("Sample sizes must agree")
#  if (! (all(is.finite(c(x, y)))))
#  stop("Data contains missing or infinite values")
#  Akl <- function(x) {
#    d <- as.matrix(dist(x))
#    m <- rowMeans(d)
#    M <- mean(d)
#    a <- sweep(d, 1, m)
#    b <- sweep(a, 2, m)
#    return(b + M)
#  }
#  A<-Akl(x)
#  B<-Akl(y)
#  dCov<-sqrt(mean(A*B))
#  dCov
#}
#ndcov2 <- function(z, ix, dims) {
##dims contains dimensions of x and y
#p <- dims[1]
#q1 <- dims[2] + 1
#d <- p + dims[2]
#x <- z[ , 1:p] #leave x as is
#y <- z[ix, q1:d] #permute rows of y
#return(nrow(z) * (dCov(x,y))^2)
#}
#size<-seq(10,205,by=5)
## sample size sequence
#for(j in 1:40){
#  n<-size[j]
#  pow<-matrix(NA,40,4)
#  p.values<-matrix(NA,M,4)
#  for(i in 1:M){
#    set.seed(i*j*321)
#    x<-mvrnorm(n,rep(0,2),sigma)
#    e<-mvrnorm(n,rep(0,2),sigma)
#    y1<-x/4+e
#    y2<-x/4*e
#    p.values[i,1]<-bcov.test(x=x,y=y1,R=10,seed=i*121)$p.value
#    # p-value of Ball correlation test of model 1
#    m1<-matrix(cbind(x,y1),n,4)
#    boot.o1<-boot(data=m1,statistic=ndcov2,R=99,sim="permutation",dims=c(2,2))
#    t.b1<-c(boot.o1$t0,boot.o1$t)
#    p.values[i,2]<-mean(t.b1>=t.b1[1])
#    # p-value of distance correlation test of model 1
#    p.values[i,3]<-bcov.test(x=x,y=y2,R=99,seed=i*123)$p.value
#    # p-value of Ball correlation test of model 2
#    m2<-matrix(cbind(x,y2),n,4)
#    boot.o2<-boot(data=m2,statistic=ndcov2,R=99,sim="permutation",dims=c(2,2))
#    t.b2<-c(boot.o2$t0,boot.o2$t)
#    p.values[i,4]<-mean(t.b2>=t.b2[1])
#    # p-value of distance correlation test of model 2
#  }
#  pow[j,] <- colMeans(p.values<alpha)
#}
#plot(size,pow[,1], type = "l", xlab = "n", ylim = c(0,1),ylab="power",main="Model 1")
#  lines(size,pow[,2],lty=2)
#  legend("topleft",1,c("Ball","dCov"),lty=c(1,2),inset=.02)
#plot(size,pow[,3], type = "l", xlab = "n", ylim = c(0,1),ylab="power",main="Model 2")
#lines(size,pow[,4],lty=2)
#legend("topleft",1,c("Ball","dCov"),lty=c(1,2),inset=.02)

## ----results='asis'-----------------------------------------------------------
rm(list=ls())
# remove variables
library(GeneralizedHyperbolic)
# to obtain function
set.seed(123)
N <- 2000
# iteration time
sigma <- c(.05, .5, 2, 16)
# different choices of sigma
x0<-25
# initial setting
myrw.Metropolis <- function(a, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], a)
    if (u[i] <= (dskewlap(y) / dskewlap(x[i-1])))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      }
  }
  return(list(x=x, k=k))
}
# function to generate the chain
r1 <- myrw.Metropolis(sigma[1], x0, N)
r2 <- myrw.Metropolis(sigma[2], x0, N)
r3 <- myrw.Metropolis(sigma[3], x0, N)
r4 <- myrw.Metropolis(sigma[4], x0, N)
c<-matrix(c(r1$k/N, r2$k/N, r3$k/N, r4$k/N),1,4)
colnames(c)<-c("sigma=0.05","sigma=0.5","sigma=2","sigma=16")
rownames(c)<-"rejection rate"
pander::pandoc.table(c)

## -----------------------------------------------------------------------------
#par(mfrow=c(2,2))
plot(r1$x,type='l',ylab="X",xlab=expression(paste(sigma,'=0.05')))
# generated sample vs the time index when the variance of proposal distribution is 0.05
plot(r2$x,type='l',ylab="X",xlab=expression(paste(sigma,'=0.5')))
# generated sample vs the time index when the variance of proposal distribution is 0.5
 abline(h=-3,col='red',lwd=.3)
 abline(h=3,col='red',lwd=.3)
 # convergence range
plot(r3$x,type='l',ylab="X",xlab=expression(paste(sigma,'=2')))
# generated sample vs the time index when the variance of proposal distribution is 4
abline(h=-3,col='red',lwd=.3)
abline(h=3,col='red',lwd=.3)
# convergence range
plot(r4$x,type='l',ylab="X",xlab=expression(paste(sigma,'=16')))
# generated sample vs the time index when the variance of proposal distribution is 16
abline(h=-3,col='red',lwd=.3)
abline(h=3,col='red',lwd=.3)# convergence range 


## -----------------------------------------------------------------------------
#par(mfrow=c(2,2))
b <- 201 ; d<-101
#discard the burnin sample
y1 <- r2$x[b:N]
y2 <- r3$x[d:N]
a <- ppoints(100)
QR <- qskewlap(1-a) 
# quantiles of Rayleigh
Q1 <- quantile(y1, a)
Q2 <- quantile(y2, a)
qqplot(QR, Q1, main="",xlab="standard Laplace distibution", ylab="Sample Quantiles")
hist(y1, breaks="scott", main="", xlab="", freq=FALSE)
lines(QR, dskewlap(QR))
qqplot(QR, Q2, main="",xlab="standard Laplace distibution", ylab="Sample Quantiles")
hist(y2, breaks="scott", main="", xlab="", freq=FALSE)
lines(QR, dskewlap(QR))

## ----results='asis'-----------------------------------------------------------
rm(list=ls())
# remove variables
x<-1e-50
y1<-log(exp(x))
y2<-exp(log(x))
c<-matrix(c(x,y1,y2),1,3)
colnames(c)<-c("x","log(exp(x))","exp(log(x))")
rownames(c)<-"value"
pander::pandoc.table(c)

## -----------------------------------------------------------------------------
isTRUE(all.equal(x,y1))

## -----------------------------------------------------------------------------
rm(list=ls())
K<-c(4,9)
#par(mfrow=c(1,2))
A1<-seq(0,sqrt(K[1]),0.01)
A2<-seq(0,sqrt(K[2]),0.01)
g1<-function(x){(1+(x^2)/(k-1))^(-k/2)}
g2<-function(x){(1+(x^2)/k)^(-(k+1)/2)}
ck<-function(k,a){
  return(sqrt(a^2*k/(k+1-a^2)))
}
y11<-y12<-numeric(length(A1))
y21<-y22<-numeric(length(A2))
for(i in 1:length(A1)){
  a<-A1[i]
  k<-K[1]
y11[i]<-2*gamma(k/2)/(sqrt(pi*(k-1))*gamma((k-1)/2))*integrate(g1, lower=0, upper=ck(k-1,a))$value
y12[i]<-2*gamma((k+1)/2)/(sqrt(pi*k)*gamma(k/2))*integrate(g2,lower=0,upper=ck(k,a))$value
}
plot(A1,y11,type='l',col='red')
lines(A1,y12,col='blue')
for(i in 1:length(A2)){
  a<-A2[i]
  k<-K[2]
y21[i]<-2*gamma(k/2)/(sqrt(pi*(k-1))*gamma((k-1)/2))*integrate(g1, lower=0, upper=ck(k-1,a))$value
y22[i]<-2*gamma((k+1)/2)/(sqrt(pi*k)*gamma(k/2))*integrate(g2,lower=0,upper=ck(k,a))$value
}
plot(A2,y21,type='l',col='red')
lines(A2,y22,col='blue')


## ----results='asis'-----------------------------------------------------------
rm(list=ls())
# remove variables
K<-c(4:25,100)
# different k
N<-length(K)
f1 <- function(x) {
  2*gamma(k/2)/(sqrt(pi*(k-1))*gamma((k-1)/2))*(1+(x^2)/(k-1))^(-k/2)
}
# integrand 1
f2<-function(x) {
  2*gamma((k+1)/2)/(sqrt(pi*k)*gamma(k/2))*(1+(x^2)/k)^(-(k+1)/2)
}
# integrand 2
ck<-function(k,a){
  return(sqrt(a^2*k/(k+1-a^2)))
}
# integral upper limit
A<-numeric(N)
for(i in 1:N){
  k<-K[i]
  ss<-uniroot(function(a){
  integrate(f1, lower=0, upper=ck(k-1,a))$value-integrate(f2,lower=0,upper=ck(k,a))$value
    },lower=0.5,upper=sqrt(k)/2+1)
  A[i]<-ss$root
}

## ----results='asis'-----------------------------------------------------------
rm(list=ls())
# remove variables
library(rootSolve)
N <- 2e2
# max. number of the iteration
nA<-28;nB<-24;nOO<-41;nAB<-70
G <- c(.5,.4)
# initial estimates 
tol <- .Machine$double.eps^0.5
G.old <- G+1
E <- numeric(N)
for(j in 1:N){
  E[j]<-2*G[1]*nA*log(G[1])/(2-G[1]-2*G[2])+2*G[2]*nB*log(G[2])/(2-G[2]-2*G[1])+2*nOO*log(1-G[1]-G[2])+nA*(2-2*G[1]-2*G[2])*log(2*G[1]*(1-G[1]-G[2]))/(2-G[1]-2*G[2])+nB*(2-2*G[1]-2*G[2])*log(2*G[2]*(1-G[1]-G[2]))/(2-G[2]-2*G[1])+nAB*log(2*G[1]*G[2])
  model<-function(x){
    F1<-2*G[1]*nA/((2-G[1]-2*G[2])*x[1])-2*nOO/(1-x[1]-x[2])+nA*(2-2*G[1]-2*G[2])*(1-2*x[1]-x[2])/((2-G[1]-2*G[2])*x[1]*(1-x[1]-x[2]))-nB*(2-2*G[1]-2*G[2])/((2-G[2]-2*G[1])*(1-x[1]-x[2]))+nAB/x[1]
    F2<-2*G[2]*nB/((2-G[2]-2*G[1])*x[2])-2*nOO/(1-x[1]-x[2])-nA*(2-2*G[1]-2*G[2])/((2-G[1]-2*G[2])*(1-x[1]-x[2]))+nB*(2-2*G[1]-2*G[2])*(1-2*x[2]-x[1])/((2-G[2]-2*G[1])*x[2]*(1-x[1]-x[2]))+nAB/x[2]
    c(F1=F1,F2=F2)
  }
  ss<-multiroot(f=model,star=c(.1,.1))
  G<-ss$root
  # update p and q
  if (sum(abs(G-G.old)/G.old)<tol) break
  G.old<-G
}
print(G.old)

## -----------------------------------------------------------------------------
plot(E[1:11],type='l',xlab="iteration time",ylab="likelihood value")

## ----results='asis'-----------------------------------------------------------
rm(list=ls())
# remove variables
attach(mtcars)
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
fit11<-vector("list",length(formulas))
cof11<-matrix(0,4,3)
for(i in seq_along(formulas)){
  fit11[[i]]<-lm(formulas[[i]])
  cof11[i,1:length(fit11[[i]]$coefficients)]<-fit11[[i]]$coefficients
}
colnames(cof11)<-c("intercept","coefficient1","coefficient2")
rownames(cof11)<-c("mpg ~ disp","mpg ~ I(1 / disp)","mpg ~ disp + wt","mpg ~ I(1 / disp) + wt")
pander::pandoc.table(cof11)

## ----results='asis'-----------------------------------------------------------
fit12<-lapply(formulas,lm)
cof12<-matrix(0,4,3)
for(i in seq_along(formulas)){
  cof12[i,1:length(fit12[[i]]$coefficients)]<-fit12[[i]]$coefficients
}
colnames(cof12)<-c("intercept","coefficient1","coefficient2")
rownames(cof12)<-c("mpg ~ disp","mpg ~ I(1 / disp)","mpg ~ disp + wt","mpg ~ I(1 / disp) + wt")
pander::pandoc.table(cof12)

## ----results='asis'-----------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})
fit21<-vector("list",length(bootstraps))
cof21<-matrix(0,10,2)
for(i in seq_along(bootstraps)){
  fit21[[i]]<-lm(bootstraps[[i]]$mpg~bootstraps[[i]]$disp)
  cof21[i,]<-fit21[[i]]$coefficients
}
colnames(cof21)<-c("intercept","coefficient")
rownames(cof21)<-c("rep1","rep2","rep3","rep4","rep5","rep6","rep7","rep8","rep9","rep10")
pander::pandoc.table(cof21)

## ----results='asis'-----------------------------------------------------------
fit22<-lapply(bootstraps,function(x){
  lm(x$mpg~x$disp)
})
cof22<-matrix(0,10,2)
for(i in seq_along(bootstraps)){
  cof22[i,]<-fit22[[i]]$coefficients
}
colnames(cof22)<-c("intercept","coefficient")
rownames(cof22)<-c("rep1","rep2","rep3","rep4","rep5","rep6","rep7","rep8","rep9","rep10")
pander::pandoc.table(cof22)

## ----results='asis'-----------------------------------------------------------
fit23<-lapply(bootstraps,lm,formula=mpg~disp)
cof23<-matrix(0,10,2)
for(i in seq_along(bootstraps)){
  cof23[i,]<-fit23[[i]]$coefficients
}
colnames(cof23)<-c("intercept","coefficient")
rownames(cof23)<-c("rep1","rep2","rep3","rep4","rep5","rep6","rep7","rep8","rep9","rep10")
pander::pandoc.table(cof23)

## ----results='asis'-----------------------------------------------------------
rsq <- function(mod) summary(mod)$r.square
rsq1 <- lapply(fit12,rsq)
r_squared1<-matrix(rsq1,4,1)
colnames(r_squared1)<-"R_squared"
rownames(r_squared1)<-c("mpg ~ disp","mpg ~ I(1 / disp)","mpg ~ disp + wt","mpg ~ I(1 / disp) + wt")
pander::pandoc.table(r_squared1)

## ----results='asis'-----------------------------------------------------------
rsq2 <- lapply(fit23,rsq)
r_squared2<-matrix(rsq2,10,1)
colnames(r_squared2)<-"R_squared"
rownames(r_squared2)<-c("rep1","rep2","rep3","rep4","rep5","rep6","rep7","rep8","rep9","rep10")
pander::pandoc.table(r_squared2)

## ----results='asis'-----------------------------------------------------------
rm(list=ls())
# remove variables
trials <- replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)
p_value<-sapply(trials,function(x) return(x$p.value))
p<-matrix(p_value[1:5],5,1)
colnames(p)<-"p_value"
rownames(p)<-c("trial1","trial2","trial3","trial4","trial5")
pander::pandoc.table(p)

## ----results='asis'-----------------------------------------------------------
p_value2<-sapply(trials,'[[',3)
p2<-matrix(p_value2[1:5],5,1)
colnames(p2)<-"p_value"
rownames(p2)<-c("trial1","trial2","trial3","trial4","trial5")
pander::pandoc.table(p2)

## -----------------------------------------------------------------------------
rm(list=ls())
# remove variables
library(parallel)
attach(mtcars)
cores <- detectCores()
cluster <- makePSOCKcluster(cores)
#  set up a local cluster
mcsapply <- function(x, f, ...) {
res <- parLapply(cluster, x, f)
simplify2array(res)
} 


## ----results='asis'-----------------------------------------------------------
boot_df <- function(x) x[sample(nrow(x), rep = T), ]
boot_lm <- function(i) {
  return(summary(lm(mpg ~ wt + disp, data = mtcars[sample(nrow(mtcars), rep = T), ]))$r.square)
}
T1<-system.time(sapply(1:500, boot_lm))
T2<-system.time(mcsapply(1:500, boot_lm))
c<-matrix(c(T1[1:3],T2[1:3]),2,3)
colnames(c)<-c("user","system","elapsed")
rownames(c)<-c("sapply","mcsapply")
pander::pandoc.table(c)

## ----results='asis'-----------------------------------------------------------
rm(list=ls())
# remove variables
library(GeneralizedHyperbolic)
# to obtain function
library(Rcpp)
# Attach R package "Rcpp"
#dir_cpp<-'Rcpp/'
# Can create source file in Rstudio
#sourceCpp(paste0(dir_cpp,"MetropolisC.cpp"))
library(microbenchmark)
set.seed(123)
N <- 2000
# iteration time
sigma <- c(.05, .5, 2, 16)
# different choices of sigma
x0<-25
# initial setting
myrw.Metropolis <- function(a, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], a)
    if (u[i] <= (dskewlap(y) / dskewlap(x[i-1])))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      }
  }
  return(list(x=x, k=k))
}
# function to generate the chain
r1 <- myrw.Metropolis(sigma[1], x0, N)
r2 <- myrw.Metropolis(sigma[2], x0, N)
r3 <- myrw.Metropolis(sigma[3], x0, N)
r4 <- myrw.Metropolis(sigma[4], x0, N)
# results by R
#r11 <-MetropolisC(sigma[1], x0, N)
#r21 <- MetropolisC(sigma[2], x0, N)
#r31 <- MetropolisC(sigma[3], x0, N)
#r41 <- MetropolisC(sigma[4], x0, N)
# results by Rcpp
c<-matrix(c(r1$k/N, r2$k/N, r3$k/N, r4$k/N),1,4)
colnames(c)<-c("sigma=0.05","sigma=0.5","sigma=2","sigma=16")
rownames(c)<-"rejection rate(R)"
pander::pandoc.table(c)
#b<-matrix(c(r11[N+1]/N, r21[N+1]/N, r31[N+1]/N, r41[N+1]/N),1,4)
#colnames(b)<-c("sigma=0.05","sigma=0.5","sigma=2","sigma=16")
#rownames(b)<-"rejection rate(Rcpp)"
#pander::pandoc.table(b)

## -----------------------------------------------------------------------------
#par(mfrow=c(2,2))
#plot(r11[1:N],type='l',ylab="X",xlab=expression(paste(sigma,'=0.05')))
# generated sample vs the time index when the variance of proposal distribution is 0.05
#plot(r21[1:N],type='l',ylab="X",xlab=expression(paste(sigma,'=0.5')))
# generated sample vs the time index when the variance of proposal distribution is 0.5
 #abline(h=-3,col='red',lwd=.3)
 #abline(h=3,col='red',lwd=.3)
 # convergence range
#plot(r31[1:N],type='l',ylab="X",xlab=expression(paste(sigma,'=2')))
# generated sample vs the time index when the variance of proposal distribution is 4
#abline(h=-3,col='red',lwd=.3)
#abline(h=3,col='red',lwd=.3)
# convergence range
#plot(r41[1:N],type='l',ylab="X",xlab=expression(paste(sigma,'=16')))
# generated sample vs the time index when the variance of proposal distribution is 16
#abline(h=-3,col='red',lwd=.3)
#abline(h=3,col='red',lwd=.3)# convergence range 


## -----------------------------------------------------------------------------
#par(mfrow=c(2,2))
b <- 201 ; d<-101
#discard the burnin sample
y1 <- r2$x[b:N]
#y11<-r21[b:N]
y2 <- r3$x[d:N]
#y21<-r31[b:N]
a <- ppoints(100)
#QR <- qskewlap(1-a) 
# quantiles of Rayleigh
Q1 <- quantile(y1, a)
#Q11<-quantile(y11, a)
Q2 <- quantile(y2, a)
#Q21<-quantile(y21, a)
#qqplot(Q1, Q11, main="",xlab="R Sample Quantiles", ylab="Rcpp Sample Quantiles")
#hist(y11, breaks="scott", main="", xlab="", freq=FALSE)
#lines(QR, dskewlap(QR))
#qqplot(Q2, Q21, main="",xlab="R Sample Quantiles", ylab="Rcpp Sample Quantiles")
#hist(y21, breaks="scott", main="", xlab="", freq=FALSE)
#lines(QR, dskewlap(QR))

## -----------------------------------------------------------------------------
#ts <- microbenchmark(MetropolisR=myrw.Metropolis(sigma[1], x0, N),
#                     MetropolisCpp=MetropolisC(sigma[1], x0, N))
#summary(ts)[,c(1,3,5,6)]

