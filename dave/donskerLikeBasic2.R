 # use the Browninan Bridge likelihood for the empirical distribution Fn
 # to motivate a likelihood for an iid sample u_1,...,u_n

 # Donsker's Thm: sqrt(n)(Fn - F) -> GP(mean=0,Cov(s,t)) with
 # Cov(s,t) = min(s,t)-s*t.
 # This GP lives over (0,1);
 # F(Y) is cumsum(1:N)/N -- ie. mass 1/N for each pop member
 # We observe u[1],...,u[n] iid ~ F(u) =d disc(Y,prob=rep(1/N,N))
 
 # For this file, we'll only look at the locations u_1,...,u_n from
 # our iid sample from disc(Y,prob=rep(1/N,N))

 # population size
N = 40
Y = 1:N
FY = cumsum(rep(1,N))/N
plot(1:N,FY,type='s')

 # sample size n
n = 6
x = 1:n  # a crummy sample
x = round(seq(1,N,length=n))
x = sort(sample(Y,n,replace=T))
lines(c(0,x,N),c(0,(1:n)/n,1),type='s',col='cyan',lty=2)


 # source in the file unifPopDraws.R to produce some prior draws of Y
 # Here the population size is N=41; we have 10000 draws from 3 different
 # priors: u1 (neff = 41); u2 (neff=400); u3 (neff=.01)
# source('~/Dave/Research/PopInf/unifPopDraws.R')
 # population size
# N = nrow(u1)
 # number of draws from the prior
# NMC = ncol(u1)


 # make a function to evaluate F(Y) (which is determined by y and emp dist Fy)
Feval <- function(Y,y,Fy=seq(1,length=length(y))/length(y)){
  # evaluate the distribution F(.) of the population at values u
  cumprobs <- c(0,Fy)
  approx(c(-9e99,y),cumprobs,xout=Y,method="constant",rule=2,ties = list("ordered", max))$y
}

 # try it out
# evaluate the empirical Fy at all values of the pop Y
Feval(1:N,c(17,22,38))
points(1:N,Feval(1:N,c(17,22,38)),pch=16,cex=.5,col='cyan')

 # evaluate the pop Y under the true dist FY 
FY = Feval(Y,Y)
 # evaluate the pop Y under the sample emp dist Fy
FxY = Feval(Y,x)

 # plot the two distributions
plot(Y,FY,type='s')
lines(Y,FxY,type='s',col='cyan')

 # plot their difference
plot(Y,sqrt(n)*(FY-FxY),type='l',ylim=c(-1.3,1.3))

 # this thing follows a brownian bridge distribution

 # now a function to make the Brownian Bridge covariance
bbCov <- function(p){
  # make the Brownian Bridge covariance at locations given by p
  # assume p is ordered
  outer(p,p,FUN="pmin") - outer(p,p,FUN="*")
}
# look at the precision matrix. It is a random walk with the ends
# tied to 0 with variance 1.  Interior increments are N(0,1); the
# ends (1st and last elements) are also N(0,1).
solve(bbCov((1:5)/6))/6

 # try it out
bbCov(c(.1,.5,.9))
image(bbCov(FY[-N]))
range(diag(bbCov(FY[-N])))
 # build the N-1 x N-1 cov matrix
CYbb = bbCov(FY[-N])

 # use mvtnorm to evaluate the multivariate normal density
library(mvtnorm)

 # evaluate the density of FxY | FY
dmvnorm(FY[-N], FxY[-N], (1/n)*CYbb, log=TRUE)

evalLogLike <- function(y,Y){
  # evaluates the log likelihood of the sample y ~ iid(Y)
  FY = Feval(Y,Y)
  FyY = Feval(Y,y)
  # browser()
  # ynew = c(6,12,18,24,32,38)-2; y=ynew  # a very even sample
  # FyY = Feval(Y,y)
  N = length(Y); n = length(y)
  # dd = FY-FyY
  # sum((FY-FyY)^2)
   # manual sum of squares
  # sum(diff(FY-FyY)^2) + (FY[1]-FyY[1])^2 + (FY[N]-FyY[N])^2
  # sum(diff(c(0,dd,0))^2)*N*n
  CYbb = bbCov((1:(N-1))/N)
  # solve(bbCov((1:5)/6))/6
  # t(dd[-N])%*%solve(CYbb)%*%dd[-N]*n # this matches the manual version above
  # dmvnorm(dd[-N], sigma=(1/n)*CYbb, log=TRUE)
  dmvnorm(FY[-N]-FyY[-N], sigma=(1/n)*CYbb, log=TRUE)
}

 # redo this function by only evaluating it at the sample values y
evalLogLike <- function(y,Y){
  # evaluates the log likelihood of the sample y ~ iid(Y), but only
  # at the obseerved locations
  FY = Feval(y,Y)
   # if we see an observation with 0 probability, return -9e99
  if(any(FY==0)) return(-9e99)
  FyY = Feval(y,y)
  N = length(Y); n = length(y)
  # dd = FY-FyY
  # sum((FY-FyY)^2)
  # manual sum of squares
  # sum(diff(FY-FyY)^2) + (FY[1]-FyY[1])^2 + (FY[N]-FyY[N])^2
  # sum(diff(c(0,dd,0))^2)*N*n
  CYbb = bbCov(FY[-n])
  # solve(bbCov((1:5)/6))/6
  # t(dd[-N])%*%solve(CYbb)%*%dd[-N]*n # this matches the manual version above
  # dmvnorm(dd[-N], sigma=(1/n)*CYbb, log=TRUE)
  dmvnorm(FY[-n]-FyY[-n], sigma=(1/n)*CYbb, log=TRUE)
}

x=c(4,5,9,19,31,32)
x=c(4:8)
x=c(8,10,12,14,16)
x=seq(round(.5/n*N),round((n-.5)/n*N),length=5)
for(k in seq(-3,3,by=1)){
  ll = evalLogLike(x+0,Y+k)
  print(ll)
}
## Wow!  Each population has exactly the same likelihood if the
## sample (or population) is shifted.  But this is only true when
## the Brownian Bridge is evaluated at all N locations of the 
## population support.  If we restrict to the n locations given
## by the sample, then the Donsker likelihood behaves more as
## expected.  The 2nd version of evalLogLike() does this.

## This is also the case if we just look at different draws, the likelihood
## depends only on n

 # make some pictures showing the likelihood for different
 # populations.

 # make a populations centered at 0 with different uniform spreads

bY = 2  # spread in population
NY = 12
Y = seq(-NY*bY,NY*bY,length=NY*2+1)

 # for a given y, what is the "best" Y of size N=41 that is 
 # uniformly spaced between -A and A?
y = c(-6,-4,-2,0,2,4,6)

Avals = seq(5,20,length=50)
ll = rep(NA,length(Avals))
for(k in 1:length(Avals)){
  A = Avals[k]
  YA = seq(-A,A,length=41)
  ll[k] = evalLogLike(y,YA)
  print(ll[k])
}
par(mfrow=c(2,1),oma=c(4,4,1,1),mar=c(0,0,0,0))
plot(Avals[ll>-9e9],ll[ll>-9e9],axes=F,xlab='',ylab='likelihood',xlim=c(5,20))
axis(2); box()
plot(range(Avals),c(-9,9),xlab='A',ylab='population',axes=F,type='n')
axis(1); axis(2); box()
 # show data
for(i in 1:length(y)) abline(a=y[i],b=0)
 # show population
for(k in 1:length(Avals)){
  A = Avals[k]
  YA = seq(-A,A,length=41)
  points(rep(A,length(YA)),YA,pch='-',col='grey70')
}



exit()

