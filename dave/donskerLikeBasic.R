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

for(k in seq(-5,5,by=1)){
  ll = evalLogLike(x+0,Y+k)
  print(ll)
}
## Wow!  Each population has exactly the same likelihood if the
## sample (or population) is shifted.

## This is also the case if we just look at different draws, the likelihood
## depends only on n
for(k in seq(-5,5,by=1)){
  x = sort(sample(Y,n,replace=F))
  ll = evalLogLike(x+0,Y+0)
  print(ll)
}

exit()

######---------------#####
 # suppose we observe n=3 u's:
n = 3;
uobs = c(.1,.11,.12)  # not very consistent with neff=400 ==> big likelihood
uobs = c(.25,.5,.75)  # very consistent with neff=400 ==> small likelihood
Fnuobs = (1:n)/n
 # evaluate F at the u values (we'll use F for pop=u2[,666])
Fuobs = Feval(uobs,u2[,666])
 # get the brownian bridge covariance
SigF = bbCov(Fuobs)
 # compute the likelihood of uobs|Y=u2[,666]
dmvnorm(Fnuobs, Fuobs, SigF, log=TRUE)

 # now, let's compute importance weights for all of the prior
 # draws u1, u2, or u3.  The likelihood is L(uobs|Y=y)

 # initialize weight vectors
w1 = rep(NA,NMC)
w2 = rep(NA,NMC)
w3 = rep(NA,NMC)

uobs = sort(runif(10))
uobs = seq(.1,.9,by=.1)  # very consistent with neff=400 ==> small likelihood
uobs = c(.1,.11,.12,.2,.22)
uobs = c(.25,.5,.75)  # very consistent with neff=400 ==> small likelihood
uobs = c(.1,.11,.12,.2,.22,.25)
n = length(uobs)
Fnuobs = (1:n)/n

 # loop through the NMC realizations and compute weights
for(i in 1:NMC){
  Fuobs = Feval(uobs,u1[,i])
  SigF = bbCov(Fuobs)
  w1[i] = dmvnorm(Fnuobs, Fuobs, SigF, log=TRUE) # L(uobs|Y=u1[,i])
}
# compute 20 draws from this importance distribution
isSamp1 = sample(1:NMC, size=1500, prob=exp(w1))

for(i in 1:NMC){
  Fuobs = Feval(uobs,u2[,i])
  SigF = bbCov(Fuobs)
  w2[i] = dmvnorm(Fnuobs, Fuobs, SigF, log=TRUE) # L(uobs|Y=u2[,i])
}
# compute 15 draws from this importance distribution
isSamp2 = sample(1:NMC, size=1500, prob=exp(w2))

for(i in 1:NMC){
  Fuobs = Feval(uobs,u3[,i])
  SigF = bbCov(Fuobs)
  w3[i] = dmvnorm(Fnuobs, Fuobs, SigF, log=TRUE) # L(uobs|Y=u3[,i])
}
 # compute 20 draws from this importance distribution
isSamp3 = sample(1:NMC, size=1500, prob=exp(w3))

 # have a look
 # draws from the prior
par(mfrow=c(3,3),oma=c(0,0,0,0),mar=c(4,4,1.5,1))
 # neff = N
plot(c(0,1),c(0,1),type='l',xlab='u',ylab='cum prob')
matlines(u1[,1:15],cumsum(rep(1/N,N)),type='s')
mtext('draws of y from prior',side=3,line=.3)
plot(c(0,1),c(0,1),type='l',xlab='u',ylab='cum prob')
matlines(u1[,isSamp1[1:15]],cumsum(rep(1/N,N)),type='s')
mtext('draws of y from posterior',side=3,line=.3)
points(uobs,Fnuobs,pch=16)
hist(apply(u1[,isSamp1],2,median),xlim=c(0,1),xlab='median(y)',main='neff=N')
 # neff = 400
plot(c(0,1),c(0,1),type='l',xlab='u',ylab='cum prob')
matlines(u2[,1:15],cumsum(rep(1/N,N)),type='s')
mtext('draws of y from prior',side=3,line=.3)
plot(c(0,1),c(0,1),type='l',xlab='u',ylab='cum prob')
matlines(u2[,isSamp2[1:15]],cumsum(rep(1/N,N)),type='s')
mtext('draws of y from posterior',side=3,line=.3)
points(uobs,Fnuobs,pch=16)
hist(apply(u2[,isSamp2],2,median),xlim=c(0,1),xlab='median(y)',main='neff=400')
# neff = 1
plot(c(0,1),c(0,1),type='l',xlab='u',ylab='cum prob')
matlines(u3[,1:15],cumsum(rep(1/N,N)),type='s')
mtext('draws of y from prior',side=3,line=.3)
plot(c(0,1),c(0,1),type='l',xlab='u',ylab='cum prob')
matlines(u3[,isSamp3[1:15]],cumsum(rep(1/N,N)),type='s')
mtext('draws of y from posterior',side=3,line=.3)
points(uobs,Fnuobs,pch=16)
hist(apply(u3[,isSamp3],2,median),xlim=c(0,1),xlab='median(y)',main='neff=1')
#hist(apply(u3[,],2,median),xlim=c(0,1),xlab='median(y)',main='neff=1')




