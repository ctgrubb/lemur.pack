# look at variants of the empirical distribution Fn
# Plots of realizations; mean(u) and sd(u)
# support is u \in [0,1] where u = F^{-1}(x); F is pop cdf
# The orderstatistics ui = F^{-1}(xi); x1,...,xn iid F

# you may need to run some other code first to load in functions
setwd("~/Dave/Research/PopInf")
FIRST = FALSE
if(FIRST) source("empVariation.R")
if(FIRST) source("dirichletLikeBasic.R")
# a useful function for computing alpha(p)
getin = function(x,y,xout){
  return(yout = approx(c(0, x), c(0, y), xout = xout)$y)
}

library(extraDistr)

## Let's do some checking to see if this Dirichlet-based description
## seems about right?

# population size
N = 4000
Y = 1:N
n = 40
Nmc = 1000
countmat = matrix(NA, nrow = n + 1, ncol = Nmc)
for(k in 1:Nmc) {
  x = sort(sample(Y, n,replace = T));
  countmat[, k] = tabulate(findInterval(Y, x, left.open = TRUE) + 1, nbin = n + 1)
}
# the distribution of these counts looks to be the same for all gaps
apply(countmat, 1, mean) / N
apply(countmat, 1, sd) / N
# expected mean
1 / (n + 1)
#expected sd
sqrt(1 * n / (n + 1) ** 2 / (n + 2))

alpha = x / N * (n + 1)  # or is it n?
xval = rep(1 / (n + 1), n + 1)

# generate Fn^dir and find the values where Fn^dir(u) = i/n
# are these the same (in distribution) as the simulated draws from F?

dsamp = t(rdirichlet(Nmc, alpha = rep((n + 1) / N, N)))  # sum(alpha) is n to match variances
psamp = apply(dsamp, 2, cumsum)
xpsamp = (1:N) / N
# plot the cumulative dirichlet realizations
par(mfrow = c(2,2), oma = c(0, 0, 0, 0), mar = c(4, 4, 1, 1), pty = 's')
matplot(xpsamp, psamp[, 1:10], type = 'l')
# add the pointwise mean to the plot
lines((1:N)/N,apply(psamp,1,mean),type='l',lwd=2)
lines((1:N)/N,apply(psamp,1,mean)+1*apply(psamp,1,sd),type='l',lwd=1,lty=2)
lines((1:N)/N,apply(psamp,1,mean)-1*apply(psamp,1,sd),type='l',lwd=1,lty=2)
abline(c(0,1),lty=3,col='cyan')
# plot the reflection about x=y (i.e. switch x and y)
matplot(psamp[,1:10],xpsamp,type='l')


# invert the realizations (i.e. equal spaced in the y values, use y as support)
usamp = apply(psamp,2,getin,y=xpsamp,xout=(0:N)/N)
matplot((0:N)/N,usamp[,1:10],type='l')
# add the pointwise mean to the plot
lines((0:N)/N,apply(usamp,1,mean),type='l',lwd=2) # slope is not quite 1
lines((0:N)/N,apply(usamp,1,mean)+1*apply(usamp,1,sd),type='l',lwd=1,lty=2)
lines((0:N)/N,apply(usamp,1,mean)-1*apply(usamp,1,sd),type='l',lwd=1,lty=2)
abline(c(0,1),lty=2)
# add in the theoretical expected value line and the (.5,.5) point
abline(c(.5*1/(n+1),n/(n+1)),lty=3,col='red')
points(.5,.5,pch=16,col='red')

# so if I take a(p) = [ alpha(p)-.5*(1/(n+1)) ]*(n+1)/n
# a(p) has mean p; variance (n+1)^2/n^2*var(alpha(p))

# Also, n/n+1*Fn !=d alpha(p), even though both give i/(n+1) at u_i
# Why? E alpha(p) = .5/(n+1) + n/(n+1)*p [see red abline above]


for(k in 1:Nmc){
  x = sort(sample(Y,n,replace=T));
  countmat[,k] = tabulate(findInterval(Y,x,left.open=TRUE)+1,nbin=n+1)
}

# check - Is Fn* essentially the Fn*n/(n+1)?
# compute usamp(1/n+1,...,n/n+1) ?=d order stats of a uniform
#                                ?=d order stats of pop Y?
# just sample the points ui,i/(n+1)
xm = matrix(runif(n*Nmc),ncol=Nmc)
xm = apply(xm,2,sort);
fstarn = matrix((1:n)/(n+1),nrow=n,ncol=Nmc)
matplot(xm,jitter(fstarn),pch='.',ylim=c(0,1))
abline(lsfit(as.vector(xm),as.vector(fstarn))$coef)

# now look at Fn*n/(n+1) over a fine grid
Ygrid = (1:N)/N
fnm = matrix(NA,nrow=N,ncol=Nmc)
for(k in 1:Nmc)  fnm[,k] = Feval(Ygrid,xm[,k])
# get the mean
Mfnm = apply(fnm*n/(n+1),1,mean)
lines(Ygrid,Mfnm,lty=1)
# line with slope of n/(n+1)
abline(c(0,n/(n+1)),lty=3,col='red')

# now try Fn*n/(n+1) but connecting the dots at 0 and 1.
fnL = matrix(NA,nrow=N,ncol=Nmc)
for(k in 1:Nmc)  fnL[,k] = approx(c(0,xm[,k],1),seq(0,1,length=n+2),xout=Ygrid)$y
# get the mean
MfnL = apply(fnL,1,mean)
lines(Ygrid,MfnL,lty=1,col='cyan')
sdfnL = apply(fnL,1,sd)

# make a separate plot
par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(4,4,1,1),pty='s')
matplot((0:N)/N,usamp[,1:10],type='l',col='grey70')
# add the pointwise mean to the plot
lines((0:N)/N,apply(usamp,1,mean),type='l',lwd=2) # slope is not quite 1
lines((0:N)/N,apply(usamp,1,mean)+1*apply(usamp,1,sd),type='l',lwd=1,lty=2)
lines((0:N)/N,apply(usamp,1,mean)-1*apply(usamp,1,sd),type='l',lwd=1,lty=2)
abline(c(0,1),lty=2)
# add in the theoretical expected value line and the (.5,.5) point
abline(c(.5*1/(n+1),n/(n+1)),lty=3,col='red')
points(.5,.5,pch=16,col='red')
# the FnLinear mean
lines(Ygrid,MfnL,lty=3,col='cyan')
lines(Ygrid,MfnL+1*sdfnL,lty=3,col='cyan')
lines(Ygrid,MfnL-1*sdfnL,lty=3,col='cyan')

# now try FnL2 = .5/(n+1) at u1; i/n at ui; 1 at 1 + connecting the dots.
fnL2 = matrix(NA,nrow=N,ncol=Nmc)
yfnL2incs = c(0,.5/(n),rep(1/n,n-1),.5/(n))
yfnL2vals = cumsum(yfnL2incs)
for(k in 1:Nmc)  fnL2[,k] = approx(c(0,xm[,k],1),yfnL2vals,xout=Ygrid)$y
# get the mean
matplot((0:N)/N,usamp[,1:10],type='l',col='grey80')
matlines(Ygrid,fnL2[,1:20],type='l',col='grey60')
MfnL2 = apply(fnL2,1,mean)
lines(Ygrid,MfnL2,lty=1,col='orange')
sdfnL2 = apply(fnL2,1,sd)
matlines(Ygrid,cbind(MfnL2+1*sdfnL2,MfnL2-1*sdfnL2),lty=3,col='orange')
thsd = sqrt(Ygrid*(1-Ygrid)/(n+1))
matlines(Ygrid,cbind(MfnL2+1*thsd,MfnL2-1*thsd),lty=3,col='black')
abline(c(0,1),lty=2)
# this suggests FnL2 is well approximated by Dir(F*n)

# so for a iid sample y from Y, we get the following likelihood
# approximation evalLogLikeDirLines:
evalLogLikeDirLines <- function(y,Y){
  # evaluates the "log likelihood" of the sample y ~ srs(Y)
  #browser()
  n=length(y); N=length(Y)
  x = tabulate(findInterval(Y,y,left.open=TRUE)+1,nbin=n+1)
  alpha = x/N*(n+1)  # or is it n?
  xval = c(.5/n,rep(1/n,n-1),.5/n)
  return(log(ddirichlet(xval,alpha)))   # just from the picture!
}

# compare this likelihood to others with the
# populations centered at 0 with different uniform spreads

bY = 2  # spread in population
NY = 12
Y = seq(-NY*bY,NY*bY,length=NY*2+1)

# for a given y, what is the "best" Y of size N=41 that is
# uniformly spaced between -A and A?
y = c(-6,-4,-2,0,2,4,6)

Avals = seq(6,20,length=50)
ll = rep(NA,length(Avals)); ll2 = ll; ll3 = ll;
for(k in 1:length(Avals)){
  A = Avals[k]
  YA = seq(-A,A,length=41)
  ll[k] = evalLogLikeDir(y,YA)
  ll2[k] = evalLogLike(y,YA)
  ll3[k] = evalLogLikeDirLines(y,YA)
  print(ll[k])
}
par(mfrow=c(2,1),oma=c(4,4,1,1),mar=c(0,0,0,0),pty='m')
plot(Avals[ll>-9e9],ll[ll>-9e9],axes=F,xlab='',ylab='likelihood',xlim=c(5,20))
points(Avals[ll2>-9e9],ll2[ll2>-9e9],col='green')
points(Avals[ll3>-9e9],ll3[ll3>-9e9],col='blue')
axis(2); box(); mtext('likelihood',side=2,line=3)
plot(range(Avals),c(-9,9),xlab='A',ylab='population',axes=F,type='n',xlim=c(5,20))
axis(1); axis(2); box(); mtext('population values',side=2,line=3)
mtext('A',side=1,line=2.5)
# show data
for(i in 1:length(y)) abline(a=y[i],b=0)
# show population
for(k in 1:length(Avals)){
  A = Avals[k]
  YA = seq(-A,A,length=41)
  if(k > 1){
    if(ll[k] != ll[k-1]) points(rep(A,length(YA)),YA,pch='-',col='grey70')
  }
}

