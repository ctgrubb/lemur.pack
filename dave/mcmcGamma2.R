# code to run mcmc for the Dave Process Neff=low Gamma sample.
Neff = 10;
N = 100
n=7
library(gtools)
set.seed(1976)
PDF = FALSE

# The prior for the neff=BIG and neff=N formulations
# priors: mean ~ N(mean=110000,sd=20000)
#           sd ~ Unif(90000-20000,90000+20000)

# note: for the gamma-neff formulation,
# the parameters are fixed so that the mean of
# the gamma is 110000 and the sd is 90000.

 # initialize parameters params = c(mean,sd) of Gamma
xsamp = rgamma(n,shape=1.44,rate=1.316e-5)
xsamp = c(47075.5,126150.4,254994.81,63954.58,16841.55,40743.77,99430.52)
beta = 1.316194e-05
alpha = 1.432562
xpopinit = c(xsamp,rgamma(N-n,shape=alpha,rate=beta))
upopinit = pgamma(xpopinit,shape=alpha,rate=beta)

logPost2=function(upop,neff=length(upop)){
  # log posterior of the result of a cumsum of a dir(c(1,...,1)*Neff/N) realization
  iord = order(upop)
  deltas = diff(c(0,upop[iord],1))
  N = length(upop)
  #browser()
  aa = rep(1,N+1)/(N+1)*(neff+1)
  logpost = (aa-1)*log(deltas)
  return(sum(logpost))
}
logPost2(upopinit,neff=10)

 # number of mcmc iterations
Niter = 40000
 # width of metropolis proposal

uout <- matrix(0,nrow=N,ncol=Niter)

# initialize x[1] with a starting value
uout[,1] <- upopinit

# compute the log of the posterior density for x[1]
logpost <- logPost2(uout[,1],neff=Neff)

# carry out the Metropolis sampling
for(i in 2:Niter){
  uout[,i] = uout[,i-1]
  for(k in (n+1):N){
    can = uout[,i]
    can[k] = runif(1)
    logpostcan = logPost2(can,neff=Neff)
    u <- runif(1)
    if(log(u) < (logpostcan - logpost)){
      # accept the proposal
      uout[k,i] <- can[k]
      logpost <- logpostcan
    } else{
      # uout[k,i] is already set to the previous value - no need to do anything
    }
  }
}
 # thin out uout to 1000
ikeep = round(seq(100,Niter,length=1000))
uout = uout[,ikeep]

par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(4,4,1,1))
plot(1:ncol(uout),uout[8,],type='l')
plot(1:ncol(uout),uout[N,],type='l')
plot(uout[8,],uout[N,],xlab='u8',ylab='uN',pch='.')
hist(uout[8,-(1:100)],main=paste('neff =',Neff,sep=''),xlab='x')

# look at the proportion of acceptances
accept <- !(diff(uout[8,])==0)
mean(accept)

apply(uout, 1, function(x) {sum(!(diff(x) == 0))})

 # a plot of the posterior sample
xout = apply(uout,2,qgamma,shape=alpha,rate=beta)
xmeans = apply(xout,2,mean)
xsds = apply(xout,2,sd)
if(PDF) pdf('postNeffLow.pdf',width=3.9,height=3.4)
par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(4,4,1,1))
plot(xmeans,xsds,xlab='mean',ylab='sd',pch=1,ylim=c(60000,180000),
     xlim=c(70000,140000))
#identify(xmeans,xsds,plot=TRUE,col='red')
iselect = c(647,649,933)
points(xmeans[iselect],xsds[iselect],col=c('blue','red','purple'),pch=16,cex=2)
if(PDF) dev.off()

 # show regular populations corresponding to post draws iselect
N=100
 # Low Neff
if(PDF) pdf('popNeffLow.pdf',width=4.0,height=1.9)
par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(4,1,1,1))
matplot(jitter(xout[,iselect],factor=10),cbind(rep(1,N),rep(2,N),rep(3,N)),pch='|',xlim=c(0,600000),ylim=c(.4,3.6),
        col=c('blue','red','purple'),type='p',axes=F,xlab='household income',ylab='.')
axis(1)
if(PDF) dev.off()

 # show some prior realizations of the cdfs
xline = seq(0,6e5,length=400)
dfline = pgamma(xline,shape=alpha,rate=beta)
par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(4,4,1,1))
plot(xline,dfline,xlim=c(0,600000),ylim=c(0,1),xlab='household income',
     ylab='',axes=T,type='l')

pHigh = rdirichlet(1,rep(1,N+1)*1000)
uHigh = cumsum(pLow)[-(N+1)]
 # make 10 realizations of the population prior
prMeanDraws = rnorm(10,mean=110000,sd=30000)
prSdDraws = runif(10,min=60000,max=180000)
prPopNeffHigh = matrix(NA,nrow=N,ncol=10)
beta = prMeanDraws/prSdDraws^2
alpha = prMeanDraws*beta
for(k in 1:10){
  prPopNeffHigh[,k] = qgamma(uHigh,shape=alpha[k],rate=beta[k])
}
matlines(rbind(0,prPopNeffHigh,1e7),c(0,seq(0,1,length=N+1)),type='s')
lines(xline,dfline,lwd=3)

 # Show the observations (and the prior central density)
if(PDF) pdf('gammaData7.pdf',width=4.0,height=1.9)
par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(4,1,1,1))
plot(xmidpoint,wmidpoint,xlim=c(0,600000),ylim=c(0,3400),xlab='household income',
     ylab='count per $5000',axes=F,type='n')
x = seq(0,7e5,length=200)
#lines(x,sum(counts)*5000*dgamma(x,shape=alpha,rate=beta),lty=1,col='green')
polygon(x,sum(counts)*5000*dgamma(x,shape=alpha,rate=beta),col='grey80',border=NA)
axis(1)
points(xsamp,0*xsamp+115,pch='|',lwd=2)
if(PDF) dev.off()


