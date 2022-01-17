
# Distribution of spacings of orderstatistics
# u_(1),...,u_(n) from n srs draws from the population Y={1,...,N}
# distribution.
# Assume Y = 1:N for this problem

# use N = 20 for now
n = 10
n = 3
N = 20
Y = 1:N

 # generate n iid draws
u <- sample(Y,size=n,replace=FALSE)
 # get the spacings
deltas <- diff(c(0,sort(u),N))


 # make a function giving P(X(1))
pminSrsDisc <- function(N,n){
  # computes the density over the outcome space {1,...,N} of the 1st
  # order statistic from a srs(n) from {1,...,N}
  i1 = 1:(N-n+1); i0 = i1-1
   # compute P(X(1) == i1)
  px1 = 1 - choose(N-i1,n)/choose(N,n)

  return(diff(c(0,px1)))
}
 # try it out
pminSrsDisc(20,3)
sum(pminSrsDisc(20,3))
 # check via MC - it looks good
NMC = 1000000; u <- matrix(NA,nrow=3,ncol=NMC)
for(k in 1:NMC) {
  u[, k] <- sample(1:20, size = 3, replace = FALSE)
}

xmin = apply(u,2,min)
#xmed = apply(u,2,median)
round(tabulate(xmin)/NMC,4)
round(pminSrsDisc(20,3),4)
 # if x(1) were 3, then p(x(2)) is given by
pminSrsDisc(20-3,3-1)
pminSrsDisc(20-3,3-1) > 0

sum(pminSrsDisc(20-3,3-1))

 # getting the counts of individuals between order stats:
 # left.open means count[1] is # of x <= vec[1]
findInterval(1:40,c(10,20),left.open=TRUE)
table(findInterval(1:40,c(10,20),left.open=TRUE))
tabulate(findInterval(1:40,c(10,20),left.open=TRUE)+1)

 # now compute the probability of the counts of the population
 # binned by the ordered iid sample

 # first, build a function that uses the overkill pminDisc() function
likeSampSrs <- function(counts){
  # computes the likelihood of counts between order statistics from
  # a sample of size n=length(counts)-1.
  # counts[1] <= x_(1), ... , counts[n] <= x_(n), counts[n+1] > x_(n).
  N = sum(counts)
  cumcounts = c(0,cumsum(counts))
  adjcounts = counts; adjcounts[1] = counts[1]-1
  n = length(counts)-1
  probs = rep(NA,n)
  #browser()
  for(k in 1:n){
    p1 = pminSrsDisc(N-cumcounts[k],n-k+1)
    probs[k] = p1[counts[k]]
  }
  print(probs)
  return(prod(probs))
}
# likeSampSrs(c(10,10,20))

 # try with previous problem
pminSrsDisc(20,3)
# if x(1) were 4, then p(x(2)) is given by
pminSrsDisc(20-4,3-1)
# if x(2) were 10 the p(x(3)|x(2),x(1)) is given by
pminSrsDisc(20-10,3-2)

 # now using the pminLike
tabulate(findInterval(1:20,c(4,10,14),left.open=TRUE)+1)
likeSampSrs(c(4,6,4,6))

# make all possible ordered samples for n=3, N=20 and
# compute the probabilities of each
n=3; N=20
x1=rep(NA,(N-2)*(N-1)*N)
x2=x1; x3=x1; prob20=rep(NA,1140)
countmat = matrix(NA,ncol=1140,nrow=4)
xmat = matrix(NA,ncol=1140,nrow=3)
ii = 1
for(i in 1:(N-2)){
  for(j in (i+1):(N-1)){
    for(k in (j+1):N){
      xmat[,ii] = c(i,j,k)
      print(xmat[,ii])
      countmat[,ii] = tabulate(findInterval(1:N,c(i,j,k),left.open=TRUE)+1,nbins=n+1)
      prob20[ii] = likeSampSrs(countmat[,ii])
      ii = ii+1
    }
  }
}

hist(prob20)
ibad=is.na(prob20)
sum(prob20,na.rm=T)  # seems not obviously wrong

 # look at the margins
# make a dplyer data frame:
xdat = data.frame(x1=xmat[1,],x2=xmat[2,],x3=xmat[3,],prob=prob20)
xdat = xdat[!ibad,]
library(dplyr)
xdat %>% group_by(x1) %>% summarize(prob=sum(prob)) -> mx1
xdat %>% group_by(x2) %>% summarize(prob=sum(prob)) -> mx2
xdat %>% group_by(x3) %>% summarize(prob=sum(prob)) -> mx3

par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(4,4,2,1))
plot(c(1,20),range(cbind(mx1$prob,mx2$prob,mx3$prob)),type='n',
     xlab='x_(k)',ylab='probability')
points(mx1$x1,mx1$prob,col='black')
points(mx2$x2,mx2$prob,col='blue')
points(mx3$x3,mx3$prob,col='red')
mtext('marginal distribuutions of srs order stats: N=20, n=3',line=.5,
      side=3,cex=.9)
legend('top',pch=1,col=c('black','blue','red'),
       legend = c('x_(1)','x_(2)','x_(3)'),cex=.8)


# looks great, but the probabilities of each draw are all the same - just like
# Chris was saying, and like theory predicts.

# Maybe we should be looking at the distribution of the empirical:
# sqrt{n}[F_n(F^-1(u)) - F(F^-1(u))] ~ Bbridg(u) approx?
# This distribution will not be indifferent to a iid sample of
# {1,2,3} from pop = {1,...,N}, right?

# If this really is different, what is the difference between
# the distribution of the empirical dist Fn and the order stats?



