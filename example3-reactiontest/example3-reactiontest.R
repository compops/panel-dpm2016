setwd("~/projects/dpm-panel2015/src-draft1/example3-reactiontest")

set.seed( 87655678 )
library("mvtnorm")

# Setup plot colors
library("RColorBrewer")
plotColors = brewer.pal(8, "Dark2");

source("samplers-mixedeffects-sleepmodel.R")

nIter         <- 1000
nBurnIn       <- 250


###############################################################################
# Load data
###############################################################################
library("lme4")
str(sleepstudy)

nAttributes   <- 2
nIndividuals  <- 18
nMaxClusters  <- nIndividuals
nObservations <- 10

y      <- matrix( sleepstudy$Reaction, nrow=18, ncol=10, byrow = TRUE)
x      <- array( 0, dim=c(nAttributes,nIndividuals,nObservations) )
x[1,,] <- matrix( 1,               nrow=18, ncol=10, byrow = TRUE)
x[2,,] <- matrix( sleepstudy$Days, nrow=18, ncol=10, byrow = TRUE)
d <- list( x=x, y=y)

###############################################################################
# Analysis for each individual
###############################################################################

resLM <- matrix(0, nrow=nIndividuals, ncol=3)

for (ii in 1:nIndividuals) {
  res           <- lm( y[ii,] ~ x[2,ii,] )
  resLM[ii,1:2] <- res$coefficients
  resLM[ii,3]   <- var( res$residuals )
}

resLMmeans <- apply(resLM,2,mean)
resLMvars  <- apply(resLM,2,var)

###############################################################################
# Define priors
###############################################################################
a0star <- c(resLMmeans[1], resLMmeans[2], rep(0, 2*nMaxClusters))
A0star <- diag( c(0.04, 0.04, 0.04*rep(1, nAttributes*nMaxClusters ) ) )

# Wishart prior on inverse covariance
c0Q    <- 10
C0Q    <- 0.05 * diag(2)
#C0Q    <- 0.5 * diag(2)
( c0Q - (2+1)/2 )^(-1) * C0Q # mode of variance matrix

c0e    <- 0
C0e    <- 0

a0     <- 0.1
nu     <- 5

prior = list( a0star=a0star, A0star=A0star, c0Q = c0Q, C0Q = C0Q,
              c0e=c0e, C0e=C0e, a0=a0, nu=nu)


###############################################################################
# Run the finite model with sparseness prior
###############################################################################

prior$eta0      <- 1/nMaxClusters
outFiniteSparse <- gibbs_finitemixture(d$y, d$x, d$x, nIter, nMaxClusters, 
                                       prior, postGrid)

###############################################################################
# Run the DPM model
###############################################################################

prior$alpha0 <- 0.01
prior$beta0  <- 0.01

outDPM       <- gibbs_dpm(d$y, d$x, d$x, nIter, nMaxClusters, prior, postGrid)

###############################################################################
# Plotting
###############################################################################

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))  
par(mar=c(4,5,1,1)) 

plot(density(outFiniteSparse$alpha[nBurnIn:nIter,1]),lwd=3, col=plotColors[1], type="l", bty="n", ylab="density", xlab=expression(alpha[0]), main="" )
lines(density(outDPM$alpha[nBurnIn:nIter,1]),lwd=3, col=plotColors[2])

plot(density(outFiniteSparse$alpha[nBurnIn:nIter,2]),lwd=3, col=plotColors[1], type="l", bty="n", ylab="density", xlab=expression(alpha[1]), main="" )
lines(density(outDPM$alpha[nBurnIn:nIter,2]),lwd=3, col=plotColors[2])

plot(density(outFiniteSparse$betaInd[nBurnIn:nIter,,1]),lwd=3, col=plotColors[1], type="l", bty="n", ylab="density", xlab=expression(beta[0]), main="" )
lines(density(outDPM$betaInd[nBurnIn:nIter,1]),lwd=3, col=plotColors[2])
plot(density(outFiniteSparse$betaInd[nBurnIn:nIter,,2]),lwd=3, col=plotColors[1], type="l", bty="n", ylab="density", xlab=expression(beta[1]), main="" )
lines(density(outDPM$betaInd[nBurnIn:nIter,2]),lwd=3, col=plotColors[2])



# lines(postGrid, colMeans(outDPM$posterior_beta[nBurnIn:nIter,]), lwd=3, col=plotColors[3])
# legend("topleft",c("Sparseness prior","DP prior"),col=plotColors[2:3],lwd=2,box.col=NA)
# rug(resLM[,2])
# 
# hist(outFiniteSparse$alpha[nBurnIn:nIter],main="",breaks=floor(sqrt(nIter)),freq=FALSE,xlab=expression(alpha),col=rgb(t(col2rgb(plotColors[5]))/256,alpha=0.25),border=NA, ylab="density",xlim=c(250,253))
# rug(resLM[,1])
# lines(density(outFiniteSparse$alpha[nBurnIn:nIter]),lwd=3,col=plotColors[5])
# hist(outFiniteSparse$nOccupiedClusters[nBurnIn:nIter],breaks=floor(sqrt(nIter)),main="",freq=FALSE,xlab="no. occupied clusters",col="darkgrey",border=NA,xlim=c(0,nMaxClusters), ylab="density")
# 
# hist(outDPM$alpha[nBurnIn:nIter],main="",breaks=floor(sqrt(nIter)),freq=FALSE,xlab=expression(alpha),col=rgb(t(col2rgb(plotColors[6]))/256,alpha=0.25),border=NA, ylab="density",xlim=c(250,253))
# rug(resLM[,1])
# lines(density(outDPM$alpha[nBurnIn:nIter]),lwd=3,col=plotColors[6])
# hist(outDPM$alphaSB[nBurnIn:nIter],main="",breaks=floor(sqrt(nIter)),freq=FALSE,xlab=expression(eta[0]),col=rgb(t(col2rgb(plotColors[7]))/256,alpha=0.25),border=NA, ylab="density",xlim=c(0,30))
# lines(density(outDPM$alphaSB[nBurnIn:nIter]),lwd=3,col=plotColors[7])
# hist(outDPM$nOccupiedClusters[nBurnIn:nIter],breaks=floor(sqrt(nIter)),main="",freq=FALSE,xlab="no. occupied clusters",col="darkgrey",border=NA,xlim=c(0,nMaxClusters), ylab="density")