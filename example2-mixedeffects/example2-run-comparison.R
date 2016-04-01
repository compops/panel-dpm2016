###############################################################################
# Script for replicating example 2 in 
# "Bayesian inference for mixed effects models with heterogeneity"
#
# Script to run Gibbs sampling for a single data set
#
# (c) Johan Dahlin 2016 ( johan.dahlin (at) liu.se )
###############################################################################

#setwd("~/projects/dpm-panel2015/src-draft1/example2-mixedeffects")

set.seed( 87655678 )
library("mvtnorm")

# Setup plot colors
library("RColorBrewer")
plotColors = brewer.pal(8, "Dark2");

source("samplers-mixedeffects-1dim.R")

nAttributes   <- 1
nMaxClusters  <- 20
nIndividuals  <- 20
nObservations <- 20
nIter         <- 1000
nBurnIn       <- 250
postGrid      <- seq(-4,4,0.1)

###############################################################################
# Generate data
###############################################################################
d <- generateData( beta=c( 2, 3, -2), Q=c( 1, 0.2, 0.2), w=c(0.4, 0.3, 0.3),
                   s2e=1, alpha=2, omega=abs(rnorm(nIndividuals)),
                   dims=c(nIndividuals,nAttributes,nObservations))

###############################################################################
# Analysis for each individual
###############################################################################

resLM <- matrix(0, nrow=nIndividuals, ncol=3)

for (ii in 1:nIndividuals) {
  res           <- lm( d$y[ii,] ~ d$x[,ii,] )
  resLM[ii,1:2] <- res$coefficients
  resLM[ii,3]   <- var( res$residuals )
}

resLMmeans <- apply(resLM,2,mean)
resLMvars  <- apply(resLM,2,var)

###############################################################################
# Define priors
###############################################################################
a0star <- c(resLMmeans[1], rep(resLMmeans[2], nMaxClusters))
A0star <- diag(rep(1, 1+nAttributes*nMaxClusters))

c0Q    <- 10
C0Q    <- 1
c(C0Q/(c0Q+1),sqrt(C0Q^2*(c0Q-1)^(-2)*(c0Q-2)^(-1)),(C0Q/(1+c0Q)))

c0e    <- 0
C0e    <- 0

a0     <- 0.1
nu     <- 5

prior = list( a0star=a0star, A0star=A0star, c0Q = c0Q, C0Q = C0Q,
              c0e=c0e, C0e=C0e, a0=a0, nu=nu)


###############################################################################
# Run the finite model with known number of clusters
###############################################################################

prior$a0star    <- c(resLMmeans[1], rep(resLMmeans[2], 3))
prior$A0star    <- diag(rep(1, 1+3))
prior$eta0      <- 1
outFiniteKnown  <- gibbs_finitemixture(d$y, d$x, nIter, 3, prior, postGrid)

###############################################################################
# Run the finite model with sparseness prior
###############################################################################
prior$a0star    <- c(resLMmeans[1], rep(resLMmeans[2], nMaxClusters))
prior$A0star    <- diag(rep(1, 1+nAttributes*nMaxClusters))
prior$eta0      <- 1/nMaxClusters

outFiniteSparse <- gibbs_finitemixture(d$y, d$x, nIter, nMaxClusters, 
                                       prior, postGrid)

###############################################################################
# Run the DPM model
###############################################################################

prior$alpha0 <- 0.01
prior$beta0  <- 0.01

outDPM       <- gibbs_dpm(d$y, d$x, nIter, nMaxClusters, prior, postGrid)

###############################################################################
# Plotting
###############################################################################

layout(matrix(c(1,1,1,2,3,4,5,6,7), 3, 3, byrow = TRUE))  
par(mar=c(4,5,1,1)) 

plot(postGrid, colMeans(outFiniteKnown$posterior_beta[nBurnIn:nIter,]), lwd=3, col=plotColors[1], type="l", bty="n", ylab="density", xlab=expression(beta), ylim=c(0,2),xlim=c(-4,4))
lines(postGrid, colMeans(outFiniteSparse$posterior_beta[nBurnIn:nIter,]), lwd=3, col=plotColors[2])
lines(postGrid, colMeans(outDPM$posterior_beta[nBurnIn:nIter,]), lwd=3, col=plotColors[3])
legend("topleft",c("M known","Sparseness prior","DP prior"),col=plotColors[1:3],lwd=2,box.col=NA)

hist(outFiniteKnown$alpha[nBurnIn:nIter],main="",breaks=floor(sqrt(nIter)),freq=FALSE,xlab=expression(alpha),col=rgb(t(col2rgb(plotColors[1]))/256,alpha=0.25),border=NA, ylab="density",xlim=c(1.7,2.2))
lines(density(outFiniteKnown$alpha[nBurnIn:nIter],from=1.7,to=2.2),lwd=3,col=plotColors[1])

hist(outFiniteSparse$alpha[nBurnIn:nIter],main="",breaks=floor(sqrt(nIter)),freq=FALSE,xlab=expression(alpha),col=rgb(t(col2rgb(plotColors[2]))/256,alpha=0.25),border=NA, ylab="density",xlim=c(1.7,2.2))
lines(density(outFiniteSparse$alpha[nBurnIn:nIter],from=1.7,to=2.2),lwd=3,col=plotColors[2])

hist(outDPM$alpha[nBurnIn:nIter],main="",breaks=floor(sqrt(nIter)),freq=FALSE,xlab=expression(alpha),col=rgb(t(col2rgb(plotColors[3]))/256,alpha=0.25),border=NA, ylab="density",xlim=c(1.7,2.2))
lines(density(outDPM$alpha[nBurnIn:nIter],from=1.7,to=2.2),lwd=3,col=plotColors[3])

hist(outFiniteSparse$nOccupiedClusters[nBurnIn:nIter],breaks=floor(sqrt(nIter)),main="",freq=FALSE,xlab="no. occupied clusters",col=plotColors[2],border=plotColors[2],xlim=c(0,nMaxClusters), ylab="density")
hist(outDPM$nOccupiedClusters[nBurnIn:nIter],breaks=floor(sqrt(nIter)),main="",freq=FALSE,xlab="no. occupied clusters",col=plotColors[3],border=plotColors[3],xlim=c(0,nMaxClusters), ylab="density")

hist(outDPM$alphaSB[nBurnIn:nIter],main="",breaks=floor(sqrt(nIter)),freq=FALSE,xlab=expression(alpha[DP]),col=rgb(t(col2rgb(plotColors[3]))/256,alpha=0.25),border=NA, ylab="density",xlim=c(0,60))
lines(density(outDPM$alphaSB[nBurnIn:nIter],from=0,to=60),lwd=3,col=plotColors[3])

###############################################################################
# End of file
###############################################################################
