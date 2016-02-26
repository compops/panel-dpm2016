setwd("~/projects/dpm-panel2015/src-draft1/example2-mixedeffects")

set.seed( 87655678 )
library("mvtnorm")

# Setup plot colors
library("RColorBrewer")
plotColors = brewer.pal(8, "Dark2");

source("samplers-mixedeffects-1dim.R")

nAttributes   <- 1
nMaxClusters  <- 20
nIndividuals  <- 50
nObservations <- 20
nIter         <- 10000
nBurnIn       <- 2500
postGrid      <- seq(-4,4,0.1)

###############################################################################
# Generate data
###############################################################################
d <- generateData( beta=c( 2, 3, -2), Q=c( 1, 0.2, 0.2), w=c(0.4, 0.3, 0.3),
                   s2e=1, alpha=2, omega=abs(rnorm(nIndividuals)),
                   dims=c(nIndividuals,nAttributes,nObservations))

###############################################################################
# Define priors
###############################################################################
a0star <- rep(0,   1 + nAttributes * nMaxClusters)
A0star <- 4 * diag(1 + nAttributes * nMaxClusters)

c0Q    <- 2
C0Q    <- 1

c0e    <- 1
C0e    <- 1

a0     <- 0.1
nu     <- 0.01

prior = list( a0star=a0star, A0star=A0star, c0Q = c0Q, C0Q = C0Q,
              c0e=c0e, C0e=C0e, a0=a0, nu=nu)


###############################################################################
# Run the finite model with known number of clusters
###############################################################################

prior$a0        <- 1
prior$a0star    <- rep(0,   1 + nAttributes * 3)
prior$A0star    <- 4 * diag(1 + nAttributes * 3)
outFiniteKnown  <- gibbs_finitemixture(d$y, d$x, nIter, 3, prior, postGrid)

###############################################################################
# Run the finite model with sparseness prior
###############################################################################

prior$a0        <- 1/nMaxClusters
prior$a0star    <- rep(0,   1 + nAttributes * nMaxClusters)
prior$A0star    <- 4 * diag(1 + nAttributes * nMaxClusters)
outFiniteSparse <- gibbs_finitemixture(d$y, d$x, nIter, nMaxClusters, 
                                       prior, postGrid)

###############################################################################
# Run the DPM model
###############################################################################

prior$alpha0 <- 0.01
prior$beta0  <- 0.01

source("samplers-mixedeffects-1dim.R")
outDPM       <- gibbs_dpm(d$y, d$x, nIter, nMaxClusters, prior, postGrid)

###############################################################################
# Plotting
###############################################################################

layout(matrix(c(1,1,1,2,2,3,4,4,5,6,7,8), 4, 3, byrow = TRUE))  
par(mar=c(4,5,1,1)) 

plot(postGrid, colMeans(outFiniteKnown$posterior_beta[nBurnIn:nIter,]), lwd=3, col=plotColors[1], type="l", bty="n", ylab="density", xlab=expression(beta), ylim=c(0,1),xlim=c(-3,3))
lines(postGrid, colMeans(outFiniteSparse$posterior_beta[nBurnIn:nIter,]), lwd=3, col=plotColors[2])
lines(postGrid, colMeans(outDPM$posterior_beta[nBurnIn:nIter,]), lwd=3, col=plotColors[3])
legend("topleft",c("M known","Sparseness prior","DP prior"),col=plotColors[1:3],lwd=2,box.col=NA)

hist(outFiniteKnown$alpha[nBurnIn:nIter],main="",breaks=floor(sqrt(nIter)),freq=FALSE,xlab=expression(alpha),col=rgb(t(col2rgb(plotColors[4]))/256,alpha=0.25),border=NA, ylab="density",xlim=c(1.4,2.2))
lines(density(outFiniteKnown$alpha[nBurnIn:nIter]),lwd=3,col=plotColors[4])
hist(outFiniteKnown$nOccupiedClusters[nBurnIn:nIter],breaks=floor(sqrt(nIter)),main="",freq=FALSE,xlab="no. occupied clusters",col="darkgrey",border=NA,xlim=c(0,nMaxClusters), ylab="density")

hist(outFiniteSparse$alpha[nBurnIn:nIter],main="",breaks=floor(sqrt(nIter)),freq=FALSE,xlab=expression(alpha),col=rgb(t(col2rgb(plotColors[5]))/256,alpha=0.25),border=NA, ylab="density",xlim=c(1.4,2.2))
lines(density(outFiniteSparse$alpha[nBurnIn:nIter]),lwd=3,col=plotColors[5])
hist(outFiniteSparse$nOccupiedClusters[nBurnIn:nIter],breaks=floor(sqrt(nIter)),main="",freq=FALSE,xlab="no. occupied clusters",col="darkgrey",border=NA,xlim=c(0,nMaxClusters), ylab="density")

hist(outDPM$alpha[nBurnIn:nIter],main="",breaks=floor(sqrt(nIter)),freq=FALSE,xlab=expression(alpha),col=rgb(t(col2rgb(plotColors[6]))/256,alpha=0.25),border=NA, ylab="density",xlim=c(1.4,2.2))
lines(density(outDPM$alpha[nBurnIn:nIter]),lwd=3,col=plotColors[6])
hist(outDPM$alphaSB[nBurnIn:nIter],main="",breaks=floor(sqrt(nIter)),freq=FALSE,xlab=expression(alpha[DP]),col=rgb(t(col2rgb(plotColors[7]))/256,alpha=0.25),border=NA, ylab="density",xlim=c(0,12))
lines(density(outDPM$alphaSB[nBurnIn:nIter]),lwd=3,col=plotColors[7])
hist(outDPM$nOccupiedClusters[nBurnIn:nIter],breaks=floor(sqrt(nIter)),main="",freq=FALSE,xlab="no. occupied clusters",col="darkgrey",border=NA,xlim=c(0,nMaxClusters), ylab="density")