setwd("~/projects/dpm-panel2015/src-draft1/example2-armodel")
set.seed( 87655678 )

# Setup plot colors
library("RColorBrewer")
library("MASS")
plotColors = brewer.pal(6, "Dark2");

###############################################################################
# Generate data from AR processes with different means
###############################################################################

source("samplers-armodel.R")

# Settings
nData         <- 50
nClusters     <- 3
nObservations <- 20

p     <- c(0.2, 0.6, 0.2)
mu    <- c(-2, 1, 3)
rho    <- 0.5
sigmae <- 0.5

# Pre-allocate data and indicators
y     <- matrix(0, nrow=nData, ncol=nObservations)
x     <- matrix(0, nrow=nData, ncol=1)
gamma <- matrix(0, nrow=nData, ncol=1)

# Generate data
for (ii in 1:nData) {
  
  x[ii]     <- sample(nClusters, 1, replace=TRUE, p[1:nClusters]/sum(p[1:nClusters]))
  y[ii,1]   <- rnorm(1);
  
  gamma[ii] <- mu[x[ii]]
  
  for (tt in 2:nObservations) {
    y[ii,tt] <- gamma[ii] + rho * y[ii,tt-1] + sigmae * rnorm(1)
  }
}



# # Build the vector of observations and the regressor matrix
# ySub     <- matrix(y[,-1],             nrow=nData * (nObservations - 1), ncol=1, byrow=TRUE)
# XSub     <- matrix(1,                  nrow=nData * (nObservations - 1), ncol=2)
# XSub[,2] <- matrix(y[,-nObservations], nrow=nData * (nObservations - 1), ncol=1, byrow=TRUE)
# 
# # Compute the posterior parameters
# m0          = rep(0,2);
# V0          = diag( rep(1,2) );
# VN = solve( t(XSub) %*% XSub + V0 )
# mN = VN %*% ( t(XSub) %*% ySub + solve(V0) %*% m0 )


###############################################################################
# MCMC
###############################################################################

nMCMC           <- 3000
nBurnIn         <- 500
nMaxClusters    <- 3

###############################################################################
# Run the finite model with known number of clusters
###############################################################################
prior           <- list(m0=c(0,0), V0=diag( c(2^2,0.4^2) ), a0=0.1, b0=0.1, a0=10)
grids           <- list(gamma=seq(-2, 2, 0.01), rho=seq(0, 1.5, 0.01), sigma=seq(0.01,1,0.01))
outFiniteKnown  <- gibbs_finitemixture(y, nMCMC, nBurnIn, nMaxClusters, prior, grids)

layout(matrix(c(1,2,3,4,4,5,6,6,6), 3, 3, byrow = TRUE))  
par(mar=c(4,5,1,1)) 

# Plot the parameter posteriors
plot( grids$gamma, outFiniteKnown$postGamma, type="l", col=plotColors[3], lwd=2, ylab="posterior density estimate", xlab=expression(gamma), bty="n" )
plot( grids$rho,   outFiniteKnown$postRho, type="l", col=plotColors[4], lwd=2, ylab="posterior density estimate", xlab=expression(rho), bty="n" )
plot( grids$sigma, outFiniteKnown$postSigma, type="l", col=plotColors[4], lwd=2, ylab="posterior density estimate", xlab=expression(sigma), bty="n" )

# Plot the number of active clusters and its distribution
plot(outFiniteKnown$nOccupiedClusters,type="l",col=plotColors[6],xlab="iteration",ylab="no. active clusters",ylim=c(0,nMaxClusters), bty="n")
hist(outFiniteKnown$nOccupiedClusters[nBurnIn:nMCMC,1],breaks=nMaxClusters,main="",freq=F,col=plotColors[6],border=NA,xlim=c(0,nMaxClusters+1))

# Plot the likelihood
plot(outFiniteKnown$ll,type="l",col=plotColors[1],xlab="iteration",ylab="log-likelihood", bty="n")
