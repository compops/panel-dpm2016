###############################################################################
# Script for replicating example 1 in 
# "Bayesian inference for mixed effects models with heterogeneity"
#
# (c) Johan Dahlin 2016 ( johan.dahlin (at) liu.se )
###############################################################################

#setwd("~/storage/owncloud/projects-work/dpm-panel2015/src-draft1/example1-mixture-normals")

set.seed( 87655678 )

# Setup plot colors
library("RColorBrewer")
plotColors = brewer.pal(6, "Dark2");

###############################################################################
# Generate data from a Gaussian mixture with three components
###############################################################################

source("samplers-normalmixture.R")

# Settings
nData     <- 100
nClusters <- 3

p     <- c(0.3, 0.5, 0.2)
mu    <- c(-1, 0, 2)
sigma <- c(0.2, 1, 0.05)

# Pre-allocate data and indicators
y <- matrix(0, nrow=nData, ncol=1)
x <- matrix(0, nrow=nData, ncol=1)

# Generate data
for (ii in 1:nData) {
  x[ii] <- sample(nClusters, 1, replace=TRUE, p)
  y[ii] <- rnorm(1, mean=mu[x[ii]], sd=sigma[x[ii]])
}

###############################################################################
# MCMC
###############################################################################

grid          <- seq( -4, 4, 0.01 )
nMCMC         <- 10000
nBurnIn       <- 2500
nMaxClusters  <- 20
prior         <- list( mu0=0, sd0=0.2, kappa0=1, nu0=1, a0=1)

###############################################################################
# Run the finite model with known number of clusters
###############################################################################

outFiniteKnown <- gibbs_finitemixture(y, nMCMC, nClusters, prior, grid)

###############################################################################
# Run the finite model with sparseness prior
###############################################################################

prior$a0        <- 1 / nMaxClusters
outFiniteSparse <- gibbs_finitemixture(y, nMCMC, nMaxClusters, prior, grid)

###############################################################################
# Run the DPM model
###############################################################################

prior$alpha0    <- 1
prior$beta0     <- 0.5
outDPM          <- gibbs_dpm(y, nMCMC, nMaxClusters, prior, grid)

###############################################################################
# Plot the output
###############################################################################

#cairo_pdf("~/projects/dpm-panel2015/paper/dpm-panel2015-draft1/figures/example1-comparison.pdf", width=10, height=8)

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE)) 
par(mar=c(4,5,1,1))

plot(grid, colMeans(outFiniteKnown$gridPost[nBurnIn:nMCMC,]), lwd=3, col=plotColors[1], type="l", bty="n", ylab="density", xlab=expression(beta[i]^{s}), ylim=c(0,0.75),xlim=c(-3,4))
lines(grid, colMeans(outFiniteSparse$gridPost[nBurnIn:nMCMC,]), lwd=3, col=plotColors[2])
lines(grid, colMeans(outDPM$gridPost[nBurnIn:nMCMC,]), lwd=3, col=plotColors[3])
rug(y)

hist(outFiniteSparse$nOccupiedClusters[nBurnIn:nMCMC], breaks=floor(sqrt(nMCMC-nBurnIn))/3, freq=FALSE, col=plotColors[2], xlab="no. occupied components", ylab="density", border=NA, main="", xlim=c(0,20), ylim=c(0,1.0))
abline(v=nClusters-0.5, lty="dotted")
text(20,1.0,"FMs",pos=2)

hist(outDPM$nOccupiedClusters[nBurnIn:nMCMC], breaks=floor(sqrt(nMCMC-nBurnIn))/2, freq=FALSE, col=plotColors[3], xlab="no. occupied components", ylab="density", border=NA, main="", xlim=c(0,20), ylim=c(0,0.4))
abline(v=nClusters-0.5, lty="dotted")
text(20,0.4,"DPM model",pos=2)

#dev.off()

###############################################################################
# End of file
###############################################################################