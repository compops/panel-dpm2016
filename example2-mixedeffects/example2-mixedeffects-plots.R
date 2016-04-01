###############################################################################
# Script for replicating example 2 in 
# "Bayesian inference for mixed effects models with heterogeneity"
#
# Script to make plots of the heterogenity after a run of
# example2-run-gamma.R
#
# (c) Johan Dahlin 2016 ( johan.dahlin (at) liu.se )
###############################################################################

#setwd("~/projects/dpm-panel2016/results-draft1/example2-mixedeffects")

# Settings for the MCMC sampler
nIter   <- 10000
nBurnIn <- 2500
nIndividualsVector <- c( 10, 20, 50, 100, 200, 500)

# Pre-allocate vector
nIndividuals  <- nIndividualsVector[6]
dFiniteKnown  <- array(0, dim=c(6,40,512))
dFiniteSparse <- array(0, dim=c(6,40,512))
dDPM          <- array(0, dim=c(6,40,512))

###############################################################################
# Load data
###############################################################################

# Loop over the six first data sets
for ( jj in 1:6 ) {
  
  # Loop over all the elements in nIndividualVector
  for ( ii in 1:length(nIndividualsVector) ) {
    nIndividuals  <- nIndividualsVector[ii]
  
    # Load data
    load(paste(paste(paste(paste("example2-mixedeffects/mixedeffects-nInd-",nIndividuals,sep=""),"-",sep=""),jj,sep=""),".RData",sep=""))
    
    # Compute estimates of the posteriors via density estimates
    dFiniteSparse[ii,jj,]    <- density(outFiniteSparse$betaInd[nBurnIn:nIter], kernel="e", from=-5, to=5)$y
    dDPM[ii,jj,]             <- density(outDPM$betaInd[nBurnIn:nIter],          kernel="e", from=-5, to=5)$y
    
  }
}

###############################################################################
# Plot data
###############################################################################

#cairo_pdf("~/projects/dpm-panel2015/paper/dpm-panel2015-draft1/figures/example2-comparison.pdf", width=10, height=14)

# Setup plot colors
library("RColorBrewer")
plotColors = brewer.pal(6, "Dark2");

grid <- seq(-5,5,length.out = 512)

layout(matrix(1:36, 6, 4, byrow = FALSE)) 
par(mar=c(4,4,0,0))

jj <- 1
for ( ii in 1:(length(nIndividualsVector)-1) ) {
  plot(grid,dFiniteSparse[ii,jj,],type="l",col=plotColors[2],xlab="",ylab="density", bty="n",lwd=2) 
  lines(grid,dDPM[ii,jj,],col=plotColors[3],lwd=2)
}

plot(grid,dFiniteSparse[ii,6,],type="l",col=plotColors[2],xlab=expression(beta),ylab="density", bty="n",lwd=2) 
lines(grid,dDPM[ii,6,],col=plotColors[3],lwd=2)

for ( jj in 2:4 ) {
  for ( ii in 1:length(nIndividualsVector) ) {
    if ( ii == 6 ){
      plot(grid,dFiniteSparse[ii,6,],type="l",col=plotColors[2],xlab=expression(beta),ylab="", bty="n",lwd=2) 
      lines(grid,dDPM[ii,6,],col=plotColors[3],lwd=2)
    } else {
      plot(grid,dFiniteSparse[ii,jj,],type="l",col=plotColors[2],xlab="",ylab="", bty="n",lwd=2) 
      lines(grid,dDPM[ii,jj,],col=plotColors[3],lwd=2)
    }
  }
}

#dev.off()

###############################################################################
# End of file
###############################################################################