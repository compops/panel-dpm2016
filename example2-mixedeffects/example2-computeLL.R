###############################################################################
# Script for replicating example 2 in 
# "Bayesian inference for mixed effects models with heterogeneity"
#
# Script to compute the IACT of the log-likelihood and cluster allocation
# after a run of example2-run-gamma.R
#
# (c) Johan Dahlin 2016 ( johan.dahlin (at) liu.se )
###############################################################################

###############################################################################
# Compute the log-likelihood for each run and iteration
###############################################################################

calcLikelihood <- function( d, out, nIter, nBurnIn ) {
  ll <- matrix(0, nrow=nIter-nBurnIn, ncol=1)
  
  for ( kk in 1:length(ll) ) {
    for ( ii in 1:nIndividuals ) {
      jj <- kk + nBurnIn
      ll[kk] <- ll[kk] + sum( dnorm(d$y[ii,], out$alpha[jj+1] + out$betaInd[jj+1,ii] * d$x[1,ii,], sqrt(out$var[jj+1,ii]), log=TRUE ) )
    }
  }
  
  ll
}

###############################################################################
# Computes the IACT values
# Code taken from the LaplacesDemon package available at:
# https://github.com/ecbrown/LaplacesDemon/blob/master/R/IAT.R
###############################################################################

iact <- function(x,xt,lag=0) {
  x = x - xt;
  xcorr = acf(x, lag.max = floor( length(x)/2 ), plot=F )
  lim   = 2/sqrt(length(x));
  
  if ( sum(is.nan(xcorr$acf)) > 0 ) {
    out = NA;
  } else {
    if (lag==0) {
      newli = min(which(xcorr$acf < lim)[1],floor( length(x) / 10),na.rm=T);
    } else {
      newli = lag
    }
    out = 1 + 2 * sum( xcorr$acf[1:newli] );
  }
  S = length(x);
  return(out);
}

###############################################################################
# Compute IACT of the log-likelihood
###############################################################################

calcIACTLL <- function(d, out, nIter, nBurnIn){
  ll <- calcLikelihood( d, out, nIter, nBurnIn )
  iact(ll, xt=mean(ll))
}

###############################################################################
# Main function
###############################################################################

#setwd("~/projects/dpm-panel2016/results-draft1/example2-mixedeffects")

# Settings for MCMC sampler and model
nIter   <- 10000
nBurnIn <- 2500
nIndividualsVector <- c( 10, 20, 50, 100, 200, 500)

IACTs <- array(0, dim=c(length(nIndividualsVector),40,3))

###############################################################################
# Loop over all 40 runs
###############################################################################
for ( jj in 1:40 ) {
  
  for ( ii in 1:length(nIndividualsVector) ) {
    nIndividuals  <- nIndividualsVector[ii]
    
    # Load output from runs and data
    load(paste(paste(paste(paste("example2-mixedeffects/mixedeffects-nInd-",nIndividuals,sep=""),"-",sep=""),jj,sep=""),".RData",sep=""))
    load(file=paste(paste(paste(paste("example2-data/mixedeffects-nInd-",nIndividuals,sep=""),"-",sep=""),jj,sep=""),".RData",sep=""))
    
    # Compute IACTs
    IACTs[ii,jj,1] <- calcIACTLL(d, outFiniteKnown,  nIter, nBurnIn)
    IACTs[ii,jj,2] <- calcIACTLL(d, outFiniteSparse, nIter, nBurnIn)
    IACTs[ii,jj,3] <- calcIACTLL(d, outDPM,          nIter, nBurnIn)
  }
  
  print(c(jj,40))
}

# Compute the mean
apply(IACTs,1,colMeans)

###############################################################################
# End of file
###############################################################################