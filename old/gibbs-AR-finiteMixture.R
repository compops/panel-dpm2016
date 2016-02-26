# Setup plot colors
library("RColorBrewer")
plotColors = brewer.pal(6, "Dark2");

nMCMC         <- 10000
nBurnIn       <- 2500
nData         <- 10
nObservations <- 100
nClusters     <- 3
nMaxClusters  <- 3

##############################################################################
# Generate data
##############################################################################

p    <- c(0.4, 0.3, 0.3)
par  <- c(0.5, 0.8, 0.2)
par1 <- c(-2,0,1)

y    <- matrix(0, nrow=nData, ncol=nObservations)
x    <- matrix(0, nrow=nData, ncol=1)

for (ii in 1:nData){
  x[ii]   <- sample(nClusters, 1, replace=TRUE, p)
  y[ii,1] <- rnorm(1);
  foo     <- par1[x[ii]]

  for (tt in 2:nObservations) {
    y[ii,tt] <- foo + par[2] * ( y[ii,tt-1] - foo ) + par[3] * rnorm(1)
  }
}

par(mfrow=c(2,5))
for (ii in 1:nData) {
  plot(y[ii,],type="l",ylim=c(-1.5,1.5))
}

##############################################################################
# Gibbs sampler
##############################################################################
a0 <- 1
# Pre-allocate matrices
I                 <- matrix(1, nrow=nMCMC+1, ncol=nData )
nMembersCluster   <- matrix(0, nrow=nMCMC,   ncol=nMaxClusters)
nOccupiedClusters <- matrix(0, nrow=nMCMC,   ncol=1)
piPost            <- matrix(1, nrow=nMCMC,   ncol=nMaxClusters)

mu                <- matrix(0, nrow=nMCMC, ncol=nMaxClusters)
rho               <- matrix(0, nrow=nMCMC, ncol=1)
sigma             <- matrix(1, nrow=nMCMC, ncol=1)

mu[1,]            <- rnorm(3)
sigma[1,]         <- abs(rnorm(1))
I[2,]             <- sample( nMaxClusters, nData, replace=T )

for ( ii in 2:nMCMC ) {

  #===========================================================================
  # Count the no. members in each cluster and sample the cluster weight
  #===========================================================================
  piCluster <- rep(0, nMaxClusters)
  
  for (jj in 1:nMaxClusters) {
    
    # Count the number of members in each cluster
    nMembersCluster[ii,jj] <- length(which(I[ii,] == jj))
    
    # Count the number of occupied clusters
    if (nMembersCluster[ii,jj] > 0) { 
      nOccupiedClusters[ii] <- nOccupiedClusters[ii] + 1 
    }
    
    # Sample probabilities p from Dirichlet (as normalised Gamma variables)
    piCluster[jj] <- rgamma(1, a0 + nMembersCluster[ii,jj], 1)
  
  }
  
  piPost[ii,] <- piCluster / sum(piCluster)
  
  #===========================================================================
  # Update the parameters (only mu is different for each cluster)
  #===========================================================================  

  #---------------------------------------------------------------------------
  # rho | mu, sigma (normal)
  #---------------------------------------------------------------------------
  idx     <- I[ii,]
  muC     <- mu[ii-1,idx] 
  rhobar  <- sum( ( y[,-1] - muC ) * ( y[,-nObservations] - muC ) ) / sum( ( y[,-1] - muC )^2 )
  rho[ii] <- rnorm(1, rhobar, sqrt( sigma[ii-1] / sum( ( y[,-1] - mu[ii-1,idx] )^2 ) ) )
  
  if (rho[ii] > 0.99){ print("Warning rho=1"); rho[ii]=0.2 }
  
  #---------------------------------------------------------------------------
  # mu | rho, sigma (normal)
  # for each cluster
  #---------------------------------------------------------------------------
  for (jj in 1:nMaxClusters) {
    
    idx        <- which(I[ii,] == jj)
    ysub       <- y[idx,]
    nTotalData <- length(idx) * (nObservations - 1)
    
    if (length(idx) > 2) {
      mubar     <- mean( ysub[,-1] - rho[ii] * ysub[,-nObservations] ) / (1 - rho[ii] )
      mu[ii,jj] <- rnorm(1, mubar, 
                         sqrt( sigma[ii-1] / (nTotalData * (1 - rho[ii] )^2 ) ) )
    } else {
      mu[ii,jj] <- 0.3 * rnorm(1)
    }
  }
  
  #---------------------------------------------------------------------------
  # sigma | mu, rho (inverse Gamma)
  #---------------------------------------------------------------------------
  idx       <- I[ii,]
  sigmahat  <- sd( y[,-1] - mu[ii,idx] - rho[ii-1] * y[,-nObservations] )
  
  aPost     <- 0.5 * ( nObservations * nData ) + 1
  bPost     <- 0.5 * sum(((y[,-1]-mu[ii,idx])-rho[ii]*(y[,-nObservations]-mu[ii,idx]))^2)
  sigma[ii] <- 1 / rgamma(1, aPost, bPost )

  #===========================================================================
  # Sample the indicators
  #===========================================================================
  piPost[ii,] <- piCluster / sum(piCluster)
  
  w <- matrix(0, nrow = nData, ncol = nMaxClusters)
  
  for (jj in 1:nMaxClusters) {
    tmp    <- dnorm( y[,-1], mu[ii,I[ii,jj]] + rho[ii] * (y[,-nObservations] - mu[ii,I[ii,jj]]), sqrt(sigma[ii]), log=TRUE )
    w[,jj] <- rowSums(tmp)
  }
  
  for ( jj in 1:nData ) {
    wmax  <- max( w[jj,] )
    v     <- exp( w[jj,] - wmax )
    
    I[ii+1,jj] <- sample(nMaxClusters, 1, replace=TRUE, v / sum(v) )
  }  
}

##############################################################################
# Plot the results
##############################################################################

layout(matrix(c(1,1,1,2,2,3,4,4,5,6,6,7), 4, 3, byrow = TRUE))  
par(mar=c(4,5,1,1)) 

plot(y[1,], type="l", bty="n", col=plotColors[1], xlab="time", ylab="observation",ylim=c(-3,3))
for (ii in 2:nData){
  lines(y[ii,], col=plotColors[1])
}


plot(mu[,1], type="l", bty="n", col=plotColors[2], xlab="time", ylab=expression(mu),ylim=c(-3,3))
lines(mu[,2], col=plotColors[5])
lines(mu[,3], col=plotColors[6])

hist(mu[nBurnIn:nMCMC],          breaks=floor(sqrt(nMCMC-nBurnIn)), main=NA, 
     xlim=c(0.3,0.7), freq=FALSE, col="darkgrey", border=NA, 
     xlab=expression(mu), ylab="posterior estimate")
lines(density(mu[nBurnIn:nMCMC],from=0.3,to=0.7),lwd=2, col=plotColors[2])

plot(rho, type="l", bty="n", col=plotColors[3], xlab="time", ylab=expression(rho),ylim=c(-1,1))

hist(rho[nBurnIn:nMCMC],         breaks=floor(sqrt(nMCMC-nBurnIn)), main=NA, 
     xlim=c(0.6,1), freq=FALSE, col="darkgrey", border=NA, 
     xlab=expression(rho), ylab="posterior estimate")
lines(density(rho[nBurnIn:nMCMC],from=0.6,to=1.0),lwd=2, col=plotColors[3])

plot(sqrt(sigma), type="l", bty="n", col=plotColors[4], xlab="time", ylab=expression(sigma),ylim=c(0,0.5))

hist(sqrt(sigma[nBurnIn:nMCMC]), breaks=floor(sqrt(nMCMC-nBurnIn)), main=NA, 
     xlim=c(0.15,0.25), freq=FALSE, col="darkgrey", border=NA, 
     xlab=expression(sigma), ylab="posterior estimate")
lines(density(sqrt(sigma[nBurnIn:nMCMC]),from=0.15,to=0.25),lwd=2, col=plotColors[4])
