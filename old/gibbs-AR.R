# Setup plot colors
library("RColorBrewer")
plotColors = brewer.pal(6, "Dark2");

nMCMC         <- 10000
nBurnIn       <- 2500
nObservations <- 100

##############################################################################
# Generate data
##############################################################################

y    <- matrix(0, nrow=nObservations, ncol=1)
y[1] <- rnorm(1);
par  <- c( 0.5, 0.8, 0.2)
  
for (tt in 2:nObservations) {
  y[tt] <- par[1] + par[2] * ( y[tt-1] - par[1] ) + par[3] * rnorm(1)
}

##############################################################################
# Gibbs sampler
##############################################################################

mu    <- matrix(0, nrow=nMCMC, ncol=1)
rho   <- matrix(0, nrow=nMCMC, ncol=1)
sigma <- matrix(1, nrow=nMCMC, ncol=1)

for ( ii in 2:nMCMC ) {
  # rho | mu, sigma (normal)
  rhobar  <- sum( ( y[-1] - mu[ii-1] ) * ( y[-nObservations] - mu[ii-1] ) ) / sum( ( y[-1] - mu[ii-1] )^2 )
  rho[ii] <- rnorm(1, rhobar, sqrt( sigma[ii-1] / sum( ( y[-1] - mu[ii-1] )^2 ) ) )
  
  # mu | rho, sigma (normal)
  mubar  <- mean( y[-1] - rho[ii] * y[-nObservations] ) / (1 - rho[ii] )
  mu[ii] <- rnorm(1, mubar, 
                  sqrt( sigma[ii-1] / (nObservations * (1 - rho[ii] )^2 ) ) )
  
  # sigma | mu, rho (inverse Gamma)
  sigmahat  <- sd( y[-1] - mu[ii] - rho[ii-1] * y[-nObservations] )
  sigma[ii] <- 1 / rgamma(1, 0.5*nObservations+1, 0.5*sum(((y[-1]-mu[ii])-rho[ii]*(y[-nObservations]-mu[ii]))^2))
}

##############################################################################
# Plot the results
##############################################################################

layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow = TRUE))  
par(mar=c(4,5,1,1)) 

plot(y, type="l", bty="n", col=plotColors[1], xlab="time", ylab="observation")

hist(mu[nBurnIn:nMCMC],          breaks=floor(sqrt(nMCMC-nBurnIn)), main=NA, 
     xlim=c(0.3,0.7), freq=FALSE, col="darkgrey", border=NA, 
     xlab=expression(mu), ylab="posterior estimate")
lines(density(mu[nBurnIn:nMCMC],from=0.3,to=0.7),lwd=2, col=plotColors[2])

hist(rho[nBurnIn:nMCMC],         breaks=floor(sqrt(nMCMC-nBurnIn)), main=NA, 
     xlim=c(0.6,1), freq=FALSE, col="darkgrey", border=NA, 
     xlab=expression(rho), ylab="posterior estimate")
lines(density(rho[nBurnIn:nMCMC],from=0.6,to=1.0),lwd=2, col=plotColors[3])

hist(sqrt(sigma[nBurnIn:nMCMC]), breaks=floor(sqrt(nMCMC-nBurnIn)), main=NA, 
     xlim=c(0.15,0.25), freq=FALSE, col="darkgrey", border=NA, 
     xlab=expression(sigma), ylab="posterior estimate")
lines(density(sqrt(sigma[nBurnIn:nMCMC]),from=0.15,to=0.25),lwd=2, col=plotColors[4])
