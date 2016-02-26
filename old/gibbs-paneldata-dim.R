library("mvtnorm")
library("MCMCpack")

nAttributes   <- 1
nMaxClusters  <- 20
nIndividuals  <- 50
nObservations <- 20
nIter         <- 1000
postGrid      <- seq(-4,4,0.1)

###############################################################################
# Generate data
###############################################################################

s2eT   <- 1
wT        <- c(0.4, 0.3, 0.3)
betaT     <- c( 2, 3, -2)
QT        <- c( 1, 0.2, 0.2)
betaIndT  <- matrix(0, nrow=nIndividuals, ncol=nAttributes)
alphaT <- 2.0
omegaT <- abs(rnorm(nIndividuals))
x     <- array( rnorm(nAttributes*nIndividuals*nObservations), dim=c(nAttributes,nIndividuals,nObservations) )
y     <- matrix(0, nrow=nIndividuals, ncol=nObservations)
IT    <- matrix(0, nrow=nIndividuals, ncol=1)

for ( ii in 1:nIndividuals ) {
  IT[ii]       <- sample(length(wT), 1, prob=wT)
  betaIndT[ii] <- betaT[IT[ii]] + sqrt(QT[IT[ii]]) * rnorm(1)
  
  for ( tt in 1:nObservations ) {
    y[ii,tt] = alphaT + sum( betaIndT[ii]*x[,ii,tt] ) + sqrt( s2eT/omegaT[ii] ) * rnorm( 1 )
  }
}

###############################################################################
# Define priors
###############################################################################
a0star <- rep(0,1+nAttributes*nMaxClusters)
A0star <- 4*diag(1+nAttributes*nMaxClusters)

c0Q    <- 2
C0Q    <- 1

c0e    <- 1
C0e    <- 1

a0     <- 0.1
nu     <- 0.01

###############################################################################
# Preallocate vectors and matrices
###############################################################################

post_alpha       <- matrix(0, nrow=nIter,   ncol=1)
post_omega       <- matrix(1, nrow=nIter+1, ncol=nIndividuals)
post_eta         <- matrix(0, nrow=nIter,   ncol=nMaxClusters)
post_s2e         <- matrix(1, nrow=nIter+1, ncol=1)

post_beta         <- array(0, dim=c(nIter,nMaxClusters,nAttributes))
Qinv              <- array(0, dim=c(nMaxClusters,nAttributes,nAttributes))

post_betaInd      <- matrix(0, nrow=nIndividuals, ncol=nAttributes)

nMembers          <- matrix(0, nrow=nIter, ncol=nMaxClusters)
nOccupiedClusters <- matrix(0, nrow=nIter, ncol=1)

I                 <- matrix(0, nrow=nIter+1, ncol=nIndividuals)

post_randomeffects <- matrix(0, nrow=nIter, ncol=length(postGrid))

###############################################################################
# Initalise quantaties
###############################################################################

post_betaInd  <- rmvnorm(nIndividuals, mean = rep(0,nAttributes))
for ( jj in 1:nMaxClusters ) { Qinv[jj,,] <- diag(nAttributes) }
I[1,]         <- sample(nMaxClusters, nIndividuals, replace=TRUE)

post_s2e[1]    <- 1
post_omega[1,] <- rep(1, nIndividuals)

###############################################################################
# Main MCMC loop
###############################################################################

for (kk in 1:nIter) {
  
  ##############################################################################
  ## Parameter estimation conditional on cluster allocation
  ##############################################################################
  
  #-----------------------------------------------------------------------------
  # Sample eta 
  # [ Dirichlet distribution ]
  #-----------------------------------------------------------------------------
  
  etaUnnorm <- rep(0, nMaxClusters)
  
  for (jj in 1:nMaxClusters) {
    
    # Count the number of members in each cluster
    nMembers[kk,jj] <- length(which(I[kk,] == jj))
    
    # Count the number of occupied clusters
    if (nMembers[kk,jj] > 0) { 
      nOccupiedClusters[kk] <- nOccupiedClusters[kk] + 1 
    }
    
    # Sample probabilities p from Dirichlet (as normalised Gamma variables)
    etaUnnorm[jj] <- rgamma(1, a0 + nMembers[kk,jj], 1)
    
  }
  
  post_eta[kk,] <- etaUnnorm / sum(etaUnnorm)
  
  #-----------------------------------------------------------------------------
  # Sample regression coefficients (fixed effects)
  # Basically multivariate Bayesian linear regression
  # [ Multivariate Gaussian ]
  #-----------------------------------------------------------------------------
  
  aNstar <- matrix( 0, nrow=1+nAttributes*nMaxClusters, ncol=1)
  ANstar <- matrix( 0, nrow=1+nAttributes*nMaxClusters, ncol=1+nAttributes*nMaxClusters)
  
  for ( ii in 1:nIndividuals ) {
    # As defined under (8.54)
    Z     <- matrix( 0, nrow=nObservations, ncol=nAttributes*nMaxClusters+1 )
    Z[,1] <- 1
    Z[,1+seq((I[kk,ii]-1)*nAttributes+1,I[kk,ii]*nAttributes)] <- t(x[,ii,])
    
    # Equation (8.50)
    #V  <- t(x[,ii,]) %*% solve(Qinv[I[kk,ii],,]) %*% x[,ii,] + post_s2e[kk] / post_omega[ kk, ii ] * diag(nObservations)
    V  <- x[,ii,] %*% t(x[,ii,]) / Qinv[I[kk,ii]] + post_s2e[kk] / post_omega[ kk, ii ] * diag(nObservations)
    
    # Update from page 267
    ANstar <- ANstar + t(Z) %*% solve(V) %*% Z + solve( A0star )
    aNstar <- aNstar + t(Z) %*% solve(V) %*% y[ii,] + solve(A0star) %*% a0star
  }
  
  ANstar <- solve(ANstar)
  aNstar <- ANstar %*% aNstar
  
  # Sample (alpha, beta_1, ..., beta_K)
  alphastar <- rmvnorm(1, mean=aNstar, sigma=ANstar )
  
  post_alpha[kk]    <- alphastar[1]
  post_beta[kk,,]   <- matrix(alphastar[-1],nrow=nMaxClusters,byrow = TRUE)
  
  #-----------------------------------------------------------------------------
  # Sample covariance matrix
  # [ Wishart ]
  #-----------------------------------------------------------------------------
  for (jj in 1:nMaxClusters) {
    
    # Check if we should sample from the posterior or prior
    if ( nMembers[kk,jj] > 1 ) {
      
      # Equations from bottom of page 268
      cQ <- c0Q + 0.5 * nMembers[kk,jj]
      CQ <- matrix(0, nrow=nAttributes, ncol=nAttributes)
      
      # Sum over all individuals in cluster jj
      for ( ii in which(I[kk,]==jj) ) {
        CQ <- CQ + ( post_betaInd[ii,] - post_beta[kk,jj,] ) %*% t( post_betaInd[ii,] - post_beta[kk,jj,] )
      }
      
      CQ <- C0Q + 0.5 * C0Q
      
      # Sampling Q_k from Wishart posterior or Gamma posterior
      #Qinv[jj] <- rWishart(1, cQ, CQ)
      Qinv[jj] <- rgamma(1, shape=cQ, rate=CQ)
    } else {
      
      # Sampling Q_k from Wishart prior or Gamma prior
      #Qinv[jj] <- rWishart(1, c0Q, C0Q)
      Qinv[jj] <- rgamma(1, shape=c0Q, rate=C0Q)
    }
    
    # Compute p( beta_i^s | y ) via (8.48) on page 264
    post_randomeffects[kk,] = post_randomeffects[kk,] + post_eta[kk,jj] * dnorm( postGrid, post_beta[kk,,][jj],sqrt(1/Qinv[jj]))
  }
  
  #-----------------------------------------------------------------------------
  # Sample variance
  # [ Inverse Gamma ]
  #-----------------------------------------------------------------------------
  
  # Compute the sufficient statistics from page 269
  cNeps <- c0e + 0.5 * nIndividuals * nObservations
  CNeps <- 0
  
  # Compute contribution from each individual
  for ( ii in 1:nIndividuals ) {
    z     <- y[ii,] - post_alpha[kk,] - post_betaInd[ii] * x[,ii,]
    CNeps <- CNeps  + post_omega[kk,ii] * t(z) %*% z
  }
  
  CNeps <- C0e + 0.5 * CNeps
  
  # Sample from inverse gamma
  post_s2e[kk+1] <- 1 / rgamma(1, shape=cNeps, scale=CNeps)
  
  ##############################################################################
  ## Classify each individual
  ## [ Multinomial ]
  ##############################################################################
  
  pi <- matrix(0, nrow=nIndividuals, ncol=nMaxClusters)
  
  # Compute and normalise the weights
  for (ii in 1:nIndividuals) {
    for (jj in 1:nMaxClusters) {
      
      # Compute the covariance matrix by (8.50)
      #V  <- t(x[,ii,]) %*% solve(Qinv[jj,,]) %*% x[,ii,] + post_s2e[kk] / post_omega[ kk, ii ] * diag(nObservations)
      V  <- x[,ii,] %*% t(x[,ii,]) / Qinv[jj,,] + post_s2e[kk] / post_omega[ kk, ii ] * diag(nObservations)
      
      # Compute the weigth from (8.55)
      #pi[ii,jj] <- log(post_eta[jj]) + dmvnorm( y[ii,], post_alpha[kk,] - colSums( post_beta[kk,jj,] * x[,ii,] ), V, log=TRUE )
      pi[ii,jj] <- log(post_eta[jj]) + dmvnorm( y[ii,], post_alpha[kk,] - post_beta[kk,jj,] * x[,ii,], V, log=TRUE )
    }
    
    # Log-weight trick to get rid of numerical problems
    wmax       <- max( pi[ii,] )
    wtra       <- exp( pi[ii,] - wmax )
    
    # Sample cluster index from multinomial
    I[kk+1,ii] <- sample( nMaxClusters, 1, prob=wtra/sum(wtra) )
  }
  
  ##############################################################################
  ## Sample the random effects
  ## [ Gaussian ]
  ## IID data from unknown mean and covariance
  ##############################################################################
  
  for (ii in 1:nIndividuals) {
    
    # Compute the posterior covariance using (8.52)
    bS <- solve( Qinv[I[kk+1,ii]] + t(x[,ii,]) %*% x[,ii,] * post_omega[kk,ii] / post_s2e[kk+1] )
    
    # Compute the posterior mean using (8.52)
    bs <- bS %*% ( Qinv[I[kk+1,ii]] %*% post_beta[kk,I[kk+1,ii],] + x[,ii,] %*% (y[ii,] - post_alpha[kk] ) * post_omega[kk,ii] / post_s2e[kk+1] )
    
    # Sample from Gaussian
    post_betaInd[ii,] <- rmvnorm(1, mean=bs, sigma=bS)
  }
  
  ##############################################################################
  ## Sample the variance heterogenity
  # [ Gamma ]
  ##############################################################################
  
  for (ii in 1:nIndividuals) {
    
    # Compute sufficient statistics using equations from middle of page 269
    z     <- y[ii,] - post_alpha[kk,] - post_betaInd[ii,] * x[,ii,]
    cw <- 0.5 * nu + 0.5 * nObservations
    Cw <- 0.5 * nu + 0.5 / post_s2e[kk+1] * t(z) %*% z
    
    # Sample the parameter from a gamma distribution
    post_omega[kk+1,ii] <- rgamma(1, shape=cw, rate=Cw)
  }
  
}

# Setup plot colors
library("RColorBrewer")
plotColors = brewer.pal(6, "Dark2");

truePosterior <- wT[1] * dnorm(postGrid,betaT[1],sqrt(QT[1])) + wT[2] * dnorm(postGrid,betaT[2],sqrt(QT[2])) + wT[3] * dnorm(postGrid,betaT[3],sqrt(QT[3]))

layout(matrix(c(1,1,1,2,3,4,5,6,6), 3, 3, byrow = TRUE))  
par(mar=c(4,5,1,1)) 

plot(postGrid,truePosterior,type="l", bty="n",ylim=c(0,0.6),lwd=2,lty="dashed")
lines(postGrid,colMeans(post_randomeffects),col=plotColors[1],lwd=2)

foo <- apply(post_beta,1,sort)
hist(foo[1,],main="",breaks=floor(sqrt(nIter)),freq=FALSE,xlab=expression(beta[1]),col="darkgrey",border=NA)
lines(density(foo[1,]),lwd=2,col=plotColors[2])
abline(v=sort(betaT)[1],lty="dotted")

hist(foo[2,],main="",breaks=floor(sqrt(nIter)),freq=FALSE,xlab=expression(beta[2]),col="darkgrey",border=NA)
lines(density(foo[2,]),lwd=2,col=plotColors[3])
abline(v=sort(betaT)[2],lty="dotted")

hist(foo[3,],main="",breaks=floor(sqrt(nIter)),freq=FALSE,xlab=expression(beta[3]),col="darkgrey",border=NA)
lines(density(foo[3,]),lwd=2,col=plotColors[4])
abline(v=sort(betaT)[3],lty="dotted")

hist(post_alpha,main="",breaks=floor(sqrt(nIter)),freq=FALSE,xlab=expression(alpha),col="darkgrey",border=NA)
lines(density(post_alpha),lwd=2,col=plotColors[6])
abline(v=alphaT,lty="dotted")

hist(nOccupiedClusters[,1],breaks=floor(sqrt(nIter)),main="",freq=FALSE,xlab="no. occupied clusters",col="darkgrey",border=NA,xlim=c(0,nMaxClusters))