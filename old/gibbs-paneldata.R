library("mvtnorm")
library("MCMCpack")

nAttributes   <- 1
nMaxClusters  <- 3
nIndividuals  <- 20
nObservations <- 10
nIter         <- 1000

###############################################################################
# Generate data
###############################################################################

s2eT   <- 1
wT        <- c(0.7, 0.1, 0.2)
betaT     <- c( 2, 3, -1)
QT        <- c( 1, 3, 1)
betaIndT  <- matrix(0, nrow=nIndividuals, ncol=nAttributes)
alphaT <- 0.3
omegaT <- abs(rnorm(nIndividuals))
x     <- array( rnorm(nAttributes*nIndividuals*nObservations), dim=c(nAttributes,nIndividuals,nObservations) )
y     <- matrix(0, nrow=nIndividuals, ncol=nObservations)
IT    <- matrix(0, nrow=nIndividuals, ncol=1)

for ( ii in 1:nIndividuals ) {
  IT[ii]       <- sample(nMaxClusters, 1, prob=wT)
  betaIndT[ii] <- betaT[IT[ii]] + sqrt(QT[IT[ii]]) * rnorm(1)
  
  for ( tt in 1:nObservations ) {
    y[ii,tt] = alphaT * sum( betaIndT[ii]*x[,ii,tt] ) + sqrt( s2eT/omegaT[ii] ) * rnorm( 1 )
  }
}

###############################################################################
# Define priors
###############################################################################
a0star <- rep(0,1+nAttributes*nMaxClusters)
A0star <- diag(1+nAttributes*nMaxClusters)

c0Q    <- 2
C0Q    <- diag(2)

c0e    <- 1
C0e    <- 1

a0     <- 1
nu     <- 1

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
nOccupiedClusters <- matrix(0, nrow=nIter, ncol=nMaxClusters)

I                 <- matrix(0, nrow=nIter+1, ncol=nIndividuals)

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
  # [ Multivariate Gaussian ]
  #-----------------------------------------------------------------------------
  
  aNstar <- matrix( 0, nrow=1+nAttributes*nMaxClusters, ncol=1)
  ANstar <- matrix( 0, nrow=1+nAttributes*nMaxClusters, ncol=1+nAttributes*nMaxClusters)
  
  for ( ii in 1:nIndividuals ) {
    Z     <- matrix( 0, nrow=nObservations, ncol=nAttributes*nMaxClusters+1 )
    Z[,1] <- 1
    Z[,1+seq((I[kk,ii]-1)*nAttributes+1,I[kk,ii]*nAttributes)] <- t(x[,ii,])
    
    V  <- t(x[,ii,]) %*% solve(Qinv[I[kk,ii],,]) %*% x[,ii,] + post_s2e[kk] / post_omega[ kk, ii ] * diag(nObservations)
  
    ANstar <- ANstar + t(Z) %*% solve(V) %*% Z + solve( A0star )
    aNstar <- aNstar + t(Z) %*% solve(V) %*% y[ii,] + solve(A0star) %*% a0star
  }
  
  aNstar <- solve(ANstar) %*% aNstar
  
  alphastar <- rmvnorm(1, mean=aNstar, sigma=ANstar )
  
  post_alpha[kk,]  <- alphaT# alphastar[1]
  post_beta[kk,,]  <- matrix(alphastar[-1],nrow=nMaxClusters,byrow = TRUE)
  
  #-----------------------------------------------------------------------------
  # Sample covariance matrix
  # [ Wishart ]
  #-----------------------------------------------------------------------------
  for (jj in 1:nMaxClusters) {
    if ( nMembers[kk,jj] > 1 ) {
      cQ <- c0Q + 0.5 * nMembers[kk,jj]
      CQ <- matrix(0, nrow=nAttributes, ncol=nAttributes)
      
      for ( ii in which(I[kk,]==jj) ) {
        CQ <- CQ + ( post_betaInd[ii,] - post_beta[kk,jj,] ) %*% t( post_betaInd[ii,] - post_beta[kk,jj,] )
      }
      
      CQ <- C0Q + 0.5 * C0Q
      
      Qinv[jj,,] <- diag(nAttributes) #rWishart(1, cQ, CQ)
    } else {
      Qinv[jj,,] <- rWishart(1, c0Q, C0Q)
    }
  }
  
  #-----------------------------------------------------------------------------
  # Sample variance
  # [ Inverse Gamma ]
  #-----------------------------------------------------------------------------
  cNeps <- c0e + 0.5 * nIndividuals * nObservations
  CNeps <- 0
  
  for ( ii in 1:nIndividuals ) {
    z     <- y[ii,] - post_alpha[kk,] - colSums( post_betaInd[ii,] * x[,ii,] )
    CNeps <- CNeps + post_omega[ii] * t(z) %*% z
  }
  
  CNeps <- C0e + 0.5 * CNeps
  
  post_s2e[kk+1] <- s2eT# rinvgamma(1, cNeps, CNeps)
  
  ##############################################################################
  ## Classify each individual
  ## [ Multinomial ]
  ##############################################################################
  
  pi <- matrix(0, nrow=nIndividuals, ncol=nMaxClusters)
  
  # Compute and normalise the weights
  for (ii in 1:nIndividuals) {
    for (jj in 1:nMaxClusters) {
      V  <- t(x[,ii,]) %*% solve(Qinv[jj,,]) %*% x[,ii,] + post_s2e[kk] / post_omega[ kk, ii ] * diag(nObservations)
      pi[ii,jj] <- log(post_eta[jj]) + dmvnorm( y[ii,], post_alpha[kk,] - colSums( post_beta[kk,jj,] * x[,ii,] ), V, log=TRUE )
    }
    
    # Sample cluster belonging
    wmax       <- max( pi[ii,] )
    wtra       <- exp( pi[ii,] - wmax )
    I[kk+1,ii] <- sample( nMaxClusters, 1, prob=wtra/sum(wtra) )
  }
  
  ##############################################################################
  ## Sample the random effects
  ## [ Gaussian ]
  ##############################################################################
  
  for (ii in 1:nIndividuals) {
    bS <- solve( Qinv[I[kk+1,ii],,] + x[,ii,] %*% t(x[,ii,]) * post_omega[kk,ii] / post_s2e[kk+1] )
    bs <- bS %*% ( Qinv[I[kk+1,ii],,] %*% post_beta[kk,I[kk+1,ii],] + x[,ii,] %*% (y[ii,] - post_alpha[kk] ) * post_omega[kk,ii] / post_s2e[kk+1] )
    
    post_betaInd[ii,] <- betaIndT[ii] #rmvnorm(1, mean=bs, sigma=bS)
  }
  
  
  ## Sample the variance heterogenity
  # [ Gamma ]
  for (ii in 1:nIndividuals) {
    z     <- y[ii,] - post_alpha[kk,] - colSums( post_betaInd[ii,] * x[,ii,] )
    cw <- 0.5 * nu + 0.5 * nObservations
    Cw <- 0.5 * nu + 0.5 / post_s2e[kk+1] * t(z) %*% z
    
    post_omega[kk+1,ii] <- omegaT[ii] #rgamma(1, cw, Cw)
  }

}