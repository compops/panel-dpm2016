library("mvtnorm")

###############################################################################
# Generate data
###############################################################################

generateData <- function( beta, Q, w, s2e, alpha, omega, dims ) {
  
  nIndividuals  <- dims[1]
  nAttributes   <- dims[2]
  nObservations <- dims[3]
  
  x        <- array( rnorm(nAttributes * nIndividuals * nObservations), 
                     dim=c(nAttributes,nIndividuals,nObservations) )
  y        <- matrix(0, nrow=nIndividuals, ncol=nObservations)
  IT       <- matrix(0, nrow=nIndividuals, ncol=1)
  betaIndT <- matrix(0, nrow=nIndividuals, ncol=nAttributes)
  
  for ( ii in 1:nIndividuals ) {
    IT[ii]       <- sample(length(w), 1, prob=w)
    betaIndT[ii] <- beta[IT[ii]] + sqrt(Q[IT[ii]]) * rnorm(1)
    
    for ( tt in 1:nObservations ) {
      y[ii,tt] <- alpha + sum( betaIndT[ii] * x[,ii,tt] ) + sqrt( s2e/omega[ii] ) * rnorm( 1 )
    }
  }
  
  list(y=y, x=x, betaInd =betaIndT, I=IT)
}

###############################################################################
# Gibbs sampling for finite mixture
###############################################################################

gibbs_finitemixture <- function(y, xf, xr, nIter, nMaxClusters, prior, postGrid)  {
  
  nIndividuals  <- dim(y)[1]
  nObservations <- dim(y)[2]
  nAttributes   <- 2
  
  ###############################################################################
  # Define priors
  ###############################################################################
  a0star <- prior$a0star
  A0star <- prior$A0star
  
  c0Q    <- prior$c0Q
  C0Q    <- prior$C0Q
  
  c0e    <- prior$c0e
  C0e    <- prior$C0e
  
  eta0   <- prior$eta0
  nu     <- prior$nu
  
  ###############################################################################
  # Preallocate vectors and matrices
  ###############################################################################
  
  post_alpha        <- matrix(0, nrow=nIter+1, ncol=2)
  post_omega        <- matrix(1, nrow=nIter+1, ncol=nIndividuals)
  post_eta          <- matrix(0, nrow=nIter+1, ncol=nMaxClusters)
  post_s2e          <- matrix(1, nrow=nIter+1, ncol=1)
  
  post_beta         <- array(0, dim=c(nIter+1,nMaxClusters,nAttributes))
  post_Qinv         <- array(0, dim=c(nIter+1,nMaxClusters,nAttributes,nAttributes))
  post_betaInd      <- array(1, dim=c(nIter+1,nIndividuals,nAttributes))
  
  nMembers          <- matrix(0, nrow=nIter+1, ncol=nMaxClusters)
  nOccupiedClusters <- matrix(0, nrow=nIter+1, ncol=1)
  
  I                 <- matrix(0, nrow=nIter+1, ncol=nIndividuals)
  
  #post_randomeffects <- array(0, dim=c(nIter,length(postGrid[1,]),2))
  
  ###############################################################################
  # Initalise quantaties
  ###############################################################################
  
  I[1,]             <- sample(nMaxClusters, nIndividuals, replace=TRUE)
  post_s2e[1]       <- 1
  post_omega[1,]    <- rep(1, nIndividuals)
  
  for (jj in 1:nMaxClusters) { post_Qinv[1,jj,,] = diag(nAttributes )}
  
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
      etaUnnorm[jj] <- rgamma(1, eta0 + nMembers[kk,jj], 1)
      
    }
    
    post_eta[kk+1,] <- etaUnnorm / sum(etaUnnorm)
    
    #-----------------------------------------------------------------------------
    # Sample regression coefficients (fixed effects)
    # Basically multivariate Bayesian linear regression
    # [ Multivariate Gaussian ]
    #-----------------------------------------------------------------------------
    
    aNstar <- matrix( 0, nrow=2+nAttributes*nMaxClusters, ncol=1)
    ANstar <- matrix( 0, nrow=2+nAttributes*nMaxClusters, ncol=2+nAttributes*nMaxClusters)
    
    for ( ii in 1:nIndividuals ) {
      # As defined under (8.54)
      Z     <- matrix( 0, nrow=nObservations, ncol=nAttributes*nMaxClusters+2 )
      Z[,1:2] <- xf[,ii,]
      Z[,2+seq((I[kk,ii]-1)*nAttributes+1,I[kk,ii]*nAttributes)] <- t(xr[,ii,])
      
      # Equation (8.50)
      V  <- t(xr[,ii,]) %*% solve(post_Qinv[kk,I[kk,ii],,]) %*% xr[,ii,] + post_s2e[kk] / post_omega[ kk, ii ] * diag(nObservations)
      
      # Update from page 267
      ANstar <- ANstar + t(Z) %*% solve(V) %*% Z + solve( A0star )
      aNstar <- aNstar + t(Z) %*% solve(V) %*% y[ii,] + solve(A0star) %*% a0star
    }
    
    ANstar <- solve(ANstar)
    aNstar <- ANstar %*% aNstar
    
    # Sample (alpha, beta_1, ..., beta_K)
    alphastar <- rmvnorm(1, mean=aNstar, sigma=ANstar )
    
    post_alpha[kk+1,]   <- alphastar[1:2]
    post_beta[kk+1,,]   <- matrix(alphastar[-(1:2)],nrow=nMaxClusters,byrow = TRUE)
    
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
          CQ <- CQ + ( post_betaInd[kk,ii,] - post_beta[kk+1,jj,] ) %*% t( post_betaInd[kk,ii,] - post_beta[kk+1,jj,] )
        }
        
        CQ <- C0Q + 0.5 * C0Q
        
        # Sampling Q_k from Wishart posterior or Gamma posterior
        post_Qinv[kk+1,jj,,] <- rWishart(1, df=cQ, Sigma=CQ)
      } else {
        
        # Sampling Q_k from Wishart prior or Gamma prior
        post_Qinv[kk+1,jj,,] <- rWishart(1, df=c0Q, Sigma=C0Q)
      }
      
      # Compute p( beta_i^s | y ) via (8.48) on page 264
      #post_randomeffects[kk,,1] = post_randomeffects[kk,,1] + post_eta[kk+1,jj] * dnorm( postGrid[1,], post_beta[kk+1,jj,1],sqrt(1/post_Qinv[kk+1,jj,1,1]))
      #post_randomeffects[kk,,2] = post_randomeffects[kk,,2] + post_eta[kk+1,jj] * dnorm( postGrid[2,], post_beta[kk+1,jj,2],sqrt(1/post_Qinv[kk+1,jj,2,2]))
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
      z     <- y[ii,] - post_alpha[kk+1,] %*% xf[,ii,] - post_betaInd[kk,ii,] %*% xr[,ii,]
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
        V  <- t(xr[,ii,]) %*% solve(post_Qinv[kk+1,jj,,]) %*% xr[,ii,] + post_s2e[kk+1] / post_omega[ kk, ii ] * diag(nObservations)
        
        # Compute the weigth from (8.55)
        pi[ii,jj] <- log(post_eta[kk+1,jj]) + dmvnorm( y[ii,], post_alpha[kk+1,] %*% xf[,ii,] + post_beta[kk+1,jj,] %*% xr[,ii,] , V, log=TRUE )
        
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
      bS <- solve( post_Qinv[kk+1,I[kk+1,ii],,] + xr[,ii,] %*% t(xr[,ii,]) * post_omega[kk,ii] / post_s2e[kk+1] )
      
      # Compute the posterior mean using (8.52)
      bs <- bS %*% ( post_Qinv[kk+1,I[kk+1,ii],,] %*% post_beta[kk,I[kk+1,ii],] + xr[,ii,] %*% t(y[ii,] - post_alpha[kk+1,] %*% xf[,ii,]) * post_omega[kk,ii] / post_s2e[kk+1] )
      
      # Sample from Gaussian
      post_betaInd[kk+1,ii,] <- rmvnorm(1, mean=bs, sigma=bS)
    }
    
    ##############################################################################
    ## Sample the variance heterogenity
    # [ Gamma ]
    ##############################################################################
    
    for (ii in 1:nIndividuals) {
      
      # Compute sufficient statistics using equations from middle of page 269
      z  <- y[ii,] - post_alpha[kk+1,] %*% xf[,ii,] - post_betaInd[kk+1,ii,] %*% xr[,ii,]
      cw <- 0.5 * nu + 0.5 * nObservations
      Cw <- 0.5 * nu + 0.5 / post_s2e[kk+1] * t(z) %*% z
      
      # Sample the parameter from a gamma distribution
      post_omega[kk+1,ii] <- rgamma(1, shape=cw, rate=Cw)
    }
    
    ##############################################################################
    # Print iteration counter
    ##############################################################################
    if (kk %% 100 == 0) {
      print(paste(paste(paste("iteration: ",kk,sep="")," of ", nIter, sep="")," complete.",sep=""))
    }
  }
  
  list(alpha=post_alpha, var=sweep(1/post_omega,1,post_s2e,'*'), eta=post_eta, 
       nOccupiedClusters=nOccupiedClusters, nMembers=nMembers,
       #posterior_beta=post_randomeffects, 
       betaInd=post_betaInd, beta=post_beta, Qinv=post_Qinv)
}


###############################################################################
# Gibbs sampling for DPM
###############################################################################

gibbs_dpm <- function(y, xf, xr, nIter, nMaxClusters, prior, postGrid)  {
  
  nIndividuals  <- dim(y)[1]
  nObservations <- dim(y)[2]
  nAttributes   <- 2
  
  ###############################################################################
  # Define priors
  ###############################################################################
  a0star <- prior$a0star
  A0star <- prior$A0star
  
  c0Q    <- prior$c0Q
  C0Q    <- prior$C0Q
  
  c0e    <- prior$c0e
  C0e    <- prior$C0e
  
  eta0   <- prior$eta0
  nu     <- prior$nu
  
  alpha0 <- prior$alpha0
  beta0  <- prior$beta0    
  
  ###############################################################################
  # Preallocate vectors and matrices
  ###############################################################################
  
  post_alphaSB      <- matrix(1, nrow=nIter+1,      ncol=1)
  cCluster          <- matrix(1, nrow=nMaxClusters, ncol=1)

  post_alpha        <- matrix(0, nrow=nIter+1, ncol=2)
  post_omega        <- matrix(1, nrow=nIter+1, ncol=nIndividuals)
  post_eta          <- matrix(0, nrow=nIter+1, ncol=nMaxClusters)
  post_s2e          <- matrix(1, nrow=nIter+1, ncol=1)
  
  post_beta         <- array(0, dim=c(nIter+1,nMaxClusters,nAttributes))
  post_Qinv         <- array(0, dim=c(nIter+1,nMaxClusters,nAttributes,nAttributes))
  post_betaInd      <- array(1, dim=c(nIter+1,nIndividuals,nAttributes))
  
  nMembers          <- matrix(0, nrow=nIter+1, ncol=nMaxClusters)
  nOccupiedClusters <- matrix(0, nrow=nIter+1, ncol=1)
  
  I                 <- matrix(0, nrow=nIter+1, ncol=nIndividuals)
  
  #post_randomeffects <- array(0, dim=c(nIter,length(postGrid[1,]),2))
  
  ###############################################################################
  # Initalise quantaties
  ###############################################################################
  
  I[1,]             <- sample(nMaxClusters, nIndividuals, replace=TRUE)
  post_s2e[1]       <- 1
  post_omega[1,]    <- rep(1, nIndividuals)
  
  for (jj in 1:nMaxClusters) { post_Qinv[1,jj,,] = diag(nAttributes )}
  
  ###############################################################################
  # Main MCMC loop
  ###############################################################################
  
  for (kk in 1:nIter) {
    
    ##############################################################################
    ## Parameter estimation conditional on cluster allocation
    ##############################################################################
    
    rStick <- 1.0
    
    #-----------------------------------------------------------------------------
    # Sample eta 
    # [ Stick-breaking ]
    #-----------------------------------------------------------------------------
    
    etaUnnorm <- rep(0, nMaxClusters)
    
    for (jj in 1:nMaxClusters) {
      
      # Count the number of members in each cluster
      nMembers[kk,jj] <- length(which(I[kk,] == jj))
      
      # Count the number of occupied clusters
      if (nMembers[kk,jj] > 0) { 
        nOccupiedClusters[kk] <- nOccupiedClusters[kk] + 1 
      }
      
      # Compute the proportion of the stick to break off
      nR           <- nIndividuals - sum(nMembers[kk,1:jj])
      cCluster[jj] <- rbeta( 1, 1 + nMembers[kk,jj], post_alphaSB[kk] + nR );  
      
      # If this is the last cluster, use the remainder of the stick
      if (jj == nMaxClusters) { 
        cCluster[jj] <- 1 
      };
      
      etaUnnorm[jj] <- rStick * cCluster[jj];
      rStick        <- rStick * ( 1 - cCluster[jj] );      
    }
    
    post_eta[kk+1,] <- etaUnnorm / sum(etaUnnorm)
    
    #-----------------------------------------------------------------------------
    # Sample regression coefficients (fixed effects)
    # Basically multivariate Bayesian linear regression
    # [ Multivariate Gaussian ]
    #-----------------------------------------------------------------------------
    
    aNstar <- matrix( 0, nrow=2+nAttributes*nMaxClusters, ncol=1)
    ANstar <- matrix( 0, nrow=2+nAttributes*nMaxClusters, ncol=2+nAttributes*nMaxClusters)
    
    for ( ii in 1:nIndividuals ) {
      # As defined under (8.54)
      Z     <- matrix( 0, nrow=nObservations, ncol=nAttributes*nMaxClusters+2 )
      Z[,1:2] <- xf[,ii,]
      Z[,2+seq((I[kk,ii]-1)*nAttributes+1,I[kk,ii]*nAttributes)] <- t(xr[,ii,])
      
      # Equation (8.50)
      V  <- t(xr[,ii,]) %*% solve(post_Qinv[kk,I[kk,ii],,]) %*% xr[,ii,] + post_s2e[kk] / post_omega[ kk, ii ] * diag(nObservations)
      
      # Update from page 267
      ANstar <- ANstar + t(Z) %*% solve(V) %*% Z + solve( A0star )
      aNstar <- aNstar + t(Z) %*% solve(V) %*% y[ii,] + solve(A0star) %*% a0star
    }
    
    ANstar <- solve(ANstar)
    aNstar <- ANstar %*% aNstar
    
    # Sample (alpha, beta_1, ..., beta_K)
    alphastar <- rmvnorm(1, mean=aNstar, sigma=ANstar )
    
    post_alpha[kk+1,]   <- alphastar[1:2]
    post_beta[kk+1,,]   <- matrix(alphastar[-(1:2)],nrow=nMaxClusters,byrow = TRUE)
    
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
          CQ <- CQ + ( post_betaInd[kk,ii,] - post_beta[kk+1,jj,] ) %*% t( post_betaInd[kk,ii,] - post_beta[kk+1,jj,] )
        }
        
        CQ <- C0Q + 0.5 * C0Q
        
        # Sampling Q_k from Wishart posterior or Gamma posterior
        post_Qinv[kk+1,jj,,] <- rWishart(1, df=cQ, Sigma=CQ)
      } else {
        
        # Sampling Q_k from Wishart prior or Gamma prior
        post_Qinv[kk+1,jj,,] <- rWishart(1, df=c0Q, Sigma=C0Q)
      }
      
      # Compute p( beta_i^s | y ) via (8.48) on page 264
      #post_randomeffects[kk,,1] = post_randomeffects[kk,,1] + post_eta[kk+1,jj] * dnorm( postGrid[1,], post_beta[kk+1,jj,1],sqrt(1/post_Qinv[kk+1,jj,1,1]))
      #post_randomeffects[kk,,2] = post_randomeffects[kk,,2] + post_eta[kk+1,jj] * dnorm( postGrid[2,], post_beta[kk+1,jj,2],sqrt(1/post_Qinv[kk+1,jj,2,2]))
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
      z     <- y[ii,] - post_alpha[kk+1,] %*% xf[,ii,] - post_betaInd[kk,ii,] %*% xr[,ii,]
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
        V  <- t(xr[,ii,]) %*% solve(post_Qinv[kk+1,jj,,]) %*% xr[,ii,] + post_s2e[kk+1] / post_omega[ kk, ii ] * diag(nObservations)
        
        # Compute the weigth from (8.55)
        pi[ii,jj] <- log(post_eta[kk+1,jj]) + dmvnorm( y[ii,], post_alpha[kk+1,] %*% xf[,ii,] + post_beta[kk+1,jj,] %*% xr[,ii,] , V, log=TRUE )
        
      }
      
      # Log-weight trick to get rid of numerical problems
      wmax       <- max( pi[ii,] )
      wtra       <- exp( pi[ii,] - wmax )
      
      # Sample cluster index from multinomial
      I[kk+1,ii] <- sample( nMaxClusters, 1, prob=wtra/sum(wtra) )
    }
    
    ##############################################################################
    # Update alpha for stick-breaking
    # [Gamma]
    ##############################################################################
    tmp         <- log(1 - cCluster[1:(nMaxClusters-1)]) 
    tmp[is.infinite(tmp)] <- 0 
    
    post_alphaSB[kk+1] <- rgamma(1, alpha0 + nMaxClusters - 1, rate = beta0 - sum( tmp ))
    
    if (post_alphaSB[kk+1] == 0) {
      print( alpha0 + nMaxClusters - 1 )
      print( cCluster[1:(nMaxClusters-1)] ) 
      readline("Press <return to continue") 
    }    
    
    ##############################################################################
    ## Sample the random effects
    ## [ Gaussian ]
    ## IID data from unknown mean and covariance
    ##############################################################################
    
    for (ii in 1:nIndividuals) {
      
      # Compute the posterior covariance using (8.52)
      bS <- solve( post_Qinv[kk+1,I[kk+1,ii],,] + xr[,ii,] %*% t(xr[,ii,]) * post_omega[kk,ii] / post_s2e[kk+1] )
      
      # Compute the posterior mean using (8.52)
      bs <- bS %*% ( post_Qinv[kk+1,I[kk+1,ii],,] %*% post_beta[kk,I[kk+1,ii],] + xr[,ii,] %*% t(y[ii,] - post_alpha[kk+1,] %*% xf[,ii,]) * post_omega[kk,ii] / post_s2e[kk+1] )
      
      # Sample from Gaussian
      post_betaInd[kk+1,ii,] <- rmvnorm(1, mean=bs, sigma=bS)
    }
    
    ##############################################################################
    ## Sample the variance heterogenity
    # [ Gamma ]
    ##############################################################################
    
    for (ii in 1:nIndividuals) {
      
      # Compute sufficient statistics using equations from middle of page 269
      z  <- y[ii,] - post_alpha[kk+1,] %*% xf[,ii,] - post_betaInd[kk+1,ii,] %*% xr[,ii,]
      cw <- 0.5 * nu + 0.5 * nObservations
      Cw <- 0.5 * nu + 0.5 / post_s2e[kk+1] * t(z) %*% z
      
      # Sample the parameter from a gamma distribution
      post_omega[kk+1,ii] <- rgamma(1, shape=cw, rate=Cw)
    }
    
    ##############################################################################
    # Print iteration counter
    ##############################################################################
    if (kk %% 100 == 0) {
      print(paste(paste(paste("iteration: ",kk,sep="")," of ", nIter, sep="")," complete.",sep=""))
    }
  }
  
  list(alpha=post_alpha, var=sweep(1/post_omega,1,post_s2e,'*'), eta=post_eta, 
       nOccupiedClusters=nOccupiedClusters, nMembers=nMembers,
       #posterior_beta=post_randomeffects, 
       betaInd=post_betaInd, beta=post_beta, Qinv=post_Qinv, alphaSB=post_alphaSB)  
  
}