repmat = function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}

rinvchisq <- function(n, df, scale=1/df)
{
  df <- rep(df, len=n); scale <- rep(scale, len=n)
  if(any(df <= 0)) stop("The df parameter must be positive.")
  if(any(scale <= 0)) stop("The scale parameter must be positive.")
  z <- rchisq(n, df=df)
  z <- ifelse(z == 0, 1e-100, z)
  x <- (df*scale) / z
  return(x)
}

dinvgamma <- function(x, a0, b0) {
  (b0^a0)/(gamma(a0)) * x^(-a0-1) * exp(-b0/x)
}


########################################################################################################
# Gibbs sampling for finite mixture
########################################################################################################

gibbs_finitemixture <- function(y, nMCMC, nBurnIn, nMaxClusters, prior, grids )  {
  
  nData         <- dim(y)[1]
  nObservations <- dim(y)[2]
  
  # Posterior grids
  postGamma    <- rep(0, length(grids$gamma) )
  postRho      <- rep(0, length(grids$rho) )
  postSigma    <- rep(0, length(grids$sigma) )
  
  # Pre-allocate matrices
  I                 <- matrix(1, nrow=nMCMC+1, ncol=nData )
  nMembersCluster   <- matrix(0, nrow=nMCMC,   ncol=nMaxClusters)
  nOccupiedClusters <- matrix(0, nrow=nMCMC,   ncol=1)
  ll                <- matrix(0, nrow=nMCMC,   ncol=1 )
  aN                <- matrix(1, nrow=nMCMC,   ncol=nMaxClusters )
  bN                <- matrix(1, nrow=nMCMC,   ncol=nMaxClusters )
  piPost            <- matrix(1, nrow=nMCMC,   ncol=nMaxClusters)
  mN                <- array( 1, dim=c(nMCMC,nMaxClusters,2) )
  VN                <- array( 1, dim=c(nMCMC,nMaxClusters,2,2) )
  sigma             <- matrix(1, nrow=nMCMC,   ncol=nMaxClusters)
  beta              <- array( 1, dim=c(nMCMC,nMaxClusters,2) )
  
    
  # Random initalisation of clusters
  I[2,] <- sample( nMaxClusters, nData, replace=T )
  
  # Make the regressor matrix
  Y     <- matrix(y[,-1],             nrow=nData * (nObservations - 1), ncol=1, byrow=TRUE)
  X     <- matrix(1,                  nrow=nData * (nObservations - 1), ncol=2)
  X[,2] <- matrix(y[,-nObservations], nrow=nData * (nObservations - 1), ncol=1, byrow=TRUE)
  
  # Hyperparameters in regression priors (shrinkage prior on beta)
  m0    <- prior$m0
  V0    <- prior$V0
  a0    <- prior$a0
  b0    <- prior$b0
  
  # Hyperparameters in mixture probability priors (smaller gives fewer clusters)
  a0    <- prior$a0;
  
  #############################################################################
  # Main MCMC loop
  #############################################################################
  
  for (ii in 2:nMCMC) {

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
      
      #========================================================================
      # Update the AR model for each cluster
      #========================================================================
      
      # Select the individuals in cluster jj
      idx <- which(I[ii,] == jj)
      
      # Check if members in cluster, otherwise sample the prior
      if (length(idx) > 1) {
        
        # Build the vector of observations and the regressor matrix
        ySub     <- matrix(y[idx,-1],             nrow=nMembersCluster[ii,jj] * (nObservations - 1), ncol=1, byrow=TRUE)
        XSub     <- matrix(1,                     nrow=nMembersCluster[ii,jj] * (nObservations - 1), ncol=2)
        XSub[,2] <- matrix(y[idx,-nObservations], nrow=nMembersCluster[ii,jj] * (nObservations - 1), ncol=1, byrow=TRUE)
        
        # Compute the posterior parameters
        VN[ii,jj,,] <- solve( t(XSub) %*% XSub + V0 )
        mN[ii,jj,]  <- VN[ii,jj,,] %*% ( t(XSub) %*% ySub + solve(V0) %*% m0 )
        aN[ii,jj]   <- a0 + 0.5 * length( idx );
        bN[ii,jj]   <- b0 + 0.5 * ( t(m0) %*% solve(V0) %*% m0 + t(ySub) %*% ySub - t(mN[ii,jj,]) %*% solve(VN[ii,jj,,]) %*% mN[ii,jj,] )
      } else {
        VN[ii,jj,,] <- V0;
        mN[ii,jj,]  <- m0;
        aN[ii,jj]   <- a0;
        bN[ii,jj]   <- b0;
      }
      
      sigma[ii,jj]  <- 1.0 / rgamma(1, aN[ii,jj], bN[ii,jj])
      beta[ii,jj,]  <- t(t(mvrnorm(1, rep(0, length(mN[ii,jj,])), bN[ii,jj]/aN[ii,jj]*VN[ii,jj,,]) / sqrt(2*aN[ii,jj] / rchisq(1, aN[ii,jj]))) + mN[ii,jj,])
    }
    
    #----------------------------------------------------------------------------------
    # Sample the indicators
    #----------------------------------------------------------------------------------
    piPost[ii,] <- piCluster / sum(piCluster)
    
    unnormw <- matrix(0, nrow = nData, ncol = nMaxClusters)
    
    for (jj in 1:nMaxClusters) {
      
      # Compute the likelihood of belonging to a certain cluster
      #tmp          <- dt( ( Y-X %*% mN[ii,jj,] ) / sqrt( bN[ii,jj] * diag( X %*% VN[ii,jj,,] %*% t(X) ) / aN[ii,jj] ), df=2*aN[ii,jj], log=TRUE )
      tmp          <- dnorm( Y, beta[ii,jj,1] + beta[ii,jj,2] * X, sqrt(sigma[ii,jj]), log=TRUE )
      unnormw[,jj] <- log(piPost[ii,jj]) + rowSums( matrix(tmp, nrow=nData, ncol=nObservations-1, byrow=TRUE) )
      
      if ( ii > nBurnIn ) {
        # Grid and evaluate the parameter posterior
        postGamma <- postGamma + ( piPost[ii,jj] * dt( (grids$gamma-mN[ii,jj,1])/(bN[ii,jj]*VN[ii,jj,1,1]/aN[ii,jj]), df=2*aN[ii,jj]) ) / ( nMCMC - nBurnIn )
        postRho   <- postRho   + ( piPost[ii,jj] * dt( (grids$rho-mN[ii,jj,2])/(bN[ii,jj]*VN[ii,jj,2,2]/aN[ii,jj]), df=2*aN[ii,jj]) ) / ( nMCMC - nBurnIn )
        postSigma <- postSigma + ( piPost[ii,jj] * dinvgamma(grids$sigma,aN[ii,jj],1/bN[ii,jj]))
      }
    }
    
    for ( jj in 1:nData ) {
      unnormwMax         <- max( unnormw[jj,] )
      unmormwTransformed <- exp( unnormw[jj,] - unnormwMax )
      
      I[ii+1,jj] <- sample(nMaxClusters, 1, replace=TRUE, 
                           unmormwTransformed / sum(unmormwTransformed) )
    }
    
    # Calculate likelihood
    ll[ii] <- sum(log(unmormwTransformed) + unnormwMax)
    
    #----------------------------------------------------------------------------------
    # Print iteration counter
    #----------------------------------------------------------------------------------    
    if (ii %% 500 == 0) {
      print(paste(paste(paste("iteration: ",ii,sep="")," of ", nMCMC, sep="")," complete.",sep=""))
    }
  }
  list(postGamma=postGamma, postRho=postRho, nOccupiedClusters=nOccupiedClusters, postSigma=postSigma, ll=ll)
}

########################################################################################################
# Gibbs sampling for finite mixture
########################################################################################################

gibbs_finitemixture2 <- function(y, nMCMC, nMaxClusters, prior, grid, xT)  {
  
  nData         <- dim(y)[1]
  nObservations <- dim(y)[2]
  
  # Hyperparameters in priors
  m0          <- prior$m0Rho
  v0          <- prior$v0Rho
  a0          <- prior$a0Sigma
  b0          <- prior$b0Sigma
  mu0Gamma    <- prior$mu0Gamma
  sd0Gamma    <- prior$sd0Gamma
  a0          <- prior$a0
  
  # Pre-allocate matrices
  I                 <- matrix(1, nrow=nMCMC+1, ncol=nData)
  nMembersCluster   <- matrix(0, nrow=nMCMC,   ncol=nMaxClusters)
  nOccupiedClusters <- matrix(0, nrow=nMCMC,   ncol=1)
  piPost            <- matrix(1, nrow=nMCMC,   ncol=nMaxClusters)
  sigma2e           <- matrix(1, nrow=nMCMC,   ncol=1)
  sigma2g           <- matrix(1, nrow=nMCMC,   ncol=nMaxClusters)
  rho               <- matrix(0, nrow=nMCMC,   ncol=1)
  gamma             <- matrix(0, nrow=nMCMC+1, ncol=nMaxClusters)
  muPost            <- matrix(1, nrow=nMCMC,   ncol=nMaxClusters)
  sdPost            <- matrix(1, nrow=nMCMC,   ncol=nMaxClusters)  
  gridPost          <- matrix(0, nrow=nMCMC,   ncol=length(grid))
  muN               <- matrix(1, nrow=nMCMC,   ncol=nMaxClusters)
  sdN               <- matrix(1, nrow=nMCMC,   ncol=nMaxClusters)    
  
  # Random initalisation
  I[1,]     <- sample(nMaxClusters, nData, replace=TRUE)
  gamma[1,] <- rnorm(nMaxClusters)
  
  #############################################################################
  # Main MCMC loop
  #############################################################################
  
  for (ii in 1:nMCMC) {
    
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
    
    #========================================================================
    # Sample the mean of the AR process (known variance)
    #========================================================================    
    for (jj in 1:nMaxClusters) {
      
      idx  <- which(I[ii,] == jj)
      nTotalObservations <- length(idx) * (nObservations - 1)
      
      xbar <- mean( y[idx,-1] - rho[ii] * y[idx,-nObservations] )
      w    <- ( nTotalObservations / sigma2e[ii] ) / ( ( nTotalObservations / sigma2e[ii] ) + ( 1 / sd0Gamma^2 ) )
      
      muN[ii,jj] <- w * xbar + ( 1 - w ) * mu0Gamma
      sdN[ii,jj] <- 1 / sqrt( ( nTotalObservations / sigma2e[ii] ) + ( 1 / sd0Gamma^2 ) )
      
      # Sample \mu | y, \sigma^2
      if (nMembersCluster[ii,jj] > 2) {
        # Compute the posterior
        gamma[ii,jj] <- rnorm(1, muN[ii,jj], sdN[ii,jj])
      } else {
        # Sample from the prior
        gamma[ii,jj] <- rnorm(1, mu0Gamma, sd0Gamma )
      }
      
    }
    
    #========================================================================
    # Sample the autocorrelation coefficient and noise variance
    #========================================================================
    
    nTotalObs <- nData * (nObservations - 1);
    
    # Compile regressors and observations
    # y_{it} - \gamma_i = rho * y_{i,t-1} + e_{it}
    #
    XLS <- matrix(y[,-nObservations],nrow=nTotalObs,ncol=1)
    yLS <- matrix(y[,-1]- repmat(matrix(gamma[ii,I[ii,]],nrow=nData,1),1,nObservations-1),nrow=nTotalObs,ncol=1)
    
    # Compute the posterior parameters
    VN <- solve(t(XLS) %*% XLS + v0)
    mN <- VN %*% (t(XLS) %*% yLS + solve(v0) %*% m0)
    aN <- a0 + 0.5 * nTotalObs;
    bN <- b0 + 0.5 * (t(m0) %*% solve(v0) %*% m0 + t(yLS) %*% yLS - t(mN) %*% solve(VN) %*% mN)
    
    sigma2e[ii+1] <- rinvchisq(1, df=aN, scale=bN/aN)
    rho[ii+1]     <- 0.5; #mN + sigma2e[ii+1] * VN * rnorm(1);
    
    #----------------------------------------------------------------------------------
    # Sample the indicators
    #----------------------------------------------------------------------------------
    piPost[ii,] <- piCluster / sum(piCluster)
    
    unnormw <- matrix(0, nrow = nData, ncol = nMaxClusters)
    
    for (jj in 1:nMaxClusters) {
      unnormw[,jj] <- log( piPost[ii,jj] ) + rowSums(dnorm(y[,-1], mean=gamma[ii,jj] + rho[ii+1] * y[,-nObservations], sd=sqrt(sigma2e[ii+1]), log=TRUE))
    }
    
    for ( jj in 1:nData ) {
      unnormwMax         <- max( unnormw[jj,] )
      unmormwTransformed <- exp( unnormw[jj,] - unnormwMax )
      
      I[ii+1,jj] <- sample(nMaxClusters, 1, replace=TRUE, 
                           unmormwTransformed / sum(unmormwTransformed) )
    }
    
    #----------------------------------------------------------------------------------
    # Compute density estimate
    #----------------------------------------------------------------------------------
    for ( jj in 1:nMaxClusters ) {
      gridPost[ii,] <- gridPost[ii,] + piPost[ii,jj] * dnorm(grid, muN[ii,jj], sdN[ii,jj])
    }
    
    #----------------------------------------------------------------------------------
    # Print iteration counter
    #----------------------------------------------------------------------------------    
    if (ii %% 500 == 0) {
      print(paste(paste(paste("iteration: ",ii,sep="")," of ", nMCMC, sep="")," complete.",sep=""))
    }
  }
  list(gridPost=gridPost, nOccupiedClusters=nOccupiedClusters, rho=rho, sigmae=sqrt(sigma2e), nMembersCluster=nMembersCluster, gamma=gamma )
}

########################################################################################################
# Collapsed Gibbs sampling for DPM process
########################################################################################################

gibbs_dpm <- function(y, nMCMC, nMaxClusters, prior, grid)  {
  
  nData <- length(y)
  
  # Hyperparameters in priors
  mu0    <- prior$mu0
  sd0    <- prior$sd0
  kappa0 <- prior$kappa0
  nu0    <- prior$nu0
  alpha0 <- prior$alpha0
  beta0  <- prior$beta0
  
  # Pre-allocate matrices
  I                 <- matrix(1, nrow=nMCMC+1, ncol=nData)
  nMembersCluster   <- matrix(1, nrow=nMCMC,   ncol=nMaxClusters)
  nOccupiedClusters <- matrix(0, nrow=nMCMC,   ncol=1)
  alpha             <- matrix(1, nrow=nMCMC,   ncol=1)
  piPost            <- matrix(1, nrow=nMCMC,   ncol=nMaxClusters)
  cCluster          <- matrix(1, nrow=nMCMC,   ncol=nMaxClusters)
  muPost            <- matrix(1, nrow=nMCMC,   ncol=nMaxClusters)
  sdPost            <- matrix(1, nrow=nMCMC,   ncol=nMaxClusters)
  gridPost          <- matrix(0, nrow=nMCMC,   ncol=length(grid))
    
  # Random initalisation
  I[1,] <- sample(nMaxClusters, nData, replace=TRUE)
  
  #############################################################################
  # Main MCMC loop
  #############################################################################
  
  for (ii in 1:nMCMC) {
    
    piCluster <- rep(0, nMaxClusters)
    rStick    <- 1;
    
    for (jj in 1:nMaxClusters) {
      
      #========================================================================
      # Update weights by stick-breaking
      #========================================================================

      # Count the number of members in each cluster
      nMembersCluster[ii,jj] <- length(which(I[ii,] == jj))
      
      # Count the number of occupied clusters
      if (nMembersCluster[ii,jj] > 0) { 
        nOccupiedClusters[ii] <- nOccupiedClusters[ii] + 1 
      }

      # Compute the proportion of the stick to break off
      nR = nData - sum(nMembersCluster[ii,1:jj])
      cCluster[jj] <- rbeta( 1, 1 + nMembersCluster[ii,jj], alpha[ii] + nR );  
      
      # If this is the last cluster, use the remainder of the stick
      if (jj == nMaxClusters) { 
        cCluster[jj] <- 1 
      };
      
      piCluster[jj] <- rStick * cCluster[jj];
      rStick        <- rStick * ( 1 - cCluster[jj] );

      #========================================================================
      # Sample the covariance for the component
      #========================================================================
      # Compute \nu_n
      nuN <- nu0 + nMembersCluster[ii,jj];
      
      # Compute \nu_n \sigma_n^2
      xbar  <- mean(y[which(I[ii,] == jj)])
      s2bar <- var(y[which(I[ii,] == jj)])
      
      nuNSigma2N <- nu0 * sd0^2 + (nMembersCluster[ii,jj] - 1) * s2bar + kappa0 * nMembersCluster[ii,jj] / (kappa0 + nMembersCluster[ii,jj]) * (xbar - mu0)^2
      
      # Compute \sigma_n^2
      sigma2N <- nuNSigma2N / nuN;
      
      # Sample \sigma^2 | y
      if (nMembersCluster[ii,jj] > 2) {
        # Compute the posterior
        sdPost[ii,jj] <- sqrt(nuN * sigma2N / rchisq(1, df=nuN))
      } else {
        # Sample from the prior
        sdPost[ii,jj] <- sqrt(nu0 * sd0 / rchisq(1, df=nu0))
      }
      
      #========================================================================
      # Sample the mean
      #========================================================================
      
      # Compute \mu_n
      muN <- kappa0 / (kappa0 + nMembersCluster[ii,jj]) * mu0 + nMembersCluster[ii,jj] / (kappa0 + nMembersCluster[ii,jj]) * mean(y[which(I[ii,] == jj)])
      
      # Copute \kappa_n
      kappaN <- kappa0 + nMembersCluster[ii,jj]
      
      # Sample \mu | y, \sigma^2
      if (nMembersCluster[ii,jj] > 2) {
        # Compute the posterior
        muPost[ii,jj] <- rnorm(1, muN, sdPost[ii,jj] / sqrt(kappaN))
      } else {
        # Sample from the prior
        muPost[ii,jj] <- rnorm(1, mu0, sqrt(sd0 / kappa0))
      }
      
    }
    
    #----------------------------------------------------------------------------------
    # Sample the indicators
    #----------------------------------------------------------------------------------
    piPost[ii,] <- piCluster / sum(piCluster)
    
    unnormw <- matrix(0, nrow = nData, ncol = nMaxClusters)
    for ( jj in 1:nMaxClusters ) {
      unnormw[,jj] <- piPost[ii,jj] * dnorm(y, mean=muPost[ii,jj], sd=sdPost[ii,jj])
    }
    
    for ( jj in 1:nData ) {
      I[ii+1,jj] <- sample(nMaxClusters, 1, replace=TRUE, 
                           unnormw[jj,] / sum(unnormw[jj,]))
    }

    #----------------------------------------------------------------------------------
    # Update alpha
    #----------------------------------------------------------------------------------    
    tmp         <- log(1 - cCluster[1:(nMaxClusters-1)]) 
    tmp[is.infinite(tmp)] <- 0 
       
    alpha[ii+1] <- rgamma(1, alpha0 + nMaxClusters - 1, rate = beta0 - sum( tmp ))
      
     if (alpha[ii+1] == 0) {
       print( alpha0 + nMaxClusters - 1 )
       print( cCluster[1:(nMaxClusters-1)] ) 
       readline("Press <return to continue") 
     }
   
    #----------------------------------------------------------------------------------
    # Compute density estimate
    #----------------------------------------------------------------------------------
    for ( jj in 1:nMaxClusters ) {
      gridPost[ii,] <- gridPost[ii,] + piPost[ii,jj] * dnorm(grid, muPost[ii,jj], sdPost[ii,jj])
    }

    #----------------------------------------------------------------------------------
    # Print iteration counter
    #----------------------------------------------------------------------------------    
    if (ii %% 500 == 0) {
      print(paste(paste(paste("iteration: ",ii,sep="")," of ", nMCMC, sep="")," complete.",sep=""))
    }
    
  }
  list(gridPost=gridPost, nOccupiedClusters=nOccupiedClusters, alpha=alpha)
}