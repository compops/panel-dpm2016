###############################################################################
# Subroutine for replicating example 1 in 
# "Bayesian inference for mixed effects models with heterogeneity"
#
# Gibbs sampling for finite mixture model and DPM model
#
# (c) Johan Dahlin 2016 ( johan.dahlin (at) liu.se )
###############################################################################

###############################################################################
# Gibbs sampling for finite mixture
###############################################################################

gibbs_finitemixture <- function(y, nMCMC, nMaxClusters, prior, grid)  {
  
  nData <- length(y)

  # Hyperparameters in priors
  mu0    <- prior$mu0
  sd0    <- prior$sd0
  kappa0 <- prior$kappa0
  nu0    <- prior$nu0
  a0     <- prior$a0

  # Pre-allocate matrices
  I                 <- matrix(1, nrow=nMCMC+1, ncol=nData)
  nMembersCluster   <- matrix(1, nrow=nMCMC,   ncol=nMaxClusters)
  nOccupiedClusters <- matrix(0, nrow=nMCMC,   ncol=1)
  piPost            <- matrix(1, nrow=nMCMC,   ncol=nMaxClusters)
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
    
    for (jj in 1:nMaxClusters) {
      
      # Count the number of members in each cluster
      nMembersCluster[ii,jj] <- length(which(I[ii,] == jj))
      
      # Count the number of occupied clusters
      if (nMembersCluster[ii,jj] > 0) { 
        nOccupiedClusters[ii] <- nOccupiedClusters[ii] + 1 
      }
      
      # Sample probabilities p from Dirichlet (as normalised Gamma variables)
      piCluster[jj] <- rgamma( 1, a0 + nMembersCluster[ii,jj], 1 )
      
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
      muN <- kappa0 / (kappa0 + nMembersCluster[ii,jj]) * mu0 + nMembersCluster[ii,jj] / (kappa0 + nMembersCluster[ii,jj]) * xbar
      
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
  list(gridPost=gridPost, nOccupiedClusters=nOccupiedClusters)
}

###############################################################################
# Collapsed Gibbs sampling for DPM process
###############################################################################

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
      nR              <- nData - sum(nMembersCluster[ii,1:jj])
      cCluster[ii,jj] <- rbeta( 1, 1 + nMembersCluster[ii,jj], alpha[ii] + nR );  
      
      # If this is the last cluster, use the remainder of the stick
      if (jj == nMaxClusters) { 
        cCluster[jj] <- 1 
      };
      
      piCluster[jj] <- rStick * cCluster[ii,jj];
      rStick        <- rStick * ( 1 - cCluster[ii,jj] );

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
      muN <- kappa0 / (kappa0 + nMembersCluster[ii,jj]) * mu0 + nMembersCluster[ii,jj] / (kappa0 + nMembersCluster[ii,jj]) * xbar
      
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
    tmp         <- log(1 - cCluster[ii,1:(nMaxClusters-1)]) 
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

###############################################################################
# End of file
###############################################################################