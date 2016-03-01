setwd("~/20160227-dpm-panel2015/example2-mixedeffects")
source("samplers-mixedeffects-1dim.R")
library("mvtnorm")

set.seed( 87655678 )

###############################################################################
# Handle input arguments
###############################################################################
args=(commandArgs(TRUE))

if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

###############################################################################
# Define settings
###############################################################################

nAttributes   <- 1
nMaxClusters  <- 20
nObservations <- 20
nIter         <- 10000
nBurnIn       <- 2500
postGrid      <- seq(-4,4,0.1)

nIndividualsVector <- c( 10, 20, 50, 100, 200, 500)

for ( ii in 1:length(nIndividualsVector) ) {
  
  nIndividuals  <- nIndividualsVector[ii]
  ###############################################################################
  # Load data
  ###############################################################################
  
  load(file=paste(paste(paste(paste("data/mixedeffects-nInd-",nIndividuals,sep=""),"-",sep=""),jj,sep=""),".RData",sep=""))
  
  ###############################################################################
  # Analysis for each individual
  ###############################################################################
  
  resLM <- matrix(0, nrow=nIndividuals, ncol=3)
  
  for (ii in 1:nIndividuals) {
    res           <- lm( d$y[ii,] ~ d$x[,ii,] )
    resLM[ii,1:2] <- res$coefficients
    resLM[ii,3]   <- var( res$residuals )
  }
  
  resLMmeans <- apply(resLM,2,mean)
  resLMvars  <- apply(resLM,2,var)
  
  ###############################################################################
  # Define priors
  ###############################################################################
  a0star <- c(resLMmeans[1], rep(resLMmeans[2], nMaxClusters))
  A0star <- diag(rep(1, 1+nAttributes*nMaxClusters))
  
  c0Q    <- 10
  C0Q    <- 1
  c(C0Q/(c0Q+1),sqrt(C0Q^2*(c0Q-1)^(-2)*(c0Q-2)^(-1)),(C0Q/(1+c0Q)))
  
  c0e    <- 0
  C0e    <- 0
  
  a0     <- 0.1
  nu     <- 5
  
  prior = list( a0star=a0star, A0star=A0star, c0Q = c0Q, C0Q = C0Q,
                c0e=c0e, C0e=C0e, a0=a0, nu=nu)
  
  
  ###############################################################################
  # Run the finite model with known number of clusters
  ###############################################################################
  
  prior$a0star    <- c(resLMmeans[1], rep(resLMmeans[2], 3))
  prior$A0star    <- diag(rep(1, 1+3))
  prior$eta0      <- 1
  outFiniteKnown  <- gibbs_finitemixture(d$y, d$x, nIter, 3, prior, postGrid)
  
  ###############################################################################
  # Run the finite model with sparseness prior
  ###############################################################################
  prior$a0star    <- c(resLMmeans[1], rep(resLMmeans[2], nMaxClusters))
  prior$A0star    <- diag(rep(1, 1+nAttributes*nMaxClusters))
  prior$eta0      <- 1/nMaxClusters
  
  outFiniteSparse <- gibbs_finitemixture(d$y, d$x, nIter, nMaxClusters, 
                                         prior, postGrid)
  
  ###############################################################################
  # Run the DPM model
  ###############################################################################
  
  prior$alpha0 <- 0.01
  prior$beta0  <- 0.01
  
  outDPM       <- gibbs_dpm(d$y, d$x, nIter, nMaxClusters, prior, postGrid)
  
  ###############################################################################
  # Save output
  ###############################################################################
  save(c(outFiniteKnown,outFiniteSparse,outDPM),file=paste(paste("results/results-mixedeffects-run-",jj,sep=""),".RData",sep=""))
  
}