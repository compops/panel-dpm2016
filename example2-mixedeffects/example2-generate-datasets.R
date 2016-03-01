setwd("~/projects/dpm-panel2015/src-draft1/example2-mixedeffects")

set.seed( 87655678 )
library("mvtnorm")

# Setup plot colors
library("RColorBrewer")
plotColors = brewer.pal(8, "Dark2");

source("samplers-mixedeffects-1dim.R")

nAttributes   <- 1
nMaxClusters  <- 20
nObservations <- 20
nIter         <- 10000
nBurnIn       <- 2500
postGrid      <- seq(-4,4,0.1)


nIndividualsVector <- c( 10, 20, 50, 100, 200, 500)
nRep               <- 40

for (ii in 1:length(nIndividualsVector)) {
  nIndividuals <- nIndividualsVector[ii]
  
  for (jj in 1:nRep) {
  
  ###############################################################################
  # Generate data
  ###############################################################################
  d <- generateData( beta=c( 2, 3, -2), Q=c( 1, 0.2, 0.2), w=c(0.4, 0.3, 0.3),
                     s2e=1, alpha=2, omega=abs(rnorm(nIndividuals)),
                     dims=c(nIndividuals,nAttributes,nObservations))
  
  
  save(d,file=paste(paste(paste(paste("data/mixedeffects-nInd-",nIndividuals,sep=""),"-",sep=""),jj,sep=""),".RData",sep=""))
  }
}