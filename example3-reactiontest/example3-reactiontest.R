setwd("~/projects/dpm-panel2015/src-draft1/example3-reactiontest")

set.seed( 87655678 )
library("mvtnorm")

# Setup plot colors
library("RColorBrewer")
plotColors = brewer.pal(8, "Dark2");

source("samplers-mixedeffects-1dim.R")

nIter         <- 10000
nBurnIn       <- 2500
postGrid      <- seq(7,14,0.05)

###############################################################################
# Load data
###############################################################################
library("lme4")
str(sleepstudy)

nAttributes   <- 1
nIndividuals  <- 18
nMaxClusters  <- nIndividuals
nObservations <- 10

y      <- matrix( sleepstudy$Reaction, nrow=18, ncol=10, byrow = TRUE)
x      <- array( 0, dim=c(nAttributes,nIndividuals,nObservations) )
x[1,,] <- matrix( sleepstudy$Days, nrow=18, ncol=10, byrow = TRUE)
d <- list( x=x, y=y)

###############################################################################
# Analysis for each individual
###############################################################################

resLM <- matrix(0, nrow=nIndividuals, ncol=3)

for (ii in 1:nIndividuals) {
  res           <- lm( y[ii,] ~ x[,ii,] )
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
C0Q    <- 0.02

c0e    <- 0
C0e    <- 0

a0     <- 0.1
nu     <- 5

prior = list( a0star=a0star, A0star=A0star, c0Q = c0Q, C0Q = C0Q,
              c0e=c0e, C0e=C0e, a0=a0, nu=nu)


###############################################################################
# Run the finite model with sparseness prior
###############################################################################
source("samplers-mixedeffects-1dim.R")
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
# Plotting
###############################################################################

layout(matrix(c(1,1,1,2,2,3,4,5,6), 3, 3, byrow = TRUE))  
par(mar=c(4,5,1,1)) 

plot(postGrid, colMeans(outFiniteSparse$posterior_beta[nBurnIn:nIter,]), lwd=3, col=plotColors[2], type="l", bty="n", ylab="density", xlab=expression(beta), ylim=c(0,2),xlim=c(0,20))
lines(postGrid, colMeans(outDPM$posterior_beta[nBurnIn:nIter,]), lwd=3, col=plotColors[3])
legend("topleft",c("Sparseness prior","DP prior"),col=plotColors[2:3],lwd=2,box.col=NA)

hist(outFiniteSparse$alpha[nBurnIn:nIter],main="",breaks=floor(sqrt(nIter)),freq=FALSE,xlab=expression(alpha),col=rgb(t(col2rgb(plotColors[5]))/256,alpha=0.25),border=NA, ylab="density",xlim=c(248,254))
lines(density(outFiniteSparse$alpha[nBurnIn:nIter]),lwd=3,col=plotColors[5])
hist(outFiniteSparse$nOccupiedClusters[nBurnIn:nIter],breaks=floor(sqrt(nIter)),main="",freq=FALSE,xlab="no. occupied clusters",col="darkgrey",border=NA,xlim=c(0,nMaxClusters), ylab="density")

hist(outDPM$alpha[nBurnIn:nIter],main="",breaks=floor(sqrt(nIter)),freq=FALSE,xlab=expression(alpha),col=rgb(t(col2rgb(plotColors[6]))/256,alpha=0.25),border=NA, ylab="density",xlim=c(248,254))
lines(density(outDPM$alpha[nBurnIn:nIter]),lwd=3,col=plotColors[6])
hist(outDPM$alphaSB[nBurnIn:nIter],main="",breaks=floor(sqrt(nIter)),freq=FALSE,xlab=expression(alpha[DP]),col=rgb(t(col2rgb(plotColors[7]))/256,alpha=0.25),border=NA, ylab="density",xlim=c(0,30))
lines(density(outDPM$alphaSB[nBurnIn:nIter]),lwd=3,col=plotColors[7])
hist(outDPM$nOccupiedClusters[nBurnIn:nIter],breaks=floor(sqrt(nIter)),main="",freq=FALSE,xlab="no. occupied clusters",col="darkgrey",border=NA,xlim=c(0,nMaxClusters), ylab="density")