###############################################################################
# Script for replicating example 3 in 
# "Bayesian inference for mixed effects models with heterogeneity"
#
# Running the inference
#
# (c) Johan Dahlin 2016 ( johan.dahlin (at) liu.se )
###############################################################################

#setwd("~/projects/dpm-panel2015/src-draft1/example3-reactiontest")

set.seed( 87655678 )
library("mvtnorm")

# Setup plot colors
library("RColorBrewer")
plotColors = brewer.pal(8, "Dark2");

source("samplers-mixedeffects-sleepmodel.R")

nIter         <- 10000
nBurnIn       <- 2500

###############################################################################
# Load data
###############################################################################
library("lme4")
str(sleepstudy)

nAttributes   <- 2
nIndividuals  <- 18
nMaxClusters  <- nIndividuals
nObservations <- 10

y      <- matrix( sleepstudy$Reaction, nrow=18, ncol=10, byrow = TRUE)
x      <- array( 0, dim=c(nAttributes,nIndividuals,nObservations) )
x[1,,] <- matrix( 1,               nrow=18, ncol=10, byrow = TRUE)
x[2,,] <- matrix( sleepstudy$Days, nrow=18, ncol=10, byrow = TRUE)
d <- list( x=x, y=y)

###############################################################################
# Analysis for each individual
###############################################################################

resLM <- matrix(0, nrow=nIndividuals, ncol=3)

for (ii in 1:nIndividuals) {
  res           <- lm( y[ii,] ~ x[2,ii,] )
  resLM[ii,1:2] <- res$coefficients
  resLM[ii,3]   <- var( res$residuals )
}

resLMmeans <- apply(resLM,2,mean)
resLMvars  <- apply(resLM,2,var)

###############################################################################
# Define priors
###############################################################################
a0star <- c(resLMmeans[1], resLMmeans[2], rep(0, 2*nMaxClusters))
A0star <- diag( c(0.01, 0.005, 0.02*rep(1, nAttributes*nMaxClusters ) ) )

# Wishart prior on inverse covariance
c0Q    <- 10
C0Q    <- 0.05 * diag(2)
#C0Q    <- 0.5 * diag(2)
( c0Q - (2+1)/2 )^(-1) * C0Q # mode of variance matrix

c0e    <- 1
C0e    <- 2

a0     <- 0.1
nu     <- 5

 prior = list( a0star=a0star, A0star=A0star, c0Q = c0Q, C0Q = C0Q,
               c0e=c0e, C0e=C0e, a0=a0, nu=nu)


###############################################################################
# Run the finite model with sparseness prior
###############################################################################

prior$eta0      <- 1/nMaxClusters
outFiniteSparse <- gibbs_finitemixture(d$y, d$x, d$x, nIter, nMaxClusters, 
                                       prior, postGrid)

###############################################################################
# Run the DPM model
###############################################################################

prior$alpha0 <- 0.01
prior$beta0  <- 0.01

outDPM       <- gibbs_dpm(d$y, d$x, d$x, nIter, nMaxClusters, prior, postGrid)

###############################################################################
# Plotting
###############################################################################

#cairo_pdf("~/projects/dpm-panel2015/paper/dpm-panel2015-draft1/figures/example3-reactiontest.pdf", width=8, height=8)

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))  
par(mar=c(4,5,1,1)) 

alpha0FM  <- density(outFiniteSparse$alpha[nBurnIn:nIter,1],from=220,to=260)
alpha1FM  <- density(outFiniteSparse$alpha[nBurnIn:nIter,2],from=-30,to=20)
alpha0DPM <- density(outDPM$alpha[nBurnIn:nIter,1],from=220,to=260)
alpha1DPM <- density(outDPM$alpha[nBurnIn:nIter,2],from=-30,to=20)

beta0FM  <- density(outFiniteSparse$betaInd[nBurnIn:nIter,,1],from=-20,to=30)
beta1FM  <- density(outFiniteSparse$betaInd[nBurnIn:nIter,,2],from=-10,to=40)
beta0DPM <- density(outDPM$betaInd[nBurnIn:nIter,,1],from=-20,to=30)
beta1DPM <- density(outDPM$betaInd[nBurnIn:nIter,,2],from=-10,to=40)

# Alpha_0
plot(alpha0FM$x,alpha0FM$y,col=plotColors[2], type="l", bty="n", ylab="density", xlab=expression(alpha[0]), main="", xlim=c(220,260), ylim=c(0,0.2), lwd=2 )
#polygon( c(alpha0FM$x,rev(alpha0FM$x)), c(alpha0FM$y, rep(0,length(alpha0FM$x))),border=NA,col=rgb(t(col2rgb(plotColors[1]))/256,alpha=0.25))
lines(alpha0DPM$x,alpha0DPM$y,col=plotColors[3], lwd=1.5)
#polygon( c(alpha0DPM$x,rev(alpha0DPM$x)), c(alpha0DPM$y, rep(0,length(alpha0DPM$x))),border=NA,col=rgb(t(col2rgb(plotColors[2]))/256,alpha=0.25))

# Alpha_1
plot(alpha1FM$x,alpha1FM$y,col=plotColors[2], type="l", bty="n", ylab="density", xlab=expression(alpha[1]), main="", xlim=c(-30,20), ylim=c(0,0.2), lwd=2 )
#polygon( c(alpha1FM$x,rev(alpha1FM$x)), c(alpha1FM$y, rep(0,length(alpha1FM$x))),border=NA,col=rgb(t(col2rgb(plotColors[1]))/256,alpha=0.25))
lines(alpha1DPM$x,alpha1DPM$y,col=plotColors[3], lwd=1.5)
#polygon( c(alpha1DPM$x,rev(alpha1DPM$x)), c(alpha1DPM$y, rep(0,length(alpha1DPM$x))),border=NA,col=rgb(t(col2rgb(plotColors[2]))/256,alpha=0.25))

# Beta_0
plot(beta0FM$x,beta0FM$y,col=plotColors[2], type="l", bty="n", ylab="density", xlab=expression(beta[i0]^{s}), main="", xlim=c(-20,30), ylim=c(0,0.1), lwd=2 )
#polygon( c(beta0FM$x,rev(beta0FM$x)), c(beta0FM$y, rep(0,length(beta0FM$x))),border=NA,col=rgb(t(col2rgb(plotColors[1]))/256,alpha=0.25))
lines(beta0DPM$x,beta0DPM$y,col=plotColors[3], lwd=1.5)
#polygon( c(beta0DPM$x,rev(beta0DPM$x)), c(beta0DPM$y, rep(0,length(beta0DPM$x))),border=NA,col=rgb(t(col2rgb(plotColors[2]))/256,alpha=0.25))

# Beta_1
plot(beta1FM$x,beta1FM$y,col=plotColors[2], type="l", bty="n", ylab="density", xlab=expression(beta[i1]^{s}), main="", xlim=c(-10,40), ylim=c(0,0.1), lwd=2 )
#polygon( c(beta1FM$x,rev(beta1FM$x)), c(beta1FM$y, rep(0,length(beta1FM$x))),border=NA,col=rgb(t(col2rgb(plotColors[1]))/256,alpha=0.25))
lines(beta1DPM$x,beta1DPM$y,col=plotColors[3], lwd=1.5)
#polygon( c(beta1DPM$x,rev(beta1DPM$x)), c(beta1DPM$y, rep(0,length(beta1DPM$x))),border=NA,col=rgb(t(col2rgb(plotColors[2]))/256,alpha=0.25))

#dev.off()

###############################################################################
# Compute posterior means
###############################################################################
colMeans(outFiniteSparse$alpha[nBurnIn:nIter,]) + colMeans(colMeans(outFiniteSparse$betaInd[nBurnIn:nIter,,]))
colMeans(outDPM$alpha[nBurnIn:nIter,]) + colMeans(colMeans(outDPM$betaInd[nBurnIn:nIter,,]))

###############################################################################
# End of file
###############################################################################