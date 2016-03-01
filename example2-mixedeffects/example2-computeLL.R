calcLikelihood <- function( d, out, nIter, nBurnIn ) {
  ll <- matrix(0, nrow=nIter-nBurnIn, ncol=1)
  
  for ( kk in 1:length(ll) ) {
    for ( ii in 1:nIndividuals ) {
      jj <- kk + nBurnIn
      ll[kk] <- ll[kk] + sum( dnorm(d$y[ii,], out$alpha[jj+1] + out$betaInd[jj+1,ii] * d$x[1,ii,], sqrt(out$var[jj+1,ii]), log=TRUE ) )
    }
  }
  
  ll
}

# Computes the IACT values
iact <- function(x,xt,lag=0) {
  x = x - xt;
  xcorr = acf(x, lag.max = floor( length(x)/2 ), plot=F )
  lim   = 2/sqrt(length(x));
  
  if ( sum(is.nan(xcorr$acf)) > 0 ) {
    out = NA;
  } else {
    if (lag==0) {
      newli = min(which(xcorr$acf < lim)[1],floor( length(x) / 10),na.rm=T);
    } else {
      newli = lag
    }
    out = 1 + 2 * sum( xcorr$acf[1:newli] );
  }
  S = length(x);
  return(out);
}

calcIACTLL <- function(d, out, nIter, nBurnIn){
  ll <- calcLikelihood( d, out, nIter, nBurnIn )
  iact(ll, xt=mean(ll))
}
