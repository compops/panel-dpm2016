###############################################################################
# Script for replicating example 3 in 
# "Bayesian inference for mixed effects models with heterogeneity"
#
# Plotting the data
#
# (c) Johan Dahlin 2016 ( johan.dahlin (at) liu.se )
###############################################################################

library("RColorBrewer")
plotColors = rep(brewer.pal(8, "Dark2"),4)
library("lme4")
str(sleepstudy)

###############################################################################
# Plot the data
###############################################################################

#cairo_pdf("~/projects/dpm-panel2015/paper/dpm-panel2015-draft1/figures/example3-data.pdf", width=8, height=10)
layout(matrix(1:18, 6, 3, byrow = FALSE)) 
par(mar=c(4,5,0,0))

for ( ii in 1:18 ) {
  id  <- as.numeric(names(table(sleepstudy$Subject))[ii])
  d   <- sleepstudy[which(sleepstudy$Subject == id ),]
  res <- lm(d$Reaction~d$Days)
  
  if ( ii < 6 ) {
    plot(d$Days, d$Reaction, type="p", pch=19, col="darkgrey",  
         ylim=c(100,500), bty="n", ylab = "reaction time", xlab="", xlim=c(0,10) )
    lines( d$Days, res$coefficients[1] + res$coefficients[2] * d$Days, 
           col="black")
  } 

  if ( ii == 6 ) {
    plot(d$Days, d$Reaction, type="p", pch=19, col="darkgrey", 
         ylim=c(100,500), bty="n", xlab = "day", xlim=c(0,10),
         ylab = "reaction time" )
    lines( d$Days, res$coefficients[1] + res$coefficients[2] * d$Days, 
           col="black") 
  }
    
  if ( ( ii > 6 ) && ( ii %% 6 != 0 ) ) {
    plot(d$Days, d$Reaction, type="p", pch=19, col="darkgrey", 
         ylim=c(100,500), bty="n", xlab = "", ylab = "", xlim=c(0,10))
  lines( d$Days, res$coefficients[1] + res$coefficients[2] * d$Days, 
         col="black") 
  }

  if ( ( ii > 6 ) && ( ii %% 6 == 0 ) ) {
    plot(d$Days, d$Reaction, type="p", pch=19, col="darkgrey", 
         ylim=c(100,500), bty="n", xlab = "day",ylab="", xlim=c(0,10))
    lines( d$Days, res$coefficients[1] + res$coefficients[2] * d$Days, 
           col="black") 
  }  
  
  text(10, 150, pos=2, labels=id)
}
#dev.off()

###############################################################################
# End of file
###############################################################################