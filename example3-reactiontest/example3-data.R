library("RColorBrewer")
plotColors = brewer.pal(8, "Dark2");

str(sleepstudy)
require(lattice)

cairo_pdf("~/projects/dpm-panel2015/paper/dpm-panel2015-draft1/figures/example3-data.pdf", width=10, height=8)
layout(matrix(1, 1, 1, byrow = TRUE)) 
par(mar=c(4,5,1,1))

xyplot(Reaction ~ Days | Subject, sleepstudy, type = c("g","p","r"),
       index = function(x,y) coef(lm(y ~ x))[1],
       xlab = "days of sleep deprivation",
       ylab = "average reaction time (ms)", aspect = "xy", 
       pch=19, cex=0.5, lwd=2, col=rep(plotColors))
dev.off()