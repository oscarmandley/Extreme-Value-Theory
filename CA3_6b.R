library(sn);library(mvtnorm); library(scatterplot3d); library(EnvStats); library(plotly)
library(mdatools); library(fitdistrplus); library(xtable); library(in2extRemes)
library(evd); library(SimCop); library(copula); library(mgpd)

source("fremantle.R.txt")
source("portpirie.R.txt")
slMerged <- merge(fremantle, portpirie, by = "Year")
sl <- data.frame(year = slMerged[,1], fremantle=slMerged[,2], Portpirie = slMerged[,4])
sl = as.matrix(sl)

plot(sl[,1], sl[,2])
u1 <- quantile(sl[,2], 0.3)
u2 <- quantile(sl[,3], 0.3)
thr = c(u1,u2)
nms = c("Fremantle","Portpirie")

potdata = mgpd_data(sl[,2:3], thr = thr)
plot( potdata[,1], potdata[,2], xlab=nms[1], ylab=nms[2], main="Exceedances" )
abline( h=0, v=0 , lty=2)

# init          = mgpd_init( potdata )
est.log       = fbgpd( c(0,1,-0.1,0,1,-0.1,1.2), dat=potdata, model="log", fixed=FALSE )
est.log$par

# est.amix = fbgpd(c(est.log$par,1), dat = potdata, model = "alog", fixed = FALSE)
# est.amix$par

pot1 <- pbgpd(1.95, 4.8, model = "log", mar1 = est.log$par[1:3],
              mar2 = est.log$par[4:6], dep = est.log$par[7])
pot1

# domain and probabilities for the pred regs - PRECISION COULD BE BETTER 0.001 or even finer
x       = seq(-0.22,1,  0.003)
y       = seq(-0.27,1,  0.003)
Q       = c(0.95,0.9,0.75)
# evaluating the density over the grid
z       = outer(  x,  y,  dbgpd,  model="log", mar1 = est.log$par[1:3], mar2 = est.log$par[4:6],
                  dep = est.log$par[7] )
# calculating the prediction region
reg     = dbgpd_region( x,  y,  z, quant = Q )

# plotting on the original scale (threshold levels are added back)
contour(reg$x+thr[1],  reg$y+thr[2],  reg$z,  levels=reg$q, drawlabels=FALSE, main="Logistic BGPD",
        col=c(1,3,4), xlab=paste(nms[1],"(m)"),ylab=paste(nms[2],"(m)"),lwd=1)
abline( h=thr[2],  v=thr[1],  lty=2 )
legend( "bottomright", c(expression(gamma==0.95),  expression(gamma==0.9), expression(gamma==0.75)),
        lty=1,  col=c(1,3,4), title="Regions", bty="n")
# Observations
points( potdata[,1]+thr[1],potdata[,2]+thr[2],cex=0.7)