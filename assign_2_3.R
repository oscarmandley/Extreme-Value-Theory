#This file gives examples of fitting bivariate GPD distributions using
#mgpd package in R (authored by Pal Rakonczai).

library(mgpd)

# The following is sea level data for Fremantle and Portpirie
fp <-
  structure(c(1.49, 1.46, 1.34, 1.74, 1.62, 1.46, 1.71, 1.74, 1.55,
              1.43, 1.62, 1.49, 1.58, 1.34, 1.37, 1.62, 1.31, 1.43, 1.49, 1.55,
              1.71, 1.49, 1.46, 1.52, 1.58, 1.65, 1.49, 1.52, 1.52, 1.49, 1.62,
              1.86, 1.58, 1.62, 1.46, 1.43, 1.46, 1.62, 1.68, 1.83, 1.62, 1.46,
              1.58, 1.77, 1.62, 1.71, 1.46, 1.6, 1.5, 1.6, 1.9, 1.7, 1.4, 1.8,
              1.37, 1.46, 1.61, 1.43, 1.67, 1.62, 1.57, 1.56, 1.46, 4.03, 3.83,
              3.65, 4.01, 4.08, 4.18, 3.8, 4.36, 3.96, 3.98, 4.69, 3.85, 3.96,
              3.85, 3.93, 3.75, 3.63, 3.57, 3.97, 4.05, 4.24, 4.22, 3.73, 4.37,
              4.06, 3.71, 3.96, 4.06, 4.55, 3.79, 3.89, 4.11, 3.85, 3.86, 3.86,
              4.21, 4.01, 4.11, 4.24, 3.96, 4.21, 3.74, 3.85, 3.88, 3.66, 4.11,
              3.71, 4.18, 3.9, 3.78, 3.91, 3.72, 4, 3.66, 3.62, 4.33, 4.55,
              3.75, 4.08, 3.9, 3.88, 3.94, 4.33), .Dim = c(63L, 2L))

# we must give the threshold levels first
thr            = apply(fp,2,quantile,prob=0.3)
nms            = c("Fremantle","Portpirie")
# here we convert the data to the desired form
potdata        = mgpd_data( fp, thr=thr )
plot( potdata[,1], potdata[,2], xlab=nms[1], ylab=nms[2], main="Exceedances" )
abline( h=0, v=0 , lty=2)


#===============================================================================
#fitting logistic model, (here you can also try other models; see the
#help page for fbgpd).
init          = mgpd_init( potdata )
# est.log       = fbgpd(c(init,1.2),dat=potdata,model="log",fixed=FALSE)
# default initial value (above) doesn't work so we set them manually:
est.log       = fbgpd( c(0,1,-0.1,0,1,-0.1,1.2), dat=potdata, model="log", fixed=FALSE )

## Note that the estimates of the parameters are as following: 1)
## estimates of GEV parameters for margin 1 (location, scale and shape
## parameters), 2) estimates of GEV parameters for margin 2 (location,
## scale and shape parameters) and 3) estimates the dependence
## parameters in the model.
est.log$par
# domain and probabilities for the pred regs - PRECISION COULD BE BETTER 0.001 or even finer
x       = seq(-0.22,1,  0.005)
y       = seq(-0.27,1,  0.005)
Q       = c(0.95,0.9,0.75)
# evaluating the density over the grid
z       = outer(  x,  y,  dbgpd,  model="log", mar1 = est.log$par[1:3], mar2 = est.log$par[4:6], dep = est.log$par[7] )
# calculating the prediction region
reg     = dbgpd_region( x,  y,  z, quant = Q )
# plotting on the original scale (threshold levels are added back)

contour(reg$x+thr[1],  reg$y+thr[2],  reg$z,  levels=reg$q, drawlabels=FALSE, main="Logistic BGPD",  col=c(1,3,4), xlab=paste(nms[1],"(m)"),ylab=paste(nms[2],"(m)"),lwd=1)
abline( h=thr[2],  v=thr[1],  lty=2 )
legend( "bottomright", c(expression(gamma==0.95),  expression(gamma==0.9), expression(gamma==0.75)), lty=1,  col=c(1,3,4), title="Regions", bty="n")
# and the obs
points( potdata[,1]+thr[1],potdata[,2]+thr[2],cex=0.7)

# Note that for the following you need to download the zip-file
# pbgpd.zip and source them to R.
library(fields)
z       = outer(  x,  y,  pbgpd,  model="log", mar1 = est.log$par[1:3], mar2 = est.log$par[4:6], dep = est.log$par[7] )
image.plot(x+thr[1],  y+thr[2],  z, main="Quantiles of Logistic BGPD", xlab=paste(nms[1],"(m)"),ylab=paste(nms[2],"(m)"),lwd=1)
contour(x+thr[1],  y+thr[2],     z, levels=Q, add=T)
abline( h=thr[2],  v=thr[1],  lty=2 )



#par(mfrow=c(1,2))
persp(x,y,z)
image(x,y,z)

## Calculating probabilities: note that we need to subtract the thresholds
> prob.log <- pbgpd(2-1.478,5-3.850,model="log",mar1=est.log$par[1:3],mar2=est.log$par[4:6],dep=est.log$par[7])

## Note that the following probability is 0
> pbgpd(0,0,model="log",mar1=est.log$par[1:3],mar2=est.log$par[4:6],dep=est.log$par[7])



#finally this one looks fancy if you have the fields package installed
#library(fields)
#image.plot(x,y,z)

#Tajvidi's generalized logistic model is also implemented, but for
#these data the extra parameter doesnt seem necessary
est.taj       = fbgpd( c(est.log$par[1:7],0.1), dat=potdata, model="taj", fixed=FALSE )
#likelihoods
est.log$value
est.taj$value
#parameters
est.log$par
est.taj$par

#===============================================================================
# TESTING OTHER DISTRIBUTION FUNCTIONS
#neglog modell
z       = outer(  x,  y,  pbgpd,  model="neglog", mar1 = est.log$par[1:3], mar2 = est.log$par[4:6], dep = 1.5 )
image.plot(x+thr[1],  y+thr[2],  z)
contour(x+thr[1],  y+thr[2],     z,add=T)
abline( h=thr[2],  v=thr[1],  lty=2 )

#psi modell
z       = outer(  x,  y,  pbgpd,  model="log", mar1 = est.log$par[1:3], mar2 = est.log$par[4:6], dep = est.log$par[7] )
image.plot(x+thr[1],  y+thr[2],  z)
contour(x+thr[1],  y+thr[2],     z,add=T)
z       = outer(  x,  y,  pbgpd,  model="psilog", mar1 = est.log$par[1:3], mar2 = est.log$par[4:6], dep = est.log$par[7], asy=0.3 )
contour(x+thr[1],  y+thr[2],     z,add=T)
abline( h=thr[2],  v=thr[1],  lty=2 )

#psi modell
z       = outer(  x,  y,  pbgpd,  model="neglog", mar1 = est.log$par[1:3], mar2 = est.log$par[4:6], dep = 1.2 )
image.plot(x+thr[1],  y+thr[2],  z)
contour(x+thr[1],  y+thr[2],     z,add=T)
z       = outer(  x,  y,  pbgpd,  model="psineglog", mar1 = est.log$par[1:3], mar2 = est.log$par[4:6], dep = 1.2, asy=0.5, p=4 )
contour(x+thr[1],  y+thr[2],     z,add=T)
abline( h=thr[2],  v=thr[1],  lty=2 )

#phi modell
z       = outer(  x,  y,  pbgpd,  model="log", mar1 = est.log$par[1:3], mar2 = est.log$par[4:6], dep = est.log$par[7] )
image.plot(x+thr[1],  y+thr[2],  z)
contour(x+thr[1],  y+thr[2],     z,add=T)
z       = outer(  x,  y,  pbgpd,  model="philog", mar1 = est.log$par[1:3], mar2 = est.log$par[4:6], dep = 1.5, asy=0.001,p=3 )
contour(x+thr[1],  y+thr[2],     z,add=T)
abline( h=thr[2],  v=thr[1],  lty=2 )

#phi modell
z       = outer(  x,  y,  pbgpd,  model="log", mar1 = est.log$par[1:3], mar2 = est.log$par[4:6], dep = 1.2 )
image.plot(x+thr[1],  y+thr[2],  z)
contour(x+thr[1],  y+thr[2],     z,add=T)
z       = outer(  x,  y,  pbgpd,  model="phineglog", mar1 = est.log$par[1:3], mar2 = est.log$par[4:6], dep = 1.2, asy=0.001,p=3 )
contour(x+thr[1],  y+thr[2],     z,add=T)
abline( h=thr[2],  v=thr[1],  lty=2 )

