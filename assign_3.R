library(in2extRemes)
in2extRemes()
par(mar=c(1,1,1,1))

#time = (fremantle$data$Year-min(fremantle$data$Year))/
#(max(fremantle$data$Year)-min(fremantle$data$Year))
#timeSq = time^2
#fremantle$data <- cbind(fremantle$data,time=time, timeSq=time^2)
#names(fremantle$data)

# fit 1 M1
# 1.4823417  0.1412723 -0.2174282 
# ste 0.01672527 0.01149706 0.06378114  (standard error)

# fit 2 M2
#        mu0        mu1       phi0       phi1      shape 
# 1.3918440  0.1707783 -1.9200454 -0.3270439 -0.1362358 

# fit 3 M2
#       mu0        mu1        mu2       phi0       phi1      shape 
# 1.3388206  0.4421351 -0.2530953 -1.9238939 -0.3657534 -0.1136636 

# fit 4 M2
#        mu0         mu1         mu2        phi0        phi1        phi2       shape 
# 1.40186193  0.17541020  0.06591733 -1.89868529 -0.41994918  0.26517960 -0.22291118 

# fit 5 M2
#       mu0        mu1        mu2       phi0       phi1      shape 
# 1.3958835  0.1808483  0.0642604 -2.1125622  0.2726217 -0.1879142 

# fit 6 M2
#         mu0         mu1         mu2       scale       shape 
# 1.38432969  0.19448974  0.05451703  0.12073287 -0.14998217 






x = 1:length(wooster$data[,2])
usin = function(x, a, b, d)
{
  a + b * sin(((x - d) * 2 * pi)/365.25)
}
wu = usin(x, -30, 25, -75)
winter = c(rep(c(rep(1, 61), rep(0, 273), rep(1, 31)), 5),
           1)
spring = c(rep(c(rep(0, 61), rep(1, 91), rep(0, 365 - 91 -
                                               61)), 5), 0)
summer = c(rep(c(rep(0, 61 + 91), rep(1, 91), rep(0, 365 -
                                                    91 - 61 - 91)), 5), 0)
fall = c(rep(c(rep(0, 61 + 91 + 91), rep(1, 91), rep(0, 365 -
                                                       91 - 61 - 91 - 91)), 5), 0)
rescale.covariate  = function(x)
{
  r.x = range(x)
  x.01 = (x-r.x[1])/diff(r.x)
  2*x.01-1
}
ydat = cbind(wu, sin((x * 2 * pi)/365.25), cos((x * 2 * pi)/365.25
), rescale.covariate(x), winter, spring, summer, fall)

wooster$data[,"value"] <- - wooster$data[,"value"]
wooster$data <- cbind(wooster$data,ydat)
colnames(wooster$data)[3:6] <- c("wu", "sin","cos","time")
wooster$data <- as.data.frame(wooster$data)




# Model 4
# -15.42712901   8.74783377  29.04479227   0.51210700  -0.19519330   0.92047547  -0.32958507  -0.09453393 
# xi2 0.15431136 

# Model 5
#mu0          mu1          mu2          mu3         phi0         phi1         phi2         phi3 
#-15.44363305   8.59670645  28.34965227  -2.26448879   0.52020674   0.03718099   0.49843247  -0.10816880 
#shape 
# -0.34335820


# Model 6
#          mu0          mu1          mu2         phi0         phi1         phi2          xi0          xi1 
# -15.61988030   9.18171968  28.77401674   0.45273339   0.04943602   0.82143284  -0.28963811   0.06582288 
# xi2          xi3          xi4 
# -0.13253971  -0.17384724  -0.04907405 