library(sn);library(mvtnorm); library(scatterplot3d); library(EnvStats); library(plotly)
library(mdatools); library(fitdistrplus); library(xtable); library(in2extRemes)
library(evd); library(SimCop); library(copula); library(mgpd); library(gplots)

source("fremantle.R.txt")
source("portpirie.R.txt")
slMerged <- merge(fremantle, portpirie, by = "Year")
sl <- data.frame(year = slMerged[,1], fremantle=slMerged[,2], Portpirie = slMerged[,4])
sl <- as.matrix(sl)

# E5.1 and E5.2 ====
# Ignore warnings
slNonparEst <- NonparEstDepFct(sl[,2:3], convex.hull = FALSE)
slNonparEst_c <- NonparEstDepFct(sl[,2:3], convex.hull = TRUE)
plot(slNonparEst, xlim = c(0,1), ylim = c(0.5,1), xlab = "w", 
     ylab = "A(w)", type = "b", cex = 0.7, main = "Pickands estimate of dependence function")
points(slNonparEst_c, type = "l", col = "blue")
polygon(c(0, 0.5, 1, 0), c(1, 0.5, 1, 1))

slSplfit <- SplineFitDepFct(slNonparEst)  #5.2 dep est func using cubic smoothing splines.
curve(slSplfit, n = 301, add = TRUE, lty ="dashed", col = 3)
legend("bottomright", legend=c("not convex", "convex", "cubic splines"), col=c("black", "blue", "green"), lty=1:3, cex= 0.8)

# E5.3 ==== # No MH algo, Use Smooshing splines and generate observations with GenerateRV
# New bivariate extreme value spine copula
slSplfitCop <- NewBEVSplineCopula(slSplfit)
slSplfitCopApprox <- GetApprox(slSplfitCop, type = 2)
plot(slSplfitCopApprox, xlab = expression(u[1]), ylab = expression(u[2]), col = heat.colors(c(-555,100)))
title("Heatmap of copula")

# E5.4 ====
n = 1000
BEVsc <- GenerateRV(slSplfitCopApprox, n)
# mar1_start <- c(1.50597, 0.113995, -0.107431)
# mar2_start <- c(3.866896, 0.197624, 0.007184)
gev_frem <- fgev(sl[,2], start = list(loc = 1.5, scale = 0.11, shape = -0.10))
gev_port <- fgev(sl[,3], start = list(loc = 3.8, scale = 0.2, shape = 0.0072))
l1 = gev_frem$e[1]; sc1 = gev_frem$e[2]; sh1 = gev_frem$e[3];
l2 = gev_port$e[1]; sc2 = gev_port$e[2]; sh2 = gev_port$e[3];

# Transform from copula to support [0,1] -> [0,inf]
BEVsc_frem <- qgev(BEVsc[,1], loc = l1, scale = sc1, shape = sh1)
BEVsc_port <- qgev(BEVsc[,2], loc = l2, scale = sc2, shape = sh2)

BEVsc_frem <- as.matrix(BEVsc_frem)
BEVsc_port <- as.matrix(BEVsc_port)

# Calculate probabilities
p1 = sum(BEVsc_frem > 1.7 & BEVsc_port > 4.2)/n
p2 = sum(BEVsc_frem > 1.8 & BEVsc_port > 4.4)/n
p3 = 1 - sum(BEVsc_frem < 1.478 & BEVsc_port < 3.85)/n
p = 1 - sum(BEVsc_frem < 1.95 & BEVsc_port < 4.8)/n
p4 = (p3 - p)/p3
#plot(BEVsc, col = "blue", cex = 0.4)
p1
p2
p3
p4