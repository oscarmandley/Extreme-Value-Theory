library(sn);library(mvtnorm); library(scatterplot3d); library(EnvStats); library(plotly)
library(mdatools); library(fitdistrplus); library(xtable); library(in2extRemes)
library(evd); library(SimCop); library(copula); library(mgpd)

source("fremantle.R.txt")
source("portpirie.R.txt")
slMerged <- merge(fremantle, portpirie, by = "Year")
sl <- data.frame(year = slMerged[,1], fremantle=slMerged[,2], Portpirie = slMerged[,4])
sl = as.matrix(sl)

# E6.1 ====
plot(sl[,1], sl[,2])
u1 <- quantile(sl[,2], 0.3)
u2 <- quantile(sl[,3], 0.3)
thr = c(u1,u2)
nms = c("Fremantle","Portpirie")

# E6.2 ====
pot1 <- fbvpot(sl[,2:3], thr, model = "log")
pot2 <- fbvpot(sl[,2:3], thr, model = "aneglog", std.err = FALSE)

# E6.3 ====

# par(mfrow = c(1,2)); par(mar = rep(2, 4))
# n = 63
# bvdata <- rbvevd(n, alpha = a, beta = b, model = "bilog", mar1 = mar1, mar2 = mar2)
# M1 <- fbvevd(bvdata, model = "bilog")
# plot(M1, which = 5, p = c(0.75, 0.9, 0.95), col = "blue", xlim = c(1.2,2), ylim = c(3.5,5))
plot(pot1, which = 3, p = c(0.75,0.9, 0.95), col = "1", xlim = c(1.2,2), ylim = c(3.5,5))
# plot(pot2, which = 3, p = c(0.75,0.9, 0.95), col = "2")

# E6.4 ====
gev_frem <- fgev(sl[,2], start = list(loc = 1.5, scale = 0.11, shape = -0.10))
gev_port <- fgev(sl[,3], start = list(loc = 3.8, scale = 0.2, shape = 0.0072))
l1 = gev_frem$e[1]; sc1 = gev_frem$e[2]; sh1 = gev_frem$e[3];
l2 = gev_port$e[1]; sc2 = gev_port$e[2]; sh2 = gev_port$e[3];

xy_tilde <- function(X, u, eta, sc, sh){
  xy_tilde = -(log(1 - eta*((1 + sh*(X-u)/sc)^(-1/sh))))^-1
}
X = 1.478
X_2 = 1.95
X_3 = 2
u = thr[1]
eta = sum(sl[,2] > thr[1])/length(sl[,2])
sc = pot1$e[1]
sh = pot1$e[2]
L1 = abs(sc/sh)
L1
X_tilde <- xy_tilde(X, u, eta, sc, sh)
X_tilde_2 <- xy_tilde(X_2, u, eta, sc, sh)
X_tilde_3 <- xy_tilde(X_3, u, eta, sc, sh)

Y = 3.85
Y_2 = 4.8
Y_3 = 5
u = thr[2]
eta = sum(sl[,3] > thr[2])/length(sl[,3])
sc = pot1$e[3]
sh = pot1$e[4]
L2 = abs(sc/sh)
L2
Y_tilde <- xy_tilde(Y, u, eta, sc, sh)
Y_tilde_2 <- xy_tilde(Y_2, u, eta, sc, sh)
Y_tilde_3 <- xy_tilde(Y_3, u, eta, sc, sh)

r = pot1$e[5]

p14_38 = 1 - exp(-((X_tilde^(-1/r) + Y_tilde^(-1/r))^r))
p14_38

p19_48 = 1 - exp(-((X_tilde_2^(-1/r) + Y_tilde_2^(-1/r))^r))
p19_48

p4 = (p14_38 - p19_48)/p14_38
p4

# E6.5 ====
p2_5 = 1 - exp(-((X_tilde_3^(-1/r) + Y_tilde_3^(-1/r))^r))
p2_5