# Library ====
library(sn);library(mvtnorm); library(scatterplot3d); library(copula); library(EnvStats); 
library(plotly); library(mdatools); library(evd); library(SimCop); library(fitdistrplus)
library(SpatialExtremes)
in_dat <- read.table("http://www.maths.lth.se/matstat/kurser/fmsn15masm23/datasets/insuranceData.txt",
                     header = T)
par(mfrow = c(2,2))
par(mar = rep(2, 4))
m1 =  min(in_dat[,1])
m2 = min(in_dat[,2])
ma1 = max(in_dat[,1])
ma2 = max(in_dat[,2])
plot(in_dat, log = "xy", xlim = c(m1, ma1), ylim = c(m2, ma2), main = "Data")

# Correlation ====
cor(in_dat, method = "pearson")
cor(in_dat, method = "kendall")
cor(in_dat, method = "spearman")

# Initiation  ====
data <- as.matrix(in_dat[c("Loss", "ALAE")])
n = dim(data)[1]


# lnorm seems good

flnorm <- fitdist(data[,1], "lnorm")

# log data sets
aa = data[,1]
fcauchy <- fitdist(aa, "cauchy")
flogis <- fitdist(aa, "logis")
funif <- fitdist(aa, "unif")
fweibull <- fitdist(aa, "weibull")
fgpd <- fitdist(data[,1], "gpd")
fnorm <- fitdist(data[,1], "norm")

g <- gofstat(list( fcauchy, flogis, funif, fweibull), fitnames = c("cauchy", "logis", "unif", "weibull"))
g$kstest

par(mfrow = c(2,2), mar = c(4, 4, 1, 1))
denscomp(flnorm, addlegend = FALSE, main = "", xlab = "data", fitcol = "orange")
qqcomp(flnorm, addlegend = FALSE, main = "", fitpch = 16, fitcol = "grey", line01lty = 2)
cdfcomp(flnorm, addlegend = FALSE, main = "", xlab = "data", fitcol = "orange", lines01 = TRUE)
ppcomp(flnorm, addlegend = FALSE, main = "", fitpch = 16, fitcol = "grey", line01lty = 2)


# plotdist(data[,2], distr = "norm", para = list(mean = 0, sd = 2))


