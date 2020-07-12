library(sn);library(mvtnorm); library(scatterplot3d); library(EnvStats); library(plotly)
library(mdatools); library(fitdistrplus); library(xtable); library(in2extRemes)
library(evd); library(SimCop); library(copula); library(mgpd)

# E2.1, nine para BEV models ====
par(mfrow = c(2,2)); par(mar = rep(2, 4))

abvevd(x = seq(0,1,by=0.01), dep = 0.99, model = "log", plot =T)
abvevd(x = seq(0,1,by=0.01), dep = 0.8, model = "log", plot =T, add = T, lty = 2)
abvevd(x = seq(0,1,by=0.01), dep = 0.4, model = "log", plot =T, add = T, lty = 3)

abvevd(x = seq(0,1,by=0.01), dep = 0.99, asy = c(0.8,0.2), model = "alog", plot =T)
abvevd(x = seq(0,1,by=0.01), dep = 0.8, asy = c(0.8,0.2), model = "alog", plot =T, add = T, lty = 2)
abvevd(x = seq(0,1,by=0.01), dep = 0.4, asy = c(0.8,0.2), model = "alog", plot =T, add = T, lty = 3)

abvevd(x = seq(0,1,by=0.01), dep = 0.99, model = "hr", plot =T)
abvevd(x = seq(0,1,by=0.01), dep = 0.8, model = "hr", plot =T, add = T, lty = 2)
abvevd(x = seq(0,1,by=0.01), dep = 0.4, model = "hr", plot =T, add = T, lty = 3)

abvevd(x = seq(0,1,by=0.01), dep = 0.99, model = "neglog", plot =T)
abvevd(x = seq(0,1,by=0.01), dep = 0.8, model = "neglog", plot =T, add = T, lty = 2)
abvevd(x = seq(0,1,by=0.01), dep = 0.4, model = "neglog", plot =T, add = T, lty = 3)

abvevd(x = seq(0,1,by=0.01), dep = 0.99, asy = c(0.8,0.2), model = "aneglog", plot =T)
abvevd(x = seq(0,1,by=0.01), dep = 0.8, asy = c(0.8,0.2), model = "aneglog", plot =T, add = T, lty = 2)
abvevd(x = seq(0,1,by=0.01), dep = 0.4, asy = c(0.8,0.2), model = "aneglog", plot =T, add = T, lty = 3)

abvevd(x = seq(0,1,by=0.01), alpha = 0.9, beta = 0.5, model = "bilog", plot =T)
abvevd(x = seq(0,1,by=0.01), alpha = 0.5, beta = 0.5, model = "bilog", plot =T, add = T, lty = 2)
abvevd(x = seq(0,1,by=0.01), alpha = 0.1, beta = 0.5, model = "bilog", plot =T, add = T, lty = 3)

abvevd(x = seq(0,1,by=0.01), alpha = 0.9, beta = 0.9, model = "bilog", plot =T)
abvevd(x = seq(0,1,by=0.01), alpha = 0.5, beta = 0.5, model = "bilog", plot =T, add = T, lty = 2)
abvevd(x = seq(0,1,by=0.01), alpha = 0.5, beta = 0.1, model = "bilog", plot =T, add = T, lty = 3)

abvevd(x = seq(0,1,by=0.01), alpha = 0.9, beta = 0.9, model = "negbilog", plot =T)
abvevd(x = seq(0,1,by=0.01), alpha = 0.5, beta = 0.5, model = "negbilog", plot =T, add = T, lty = 2)
abvevd(x = seq(0,1,by=0.01), alpha = 0.5, beta = 0.1, model = "negbilog", plot =T, add = T, lty = 3)

abvevd(x = seq(0,1,by=0.01), alpha = 0.9, beta = 0.9, model = "ct", plot =T)
abvevd(x = seq(0,1,by=0.01), alpha = 0.5, beta = 0.5, model = "ct", plot =T, add = T, lty = 2)
abvevd(x = seq(0,1,by=0.01), alpha = 0.5, beta = 0.1, model = "ct", plot =T, add = T, lty = 3)

abvevd(x = seq(0,1,by=0.01), alpha = 0.9, beta = 0.05, model = "amix", plot =T)
abvevd(x = seq(0,1,by=0.01), alpha = 0.5, beta = 0.2, model = "amix", plot =T, add = T, lty = 2)
abvevd(x = seq(0,1,by=0.01), alpha = 0.5, beta = 0.1, model = "amix", plot =T, add = T, lty = 3)

# dep
# log, alog - asy, hr, neglog, aneglog - asy
# alpha beta
# bilog, negbilog, ct, amix


# Coles-Tawn and Aneglog seems like interesting models so they are choosen.
abvevd(x = seq(0,1,by=0.01), alpha = 0.9, beta = 0.9, model = "ct", plot =T)
abvevd(x = seq(0,1,by=0.01), alpha = 0.5, beta = 0.5, model = "ct", plot =T, add = T, lty = 2)
abvevd(x = seq(0,1,by=0.01), alpha = 0.5, beta = 0.1, model = "ct", plot =T, add = T, lty = 3)

abvevd(x = seq(0,1,by=0.01), dep = 0.99, asy = c(0.8,0.2), model = "aneglog", plot =T)
abvevd(x = seq(0,1,by=0.01), dep = 0.8, asy = c(0.8,0.2), model = "aneglog", plot =T, add = T, lty = 2)
abvevd(x = seq(0,1,by=0.01), dep = 0.4, asy = c(0.8,0.2), model = "aneglog", plot =T, add = T, lty = 3)

# E2.2 Choose 2 and generate 200 observation ====
# CT
par(mfrow = c(1,2)); par(mar = rep(2, 4))
# Varying beta
abvevd(x = seq(0,1,by=0.01), alpha = 0.8, beta = 20, model = "ct", plot =T, main = "A(w) function of beta, alpha = 0.8", col = "blue")
abvevd(x = seq(0,1,by=0.01), alpha = 0.8, beta = 5, model = "ct", plot =T, add = T, lty = 2, col = "black")
abvevd(x = seq(0,1,by=0.01), alpha = 0.8, beta = 0.7, model = "ct", plot =T, add = T, lty = 3, col = "red")
abvevd(x = seq(0,1,by=0.01), alpha = 0.8, beta = 0.05, model = "ct", plot =T, add = T, lty = 4, col = "green")
legend("bottomright", legend=c("beta = 20", "beta = 5", "beta = 0.7", "beta = 0.5"), col=c("blue", "black", "red", "green"), lty=1:2, cex=0.8)
# Varying alpha
abvevd(x = seq(0,1,by=0.01), alpha = 20, beta = 0.8, model = "ct", plot =T, main = "A(w) function of alpha, beta = 0.8", col = "blue")
abvevd(x = seq(0,1,by=0.01), alpha = 5, beta = 0.8, model = "ct", plot =T, add = T, lty = 2, col = "black")
abvevd(x = seq(0,1,by=0.01), alpha = 0.7, beta = 0.8, model = "ct", plot =T, add = T, lty = 3, col = "red")
abvevd(x = seq(0,1,by=0.01), alpha = 0.05, beta = 0.8, model = "ct", plot =T, add = T, lty = 4, col = "green")
legend("bottomright", legend=c("alpha = 20", "alpha = 5", "alpha = 0.7", "alpha = 0.5"), col=c("blue", "black", "red", "green"), lty=1:2, cex=0.8)

# Aneglog
par(mfrow = c(1,3)); par(mar = rep(2, 4))
# Varying dep
abvevd(x = seq(0,1,by=0.01), dep = 1, asy = c(0.8,0.8), mod = "aneglog", plot =T, main = "A(w) function of dep, asy1, asy2 = 0.8", col = "blue")
abvevd(x = seq(0,1,by=0.01),  dep = 0.66, asy = c(0.8,0.8), mod = "aneglog", plot =T, add = T, lty = 2, col = "black")
abvevd(x = seq(0,1,by=0.01),  dep = 0.33, asy = c(0.8,0.8), mod = "aneglog", plot =T, add = T, lty = 3, col = "red")
abvevd(x = seq(0,1,by=0.01),  dep = 0.01, asy = c(0.8,0.8), mod = "aneglog", plot =T, add = T, lty = 4, col = "green")
legend("bottomright", legend=c("dep = 1", "dep = 0.66", "dep = 0.33", "dep = 0.01"), col=c("blue", "black", "red", "green"), lty=1:2, cex=1)
# Varying asy2
abvevd(x = seq(0,1,by=0.01), dep = 0.5, asy = c(0.8,1), mod = "aneglog", plot =T, main = "A(w) function of asy2, dep = 0.5 asy1 = 0.8", col = "blue")
abvevd(x = seq(0,1,by=0.01), dep = 0.5, asy = c(0.8,0.66), mod = "aneglog", plot =T, add = T, lty = 2, col = "black")
abvevd(x = seq(0,1,by=0.01), dep = 0.5, asy = c(0.8,0.33), mod = "aneglog", plot =T, add = T, lty = 3, col = "red")
abvevd(x = seq(0,1,by=0.01), dep = 0.5, asy = c(0.8,0.01), mod = "aneglog", plot =T, add = T, lty = 4, col = "green")
legend("bottomright", legend=c("asy2 = 1", "asy2 = 0.66", "asy2 = 0.33", "asy2 = 0.01"), col=c("blue", "black", "red", "green"), lty=1:2, cex=1)
# Varying asy1
abvevd(x = seq(0,1,by=0.01),  dep = 0.5, asy = c(1,0.8), mod = "aneglog", plot =T, main = "A(w) function of asy1, dep = 0.5 asy2 = 0.8", col = "blue")
abvevd(x = seq(0,1,by=0.01),  dep = 0.5, asy = c(0.66,0.8), mod = "aneglog", add = T, lty = 2, col = "black")
abvevd(x = seq(0,1,by=0.01), dep = 0.5, asy = c(0.33,0.8), mod = "aneglog", plot =T, add = T, lty = 3, col = "red")
abvevd(x = seq(0,1,by=0.01), dep = 0.5, asy = c(0.01,0.8), mod = "aneglog", plot =T, add = T, lty = 4, col = "green")
legend("bottomright", legend=c("asy1 = 1", "asy1 = 0.66", "asy1 = 0.33", "asy1 = 0.01"), col=c("blue", "black", "red", "green"), lty=1:2, cex=1)

n = 200
r.ct_1 <- rbvevd(n, alpha = 0.8, beta = 20, mod = "ct")
r.ct_2 <- rbvevd(n, alpha = 0.8, beta = 5, mod = "ct")
r.ct_3 <- rbvevd(n, alpha = 0.8, beta = 0.7, mod = "ct")
r.ct_4 <- rbvevd(n, alpha = 0.8, beta = 0.05, mod = "ct")
par(mfrow = c(2,2)); par(mar = rep(2, 4))
plot(r.ct_1, main =  "Coles-Tawn model, a = 0.8, b = 20", col = "black", cex = 0.5)
plot(r.ct_2, main =  "CT, a = 0.8, b = 5", col = "black", cex = 0.5)
plot(r.ct_3, main =  "CT, a = 0.8, b = 0.7", col = "black", cex = 0.5)
plot(r.ct_4, main =  "CT, a = 0.8, b = 0.05", col = "black", cex = 0.5)

# Aneglog
abvevd(x = seq(0,1,by=0.01), dep = 0.99, asy = c(0.8,0.2), model = "aneglog", plot =T)

r.anlo_1 <- rbvevd(n, dep = 0.5, asy = c(0.8,0.2), mod = "aneglog")
r.anlo_2 <- rbvevd(n, dep = 5, asy = c(0.8,0.2), mod = "aneglog")
r.anlo_3 <- rbvevd(n, dep = 50, asy = c(0.8,0.2), mod = "aneglog")
r.anlo_4 <- rbvevd(n, dep = 500, asy = c(0.8,0.2), mod = "aneglog")
par(mfrow = c(2,2)); par(mar = rep(2, 4))
plot(r.anlo_1, main =  " dep = 0.5, asy = c(0.8,0.2)", col = "black", cex = 0.5)
plot(r.anlo_2, main =  " dep = 5, asy = c(0.8,0.2)", col = "black", cex = 0.5)
plot(r.anlo_3, main =  " dep = 50, asy = c(0.8,0.2)", col = "black", cex = 0.5)
plot(r.anlo_4, main =  " dep = 500, asy = c(0.8,0.2)", col = "black", cex = 0.5)

cor(r.ct_1, method = "pearson")
cor(r.ct_2, method = "pearson")
cor(r.ct_3, method = "pearson")
cor(r.ct_4, method = "pearson")

cor(r.ct_1, method = "kendall")
cor(r.ct_2, method = "kendall")
cor(r.ct_3, method = "kendall")
cor(r.ct_4, method = "kendall")

cor(r.ct_1, method = "spearman")
cor(r.ct_2, method = "spearman")
cor(r.ct_3, method = "spearman")
cor(r.ct_4, method = "spearman")

cor(r.anlo_1, method = "pearson")
cor(r.anlo_2, method = "pearson")
cor(r.anlo_3, method = "pearson")
cor(r.anlo_4, method = "pearson")

cor(r.anlo_1, method = "kendall")
cor(r.anlo_2, method = "kendall")
cor(r.anlo_3, method = "kendall")
cor(r.anlo_4, method = "kendall")

cor(r.anlo_1, method = "spearman")
cor(r.anlo_2, method = "spearman")
cor(r.anlo_3, method = "spearman")
cor(r.anlo_4, method = "spearman")

# E2.4 ==== Generate from bivariate asymmetric logistic copula
# r > 1, th > 0, 1 > phi > 0
n = 200
cop_1 <- NewBEVAsyLogisticCopula(r = 2,theta=1, phi= 0.5)
cop_2 <- NewBEVAsyLogisticCopula(r = 20,theta=1, phi= 0.5)
cop_3 <- NewBEVAsyLogisticCopula(r = 2000,theta=1, phi= 0.5)
cop_4 <- NewBEVAsyLogisticCopula(r = 2,theta=0.5, phi= 0.5)
cop_5 <- NewBEVAsyLogisticCopula(r = 2,theta=5, phi= 0.5)
cop_6 <- NewBEVAsyLogisticCopula(r = 2,theta=500, phi= 0.5)
cop_7 <- NewBEVAsyLogisticCopula(r = 2,theta=1, phi= 0)
cop_8 <- NewBEVAsyLogisticCopula(r = 2,theta=1, phi= 0.5)
cop_9 <- NewBEVAsyLogisticCopula(r = 2,theta=1, phi= 1)

approx_1 <- GetApprox(cop_1)
approx_2 <- GetApprox(cop_2)
approx_3 <- GetApprox(cop_3)
approx_4 <- GetApprox(cop_4)
approx_5 <- GetApprox(cop_5)
approx_6 <- GetApprox(cop_6)
approx_7 <- GetApprox(cop_7)
approx_8 <- GetApprox(cop_8)
approx_9 <- GetApprox(cop_9)

sample_1 <- GenerateRV(approx_1, n)
sample_2 <- GenerateRV(approx_2, n)
sample_3 <- GenerateRV(approx_3, n)
sample_4 <- GenerateRV(approx_4, n)
sample_5 <- GenerateRV(approx_5, n)
sample_6 <- GenerateRV(approx_6, n)
sample_7 <- GenerateRV(approx_7, n)
sample_8 <- GenerateRV(approx_8, n)
sample_9 <- GenerateRV(approx_9, n)

par(mfrow = c(3,3)); par(mar = rep(2, 4))
plot(sample_1)
plot(sample_2)
plot(sample_3)
plot(sample_4)
plot(sample_5)
plot(sample_6)
plot(sample_7)
plot(sample_8)
plot(sample_9)

# E3 Preprocessing ====
source("fremantle.R.txt")
source("portpirie.R.txt")
max(fremantle[,2])
max(portpirie[,2])
slMerged <- merge(fremantle, portpirie, by = "Year")
sl <- data.frame(year = slMerged[,1], fremantle=slMerged[,2], Portpirie = slMerged[,4])
max(sl[,2]); max(sl[,3])
# Scatter plots
par(mfrow = c(1,2))#; par(mar = rep(2, 4))
plot(sl[,1], sl[,2], xlab = "year", ylab = "Sea level [M]") # Fremantle
abline(h = 1.7)
plot(sl[,1], sl[,3], xlab = "year", ylab = "Sea level [M]") # Portpirie
abline(h = 4.2)

# FML
M1 <- fbvevd(sl[,2:3], model = "ct")
M2 <- fbvevd(sl[,2:3], model = "aneglog")
M3 <- fbvevd(sl[,2:3], model = "bilog")

round(cbind(Estimates = fitted(M1), StandardErrors = std.errors(M1)),3)
round(cbind(Estimates = fitted(M2), StandardErrors = std.errors(M2)),3)
round(cbind(Estimates = fitted(M3), StandardErrors = std.errors(M3)),3)

# IFM - Estimate parameters Marginal distributions.
f1 <- fgev(sl[,2])
f2 <- fgev(sl[,3])
f1$est
f2$est

# IFM - Estimate parameters in dependence function.
M1.IFM <- fbvevd(sl[,2:3], model = "ct", loc1 = f1$est[1], scale1 = f1$est[2], 
                 shape1 = f1$est[3], loc2 = f2$est[1], scale2 = f2$est[2], shape2 = f2$est[3])
M2.IFM <- fbvevd(sl[,2:3], model = "aneglog", loc1 = f1$est[1], scale1 = f1$est[2],
                 shape1 = f1$est[3], loc2 = f2$est[1], scale2 = f2$est[2], shape2 = f2$est[3])
M3.IFM <- fbvevd(sl[,2:3], model = "bilog", loc1 = f1$est[1], scale1 = f1$est[2],
                 shape1 = f1$est[3], loc2 = f2$est[1], scale2 = f2$est[2], shape2 = f2$est[3])

# Evaluation of parametric models.
all.AIC <- AIC(M1,M2,M3, M1.IFM, M2.IFM,M3.IFM)
all.AIC[min(all.AIC[,2]) == all.AIC[,2],]

# Non-parametric fit
N1 <- abvnonpar(data = sl[,2:3], epmar = F, method = "cfg", convex = F, plot =T, col = "blue",
                main = "Dependence function estimation using cfg")
N2 <- abvnonpar(data = sl[,2:3], epmar = F, method = "cfg", convex = T, add = T, col = "black")
N3 <- abvnonpar(data = sl[,2:3], epmar = T, method = "cfg", convex = F, add = T, col = "red")
N4 <- abvnonpar(data = sl[,2:3], epmar = T, method = "cfg", convex = T, add = T, col = "green")
legend("bottomright", legend=c("Parametric transformation","Convex hull of parametric transformation", "Empiric transformation",  "Convex hull of Empiric transformation"), col=c("blue", "black", "red", "green"), lty=1:4, cex=1.0)

# E3.6 and 3.7 Est prob ====
# Parametric estimation of probabilities
sl_am <- as.matrix(sl[,2:3])
# M3 <- fbvevd(sl_am, model = "bilog")
# M3$estimate

mar1 = c(f1$e[1],f1$e[2],f1$e[3])
mar2 = c(f2$e[1],f2$e[2],f2$e[3])
a <- M3.IFM$e[1]
b <- M3.IFM$e[2]
# a <- M1.IFM$e[1]
# b <- M1.IFM$e[2]

# P1
x <- c(30,4.2)
p_CD <- pbvevd(q = x, alpha = a, beta = b, model = "bilog",
               mar1 = mar1, mar2 = mar2)
x <- c(1.7,50)
p_BD <- pbvevd(q = x, alpha = a, beta = b, model = "bilog",
               mar1 = mar1, mar2 = mar2)
x <- c(1.7,4.2)
p_D <- pbvevd(q = x, alpha = a, beta = b, model = "bilog",
              mar1 = mar1, mar2 = mar2)

p_A = 1 - p_CD - p_BD + p_D
round(p_A,4)

# P2
x <- c(30,4.4)
p_CD <- pbvevd(q = x, alpha = a, beta = b, model = "bilog",
               mar1 = mar1, mar2 = mar2)
x <- c(1.8,50)
p_BD <- pbvevd(q = x, alpha = a, beta = b, model = "bilog",
               mar1 = mar1, mar2 = mar2)
x <- c(1.8,4.4)
p_D <- pbvevd(q = x, alpha = a, beta = b, model = "bilog",
              mar1 = mar1, mar2 = mar2)

p_A = 1 - p_CD - p_BD + p_D
round(p_A,5)

# P3
x <- c(1.478,3.85)
p_14_38 <- pbvevd(q = x, alpha = a, beta = b, model = "bilog",
                  mar1 = mar1, mar2 = mar2)
p_14_38 = 1 - p_14_38
round(p_14_38,4)

# Parametric estimation of p4
mar1 = c(f1$e[1],f1$e[2],f1$e[3])
mar2 = c(f2$e[1],f2$e[2],f2$e[3])
a <- M3.IFM$e[1]
b <- M3.IFM$e[2]
# a <- M1.IFM$e[1]
# b <- M1.IFM$e[2]

x <- c(1.478,3.850)
p_D_1 <- pbvevd(q = x, alpha = a, beta = b, model = "bilog",
                mar1 = mar1, mar2 = mar2)
p_D_1
x <- c(1.95,4.8)
p_D_2 <- pbvevd(q = x, alpha = b, beta = b, model = "bilog",
                mar1 = mar1, mar2 = mar2)
p_D_2

p = ((1-p_D_1) - (1-p_D_2))/(1-p_D_1)
p

# 3.6 and 3.7 Non-Parametric estimation.
# Parametric transformation
# Estimatation of p1
y1 = c(3, 1.7, 1.7)
y2 = c(4.2, 5, 4.2)
G1 <- pgev(y1, mar1[1], mar1[2], mar1[3])
G2 <- pgev(y2, mar2[1], mar2[2], mar2[3])
p <- exp(log(G1*G2)*N2$y[ceiling(log(G2)/log(G1*G2)*100)])
p1 = 1 + p[3] - p[1] - p[2]
round(p1,4)

# Estimation of p2
y1 = c(3, 1.8, 1.8)
y2 = c(4.4, 5, 4.4)
G1 <- pgev(y1, mar1[1], mar1[2], mar1[3])
G2 <- pgev(y2, mar2[1], mar2[2], mar2[3])
p <- exp(log(G1*G2)*N2$y[ceiling(log(G2)/log(G1*G2)*100)])
p2 = 1 + p[3] - p[1] - p[2]
round(p2,4)

# Estimation of p3
y1 = 1.478
y2 = 3.85
G1 <- pgev(y1, mar1[1], mar1[2], mar1[3])
G2 <- pgev(y2, mar2[1], mar2[2], mar2[3])
p <- exp(log(G1*G2)*N2$y[ceiling(log(G2)/log(G1*G2)*100)])
p3 <- 1-p
round(p3,3)

# Estimation of p4
y1 = 1.478
y2 = 3.85
G1 <- pgev(y1, mar1[1], mar1[2], mar1[3])
G2 <- pgev(y2, mar2[1], mar2[2], mar2[3])
p_D_1 <- exp(log(G1*G2)*N2$y[ceiling(log(G2)/log(G1*G2)*100)])

y1 = 1.95
y2 = 4.8
G1 <- pgev(y1, mar1[1], mar1[2], mar1[3])
G2 <- pgev(y2, mar2[1], mar2[2], mar2[3])
p_D_2 <- exp(log(G1*G2)*N2$y[ceiling(log(G2)/log(G1*G2)*100)])

p4 = ((1-p_D_1) - (1-p_D_2))/(1-p_D_1)
p4

# Empirical transformation
# Estimatation of p1
y1 = c(3, 1.7, 1.7)
y2 = c(4.2, 5, 4.2)
G1 <- pgev(y1, mar1[1], mar1[2], mar1[3])
G2 <- pgev(y2, mar2[1], mar2[2], mar2[3])
p <- exp(log(G1*G2)*N4$y[ceiling(log(G2)/log(G1*G2)*100)])
p1 = 1 + p[3] - p[1] - p[2]
round(p1,4)

# Estimation of p2
y1 = c(3, 1.8, 1.8)
y2 = c(4.4, 5, 4.4)
G1 <- pgev(y1, mar1[1], mar1[2], mar1[3])
G2 <- pgev(y2, mar2[1], mar2[2], mar2[3])
p <- exp(log(G1*G2)*N4$y[ceiling(log(G2)/log(G1*G2)*100)])
p2 = 1 + p[3] - p[1] - p[2]
round(p2,4)

# Estimation of p3
y1 = 1.478
y2 = 3.85
G1 <- pgev(y1, mar1[1], mar1[2], mar1[3])
G2 <- pgev(y2, mar2[1], mar2[2], mar2[3])
p <- exp(log(G1*G2)*N4$y[ceiling(log(G2)/log(G1*G2)*100)])
p3 <- 1-p
round(p3,3)

# Estimation of p4
y1 = 1.478
y2 = 3.85
G1 <- pgev(y1, mar1[1], mar1[2], mar1[3])
G2 <- pgev(y2, mar2[1], mar2[2], mar2[3])
p_D_1 <- exp(log(G1*G2)*N4$y[ceiling(log(G2)/log(G1*G2)*100)])

y1 = 1.95
y2 = 4.8
G1 <- pgev(y1, mar1[1], mar1[2], mar1[3])
G2 <- pgev(y2, mar2[1], mar2[2], mar2[3])
p_D_2 <- exp(log(G1*G2)*N4$y[ceiling(log(G2)/log(G1*G2)*100)])

p4 = ((1-p_D_1) - (1-p_D_2))/(1-p_D_1)
p4

# E3.8 plotting ====
# Creat rbvevd object
par(mfrow = c(1,2)); par(mar = rep(2, 4))
n = 63
bvdata <- rbvevd(n, alpha = a, beta = b, model = "bilog",
                 mar1 = mar1, mar2 = mar2)
M1 <- fbvevd(bvdata, model = "bilog")
plot(M1, which = 5, p = c(0.75, 0.9, 0.95), col = "blue", xlim = c(1.2,2), ylim = c(3.5,5))
plot(M2, which = 5, p = c(0.75, 0.9, 0.95), col ="green", xlim = c(1.2,2), ylim = c(3.5,5))

# E3.9 Linear trend ====

slMerged <- merge(fremantle, portpirie, by = "Year")
sl_all <- data.frame(year = slMerged[,1], fremantle=slMerged[,2], SOI = slMerged[,3], Portpirie = slMerged[,4])

# Normalize year
sl_all[,1] = 2*(sl_all[,1]-sl_all[1,1])/(sl_all[63,1]-sl_all[1,1]) -1

M_trend.IFM <- fbvevd(sl[,2:3], model = "bilog", nsloc1 = sl_all[,1], nsloc2 = sl_all[,3], loc1 = f1$est[1], scale1 = f1$est[2],
                      shape1 = f1$est[3], loc2 = f2$est[1], scale2 = f2$est[2], shape2 = f2$est[3])

M_trend_FML <- fbvevd(sl[,2:3], nsloc1 = sl_all[,1], nsloc2 = sl_all[,3], model = "bilog")

anova(M_trend_FML, M_trend.IFM)

# To be compared to M3.IFM and M3
x = logLik(M_trend_FML)
y = logLik(M3)
lr.test(x = x, y = y, alpha = 0.05, df = 2)

y = logLik(M_trend.IFM)
x = logLik(M3.IFM)
lr.test(x = x, y = y, alpha = 0.05, df = 2)