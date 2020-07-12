library(sn);library(mvtnorm); library(scatterplot3d); library(copula); library(EnvStats); library(plotly)
library(mdatools); library(evd); library(SimCop); library(fitdistrplus)
n = 500

# E5 part 1 ====
myMvd_normal <- mvdc(copula = ellipCopula(family = "normal", param = 0.75), margins = c("norm", "exp"), 
                     paramMargins = list(list(mean = 0, sd = 2), list(rate = 2)))
myMvd_gumbel <- mvdc(copula = archmCopula(family = "gumbel", param = 2), margins = c("exp", "exp"), 
                     paramMargins = list(list(rate = 4), list(rate = 2)))
dat.norm <- rMvdc(n, myMvd_normal)
dat.gumb <- rMvdc(n, myMvd_gumbel)

norm_rho_p = cor(dat.norm, method = "pearson")
gumb_rho_p = cor(dat.gumb, method = "pearson")

norm_rho_s = cor(dat.norm, method = "spearman")
gumb_rho_s = cor(dat.gumb, method = "spearman")

n = 200
# E5 part 2 ====
# shape, scale for gamma > 0
myMvd_normal_1 <- mvdc(copula = ellipCopula(family = "normal", param = 0.75), margins = c("norm", "exp"), 
                       paramMargins = list(list(mean = 0, sd = 2), list(rate = 2)))
myMvd_normal_2 <- mvdc(copula = ellipCopula(family = "normal", param = 0.75), margins = c("gamma", "gamma"), 
                       paramMargins = list(list(shape = 2, scale = 2), list(shape = 1, scale = 3)))
myMvd_normal_3 <- mvdc(copula = ellipCopula(family = "normal", param = 0.75), margins = c("gamma", "exp"), 
                       paramMargins = list(list(shape = 3, scale = 1), list(rate = 2)))

set.seed(1)
dat.norm_1 <- rMvdc(n, myMvd_normal_1)
set.seed(1)
dat.norm_2 <- rMvdc(n, myMvd_normal_2)
set.seed(1)
dat.norm_3 <- rMvdc(n, myMvd_normal_3)

r.p.1 = cor(dat.norm_1, method = "pearson")
r.p.2 = cor(dat.norm_2, method = "pearson")
r.p.3 = cor(dat.norm_3, method = "pearson")

r.s.1 = cor(dat.norm_1, method = "spearman")
r.s.2 = cor(dat.norm_2, method = "spearman")
r.s.3 = cor(dat.norm_3, method = "spearman")

# E6: Graphics for md ====
par(mfrow = c(2,2)); par(mar = rep(2, 4))
plot(dat.norm, main = "normal, marg = norm, exp")
plot(dat.gumb, main = "gumbel, marg = exp, exp")
plot(dat.norm_1, main = "normal, marg = norm, exp")
plot(dat.norm_2, main = "normal, marg = gamma, gamma")

par(mfrow = c(2,2)); par(mar = rep(2, 4))
contour(myMvd_normal, dMvdc, xlim = c(-5, 3), ylim = c(-1, 2), main = "normal, marg = norm, exp")
contour(myMvd_gumbel, dMvdc, xlim = c(-0.1, 0.4), ylim = c(-0.2, 0.6), main = "gumbel, marg = exp, exp")
contour(myMvd_normal_1, dMvdc, xlim = c(-5, 3.2), ylim = c(-0.3, 1.4), main = "normal, marg = norm, exp")
contour(myMvd_normal_2, dMvdc, xlim = c(0, 5.5), ylim = c(-1, 4), main = "normal, marg = gamma, gamma")