library(sn);library(mvtnorm); library(scatterplot3d); library(copula); library(EnvStats);
library(plotly); library(mdatools); library(evd); library(SimCop); library(fitdistrplus);
library(xtable)
# E9 ====
in_dat <- read.table(
  "http://www.maths.lth.se/matstat/kurser/fmsn15masm23/datasets/insuranceData.txt",
                     header = T)
#par(mfrow = c(2,2))
#par(mar = rep(2, 4))
plot(in_dat, log = "xy", main = "Insurance costs data")

cor(in_dat, method = "pearson")
tau_k <- cor(in_dat, method = "kendall")
cor(in_dat, method = "spearman")
data <- as.matrix(in_dat[c("Loss", "ALAE")])
n = dim(data)[1]
xtable(cbind(cor(in_dat, method = "pearson"), cor(in_dat, method = "kendall"),
             cor(in_dat, method = "spearman")))

# Marginal fits ====
# Data is transformed to avoid numerical instability.
# This is only used for the Generalised Extreme Value distribution.

par(mfrow = c(2,2)); par(mar = rep(2, 4))
hist(data[,1],100, main = "Histogram of Loss")
hist(log(data[,1]),100, main = "Hist of log Loss")
hist(data[,2],100, main = "Hist of ALAE")
hist(log(data[,2]),100, main = "Hist of log ALAE")

data_100 = data/100000
# 1
fcauchy_1 <- fitdist(data[,1], "cauchy")
flnorm_1 <- fitdist(data[,1], "lnorm")
ft_1 <- fitdist(data[,1], "t", start = list(df = 2))
fexp_1 <- fitdist(data[,1], "exp", start = list(rate = 1/100), method ="mme")
fnorm_1 <- fitdist(data[,1], "norm")

fgev_1 <- fgev(data_100[,1], std.err = FALSE)

g <- gofstat(list(fcauchy_1, flnorm_1, ft_1, fexp_1), fitnames = c("cauchy", "lnorm", "t",
                                                                   "exp"))
g$kstest

# 2
fcauchy_2 <- fitdist(data[,2], "cauchy")
flnorm_2 <- fitdist(data[,2], "lnorm")
ft_2 <- fitdist(data[,2], "t", start = list(df = 2))
fexp_2 <- fitdist(data[,2], "exp", start = list(rate = 1/100), method ="mme")

par(mfrow = c(2,2)); par(mar = rep(2, 4))
qqcomp(fcauchy_2, addlegend = F, main = "Cauchy", fitpch = 16, fitcol = "grey", line01lty = 2)
qqcomp(ft_2, addlegend = F, main = "t", fitpch = 16, fitcol = "grey", line01lty = 2)
qqcomp(fexp_2, addlegend = F, main = "Exponential", fitpch = 16, fitcol = "grey", line01lty = 2)
qqcomp(fnorm_1, addlegend = F, main = "Normal", fitpch = 16, fitcol = "grey", line01lty = 2)

fgev_2 <- fgev(data_100[,2], std.err = FALSE)

g <- gofstat(list(fcauchy_2, flnorm_2, ft_2, fexp_2), fitnames = c("cauchy", "lnorm", "t",
                                                                   "exp"))
g$kstest

# denscomp(flnorm_2, addlegend = F, main = "", xlab = "data", fitcol = "orange")
# qqcomp(flnorm_2, addlegend = F, main = "", fitpch = 16, fitcol = "grey", line01lty = 2)
# cdfcomp(flnorm_2, addlegend = F, main = "", xlab = "data", fitcol = "orange", lines01 = TRUE)
# ppcomp(flnorm_2, addlegend = F, main = "", fitpch = 16, fitcol = "grey", line01lty = 2)

# Empirical quantiles
data_gev_1 = rgev(1500, 0.06455, 0.08859, 1.17632)
data_gev_2 = rgev(1500, 0.03445, 0.03969, 0.79596)
gev_1_sort = sort(100000*data_gev_1)
gev_2_sort = sort(100000*data_gev_2)


par(mfrow = c(2,2), mar = c(4, 4, 1, 1))
qqcomp(flnorm_1, addlegend = FALSE, main = "", fitpch = 16, fitcol = "grey", line01lty = 2)
qqcomp(flnorm_2, addlegend = FALSE, main = "", fitpch = 16, fitcol = "grey", line01lty = 2)
plot(data[,1], log = "xy", ylab = "Loss")
lines(gev_1_sort, col = "red")
plot(sort(data[,2]), log = "xy", ylab = "ALAE")
lines(gev_2_sort, col = "red")


# Conclusion both margins are Lognormal or GEV.

# Copula ====
myMvd_normal <- mvdc(copula = ellipCopula(family = "normal", param = 0.5), 
              margins = c("lnorm", "lnorm"), paramMargins = list(
                list(meanlog = flnorm_1$estimate[1],
              sdlog = flnorm_1$estimate[2]), list(meanlog = flnorm_2$estimate[1], 
              sdlog = flnorm_2$estimate[2])))
myMvd_t <- mvdc(copula = ellipCopula(family = "t", param = 0.5), 
              margins = c("lnorm", "lnorm"), paramMargins = list(
                list(meanlog = flnorm_1$estimate[1],
              sdlog = flnorm_1$estimate[2]), list(meanlog = flnorm_2$estimate[1], 
              sdlog = flnorm_2$estimate[2])))
myMvd_gumbel <- mvdc(copula = archmCopula(family = "gumbel", param = 2), 
              margins = c("lnorm", "lnorm"), paramMargins = list(
                list(meanlog = flnorm_1$estimate[1],
              sdlog = flnorm_1$estimate[2]), list(meanlog = flnorm_2$estimate[1], 
              sdlog = flnorm_2$estimate[2])))
                                                                          

# Start value by margin and correlation test.
c1 = flnorm_1$estimate[1]
c2 <- flnorm_1$estimate[2]
c3 <- flnorm_2$estimate[1]
c4 <- flnorm_2$estimate[2]
c5 <- tau_k
start_normal <- c(9.373159 , 1.637848 , 8.521815, 1.4293 ,0.3112)
start_gumbel <- c(9.373159 , 1.637848 , 8.521815, 1.4293 ,2)
# para_gev_1 = c(0.06455, 0.08859, 1.17632)
# para_gev_2 = c(0.03445, 0.03969, 0.79596)

# ML maxium liklihood ====
fit_normal <- fitMvdc(data, myMvd_normal, start = start_normal, 
                      optim.control = list(trace = TRUE, maxit = 2000))
fit_normal

fit_gumbel <- fitMvdc(data, myMvd_gumbel, start = start_gumbel,
                      optim.control = list(trace = TRUE, maxit = 2000))
fit_gumbel


gofCopula(copula = archmCopula(family = "gumbel", param = 1.468), data, N = 5)
# Resulting statistic = 0.025239, parameter = 1.4394, p-value = 0.1863
gofCopula(copula = ellipCopula(family = "norm", param = 0.4398), data, N = 5)
# Resulting statistic = 0.10393, parameter = 0.46563, p-value = 0.009804

myMvd_gumbel_2 <- mvdc(copula = archmCopula(family = "gumbel", param = 1.468), 
                     margins = c("lnorm", "lnorm"), paramMargins = list(list(meanlog = 9.373,
                     sdlog = 1.670), list(meanlog = 8.523, 
                     sdlog = 1.428)))

n = 1500
# Generating data for explected value calculations ====

myMvd_gumbel_gev <- mvdc(copula = archmCopula(family = "gumbel", param = 1.4), 
                         margins = c("gev", "gev"), paramMargins = list(list(loc = 0.06455,
                         scale = 0.08859, shape = 1.17632),
                         list(loc = 0.03445, scale = 0.03969, shape = 0.79596)))
data.gumbel_gev <- rMvdc(200000, myMvd_gumbel_gev)
data.gumbel_gev = data.gumbel_gev*100000 # Transform back
data_gev = data.gumbel_gev

data.gumbel_2 <- rMvdc(n, myMvd_gumbel_2)
data_gen <- rMvdc(200000, myMvd_gumbel_2)

data_ind <- cbind(rlnorm(200000, meanlog = 9.37316, sdlog = 1.63785), rlnorm(10000,
                  meanlog = 8.52182, sdlog = 1.4293))

# plots with xlim
par(mfrow = c(1,2))
m1 = min(data[,1]); m2 = min(data[,2]); ma1 = max(data[,1]); ma2 = max(data[,2])
plot(in_dat, log = "xy", xlim = c(m1,ma1), ylim = c(m2, ma2), col = "blue", pch = 1,
     main = "Lognorm margins")
points(data.gumbel_2, xlim = c(m1,ma1), ylim = c(m2, ma2), col = "red", pch = 3)
plot(in_dat, log = "xy", xlim = c(m1,ma1), ylim = c(m2, ma2), col = "blue", pch = 1,
     main = "Gev margins")
#points(data_gev, xlim = c(m1,ma1), ylim = c(m2, ma2), col = "red", pch = 3)

L = c(10^4, 5*10^5, 10^6)
RL = c(0, 0.25, 0.75, 0.95)

E_data <- matrix(0, nrow = 3, ncol = 4)
for (i in 1:3){
  for (j in 1:4){
    for (k in 1:1500){
      E_data[i,j] = E_data[i,j] + RP_function(data[k,1],data[k,2],L[i],RL[j])
    }
    E_data[i,j] = E_data[i,j]/1500
  }
}
E_data = round(E_data,-1)

E_gen <- matrix(0, nrow = 3, ncol = 4)
n_gen = dim(data_gen)[1]
for (i in 1:3){
  for (j in 1:4){
    for (k in 1:n_gen){
      E_gen[i,j] = E_gen[i,j] + RP_function(data_gen[k,1],data_gen[k,2],L[i],RL[j])
    }
    E_gen[i,j] = E_gen[i,j]/n_gen
  }
}
E_gen = round(E_gen,-1)

E_gev <- matrix(0, nrow = 3, ncol = 4)
n_gev = dim(data_gev)[1]
for (i in 1:3){
  for (j in 1:4){
    for (k in 1:n_gev){
      E_gev[i,j] = E_gev[i,j] + RP_function(data_gev[k,1],data_gev[k,2],L[i],RL[j])
    }
    E_gev[i,j] = E_gev[i,j]/n_gev
  }
}
E_gev = round(E_gev,-1)

n_ind = dim(data_ind)[1]
E_ind <- matrix(0, nrow = 3, ncol = 4)
for (i in 1:3){
  for (j in 1:4){
    for (k in 1:n_ind){
      E_ind[i,j] = E_ind[i,j] + RP_function(data_ind[k,1],data_ind[k,2],L[i],RL[j])
    }
    E_ind[i,j] = E_ind[i,j]/n_ind
  }
}
E_ind = round(E_ind,-1)

# Reinsurer´s payment function ====
RP_function <- function(X1, X2, L, RL){
  R = RL*L
  if (X1 < R){RP <- 0}
  if ((R < X1) && (X1 < L)){RP <- X1 - R + X2*(X1 - R)/X1 + 4}
  if (X1 > L){RP <- L - R + X2*(L-R)/L}
  return(RP)
}

