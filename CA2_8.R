library(sn);library(mvtnorm); library(scatterplot3d); library(copula); library(EnvStats); 
library(plotly); library(mdatools); library(evd); library(SimCop); library(fitdistrplus)
# E8 ====
n = 5000000

myMvd_normal <- mvdc(copula = ellipCopula(family = "normal", param = 0.7),
                     margins = c("lnorm", "lnorm"), paramMargins = list(list(meanlog  = 0,
                    sdlog  = 1), list(meanlog  = 0, sdlog  = 1)))
myMvd_t <- mvdc(copula = ellipCopula(family = "t", param = 0.7, df = 2),
                margins = c("lnorm", "lnorm"), paramMargins = list(list(meanlog  = 0,
                sdlog  = 1), list(meanlog  = 0, sdlog  = 1)))
myMvd_gumbel <- mvdc(copula = archmCopula(family = "gumbel", param = 2),
                     margins = c("lnorm", "lnorm"), paramMargins = list(
                    list(meanlog  = 0, sdlog  = 1), list(meanlog  = 0, sdlog  = 1)))

# Threshold calculation
alpha = 0.01
u_alpha = exp(sqrt(2)*erfinv(1 - 2*alpha))   
u_95 = exp(sqrt(2)*erfinv(1 - 2*0.05))   
u_99 = exp(sqrt(2)*erfinv(1 - 2*0.01))   
print(u_95)
print(u_99)

# Generating observations
set.seed(8)
data.norm <- rMvdc(n, myMvd_normal)
set.seed(8)
data.t <- rMvdc(n, myMvd_t)
set.seed(8)
data.gumbel <- rMvdc(n, myMvd_gumbel)

# Checking Kendall´s Tau
#data.norm.tau <- cor(data.norm, method = "kendall")
#data.t.tau <- cor(data.t, method = "kendall")
#data.gumbel.tau <- cor(data.gumbel, method = "kendall")
#print(data.norm.tau)                                    
#print(data.t.tau)                                       
#print(data.gumbel.tau)         
# Gives more or less tau = 0.5 for all copulas, regardless of seed.
data.norm.spea <- cor(data.norm, method = "spearman")
data.t.spea <- cor(data.t, method = "spearman")

# calculating number of exceedences.
data.norm.bool <- (data.norm > u_99)
temp1 = matrix(0, nrow = n, ncol = 3)
temp1[,1:2] = data.norm.bool
temp1[,3] = temp1[,1]*temp1[,2]
data.norm.bool <- temp1

data.t.bool <- (data.t > u_99)
temp2 = matrix(0, nrow = n, ncol = 3)
temp2[,1:2] = data.t.bool
temp2[,3] = temp2[,1]*temp2[,2]
data.t.bool <- temp2

data.gumbel.bool <- (data.gumbel > u_99)
temp3 = matrix(0, nrow = n, ncol = 3)
temp3[,1:2] = data.gumbel.bool
temp3[,3] = temp3[,1]*temp3[,2]
data.gumbel.bool <- temp3

colSums(data.norm.bool)/n
colSums(data.t.bool)/n
colSums(data.gumbel.bool)/n

# 4 Calculating expected loss

ind_1 <- which(temp1[,3] > 0)
ind_2 <- which(temp2[,3] > 0)
ind_3 <- which(temp3[,3] > 0)

loss_s_norm <- sum(data.norm[which(temp1[,3] > 0), 1:2])
loss_s_t <- sum(data.t[which(temp2[,3] > 0), 1:2])
loss_s_gumbel <- sum(data.gumbel[which(temp3[,3] > 0), 1:2])

print(loss_s_norm/n)
print(loss_s_t/n)
print(loss_s_gumbel/n)
