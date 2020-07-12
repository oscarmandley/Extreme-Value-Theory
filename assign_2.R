library(in2extRemes)
in2extRemes()
par(mar=c(1,1,1,1))
# File -> Read data -> choose rain.R -> R source -> name "rain"

names(rain) 
# Make fit1 by Analyze -> ext value dist, GP, 30 threshold
names(rain$models$fit1$results)

# 5

fevd(x = value, data = rain$data, threshold = 30, type = 'GP')
v <- c(0.9187666, -0.0655054, -0.0655054,0.01024179)


parcov.fevd(rain$models$fit1)
solve(rain$models$fit1$results$hessian)

ci(rain$models$fit1,type="parameter")
ci(rain$models$fit1, type = "parameter", method="proflik",which.par=2,xrange=c(-0.5,1),nint=100)

# 12. Use %*% for matrix multiplication. Dont transpose
V <- matrix(0:0, nrow = 3, ncol = 3)
V[2,2:3] <- v[1:2]
V[3,2:3] <- v[3:4]

z2 <- sum(rain$data[1:nrow(rain$data),2] > 30, na.rm = TRUE) 

u = 30
sc = 7.4403
sh = 0.1845

n = nrow(rain$data)
zeta <- z2/n
V[1,1] <- zeta*(1-zeta)/n
N = 10
m = nrow(rain$data)*N/48

zN = u + sc/sh*((m*z)^sh - 1)
dz = c(sc*m^sh*z^(sh-1), sh^-1*((m*z)^sh -1), -sc*sh^-2*((m*z)^sh - 1) + sc*sh^-1*(m*z)^sh*log(m*z))
v_z = dz%*%V%*%dz
# 12
# Answer: 10y: 65.96, var 27.6, 100y: 106.3, var 435
# CI_10 = (66 - 1.96*sqrt(27.6), 66 + 1.96*sqrt(27.6))
# Same for CI_100

# 13 10y: 65.96 CI (55.893, 76.0299), 100y: 106.3 CI (65.6224, 147.0622)

# 14 profile likelihood 10y: 65.96 CI (58.6417, 81.1786), 100y 106.34 CI (80.973, 184.582)

# 15 They seem to fit reasonably well but in the extreme, they deviate a bit.

# 16 access the whole dataset by rain$models$fit1$x

# 16 ---> wisdom, the sigma absorb sigma, so dont normalize, use n number of data points
x = rain$models$fit1$x[rain$models$fit1$x-30>0]
x = x - 30
s = sum(x)
n_exc = length(x)

l_lh = -n_exc*(1 + log(s/n_exc))

# 17 
fevd(x = value, data = rain$data, threshold = 30, type = "Exponential")
# 17 gives same results.

# 18
v_sc = 0.5429136
sc = 9.084211
CI = c(sc - 1.96*sqrt(v_sc)/sqrt(n_exc), sc + 1.96*sqrt(v_sc)/sqrt(n_exc))
# SHOULD WE USE 1/SQRT(n_exc)???

# compare 485.09 with 487.4
p = 2*(487.4 - 485.9)
xi_1 = 3.84
# => stay in H0

x_sorted = sort(x)
grid = (1:n_exc)/(n_exc + 1)
G = qexp(grid, rate = 1, lower.tail = TRUE, log.p = FALSE)
# G_inv = 1/G
plot(x_sorted,G)










