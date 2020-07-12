library(sn);library(mvtnorm); library(scatterplot3d); library(copula); library(EnvStats);
library(plotly); library(mdatools); library(evd); library(SimCop); library(fitdistrplus);
library(xtable)
set.seed(1)
n = 500

# E2.1 ====
rho = 0.4
myCop.norm <- ellipCopula(family = "normal", dim = 3, dispstr = "ex", param = rho)
n.u <- rCopula(n, myCop.norm)
scatterplot3d(n.u, xlab = "x", ylab = "y", zlab = "z", main="Observations from normal copula")
# Empirical
n.rho_p <- cor(n.u)
n.tau <- cor(n.u, method = "kendall")
n.rho_s <- cor(n.u, method = "spearman")
# Theoretical
n.rho_p_teo = cbind(c(1,rho,rho),c(rho,1,rho),c(rho,rho,1))
n.tau_teo = 2/pi*asin(n.rho_p_teo)
n.rho_s_teo = 6/pi*asin(n.rho_p_teo/2)
# Agreement? Error
n.rho_p_err = n.rho_p - n.rho_p_teo
n.tau_err = n.tau - n.tau_teo
n.rho_s_err = n.rho_s - n.rho_s_teo
# Rank calculations
rr = matrix(0, nrow = 500, ncol = 3)
rr[,1] = rank(n.u[,1])
rr[,2] = rank(n.u[,2])
rr[,3] = rank(n.u[,3])
cor(rr)

# E2.2 ====
th1 = -0.4
th2 = 0
th3 = 0.4
myCop.t1 <- ellipCopula(family = "t", dim = 2, dispstr = "toep", param = c(th1), df = 8)
myCop.t2 <- ellipCopula(family = "t", dim = 2, dispstr = "toep", param = c(th2), df = 8)
myCop.t3 <- ellipCopula(family = "t", dim = 2, dispstr = "toep", param = c(th3), df = 8)
t.u1 <- rCopula(n, myCop.t1)
t.u2 <- rCopula(n, myCop.t2)
t.u3 <- rCopula(n, myCop.t3)

par(mfrow=c(2,2))
plot(t.u1, xlab = "x", ylab = "y", main="obs from t copula with th1 = -0.4")
plot(t.u2, xlab = "x", ylab = "y", main="t copula with th2 = 0")
plot(t.u3, xlab = "x", ylab = "y",  main="t copula with th3 = 0.4")

# th1
t.1_rho_p <- signif(cor(t.u1) ,digits=3) 
t.1_tau <- signif(cor(t.u1, method = "kendall") ,digits=3) 
t.1_rho_s <- signif( cor(t.u1, method = "spearman"),digits=3)
t.1_rho_p_teo = signif(cbind(c(1,th1),c(th1,1)),digits=3) 
t.1_tau_teo = signif(2/pi*asin(t.1_rho_p),digits=3) 
t.1_rho_s_teo =signif(6/pi*asin(t.1_rho_p/2) ,digits=3) 

# th2
t.2_rho_p <- signif(cor(t.u2) ,digits=3)
t.2_tau <- signif(cor(t.u2, method = "kendall") ,digits=3) 
t.2_rho_s <- signif(cor(t.u2, method = "spearman") ,digits=3)
t.2_rho_p_teo = signif(cbind(c(1,th2),c(th2,1)),digits=3) 
t.2_tau_teo = signif(2/pi*asin(t.2_rho_p) ,digits=3) 
t.2_rho_s_teo = signif(6/pi*asin(t.2_rho_p/2) ,digits=3)

# th3
t.3_rho_p <- signif(cor(t.u3) ,digits=3) 
t.3_tau <- signif(cor(t.u3, method = "kendall") ,digits=3) 
t.3_rho_s <- signif(cor(t.u3, method = "spearman") ,digits=3) 
t.3_rho_p_teo = signif(cbind(c(1,th3),c(th3,1)),digits=3) 
t.3_tau_teo = signif(2/pi*asin(t.3_rho_p_teo) ,digits=3) 
t.3_rho_s_teo = signif(6/pi*asin(t.3_rho_p_teo/2) ,digits=3) 
t.3_rho_s_teo = signif(t.3_rho_s_teo,digits=3)

# Pearson - exporting to latex
rho_p <- cbind(t.1_rho_p, t.1_rho_p_teo, t.2_rho_p, t.2_rho_p_teo, t.3_rho_p, t.3_rho_p_teo)
write.table(rho_p, "rho_p_table.txt", quote=FALSE, eol="\\\\\n", sep=" & ")
# Kendall - exporting to latex
tau_matrix <- cbind(t.1_tau, t.1_tau_teo, t.2_tau, t.2_tau_teo, t.3_tau, t.3_tau_teo)
write.table(tau_matrix, "tau_table.txt", quote=FALSE, eol="\\\\\n", sep=" & ")
# Spearman - exporting to latex
rho_s <- cbind(t.1_rho_s, t.1_rho_s_teo, t.2_rho_s, t.2_rho_s_teo, t.3_rho_s, t.3_rho_s_teo)
write.table(rho_s, "rho_s_table.txt", quote=FALSE, eol="\\\\\n", sep=" & ")

# E2.3 ====
frank.cop <- frankCopula(param = 3, dim = 3)
f.u <- rCopula(n, frank.cop)
scatterplot3d(f.u)
f.rho_p <- cor(f.u)
f.tau <- cor(f.u, method = "kendall")
f.rho_s <- cor(f.u, method = "spearman")

data_x = round(cbind(f.rho_p, f.tau, f.rho_s), 3)
xtable(data_x)
t <- round(tau(frank.cop), 3)

f.tau_teo = 2/pi*asin(f.rho_p)
f.rho_s_teo = 6/pi*asin(f.rho_p/2)
f.rho_s_error = f.rho_s_teo - f.rho_s
f.tau_error = f.tau_teo - f.tau

# E2.4 Create funct ====
gen_ex4 <- function(n_samples, dim, theta){
  #  l = 1/((1+x)^(1/theta))    #  l_inv = (1/x)^theta - 1
  E <- rexp(n = dim*n_samples, rate = 1)
  E <- matrix(E, ncol = dim, byrow = TRUE)
  M <- rgamma(n = n_samples, shape = 1/theta, scale = 1)
  U = 1/((1+(E/M))^(1/theta))
}

n_samples = 500; dim = 3
theta = 0.5
v1 <- gen_ex4(n_samples, dim, theta)
theta = 1
v2 <- gen_ex4(n_samples, dim, theta)
theta = 5
v3 <- gen_ex4(n_samples, dim, theta)
theta = 50
v4 <- gen_ex4(n_samples, dim, theta)

par(mfrow=c(2,2)) 
par(mar = rep(2, 4))
plot(v1, main = "th = 0.5")
plot(v2, main = "th = 1")
plot(v3, main = "th = 5")
plot(v4, main = "th = 50")

# E2.5 ==== 
x = runif(1000, -1, 1)
y = atanh(x)
z = -y

cor(x,y, method = "pearson")
cor(x,y, method = "kendall")
cor(x,y, method = "spearman")

cor(x,z, method = "pearson")
cor(x,z, method = "kendall")
cor(x,z, method = "spearman")

par(mfrow=c(1,2))
plot(y, main = "Scatterplot of y")
plot(z, main = "Scatterplot of z")


par(mfrow=c(1,2))
a <- plot(x, y, type = "p", ylab = "y = T(x)", col="black", main="T(X) = arctanh(X)")
b <- plot(x, z, type = "p", ylab = "z = -T(x)", col="black", main="T(X) = -arctanh(X)")

# E3 d den, p dist ====
par(mfrow=c(2,2)); par(mar = rep(2, 4))
persp(myCop.t1, dCopula, main = "t1, density", xlab = "u1", ylab = "u2")
persp(myCop.t2, dCopula, main = "t2, density", xlab = "u1", ylab = "u2")
contour(myCop.t1, pCopula, main = "t1, distribution", xlab = "u1", ylab = "u2")
contour(myCop.t2, pCopula, main = "t2, distribution", xlab = "u1", ylab = "u2")

# E4 ====
n = 500
set.seed(1)

# Gaussian
myCop.norm <- ellipCopula(family = "normal", dim = 2, dispstr = "ex", param = 0.7)
n.u <- rCopula(n, myCop.norm)
# Gumbel
myCop.gumbel <- archmCopula(family = "gumbel", dim = 2, param = 2)
g.u <- rCopula(n, myCop.gumbel)
# Clayton
myCop.clayton <- archmCopula(family = "clayton", dim = 2, param = 2.2)
c.u <- rCopula(n, myCop.clayton)
# t Copula
myCop.t <- ellipCopula(family = "t", dim = 2, dispstr = "toep", param = 0.71, df = 4)
t.u <- rCopula(n, myCop.t)

t.u_evaluate <- cbind(dCopula(t.u, myCop.t), pCopula(t.u, myCop.t))

par(mfrow=c(2,2)) 
par(mar = rep(2, 4))
a <- plot(n.u[,1], n.u[,2], type = "p", main = "Normal")
b <- plot(g.u[,1], g.u[,2], type = "p", main = "Gumbel") 
c <- plot(c.u[,1], c.u[,2], type = "p", main = "Clayton") 
d <- plot(t.u[,1], t.u[,2], type = "p", main = "t") 

# inverse_normal = qnorm(n.u)
par(mfrow=c(2,2)) 
par(mar = rep(2, 4))
plot(qnorm(n.u), main = "Normal")
plot(qnorm(g.u), main = "Gumbel")
plot(qnorm(c.u), main = "Clayton")
plot(qnorm(t.u), main = "t")
