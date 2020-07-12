library(sn);library(mvtnorm); library(scatterplot3d); library(EnvStats); library(plotly)
library(mdatools); library(fitdistrplus); library(xtable); library(in2extRemes)
library(evd); library(SimCop); library(copula); library(mgpd)

# E4 Preprocessing ====
source("fremantle.R.txt")
source("portpirie.R.txt")
max(fremantle[,2])
max(portpirie[,2])
slMerged <- merge(fremantle, portpirie, by = "Year")
sl <- data.frame(year = slMerged[,1], fremantle=slMerged[,2], Portpirie = slMerged[,4])

# E4 ====
# Fit margins and evd models
f1 <- fgev(sl[,2]); f2 <- fgev(sl[,3])
M1 <- fbvevd(sl[,2:3], model = "hr")
M2 <- fbvevd(sl[,2:3], model = "log")
M3 <- fbvevd(sl[,2:3], model = "neglog")
sl <- as.matrix(sl)

# Initiate copulas
tau <- 1/3
hr.cop          <- huslerReissCopula(iTau(huslerReissCopula(), tau))
gumbel.cop      <- gumbelCopula     (iTau(gumbelCopula(),      tau))
galambos.cop    <- galambosCopula   (iTau(galambosCopula(),    tau))

# Create Mbdc objects
myBvd_hr <- mvdc(copula = hr.cop,
                 margins = c("gev", "gev"), paramMargins = list(list(loc = 0.3, scale = 1, shape = 0.3),
                                                                list(loc = 0.3, scale = 1, shape = 0.3)))
myBvd_gumbel <- mvdc(copula = gumbel.cop,
                     margins = c("gev", "gev"), paramMargins = list(list(loc = 0.3, scale = 1, shape = 0.3),
                                                                    list(loc = 0.3, scale = 1, shape = 0.3)))
myBvd_galambos <- mvdc(copula = galambos.cop,
                       margins = c("gev", "gev"), paramMargins = list(list(loc = 0.3, scale = 1, shape = 0.3),
                                                                      list(loc = 0.3, scale = 1, shape = 0.3)))
# Give start values
start_hr <- c(f1$est, f2$est,2)
start_gumbel <- c(f1$est, f2$est,2)
start_galambos <- c(f1$est, f2$est,0.1)
start_gumbel
start_galambos

# Fit distribution using copulas
fit.hr <- fitMvdc(sl[,2:3], myBvd_hr, start = start_hr,
                  optim.control = list(trace = TRUE, maxit = 2000))
fit.gumbel <- fitMvdc(sl[,2:3], myBvd_gumbel, start = start_gumbel,
                      optim.control = list(trace = TRUE, maxit = 2000))
fit.galambos <- fitMvdc(sl[,2:3], myBvd_galambos, start = start_galambos,
                        optim.control = list(trace = TRUE, maxit = 2000))

# Compare estimates.
fit.hr@estimate
M1$e

fit.gumbel@estimate
M2$e

fit.galambos@estimate
M3$e

# E4.2 - Inference For Margins IFM
fgev_frem <- fgev(sl[,2], std.err = FALSE)
fgev_port <- fgev(sl[,3], std.err = FALSE)

n = 300
frem_pseudo <- rgev(n, loc = fgev_frem$e[1], scale = fgev_frem$e[2], shape = fgev_frem$e[3])
port_pseudo <- rgev(n, loc = fgev_port$e[1], scale = fgev_port$e[2], shape = fgev_port$e[3])
data_pseudo <- cbind(frem_pseudo, port_pseudo)

hr_fit_IFM <- fitCopula(copula = hr.cop, data = pobs(data_pseudo), method = "ml")
gumbel_fit_IFM <- fitCopula(copula = gumbel.cop, data = pobs(data_pseudo), method = "ml", optim.method="Brent", lower = 0, upper = 3)
galambos_fit_IFM <- fitCopula(copula = galambos.cop, data = pobs(data_pseudo), method = "ml")

# E4.2 - Canonical Maximum Likelihood (CML)
hr_fit_CML <- fitCopula(pobs(sl_am), copula = hr.cop) 
gumbel_fit_CML <- fitCopula(pobs(sl_am), copula = gumbel.cop) 
galambos_fit_CML <- fitCopula(pobs(sl_am), copula = galambos.cop) 

# E4.3 Goodness of fitt
gof1 <- gofEVCopula(hr_fit_IFM@copula, sl_am, N = 70, optim.method="Brent", lower = 0, upper = 3)
gof2 <- gofEVCopula(gumbel_fit_IFM@copula, sl_am, N = 70, optim.method="Brent", lower = 0, upper = 3)
gof3 <- gofEVCopula(galambos_fit_IFM@copula, sl_am, N = 70, optim.method="Brent", lower = 0, upper = 3)
gof4 <- gofEVCopula(hr_fit_CML@copula, sl_am, N = 70, optim.method="Brent", lower = 0, upper = 3)
gof5 <- gofEVCopula(gumbel_fit_CML@copula, sl_am, N = 70, optim.method="Brent", lower = 0, upper = 3)
gof6 <- gofEVCopula(galambos_fit_CML@copula, sl_am, N = 70, optim.method="Brent", lower = 0, upper = 3)

# E4.4 Non-parametric estimator of dependence function
w = (0:99)/99
DepEstCFG <- An.biv(sl_am, w, estimator = c("CFG"))
DepEstPick <- An.biv(sl_am, w, estimator = c("Pickands"))

plot(w, DepEstCFG, xlim = c(0,1), ylim = c(0.5,1), type = "b", col = "blue", main = "Non-para est of Dependence function", ylab = "Dependence function")
lines(w, DepEstPick, type = "b")
legend("bottomright", legend=c("CFG", "Pickands"), col=c("blue", "black"), lty=1:2, cex=0.8)
polygon(c(0, 0.5, 1, 0), c(1, 0.5, 1, 1))

# E3.6 and 3.7 Est prob ====
# Parametric estimation of probabilities
sl_am <- as.matrix(sl[,2:3])
# M3 <- fbvevd(sl_am, model = "bilog")
# M3$estimate

# mar1 = c(f1$e[1],f1$e[2],f1$e[3])
# mar2 = c(f2$e[1],f2$e[2],f2$e[3])
# a <- M3.IFM$e[1]
# b <- M3.IFM$e[2]

fhr = fit.hr@estimate
fhr
mar1 = c(fhr[1],fhr[2],fhr[3])
mar2 = c(fhr[4],fhr[5],fhr[6])
dep_p = fhr[7]

# P1
x <- c(30,4.2)
p_CD <- pbvevd(q = x, dep = dep_p, model = "hr",
               mar1 = mar1, mar2 = mar2)
x <- c(1.7,50)
p_BD <- pbvevd(q = x, dep = dep_p, model = "hr",
               mar1 = mar1, mar2 = mar2)
x <- c(1.7,4.2)
p_D <- pbvevd(q = x, dep = dep_p, model = "hr",
              mar1 = mar1, mar2 = mar2)

p_A = 1 - p_CD - p_BD + p_D
round(p_A,4)

# P2
x <- c(30,4.4)
p_CD <- pbvevd(q = x, dep = dep_p, model = "hr",
               mar1 = mar1, mar2 = mar2)
x <- c(1.8,50)
p_BD <- pbvevd(q = x, dep = dep_p, model = "hr",
               mar1 = mar1, mar2 = mar2)
x <- c(1.8,4.4)
p_D <- pbvevd(q = x, dep = dep_p, model = "hr",
              mar1 = mar1, mar2 = mar2)

p_A = 1 - p_CD - p_BD + p_D
round(p_A,5)

# P3
x <- c(1.478,3.85)
p_14_38 <- pbvevd(q = x, dep = dep_p, model = "hr",
                  mar1 = mar1, mar2 = mar2)
p_14_38 = 1 - p_14_38
round(p_14_38,4)

# Parametric estimation of p4
x <- c(1.478,3.850)
p_D_1 <- pbvevd(q = x, dep = dep_p, model = "hr",
                mar1 = mar1, mar2 = mar2)
x <- c(1.95,4.8)
p_D_2 <- pbvevd(q = x, dep = dep_p, model = "hr",
                mar1 = mar1, mar2 = mar2)
p = ((1-p_D_1) - (1-p_D_2))/(1-p_D_1)
p

# Gumbel probabilties
fgu = fit.gumbel@estimate
fgu
mar1 = c(fgu[1],fgu[2],fgu[3])
mar2 = c(fgu[4],fgu[5],fgu[6])
dep_p = fgu[7]

# P1
x <- c(30,4.2)
p_CD <- pbvevd(q = x, dep = 1/dep_p, model = "log",
               mar1 = mar1, mar2 = mar2)
x <- c(1.7,50)
p_BD <- pbvevd(q = x, dep = 1/dep_p, model = "log",
               mar1 = mar1, mar2 = mar2)
x <- c(1.7,4.2)
p_D <- pbvevd(q = x, dep = 1/dep_p, model = "log",
              mar1 = mar1, mar2 = mar2)

p_A = 1 - p_CD - p_BD + p_D
round(p_A,4)

# P2
x <- c(30,4.4)
p_CD <- pbvevd(q = x, dep = 1/dep_p, model = "log",
               mar1 = mar1, mar2 = mar2)
x <- c(1.8,50)
p_BD <- pbvevd(q = x, dep = 1/dep_p, model = "log",
               mar1 = mar1, mar2 = mar2)
x <- c(1.8,4.4)
p_D <- pbvevd(q = x, dep = 1/dep_p, model = "log",
              mar1 = mar1, mar2 = mar2)

p_A = 1 - p_CD - p_BD + p_D
round(p_A,5)

# P3
x <- c(1.478,3.85)
p_14_38 <- pbvevd(q = x, dep = 1/dep_p, model = "log",
                  mar1 = mar1, mar2 = mar2)
p_14_38 = 1 - p_14_38
round(p_14_38,4)

# Parametric estimation of p4

x <- c(1.478,3.850)
p_D_1 <- pbvevd(q = x, dep = 1/dep_p, model = "log",
                mar1 = mar1, mar2 = mar2)
x <- c(1.95,4.8)
p_D_2 <- pbvevd(q = x, dep = 1/dep_p, model = "log",
                mar1 = mar1, mar2 = mar2)
p = ((1-p_D_1) - (1-p_D_2))/(1-p_D_1)
p

# Nonparametric dependence function using CFG
fhr = fit.hr@estimate
fhr
mar1 = c(fhr[1],fhr[2],fhr[3])
mar2 = c(fhr[4],fhr[5],fhr[6])

# 3.6 and 3.7 Non-Parametric estimation.
# Parametric transformation
# Estimatation of p1
y1 = c(3, 1.7, 1.7)
y2 = c(4.2, 5, 4.2)
G1 <- pgev(y1, mar1[1], mar1[2], mar1[3])
G2 <- pgev(y2, mar2[1], mar2[2], mar2[3])
p <- exp(log(G1*G2)*DepEstCFG[ceiling(log(G2)/log(G1*G2)*100)])
p1 = 1 + p[3] - p[1] - p[2]
round(p1,4)

# Estimation of p2
y1 = c(3, 1.8, 1.8)
y2 = c(4.4, 5, 4.4)
G1 <- pgev(y1, mar1[1], mar1[2], mar1[3])
G2 <- pgev(y2, mar2[1], mar2[2], mar2[3])
p <- exp(log(G1*G2)*DepEstCFG[ceiling(log(G2)/log(G1*G2)*100)])
p2 = 1 + p[3] - p[1] - p[2]
round(p2,4)

# Estimation of p3
y1 = 1.478
y2 = 3.85
G1 <- pgev(y1, mar1[1], mar1[2], mar1[3])
G2 <- pgev(y2, mar2[1], mar2[2], mar2[3])
p <- exp(log(G1*G2)*DepEstCFG[ceiling(log(G2)/log(G1*G2)*100)])
p3 <- 1-p
round(p3,3)

# Estimation of p4
y1 = 1.478
y2 = 3.85
G1 <- pgev(y1, mar1[1], mar1[2], mar1[3])
G2 <- pgev(y2, mar2[1], mar2[2], mar2[3])
p_D_1 <- exp(log(G1*G2)*DepEstCFG[ceiling(log(G2)/log(G1*G2)*100)])

y1 = 1.95
y2 = 4.8
G1 <- pgev(y1, mar1[1], mar1[2], mar1[3])
G2 <- pgev(y2, mar2[1], mar2[2], mar2[3])
p_D_2 <- exp(log(G1*G2)*DepEstCFG[ceiling(log(G2)/log(G1*G2)*100)])

p4 = ((1-p_D_1) - (1-p_D_2))/(1-p_D_1)
p4

# Gumbel probabilties -----------------------
fgu = fit.gumbel@estimate
fgu
mar1 = c(fgu[1],fgu[2],fgu[3])
mar2 = c(fgu[4],fgu[5],fgu[6])
dep_p = fgu[7]

# Empirical transformation
# Estimatation of p1
y1 = c(3, 1.7, 1.7)
y2 = c(4.2, 5, 4.2)
G1 <- pgev(y1, mar1[1], mar1[2], mar1[3])
G2 <- pgev(y2, mar2[1], mar2[2], mar2[3])
p <- exp(log(G1*G2)*DepEstCFG[ceiling(log(G2)/log(G1*G2)*100)])
p1 = 1 + p[3] - p[1] - p[2]
round(p1,4)

# Estimation of p2
y1 = c(3, 1.8, 1.8)
y2 = c(4.4, 5, 4.4)
G1 <- pgev(y1, mar1[1], mar1[2], mar1[3])
G2 <- pgev(y2, mar2[1], mar2[2], mar2[3])
p <- exp(log(G1*G2)*DepEstCFG[ceiling(log(G2)/log(G1*G2)*100)])
p2 = 1 + p[3] - p[1] - p[2]
round(p2,4)

# Estimation of p3
y1 = 1.478
y2 = 3.85
G1 <- pgev(y1, mar1[1], mar1[2], mar1[3])
G2 <- pgev(y2, mar2[1], mar2[2], mar2[3])
p <- exp(log(G1*G2)*DepEstCFG[ceiling(log(G2)/log(G1*G2)*100)])
p3 <- 1-p
round(p3,3)

# Estimation of p4
y1 = 1.478
y2 = 3.85
G1 <- pgev(y1, mar1[1], mar1[2], mar1[3])
G2 <- pgev(y2, mar2[1], mar2[2], mar2[3])
p_D_1 <- exp(log(G1*G2)*DepEstCFG[ceiling(log(G2)/log(G1*G2)*100)])

y1 = 1.95
y2 = 4.8
G1 <- pgev(y1, mar1[1], mar1[2], mar1[3])
G2 <- pgev(y2, mar2[1], mar2[2], mar2[3])
p_D_2 <- exp(log(G1*G2)*DepEstCFG[ceiling(log(G2)/log(G1*G2)*100)])

p4 = ((1-p_D_1) - (1-p_D_2))/(1-p_D_1)
p4