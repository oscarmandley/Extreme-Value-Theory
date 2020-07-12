# Ignore warnings
slNonparEst <- NonparEstDepFct(sl_am, convex = FALSE)
slNonparEst_c <- NonparEstDepFct(sl[,2:3], convex = TRUE)

# Non-parametric fit!!!
slSplfit <- SplineFitDepFct(slNonparEst)
lines(slNonparEst_c, col = 2)
curve(slSplfit, n = 301, add = TRUE, lty ="dashed", col = 3)
slSplfitCop <- NewBEVSplineCopula(slSplfit)
# New bivariate extreme value spine copula
slSplfitCopApprox1 <- GetApprox(slSplfitCop, type = 1)
slSplfitCopApprox2 <- GetApprox(slSplfitCop, type = 2)
plot(slSplfitCopApprox1, xlab = expression(u[1]), ylab = expression(u[2]), col = heat.colors(100))
plot(slSplfitCopApprox2, xlab = expression(u[1]), ylab = expression(u[2]), col = heat.colors(100))



## How to generate observations # GenerateRV
# dont bother with Metropolis hastings, use smoosmning splines and generate observations.
# data must be matrix not data frame.
Sample_np_1 <- GenerateRV(slSplfitCopApprox1, 100, type = 1)
Sample_np_2 <- GenerateRV(slSplfitCopApprox2, 100, type = 2)

Sample_rv_1 <- GenerateRV.CopApprox(slSplfitCopApprox1, 100)
Sample_rv_2 <- GenerateRV.CopApprox(slSplfitCopApprox2, 100)