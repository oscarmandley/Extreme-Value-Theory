library(sn);library(mvtnorm); library(scatterplot3d); library(copula); library(EnvStats); library(plotly)
library(mdatools); library(evd); library(SimCop); library(fitdistrplus)
# E7 ====
dat_cars <- read.csv("http://www.maths.lth.se/matstat/kurser/fmsn15masm23/datasets/cardata.txt")
pairs(dat_cars)
round(cor(dat_cars, method ="spearman"),2)

dat_horse = dat_cars[,3]
dat_price = dat_cars[,1]
dat_MPG = dat_cars[,2]

cor.test(dat_horse, dat_price,
         alternative = "less",
         method = "spearman",
         exact = NULL, conf.level = 0.95)

cor.test(dat_MPG, dat_price,
         alternative = "greater",
         method = "spearman",
         exact = NULL, conf.level = 0.95)