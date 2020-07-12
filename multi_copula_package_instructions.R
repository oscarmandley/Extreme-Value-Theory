library("copula")
> myCop.norm <- ellipCopula(family = "norma", dim = 3, dispstr = "ex",
                                             + param = 0.4)
> set.seed(1)
> myCop.t <- ellipCopula(family = "t", dim = 3, dispstr = "toep",
                         + param = c(0.8, 0.5), df = 8)
> myCop.clayton <- achmCopula(family = "clayton", dim = 3, param = 2)