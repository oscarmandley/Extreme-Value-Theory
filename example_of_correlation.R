time = c(1,2,3,4,5,6)
income = c(2,5,2,4,7,3)

cor(time, income)
cor.test(time,income, alternative = "less")
