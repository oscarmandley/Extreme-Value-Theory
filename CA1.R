a1 = 98765 - 43210
a2 = 12.3456789^2
a3 = sqrt((1.23+4.56)*(7.8+9.0))

b1 = 3+4*5
b2 = 3/4/5
b3 = 1:20/2
b4 = (1:20)/2
b5 = 1:(20/2)
b6 = c(4,2,4)^2 - c(1,3)
b7 = (1:5)^(1:2)

c1 = (1:50)*2
c2 = c((1:5),(4:1))
c3 = rep((10*(1:4)),3)

#d
y <- rnorm(100)
d1 = which(y > 1.96)
d2 = which(y < 0)
d3 = y[(1:50)+50]/y[(1:50)]

ibm <- scan("http://www.maths.lth.se/matstat/kurser/fmsn15masm23/datasets/ibm.txt")

e1 = mean(ibm)
e2 = sd(ibm)
e3 = max(ibm)
e4 = min(ibm)
e5 = which(ibm == e3)
ibm[e5] = 65.75
# recalculate
ibm_m = mean(ibm)
ibm_sd = sd(ibm)
ibm_max = max(ibm)
ibm_min = min(ibm)

ibm_return = ibm[2:length(ibm)]/ibm[((1:length(ibm)-1))]
lreturn = log(ibm_return)
hist(lreturn, breaks=30)
t.test(lreturn, mu = 0)
t.test(lreturn, mu = 0, conf.level = 0.99)

# Exercise 5
f_grid = (0:100)/10 - 5
f = f_grid^2
plot(f_grid, f, type = "l")

x <- seq(-5, 5, length= 201)
fx <- x^2
plot(x, fx, type = "l", ylab = "f(x)", col="red", main="Graph of f(x) = x^2")

# Exercise 6
x <- seq(-3, 3, length=51)
y <- seq(-3,3, length=51)
fxy <- exp(-outer(x^2, y^2, "+"))

persp(x,y,fxy, col = "cyan", theta = 45, main="Perspective of f(x,y) = exp(-(x^2 + y^2))")
contour(x,y,fxy, xlab = "x", ylab = "y", nlevels = 25, main="Contour of f(x,y) = exp(-(x^2 + y^2))")
image(x, y, fxy, col = terrain.colors(50), main = "Image of f(x,y) = exp(-(x^2 + y^2))")


# %%%%%%%%%%%%%%%%%%%%%%%%%% EXERCISE 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set.seed(1)
n = 0
for (k in 1:3){
  draw <-matrix(0, ncol=6, nrow = 3016)
  for (i in 1:3016){
    for (j in 1:6){
      draw[i,1:6] <- sample(1:49, 6, replace = F)
    }
  }
  draw_hash <-matrix(0, ncol=1, nrow = 3016)
  for (i in 1:3016){
    draw_hash[i] = sum(2^draw[i,1:6])
  }
  if (length(unique(draw_hash)) == 3016){
    n = n + 1
  }
}
p_hat = n/300

n = choose(49,6)
p_analyt = 1
for (i in 1:3015){
  p_analyt = p_analyt/n*(n-i)
}

# %%%%%%%%%%%%%%%%%%%%%%%%%% EXERCISE 8 %%%%%%%%%%%%%%%%%%%%%%%%%%

nr_algo <- function(y, x0, eps){
  if (nargs() < 2) {
    x0 = 1
  }
  if (nargs() < 3) {
    eps = 10^-4
  }
  x = x0
  print(nargs())
  while (abs(x^3 - y) >= eps)
  {
    x = (2*x + y/x/x)/3
    options(digits=10)
    print(x)
  }
}

# Example
q = 6
#xx = nr_algo(q)
q^(1/3)

























