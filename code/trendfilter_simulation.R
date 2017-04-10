library(genlasso)
library(glmgen)
source("code/trendfilter.admm.R")
source("code/trendfilter.fixpt.R")

## fused lasso
set.seed(777)
n <- 100
i <- 1:n
truth <- (i > 0 & i <= 20) * 1 +
  (i > 20 & i <= 40) * -5 +
  (i > 40 & i <= 60) * 0 +
  (i > 60 & i <= 80) * 5 +
  (i > 80 & i <= 100) * -1
z = truth + rnorm(n, sd = 0.3)
plot(truth, col = "red", type = "l", bty = "n", ylim = range(c(z, truth)))
points(truth, col = "red", pch = 19, cex = 0.5)
points(z, pch=19, cex = 0.2)

w = rep(1, n)
D = getD1d(n)


fit.genlasso.init = genlasso::trendfilter(y = z, ord = 0)
beta.init = coef(fit.genlasso.init, lambda = lambda, type = "both")$beta
rho.init = coef(fit.genlasso.init, lambda = lambda, type = "both")$u
cv.tf = cv.trendfilter(fit.genlasso.init)
lambda = cv.tf$lambda.min
fit.genlasso = glmgen::trendfilter(z, k = 0, family='gaussian', lambda=lambda)
fit.prox.l1 = tfprox(z, w, D, penalty = "l1", lambda = lambda, objtol = 1e-10, max.iter = 1e4)
fit.genlasso.init = genlasso::trendfilter(y = z, ord = 0)
beta.init = coef(fit.genlasso.init, lambda = lambda, type = "both")$beta
rho.init = coef(fit.genlasso.init, lambda = lambda, type = "both")$u
fit.prox.l1.init = tfprox(z, w, D, penalty = "l1", lambda = lambda, beta.init = beta.init, rho.init = rho.init, objtol = 1e-10, max.iter = 1e4)
fit.prox.dp = tfprox(z, w, D, penalty = "double-pareto", lambda = lambda, objtol = 1e-10, max.iter = 1e4)
fit.prox.dp.init = tfprox(z, w, D, penalty = "double-pareto", lambda = lambda, beta.init = beta.init, rho.init = rho.init, objtol = 1e-10, max.iter = 1e4)
points(fit.genlasso$beta, pch = 19, cex = 0.5)
points(fit.prox.l1$beta, pch = 21, cex = 0.5, col = "green")
points(fit.prox.dp$beta, pch = 19, cex = 0.5, col = "blue")
points(fit.prox.dp$beta, pch = 19, cex = 0.5, col = "blue")

mean((fit.genlasso$beta - truth)^2)
mean((fit.prox.l1$beta - truth)^2)
mean((fit.prox.dp$beta - truth)^2)


## linear trend filter
set.seed(777)
n <- 100
i <- 1:n
truth <- (i > 0 & i <= 20) * (2 * i) +
  (i > 20 & i <= 40) * (-1 * i) +
  (i > 40 & i <= 60) * 0 +
  (i > 60 & i <= 80) * (3 * i - 220) +
  (i > 80 & i <= 100) * (-4 * i + 400)
z = truth + rnorm(n, sd = 3)
plot(truth, col = "red", type = "l", bty = "n", ylim = range(c(z, truth)))
points(truth, col = "red", pch = 19, cex = 0.5)
points(z, pch=19, cex = 0.2)


w = rep(1, n)
D = getDtf(n, ord = 1)

fit.genlasso.init = genlasso::trendfilter(y = z, ord = 1)
cv.tf = cv.trendfilter(fit.genlasso.init)
lambda = cv.tf$lambda.min

fit.genlasso = glmgen::trendfilter(z, k = 1, family='gaussian', lambda=lambda)
fit.prox.l1 = tfprox(z, w, D, penalty = "l1", lambda = lambda, objtol = 1e-15, max.iter = 1e4)



fit.prox.dp = tfprox(z, w, D, penalty = "double-pareto", lambda = lambda * 10, objtol = 1e-10, max.iter = 1e4)

beta.init = coef(fit.genlasso.init, lambda = lambda * 10, type = "both")$beta
rho.init = coef(fit.genlasso.init, lambda = lambda * 10, type = "both")$u

fit.prox.dp.init = tfprox(z, w, D, penalty = "double-pareto", lambda = lambda * 10, beta.init = beta.init, rho.init = rho.init, objtol = 1e-15, max.iter = 2e4)



points(fit.genlasso$beta, pch = 19, cex = 0.5, col = "green")

points(fit.prox.l1$beta, pch = 21, cex = 0.5, col = "green")

points(fit.prox.dp$beta, pch = 19, cex = 0.5, col = "blue")
points(fit.prox.dp.init$beta, pch = 19, cex = 0.5, col = "purple")

mean((fit.genlasso$beta - truth)^2)
mean((fit.prox.dp.init$beta - truth)^2)


## piecewise linear

set.seed(777)
n = 100
i = 1:n
true = (i < 20) * (i * 0.3 + 0.5) + (i >= 20 & i < 30) * (-.1 * i + .7) + (i >= 30 & i < 50) * 3 + (i >= 50 & i < 70) * (-.5 * i + 40) + (i >= 70) * (.1 * i - 5)
plot(true)

D = getDtf(n, 1)
w = rep(1, n)

z = true  +   rnorm(n, sd=0.25)

fit.genlasso.init = genlasso::trendfilter(y = z, ord = 1)
cv.tf = cv.trendfilter(fit.genlasso.init)
lambda = cv.tf$lambda.min



lambda = 1
plot(z, col=rgb(0.2,0.2,0.2,0.1), pch=19, bty = "n")
fit.genlasso = glmgen::trendfilter(z, k = 1, family='gaussian', lambda=lambda)
points(fit.genlasso$beta, pch = 19, cex = 0.5)
fit.prox.l1 = tfprox(z, w, D, penalty = "l1", lambda = lambda, objtol = 1e-10, max.iter = 1e4)
points(fit.prox.l1$beta, pch = 21, cex = 0.5, col = "green")


fit.prox.dp = tfprox(z, w, D, penalty = "double-pareto", lambda = lambda * 10, objtol = 1e-14, max.iter = 1e3)
beta.init = coef(fit.genlasso.init, lambda = lambda * 10, type = "both")$beta
rho.init = coef(fit.genlasso.init, lambda = lambda * 10, type = "both")$u
fit.prox.dp.init = tfprox(z, w, D, penalty = "double-pareto", lambda = lambda * 10, beta.init = beta.init, rho.init = rho.init, objtol = 1e-15, max.iter = 1e3)



points(fit.prox.dp$beta, pch = 21, cex = 0.5, col = "blue")



fit.genlasso.init = genlasso::trendfilter(y = z, ord = 1)
beta.init = coef(fit.genlasso.init, lambda = lambda, type = "both")$beta
rho.init = coef(fit.genlasso.init, lambda = lambda, type = "both")$u
fit.prox.dp.init = tfprox(z, w, D, penalty = "double-pareto", lambda = lambda, beta.init = beta.init, rho.init = rho.init, objtol = 1e-14, max.iter = 2e4)
points(fit.prox.dp$beta, pch = 19, cex = 0.5, col = "blue")
points(true, cex = 0.1, col = "red")

points(fit.prox.dp.init$beta, pch = 19, cex = 0.5, col = "purple")

mse.genlasso = mean((fit.genlasso$beta - true)^2)
mse.prox.l1 = mean((fit.prox.l1$beta - true)^2)
mse.prox.dp = mean((fit.prox.dp$beta - true)^2)
mse.prox.dp.init = mean((fit.prox.dp.init$beta - true)^2)
mse.genlasso
mse.prox.l1
mse.prox.dp


## quantile regression














## piecewise linear example

set.seed(777)
n = 100
k = 4
pt = (1:k) * round(n / (k + 1))
bpt = c()
for (i in 1:k) {
  bpt[i] = sample((pt[i] - round(n / (3 * (k + 1)))) : (pt[i] + round(n / (3 * (k + 1)))), 1)
}
bpt = c(1, bpt, n)
b = 0.5
true = c()
for (j in 1:(k + 1)) {
  slope = runif(1, -b, b)
  a = runif(1, -5, 5)
  for (i in 1:n) {


    if (i >= bpt[j] & i <= bpt[j + 1]) {
      true[i] = slope * (i - bpt[j]) + a
    }
  }
}
plot(true)
z = true + rnorm(n, 0, sd = 0.1)

n = 1000
p = 0.99
sigma = 20
b = 0.5

x = 0
z = c()
v = runif(1, -b, b)
true = c()

for (t in 1:n) {
  z[t] = x + rnorm(1, 0, sigma)
  true[t] = x
  I = runif(1) < p
  x = x + v
  v = I * v + (1 - I) * runif(1, -b, b)
}
plot(true, col = "red", type = "l")
D = getDtf(n, 1)
w = rep(1, n)
points(z, col=rgb(0.2,0.2,0.2,0.1), pch=19, bty = "n")


lambda = 10000
fit.genlasso = glmgen::trendfilter(z, k = 1, family='gaussian', lambda=lambda)
lines(fit.genlasso$beta, pch = 19, cex = 0.5)



fit.genlasso.init = genlasso::trendfilter(y = z, ord = 1)
beta.init = coef(fit.genlasso.init, lambda = lambda, type = "both")$beta
rho.init = coef(fit.genlasso.init, lambda = lambda, type = "both")$u
fit.prox.dp.init = tfprox(z, w, D, penalty = "double-pareto", lambda = lambda, beta.init = beta.init, rho.init = rho.init, objtol = 1e-10, max.iter = 1e4)


fit.prox.l1 = tfprox(z, w, D, penalty = "l1", lambda = lambda, objtol = 1e-10, max.iter = 1e4)
fit.prox.l1.init = tfprox(z, w, D, penalty = "l1", lambda = lambda, beta.init = beta.init, rho.init = rho.init, objtol = 1e-10, max.iter = 1e4)




points(fit.prox.l1$beta, pch = 21, cex = 0.5, col = "green")
fit.prox.dp = tfprox(z, w, D, penalty = "double-pareto", lambda = lambda, objtol = 1e-14, max.iter = 2e4)
points(fit.prox.dp$beta, pch = 19, cex = 0.5, col = "blue")
points(true, cex = 0.1, col = "red")

mse.genlasso = mean((fit.genlasso$beta - true)^2)
mse.prox.l1 = mean((fit.prox.l1$beta - true)^2)
mse.prox.dp = mean((fit.prox.dp$beta - true)^2)
mse.genlasso
mse.prox.l1
mse.prox.dp



## piecewise second-order polynomial

set.seed(777)
n = 100
i = 1:n
true = ((i < 20) * (rnorm(10) * i^2 + rnorm(10) * i + rnorm(10))
        + (i >= 20 & i < 30) * (rnorm(1) * i^2 + rnorm(1) * i + rnorm(1))
        + (i >= 30 & i < 50) * rnorm(1)
        + (i >= 50 & i < 70) * (rnorm(1) * i^2 + rnorm(1) * i + rnorm(1))
        + (i >= 70) * (rnorm(1) * i^2 + rnorm(1) * i + rnorm(1)))
plot(true)

D = getDtf(n, 1)
w = rep(1, n)

z = true  +   rnorm(n, sd=0.1)

lambda = 1
plot(z, col=rgb(0.2,0.2,0.2,0.1), pch=19, bty = "n")
fit.genlasso = genlasso::trendfilter(y = z, ord = 1)
points(coef(fit.genlasso, lambda = lambda)$beta, pch = 19, cex = 0.5)




fit.prox.l1 = tfprox(z, w, D, penalty = "l1", lambda = lambda, objtol = 1e-10, max.iter = 1e4)
points(fit.prox.l1$beta, pch = 21, cex = 0.5, col = "green")
fit.prox.dp = tfprox(z, w, D, penalty = "double-pareto", lambda = lambda, objtol = 1e-14, max.iter = 2e4)
points(fit.prox.dp$beta, pch = 19, cex = 0.5, col = "blue")
points(true, cex = 0.1, col = "red")

mse.genlasso = mean((coef(fit.genlasso, lambda = lambda)$beta - true)^2)
mse.prox.l1 = mean((fit.prox.l1$beta - true)^2)
mse.prox.dp = mean((fit.prox.dp$beta - true)^2)
mse.genlasso
mse.prox.l1
mse.prox.dp
