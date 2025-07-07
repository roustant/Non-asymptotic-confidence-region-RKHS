kH1Fun <- function(x, y){
  cosh(pmin(x, y)) * cosh(1 - pmax(x, y)) / sinh(1)
}

kH1 <- function(x, y){
  outer(x, y, kH1Fun)
}

lambda <- function(N){ 
  pi^2 * (0:N)^2 
}

basisFun <- function(x, N){
  res <- 1
  if (N >= 1) res <- c(res, sqrt(2) * cos(pi * (1:N) * x))
  return(matrix(res, ncol = 1))
}
basisFunDer <- function(x, N){
  res <- 0
  if (N >= 1) res <- c(res, - sqrt(2) * sin(pi * (1:N) * x) * pi * (1:N))
  return(matrix(res, ncol = 1))
}

L <- 20-1
set.seed(1)
w <- matrix(rnorm(L+1) * (1:(L+1))^{-2}, ncol = 1)

fFun <- function(x){
  t(basisFun(x, L)) %*% w
}
fFunDer <- function(x){
  t(basisFunDer(x, L)) %*% w
}


beta <- 1 + lambda(L)
fNorm <- sqrt( sum( beta * w^2) ) 

# modify fFun and fFunDer such that it accepts vectors as arguments
f <- function(x) sapply(x, fFun) 
fDer <- function(x) sapply(x, fFunDer) 

# graphical representation of the test function and its derivative
t <- seq(0, 1, length = 200)
par(mfrow = c(1, 1))
plot(t, f(t), lty = 1, type = "l")
plot(t, fDer(t), lty = 1, type = "l")

fSq <- function(x) {
  f(x)^2 + fDer(x)^2
}


# evaluation de b = sup de f^2 + f'^2
library(rgenoud)
resOpt <- genoud(fSq, nvars = 1, max = TRUE, max.generations = 5, 
                 Domains = matrix(c(0, 1), nrow = 1, ncol = 2))

par(mfrow = c(1, 1))
plot(t, fSq(t), lty = 1, type = "l")
points(resOpt$par, resOpt$value, pch = 19, col = "blue")
b <- resOpt$value
