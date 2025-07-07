rm(list = ls())

library(rgenoud)

source("confRegionRKHS_Commons.R")

pdfPlot <- FALSE

sinc <- function(u) {
  # sinus cardinal normalisee. Son integrale vaut 1 sur R.
  ifelse (u == 0, 1, sin(u) / (u))
}

kPWfun <- function(x, y){
  eta <- 30
  sinc(eta * (x-y) ) * eta / pi
  # the normalizing factor eta / pi guarantees 
  # that the RKHS norm is the L2 norm
}

kPW <- function(x, y){
  # returns the matrix (kPW(x_i, y_j))
  outer(x, y, kPWfun) 
}

# creation of the test function, as a linear combination of k(z_j, .)
M <- 20
set.seed(0)
w <- matrix(runif(M, min = -0.1, max = 0.1), ncol = 1)
z <- ((1:M) - 0.5)/M
f <- function(x){
  kPW(x, z) %*% w
}
t <- seq(0, 1, length = 200)
plot(t, f(t), lty = 1, type = "l")

# norm of the test function
fNorm <- sqrt(t(w) %*% kPW(z, z) %*% w)

# squared test function
fSq <- function(x, eta) {
  f(x)^2
}

# numerical computation of the infinite norm of f
resOpt <- genoud(fSq, nvars = 1, max = TRUE, max.generations = 5, 
                 Domains = matrix(c(0, 1), nrow = 1, ncol = 2),
                 print.level = 0)
plot(t, fSq(t), lty = 1, type = "l")
points(resOpt$par, resOpt$value, pch = 19, col = "blue")
C1 <- sqrt(resOpt$value) # the infinite norm of f

# computation of delta_0 = \int_{R \ [0,1]} f(t)^2 dt 
resInt1 <- integrate(fSq, lower = -Inf, upper = Inf)
resInt2 <- integrate(fSq, lower = 0, upper = 1)
delta0 <- resInt1$value - resInt2$value +
  resInt1$abs.error + resInt2$abs.error

# confidence interval for the RKHS norm, as an upper bound
PWnormCI <- function(X, Y, Yder = NULL, infNorm, delta, alpha){
  n <- length(X)
  res <- mean(Y^2) + infNorm * sqrt(- log(alpha) / (2*n)) + delta
  return(sqrt(res))
}


# computation
n <- 10
seed <- 1 # about X sampling
set.seed(seed)
X <- runif(n)
Y <- f(X)

alpha <- 0.25
fNormUB <- PWnormCI(X = X, Y = Y, alpha = alpha, 
                    infNorm = C1, delta = delta0)
cat("\n RKHS norm of f : ", fNorm)
cat("\n Upper bound at level", alpha, ":", fNormUB)
cat("\n")

# coverage
fileName <- paste("PaleyWienerRegionCoverage_alpha", 
                  floor(100 * alpha), ".pdf", sep = "")
if (pdfPlot) pdf(file = fileName, width = 7, height = 4)
coverPW <- coverage(fNorm = fNorm, f = f, n = n, 
                    RKHSnormCI = PWnormCI, alpha = alpha, nMC = 200, 
                    infNorm = C1, delta = delta0)
if (pdfPlot) dev.off()

cat("\nThe empirical coverage at confidence level ", 
    100*(1-alpha), "% is: ", coverPW$coverage, "\n", sep = "")



# plot
alpha <- 0.1
set.seed(seed)
X <- runif(n)
Y <- f(X)

zalphaSmall <- PWnormCI(X = X, Y = Y, alpha = 0.5, 
                      infNorm = C1, delta = delta0)
fHatRegSmall <- fHatRegion(t, X, Y, kern = kPW, zalpha = zalphaSmall)
zalpha <- PWnormCI(X = X, Y = Y, alpha = alpha, 
                      infNorm = C1, delta = delta0)
fHatReg <- fHatRegion(t, X, Y, kern = kPW, zalpha = zalpha)


fileName <- paste("PaleyWienerRegion_alpha", 
                  floor(100 * alpha), "seed", seed, ".pdf", sep = "")
if (pdfPlot) pdf(file = fileName, width = 7, height = 5)
regionPlot(f = f, fDer = NULL, X = X, ylim = c(-3.5, 3.5), 
           regSmall = fHatRegSmall, reg = fHatReg, 
           alpha)
if (pdfPlot) dev.off()


# coverage of the region
if (0 == 1){

nMC <- 20
alpha <- 0.5

count <- 0
for (i in 1:nMC){
  X <- runif(n)
  Y <- f(X)
  zalpha <- PWnormCI(X = X, Y = Y, alpha = alpha, 
                     infNorm = C1, delta = delta0)
  myFun <- function(x){
    eps <- 1
    fHatReg <- fHatRegion(x, X, Y, kern = kPW, zalpha = zalpha)
    res <- (f(x) - fHatReg$center)^2 - (fHatReg$up - fHatReg$down)^2 / 4
    return(res)
  }
  # curve(myFun, from = 0, to = 1)
  resOpt <- genoud(myFun, nvars = 1, max = TRUE, max.generations = 5,
                   print.level = 0,
                   Domains = matrix(c(0, 1), nrow = 1, ncol = 2))
  # points(resOpt$par, resOpt$value, pch = 19, col = "blue")
  print(resOpt$value)
  tol <- 1e-6
  count <- count + (resOpt$value < tol)
}
cat("\nRegion coverage at level ", alpha, ": ", count/nMC, "\n", sep = "")
}
