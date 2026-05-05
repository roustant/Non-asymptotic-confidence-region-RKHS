library(sensitivity)
library(evd)

## ---
## functions for Hone(mu), 
## with mu = Gumbel(1013, 558) truncated on [500, 3000]
## ---

pdfPlot <- FALSE

# original specifications of the truncated Gumbel
loc <- 1013
scale <- 558
Qmin <- 500; Qmax <- 3000
# modification for the support to be [0, 1]
# pushforward by x -> (x - xmin)/(xmax - xmin)
scale01 <- scale / (Qmax - Qmin)
loc01 <- (loc - Qmin) / (Qmax - Qmin)

simFlood <- function(n){
  # Z <- rgumbel.trunc(n, loc = 1013, scale = 558, 
  #                    min = Qmin, max = Qmax)
  # X <- (Z - Qmin)/(Qmax - Qmin)
  X <- rgumbel.trunc(n, loc = loc01, scale = scale01, min = 0, max = 1)
  return(X)
}


# eigenfunctions, eigenvalues, kernel
# ---

out <- PoincareOptimal(
  distr = list("gumbel", loc = loc01, scale = scale01), 
  min = 0, max = 1, 
  method = "integral", only.values = FALSE, 
  der = TRUE, plot = FALSE, n = 500)


# eigenvalues
lambda <- function(N){
  res <- 0
  res <- c(res, out$values[2:(N+1)])
  return(res)
}

# eigenfunctions
eigenFunList <- function(N){
  # returns eigenfunctions phi0, ..., phiN
  # as a list of functions
  res <- list(function(x) 1)
  for (i in 2:(N+1)){
    res <- c(res, approxfun(x = out$knots, 
                            y = out$vectors[, i]))  
  }
  return(res)
}

par(mfrow = c(1, 1))
t <- seq(0, 1, length = 200)

plot(t, dgumbel.trunc(t, loc = loc01, scale = scale01,
                      min = 0, max = 1),
     type = "l", ylab = "Density", xlab = "")
hist(simFlood(1e5), freq = FALSE, xlab = "Q", main = "",
     add = TRUE)


if (pdfPlot) pdf(file = "QlawBasis.pdf", width = 8, height = 5)
N <- 5
funValues <- eigenFunList(N)
plot(t, rep(1, length(t)), type = "l", 
     ylim = c(min(funValues[[N+1]](t)), max(funValues[[N+1]](t))),
     xlab = "", ylab = "")
yMin <- min(funValues[[N+1]](t))
ypdf <- dgumbel.trunc(t, loc = loc01, scale = scale01, 
                      min = 0, max = 1)
polygon(x = c(t, 0, 0),
        y = c(ypdf, ypdf[length(t)], ypdf[1]) + yMin,
        col = "lightgray", border = FALSE)
for (i in 1:N){
  lines(t, funValues[[i+1]](t), type = "l", col = i+1, lty = i+1)
}

legend("topleft", legend = paste("degree", 0:N),
         col = 1:(N+1), lty = 1:(N+1), cex = 0.8, bg = "white")
if (pdfPlot) dev.off()

basisFun <- function(x, N){
  # returns the values at x of the first N+1 eigenfunctions
  # (phi0(x), phi1(x), ..., phiN(x))
  
  evalFun <- function(f, x) f(x)
  res <- lapply(eigenFunList(N), evalFun, x)
  return(unlist(res))
}

# kernel, constructed with the Mercer expansion
kH1MuFun <- function(x, y){
  # accepts vectors for x, y
  M <- 50
  S <- 0
  beta <- 1/(1 + lambda(M))
   for (i in 1:M){
     S <- S + beta[i] * eigenFunList(M)[[i]](x) * eigenFunList(M)[[i]](y)  
   }
  return(S)
} 

kH1Mu <- function(x, y){
  outer(x, y, kH1MuFun)
}

# x <- runif(1); y <- runif(1)
# c(kH1MuFun(x, y), kH1Fun(x, y))


## ---
## The function of interest: flood model
## ---
# Input names
floodInputNames <- c("Q", "Ks", "Zv", "Zm", "Hd", "Cb", "L", "B")
floodOutputNames <- c("Subverse", "Cost")

# Flood model
flood <- function(X){ 
  mat <- as.matrix(X, ncol = 8)
  output <- rep(NA, nrow(mat))
  for (i in 1:nrow(mat)){
    H <- (mat[i, 1] / (mat[i, 2] * mat[i, 8] * sqrt((mat[i, 4] - mat[i, 3]) / mat[i, 7])))^0.6
    S <- mat[i, 3] + H - mat[i, 5] - mat[i, 6] 
    output[i] <- S
  }
  return(output)
}

# Flood model derivatives (by finite-differences)
floodDer <- function(X, i, eps = 1e-7){
  der <- X
  X1 <- X
  X1[, i] <- X[, i] + eps
  der <- (flood(X1) - flood(X)) / (eps)
  return(der)
}


# Function for flood model inputs sampling
floodSample <- function(size){
  X <- matrix(NA, size, 8)
  X[, 1] <- rgumbel.trunc(size, loc = 1013, scale = 558, min = Qmin, max = Qmax)
  X[, 2] <- rnorm.trunc(size, mean = 30, sd = 8, min = 15)
  X[, 3] <- rtriangle(size, a = 49, b = 51)
  X[, 4] <- rtriangle(size, a = 54, b = 56)
  X[, 5] <- runif(size, min = 7, max = 9)
  X[, 6] <- rtriangle(size, a = 55, b = 56)
  X[, 7] <- rtriangle(size, a = 4990, b = 5010)
  X[, 8] <- rtriangle(size, a = 295, b = 305)
  return(X)
}


# nominal value for all variables
xRef <- c(NA, # truncated Gumbel : not required
        qnorm.trunc(p = 0.5, mean = 30, sd = 8, min = 15), # median normal
        50, 55, 8, 55.5, 5000, 300) # median other distributions

fFun <- function(x){
  # x : real number between 0 and 1, corr. to normalized Q
  z <- c(Qmin + (Qmax - Qmin) * x, xRef[-1])
  z <- matrix(z, nrow = 1)
  res <- flood(z) + 11
  return(res)
}

fFunDer <- function(x){
  z <- c(Qmin + (Qmax - Qmin) * x, xRef[-1])
  z <- matrix(z, nrow = 1)
  # derivative : f'(x) = dflood/dz1(z) * (Qmax - Qmin) 
  res <- floodDer(z, i = 1) * (Qmax - Qmin)
  return(res)
}

# modify fFun and fFunDer such that it accepts vectors as arguments
f <- function(x) sapply(x, fFun) 
fDer <- function(x) sapply(x, fFunDer) 

# graphical representation of the test function and its derivative
t <- seq(0, 1, length = 200)
par(mfrow = c(1, 2))
plot(t, f(t), lty = 1, type = "l")
plot(t, fDer(t), lty = 1, type = "l")

fSq <- function(x) {
  f(x)^2 + fDer(x)^2
}

X <- simFlood(1e3)
fNorm <- sqrt(mean(sapply(X, fSq)))



