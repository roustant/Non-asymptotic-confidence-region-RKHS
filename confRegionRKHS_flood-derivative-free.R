rm(list = ls())


source("confRegionRKHS_flood_DefToyCase.R")


NormSqHat <- function(N, X, Y){
  # Y = f(X)
  n <- length(X)
  phiMat <- sapply(X, basisFun, N = N) # matrix of size (N+1) * n
  phiMat <- matrix(phiMat, ncol = n)
  beta <- 1 + lambda(N)
  res <- 0
  for (l in 1:(N+1)){
    hl <- phiMat[l, ] * Y # vecteur contenant f(X)phi_l(X) 
    Tl <- sum(hl)^2 - sum(hl^2) # la U-stat (au facteur de norm. pres)
    res <- res + beta[l] * Tl   
  }
  return(res * 1 / (n * (n-1)) )
} 


# Infinite norm of eigenvalues and their derivatives

p <- 5/2
nTrunc <- 30
eigenFunNormInf <- apply(abs(out$vectors[, 1:nTrunc]), 2, max)
eigenDerNormInf <- apply(abs(out$der[, 1:nTrunc]), 2, max)
partialS <- cumsum((out$values[2:nTrunc])^{-p} * eigenFunNormInf[2:nTrunc])
partialDerS <- cumsum((out$values[2:nTrunc])^{-p} * eigenDerNormInf[2:nTrunc])
par(mfrow = c(1, 2))
plot(partialS, type = "l")
plot(partialDerS, type = "l")


phiInf <- function(N){
  # N should be less than 500
  res <- apply(abs(out$vectors[, 1:(N+1)]), 2, max)
  res[1] <- 1 # we use the value for phi0
  return(res)
}

phiDerInf <- function(N){
  # N should be less than 500
  res <- apply(abs(out$der[, 1:(N+1)]), 2, max)
  res[1] <- 0 # we use the true value for phi0'
  return(res)
}

phiInf(10)
phiDerInf(10)


n <- 1000
q <- 4
p <- 5/2

fInf <- max(abs(f(0)), abs(f(1))) # sup norm of f here
partialS <- sum(lambda(nTrunc)[2:nTrunc]^{-p} * 
                  phiInf(nTrunc)[2:nTrunc])
A <- 2 * fInf / partialS

partialDerS <- sum(lambda(nTrunc)[2:nTrunc]^{-p} * 
                     phiDerInf(nTrunc)[2:nTrunc])
fDerInf <- A * partialDerS 


# bound of | c_l |
cellUB <- function(N, fInf = NA, fDerInf = NA){
  UB <- A / lambda(N)^p # by assumption
  UB[1] <- fInf # upper bound of c0
  UB
}

hCentUB <- function(N, fInf, fDerInf){
  # upper bound of the infinite norm of h_l - c_l
  # with bounds of h_l and c_l
  bound1 <- fInf * phiInf(N) + cellUB(N, fInf, fDerInf)
  # with bounds of h'_l
  bound2 <- fDerInf * phiInf(N) + fInf * phiDerInf(N)
  # minimum of the two
  pmin(bound1, bound2)
}


### gammaN (upper bound of the support of the i.i.d r.v. in Hoeffding inequality)
sUB <- function(n, N, fInf, fDerInf){
  beta <- 1 + lambda(N)
  aux1 <- cellUB(N, fInf, fDerInf) * hCentUB(N, fInf, fDerInf)
  term1 <- 4 * sum(beta * aux1) 
#  aux2 <- fInf^2 + hCentUB(N, fInf, fDerInf)^2
  aux2 <- hCentUB(N, fInf, fDerInf)^2
  term2 <- 2 * sum(beta * aux2) / (n-1)
  return(term1 + term2)
}


### t
margin <- function(n, N, alpha = 0.05, fInf, fDerInf){
  beta <- 1 + lambda(N)
  deltaUB <- 2 * fInf^2 / (n-1) * sum(beta)
  t <- sUB(n, N, fInf, fDerInf) * sqrt(- log( alpha ) / (2 * n ))
  return(list(deltaUB = deltaUB, t = t))
}


### zeta: Rn(f) <= zeta(N) * || f ||_{H1}^2
# should decrease to 0 when N tends to infty

xi <- function(N){
  (N+1)^(-q) * 0.5 * 2^{q}  # for N = 1, we obtain 0.5
}

### upperBound
UB <- function(fNormSqHat, n, N, fInf, fDerInf, alpha = 0.1){
  statVal <- fNormSqHat
  aux <- margin(n = n, N = N, alpha = alpha, 
                fInf = fInf, fDerInf = fDerInf)
  xiBarVal <- 1 - xi(N) 
  bound <- sqrt(max(statVal + aux$t + aux$deltaUB, 0) / xiBarVal)
  return(c(bound = bound, stat = statVal, 
           t = aux$t, deltaUB = aux$deltaUB, bias = xiBarVal))
}

## computation of the optimal bounds in function of N

set.seed(10)
cat("n = ", n, "; q = ", q, "\n")
X <- simFlood(n)
Y <- f(X)

Nmax <- 10
UBres <- matrix(NA, nrow = Nmax, ncol = 5)
for (N in 1:Nmax){
  #print(n)
  UBcurrent <- UB(fNormSqHat = NormSqHat(N = N, X = X, Y = Y), 
                  n = n, N = N, fInf, fDerInf)
  UBres[N, ] <- UBcurrent
}
colnames(UBres) <- names(UBcurrent)
print(UBres)
cat("\n RKHS norm of f : ", fNorm)
cat("\n Upper bound with Inf norm : ", sqrt(fInf^2 +fDerInf^2))
cat("\n")


    
