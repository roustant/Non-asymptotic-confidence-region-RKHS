rm(list = ls())

source("confRegionRKHS_Commons.R")
source("confRegionRKHS_Sobolev_DefToyCase.R")

pdfPlot <- FALSE
method <- "hoeffding"
method <- "subGaussian"


# computation
n <- 2; seed <- 6
n <- 10; seed <- 1 # about X sampling
set.seed(seed)
X <- runif(n)
Y <- f(X)
Yder <- fDer(X)

if (method == "hoeffding"){
  # evaluation of b = sup de f^2 + f'^2
  library(rgenoud)
  resOpt <- genoud(fSq, nvars = 1, max = TRUE, max.generations = 5, 
                   Domains = matrix(c(0, 1), nrow = 1, ncol = 2))
  
  par(mfrow = c(1, 1))
  plot(t, fSq(t), lty = 1, type = "l")
  points(resOpt$par, resOpt$value, pch = 19, col = "blue")
  b <- resOpt$value
  lipSlope <- sqrt(b)
} else if (method == "subGaussian") {
  b <- subGaussianParam(sample = runif(1e5), fSq, min = -5, max = 5)
  lipSlope <- NULL
}


alpha <- 0.25
fNormUB <- H1normCI(X = X, Y = Y, Yder = Yder, 
                    b = b, alpha = alpha, method = method)
cat("\n RKHS norm of f : ", fNorm)
cat("\n Upper bound at level", alpha, ":", fNormUB)
cat("\n")


# coverage
alpha <- 0.25 #25
fileName <- paste("H1RegionCoverage_alpha", 
                  floor(100 * alpha), method, "_n", n, ".pdf", sep = "")
if (pdfPlot) pdf(file = fileName, width = 4, height = 4)
coverPW <- coverage(fNorm = fNorm, f = f, fDer = fDer, n = n, 
                    RKHSnormCI = H1normCI, alpha = alpha, nMC = 200, 
                    b = b, method = method)
if (pdfPlot) dev.off()




# plot
alpha <- 0.1
set.seed(seed)
X <- runif(n)
Y <- f(X)
Yder <- fDer(X)

t <- seq(0, 1, length = 501)

zalphaSmall <- H1normCI(X = X, Y = Y, Yder = Yder, 
                        b = b, alpha = 0.5, method = method)
fHatRegSmall <- fHatRegion(t, X, Y, kern = kH1, zalpha = zalphaSmall)
zalpha <- H1normCI(X = X, Y = Y, Yder = Yder, 
                   b = b, alpha = alpha, method = method)
fHatReg <- fHatRegion(t, X, Y, kern = kH1, zalpha = zalpha)


fileName <- paste("H1Region_alpha", 
                  floor(100 * alpha), "seed", seed, method, 
                  "_n", n, ".pdf", sep = "")
if (pdfPlot) pdf(file = fileName, width = 7, height = 5)
regionPlot(f = f, X = X, fNorm = fNorm, 
           fDer = fDer, lipSlope = lipSlope,
           regSmall = fHatRegSmall, 
           reg = fHatReg, 
           alpha = alpha
           )
if (pdfPlot) dev.off()

