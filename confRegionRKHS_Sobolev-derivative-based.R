rm(list = ls())

source("confRegionRKHS_Commons.R")
source("confRegionRKHS_Sobolev_DefToyCase.R")

pdfPlot <- FALSE



H1normCI <- function(X, Y, Yder, b, alpha){
  n <- length(X)
  sqNormHat <- mean(Y^2) + mean(Yder^2)
  res <- sqNormHat + b * sqrt(- log(alpha) / (2*n))
  return(sqrt(res))
}

# computation
n <- 10
seed <- 1 # about X sampling
set.seed(seed)
X <- runif(n)
Y <- f(X)
Yder <- fDer(X)

alpha <- 0.1
fNormUB <- H1normCI(X = X, Y = Y, Yder = Yder, 
                    b = b, alpha = alpha)
cat("\n RKHS norm of f : ", fNorm)
cat("\n Upper bound at level", alpha, ":", fNormUB)
cat("\n")


# coverage
alpha <- 0.25
fileName <- paste("H1RegionCoverage_alpha", 
                  floor(100 * alpha), ".pdf", sep = "")
if (pdfPlot) pdf(file = fileName, width = 7, height = 4)
coverPW <- coverage(fNorm = fNorm, f = f, fDer = fDer, n = n, 
                    RKHSnormCI = H1normCI, alpha = alpha, nMC = 200, 
                    b = b)
if (pdfPlot) dev.off()




# plot
alpha <- 0.1
set.seed(seed)
X <- runif(n)
Y <- f(X)
Yder <- fDer(X)

zalphaSmall <- H1normCI(X = X, Y = Y, Yder = Yder, 
                        b = b, alpha = 0.5)
fHatRegSmall <- fHatRegion(t, X, Y, kern = kH1, zalpha = zalphaSmall)
zalpha <- H1normCI(X = X, Y = Y, Yder = Yder, 
                        b = b, alpha = alpha)
fHatReg <- fHatRegion(t, X, Y, kern = kH1, zalpha = zalpha)


fileName <- paste("H1Region_alpha", 
                  floor(100 * alpha), "seed", seed, ".pdf", sep = "")
if (pdfPlot) pdf(file = fileName, width = 7, height = 5)
regionPlot(f = f, fDer = fDer, X = X, 
           regSmall = fHatRegSmall, reg = fHatReg, 
           alpha = alpha, ylim = c(-1.8, 0.3))
if (pdfPlot) dev.off()

