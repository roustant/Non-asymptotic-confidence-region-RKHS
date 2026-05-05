rm(list = ls())

source("confRegionRKHS_Commons.R")
source("confRegionRKHS_flood_DefToyCase.R")


pdfPlot <- TRUE

method <- "hoeffding"
method <- "subGaussian"



# computation
n <- 2; seed <- 6
n <- 10; seed <- 1 # about X sampling
set.seed(seed)

X <- simFlood(n)
Y <- f(X)
Yder <- fDer(X)

if (method == "hoeffding"){
  # evaluation of b = sup de f^2 + f'^2
  resOpt <- optimize(fSq, interval = c(0, 1), maximum = TRUE)
  par(mfrow = c(1, 1))
  plot(t, fSq(t), lty = 1, type = "l")
  points(resOpt$maximum, resOpt$objective, pch = 19, col = "blue")
  b <- resOpt$objective
  lipSlope <- fDer(0)
} else if (method == "subGaussian") {
  b <- subGaussianParam(simFlood(1e5), fSq, min = -5, max = 5)
  lipSlope <- NULL
}

alpha <- 0.1
fNormUB <- H1normCI(X = X, Y = Y, Yder = Yder, 
                    b = b, alpha = alpha, method = method)
cat("\n RKHS norm of f : ", fNorm)
cat("\n Upper bound at level", alpha, ":", fNormUB)
cat("\n")


# coverage
alpha <- 0.25
fileName <- paste("H1flood_RegionCoverage_alpha", 
                  floor(100 * alpha), ".pdf", sep = "")
if (pdfPlot) pdf(file = fileName, width = 7, height = 4)
coverPW <- coverage(fNorm = fNorm, f = f, fDer = fDer, n = n, 
                    RKHSnormCI = H1normCI, alpha = alpha, method = method,
                    nMC = 200, simFun = simFlood, b = b)
if (pdfPlot) dev.off()




# plot
alpha <- 0.1

zalphaSmall <- H1normCI(X = X, Y = Y, Yder = Yder, 
                        b = b, alpha = 0.5, method = method)
fHatRegSmall <- fHatRegion(t, X, Y, kern = kH1Mu, 
                           zalpha = zalphaSmall)
zalpha <- H1normCI(X = X, Y = Y, Yder = Yder, 
                        b = b, alpha = alpha, method = method)
fHatReg <- fHatRegion(t, X, Y, kern = kH1Mu, zalpha = zalpha)


fileName <- paste("H1flood_Region_alpha", 
                  floor(100 * alpha), "seed", seed, ".pdf", sep = "")
if (pdfPlot) pdf(file = fileName, width = 7, height = 5)
regionPlot(f = f, fNorm = fNorm, 
           fDer = fDer, lipSlope = lipSlope,
           X = X, 
           regSmall = fHatRegSmall, reg = fHatReg,
           alpha = alpha)
if (pdfPlot) dev.off()

