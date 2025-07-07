# predictor function (kriging mean) and power function
krigExpr <- function(x, X, Y, kern){
  kxX <- kern(x, X)
  nugget <- 1e-8
  kXXInv <- solve(kern(X, X) + nugget * diag(1, nrow = n, ncol = n))
  krigMean <- as.numeric(kxX %*% kXXInv %*% Y)
  krigCov <- kern(x, x) -  kxX %*% kXXInv %*% t(kxX) 
  krigVar <- pmax(diag(krigCov), 0)
  return(list(fHat = krigMean, powerFun = krigVar))
}


# predictor region
fHatRegion <- function(x, X, Y, kern, zalpha){
  krig <- krigExpr(x, X, Y, kern = kern)
  halfWidth <- zalpha * sqrt(krig$powerFun)
  return(list(center = krig$fHat, 
              up = krig$fHat + halfWidth,
              down = krig$fHat - halfWidth))
  # the region is [down, up]
}

# computation of the empirical coverage

coverage <- function(fNorm, f, fDer = NULL, n, 
                     RKHSnormCI, alpha, nMC = 200, ...){
  # f : test function
  # fDer : if available the derivative of f
  # fNorm : its norm
  # n : sample size
  # alpha : conf level (typically 0.1)
  # nMC : number of Monte Carlo replicates
  # RKHSnormCI: a function that computes a conf int of the RKHS norm
  #             given X (arguments), Y (=f(X)), alpha and  
  #             and possible other arguments in "..."
  # ... : the other arguments required by RKHSnormCI
  fNormUB <- vector(length = nMC)
  for (i in 1:nMC){
    X <- runif(n)
    Y <- f(X)
    Yder <- NULL
    if (!is.null(fDer)) {
      Yder <- fDer(X)
    } 
    fNormUB[i] <- RKHSnormCI(X = X, Y = Y, Yder = Yder, 
                             alpha = alpha, ...)
  }
  coverEmp <- mean(rep(fNorm, nMC) < fNormUB)
  barplot(fNormUB, col = gray(0.9), lwd = 0.1, 
          xlab = paste("Empirical coverage for 1 - alpha = ", 
                       1 - alpha, ": ", coverEmp, sep = ""))
  abline(h = fNorm, col = "blue", lwd = 2)
  legend(# x = nMC - 55, y = 0.4,
         "bottomright", 
         cex = 0.8,
         legend = c("True value", 
                    paste(floor(100*(1-alpha)), "% confidence intervals", sep = "")), 
         pch = c(NA, 15),
         border = c(NA, "darkgray"),
         col = c("blue", "darkgray"), 
         lty = c(1, NA), bg = "white")
  return(list(fNormUB = fNormUB, coverage = coverEmp))
}

regionPlot <- function(f, fDer = NULL, X, regSmall, reg, alpha, ...){
  par(mfrow = c(1, 1))
  t <- seq(0, 1, length = 200)
  #yMax <- max(c(regSmall$up, reg$up))
  #yMin <- min(c(regSmall$down, reg$down))
  plot(t, f(t), lty = "dotted", type = "l", 
       xlab = "", ylab = "", ...) #ylim = c(yMin, yMax))
  tt <- c(t, rev(t))
  yy <- c(fHatRegSmall$up, rev(fHatRegSmall$down))
  yyPlus <- c(fHatReg$up, rev(fHatReg$down))
  polygon(tt, yyPlus, col = gray(0.5), border = NA)
  polygon(tt, yy, col = gray(0.8), border = NA)
  lines(t, f(t), lty = "dashed", lwd = 2)
  lines(t, fHatReg$center, col = "violet", lwd = 2)
  if (!is.null(fDer)) {
    Yder <- fDer(X) 
    h <- 0.04
    arrows(x0 = X, y0 = Y, x1 = X + h, y1 = Y + h * Yder, 
           length = 0.08, col = "blue", lwd = 1)
    arrows(x0 = X, y0 = Y, x1 = X - h, y1 = Y - h * Yder, 
           length = 0.08, col = "blue", lwd = 1)
  }
  points(X, f(X), col = "blue", pch = 19, cex = 1)
  legend("bottomright", 
         c(paste(100 * (1 - alpha), "% confidence region", sep = ""),
           paste(50, "% confidence region", sep = "")),
         fill = c(gray(0.5), gray(0.9)), cex = 0.8, bg = "white")
}

