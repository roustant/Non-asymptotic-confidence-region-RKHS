H1normCI <- function(X, Y, Yder, b, alpha, method){
  n <- length(X)
  sqNormHat <- mean(Y^2) + mean(Yder^2)
  if (method == "hoeffding"){
    aux <- 1/2
  } else if (method == "subGaussian"){
    aux <- 2
  }
  res <- sqNormHat + b * sqrt(- aux * log(alpha) / n)
  return(sqrt(res))
}

# evaluation of the parameter s of f(X)^2 + f'(X)^2
# when assumed to be sub-gaussian, i.e. 
# E(e^{l (X - mu)}) <= e^{s^2 l^2/2}
# we replace the expectations by empirical means
# and optimize with respect to l

subGaussianParam <- function(sample, f, plot = TRUE, min = -1, max = 1){
  # maximizes g: l -> ln(E(e^{l (f(X) - E(f(X)))})) / (l^2/2) on support
  # min, max : domain on which g is maximized
  # f : a function that accepts vectors as entries
  # sample : a sample drawn from mu
  # plot : if TRUE plot the function l -> ln(E(e^{l (X- mu)}))/(l^2/2)
  Z <- f(sample)
  Zmean <- mean(Z)
  g <- function(l){
    res <- mean(exp(l * (Z - Zmean)))
    res <- log(res) * 2 / (l^2)
    return(res)
  }
  gFun <- function(x){
    sapply(x, g)
  }
  lval <- seq(from = min, to = max, length = 200)
  
  resOpt <- optimize(gFun, interval = c(min, max),
                     maximum = TRUE)
  
  if (plot){
    par(mfrow = c(1, 1))
    plot(lval, gFun(lval), type = "l", col = "blue",
         xlab = expression(lambda), ylab = expression(g(lambda)))
    points(resOpt$maximum, resOpt$objective, pch = 19, col = "blue")
  }
  res <- resOpt$objective
  return(sqrt(res))
}


#  lipschitz upper/lower bound
lipschBound <- function(x, y, slope, xmin = 0, xmax = 1){
  # Computes a piecewise affine function, which bounds
  # an interpolator f such that |f'| <= |slope|
  # If slope > 0 (resp <0) gives an upper (resp. lower) bound
  # ---------
  # Arguments
  # ---------
  # x, xmin, xmax : vector of input values, 
  #                 assumed to lie in [xmin, xmax]
  # y : vector of output values
  # slope : |slope| is an upper bound of the sup norm of f' 
  #         where f interpolates y at x
  # -----
  # Value
  # -----
  # Returns a list defining the piecewise affine function
  #   $x : knots 
  #   $y : values at knots
  
  # order x, y according to x, in ascending order
  ii <- order(x, y)
  x <- x[ii]
  y <- y[ii]

  # verification that the slope is larger than the maximal slope
  # between consecutive data points
  maxSlopeData <- max(abs(diff(y)/diff(x)))
  if (maxSlopeData > abs(slope)){
    stop(paste("The absolute value of slope", slope, "must be greater than the 
         maximal slope between consecutive points", maxSlopeData))
  }
  n <- length(x)
  xx <- xmin 
  # compute the value at xmin of the line 
  # containing (x[1], y[1]) with slope equal to -slope
  yy <- y[1] - slope * (xmin - x[1])
  for (i in 1:(n-1)){
    xx <- c(xx, x[i])
    yy <- c(yy, y[i])
    # compute the intersection between the lines
    #  left : contains (x[i], y[i]), slope equal to 'slope'
    #  right: contains (x[i+1], y[i+1]), slope eq. to -'slope'
    xxInter <- (y[i+1] - y[i]) / (2 * slope) + (x[i+1] + x[i]) / 2
    yyInter <- (y[i+1] + y[i]) / 2 + slope * (x[i+1] - x[i]) / 2
    xx <- c(xx, xxInter)
    yy <- c(yy, yyInter)
  }
  xx <- c(xx, x[n])
  yy <- c(yy, y[n])
  xx <- c(xx, xmax)
  # compute the value at xmax of the line 
  # containing (x[n], y[n]) with slope equal to 'slope'
  yy <- c(yy, y[n] + slope * (xmax - x[n]))
  return(list(x = xx, y = yy))
}

lipschPolygon <- function(x, y, slope, ...){
  z1 <- lipschBound(x, y, slope = slope)
  z2 <- lipschBound(x, y, slope = - slope)
  xx <- c(z1$x, rev(z2$x))
  yy <- c(z1$y, rev(z2$y))
  polygon(xx, yy, ...)
}

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
              powerFun = krig$powerFun, 
              up = krig$fHat + halfWidth,
              down = krig$fHat - halfWidth))
  # the region is [down, up]
}

# computation of the empirical coverage

coverage <- function(fNorm, f, fDer = NULL, n, 
                     RKHSnormCI, alpha, nMC = 200, 
                     simFun = runif, ...){
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
    X <- simFun(n)
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

regionPlot <- function(f, fNorm = NULL, 
                       fDer = NULL, lipSlope = NULL,
                       X, regSmall, reg, alpha, ...){
  par(mfrow = c(1, 1))
  t <- seq(0, 1, length = length(reg$center))
  yMax <- max(c(regSmall$up, reg$up))
  yMin <- min(c(regSmall$down, reg$down))
  ylim <- 1.2 * c(yMin, yMax)    
  plot(t, f(t), lty = "dotted", type = "l", 
       xlab = "", ylab = "", ylim = ylim, ...) #
  tt <- c(t, rev(t))
  yy <- c(fHatRegSmall$up, rev(fHatRegSmall$down))
  yyPlus <- c(fHatReg$up, rev(fHatReg$down))
  polygon(tt, yyPlus, col = gray(0.5), border = NA)
  polygon(tt, yy, col = gray(0.9), border = NA)
  lines(t, f(t), lty = "dashed", lwd = 2)
  lines(t, fHatReg$center, col = "violet", lwd = 2)
  
  leg.txt <- c(paste(100 * (1 - alpha), "% confidence region", sep = ""),
               paste(50, "% confidence region", sep = ""))
  leg.pch <- c(22, 22)
  leg.col <- c("black", "black")
  leg.pt.bg <- c(gray(0.5), gray(0.9))
  leg.lty <- c(NA, NA)
  
  if (!is.null(fDer)) {
    Yder <- fDer(X) 
    h <- 0.04
    arrows(x0 = X, y0 = Y, x1 = X + h, y1 = Y + h * Yder, 
           length = 0.08, col = "blue", lwd = 1)
    arrows(x0 = X, y0 = Y, x1 = X - h, y1 = Y - h * Yder, 
           length = 0.08, col = "blue", lwd = 1)
    if (!is.null(lipSlope)) {
      lipschPolygon(X, Y, slope = lipSlope, 
                    border = "orange", lty = "dashed")
      leg.txt <- c(leg.txt, "Lipschitz region")
      leg.pch <- c(leg.pch, NA)
      leg.col <- c(leg.col, "orange")
      leg.pt.bg <- c(leg.pt.bg, "NA")
      leg.lty <- c(leg.lty, "dashed")
    }
  }
  points(X, f(X), col = "blue", pch = 19, cex = 1)
  if (!is.null(fNorm)){
    lines(t, fHatReg$center + fNorm * sqrt(fHatReg$powerFun), 
          lty = "dotted")
    lines(t, fHatReg$center - fNorm * sqrt(fHatReg$powerFun), 
          lty = "dotted")
    leg.txt <- c(leg.txt, "Oracle region")
    leg.pch <- c(leg.pch, NA)
    leg.col <- c(leg.col, "black")
    leg.pt.bg <- c(leg.pt.bg, "NA")
    leg.lty <- c(leg.lty, "dotted")
  }
  
  legend("bottomright", legend = leg.txt,
         pch = leg.pch, col = leg.col, 
         pt.bg = leg.pt.bg, lty = leg.lty, lwd = 1,
         bg = "white") #, cex = 0.8) 
  
}

