library(sensitivity)

library(evd)
n <- 500
out <- PoincareOptimal(distr = list("gumbel", 0, 1), min = -0.92, max = 3.56, 
                       method = "integral", 
                       only.values = FALSE, 
                       der = TRUE,
                       plot = TRUE,
                       n = n)
print(out$opt)
str(out)

p <- 5/2
eigenFunNormInf <- apply(abs(out$vectors[, 2:n]), 2, max)
eigenDerNormInf <- apply(abs(out$der[, 2:n]), 2, max)
partialS <- cumsum((out$values[2:n])^{-p} * eigenFunNormInf)
partialDerS <- cumsum((out$values[2:n])^{-p} * eigenDerNormInf)
par(mfrow = c(1, 2))
plot(partialS, type = "l")
plot(partialDerS, type = "l")
