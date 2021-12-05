# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------
# Code from final MS.


# 3.1 Introduction
# ================

# ~~~~ Create Figure 3.1 ~~~~
T <- 9                              # Number of years
N <- matrix(NA, nrow=T+1, ncol=3)   # Matrix to hold size of 3 populations
lambda <- c(0.8, 1, 1.2)            # Pop. growth rate of 3 pops
N[1,] <- 10                         # Initial population size is 10
for (t in 1:T){                     # Loop over years
  for (k in 1:3){                   # Loop over 3 populations
    N[t+1,k] <- lambda[k] * N[t,k]  # Project pop. one time step forward
  } #k
} #t

op <- par(cex=1.4)
plot(N[,1], type="l", ylim=range(N), ylab="Population size", xlab="Year", axes=FALSE, lwd=2.5)
lines(N[,2], lwd=2.5)
lines(N[,3], lwd=2.5)
axis(1, at=1:(T+1), labels=1:(T+1))
axis(2, las=1)
text(x=T, y=N[5,1], expression(paste(lambda, "=0.8")))
text(x=T, y=N[10,2] + 3, expression(paste(lambda, "=1.0")))
text(x=T, y=N[8,3], expression(paste(lambda, "=1.2")))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
