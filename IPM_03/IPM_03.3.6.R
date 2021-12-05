# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

# 3.3 Classical analysis of a matrix population model
# ===================================================

# 3.3.6 Matrix population models with density-dependence and
#       demographic stochasticity
# -----------------------------------------------------------

# Define means of the demographic rates and strength of density dependence
mean.sj <- 0.3
mean.sa <- 0.55
f1.int <- 2.3         # Productivity of 1y females when population size is 0
f1.beta <- -0.02      # Strength of density dependence on productivity of 1y
                      # This is simply the slope of the regression of f1 on N
fa.int <- 2.3         # Productivity of adult females when population size is 0
fa.beta <- -0.01      # Strength of density dependence on productivity of adults

# Define the number of years with predictions
T <- 200
nsim <- 1000

# Define population matrix and initial stage-specific population sizes
N <- array(NA, dim=c(2, T+1, nsim))
N[,1,] <- c(10, 10)
r <- f1 <- fa <- matrix(NA, nrow=T, ncol=nsim)

# Project population
for (s in 1:nsim){
  if(s %% 100 == 0) {cat(paste("*** Simrep", s, "***\n")) } # Counter
  for (t in 1:T){
    f1[t,s] <- f1.int + f1.beta * sum(N[,t,s])
    fa[t,s] <- fa.int + fa.beta * sum(N[,t,s])
    N[1,t+1,s] <- rpois(1, mean.sj * (f1[t,s] * N[1,t,s] + fa[t,s] * N[2,t,s]))
    N[2,t+1,s] <- rbinom(1, sum(N[,t,s]), mean.sa)
  } #t
} #s

# ~~~~ code for Fig. 3.16 ~~~~
op <- par(mar=c(2, 4, 2, 1), las=1, cex=1.1, "mfrow")
layout(matrix(1:3, 3, 1, byrow=TRUE), widths=2, heights=c(1, 1, 1.1), TRUE)
plot(colSums(N[,,1]), type="l", col="lightgrey", ylim=c(0, 120),
    ylab="Population size", xlab="", axes=FALSE)
axis(1); axis(2)
for (s in 2:nsim){
  lines(colSums(N[,,s]), col="lightgrey")
}
Nall <- apply(N, c(2,3), sum)
lines(apply(Nall, 1, mean), col="black", lwd=1.5)
mtext("A", at=1, cex=1.5)

f1[f1==f1.int] <- NA
plot(f1[,1], type="l", col="lightgrey", ylim=c(0.3, 2.3),
    ylab=expression(italic(f)[1]), xlab=NA, axes=FALSE)
axis(1); axis(2)
for (s in 2:nsim){
  lines(f1[,s], col="lightgrey")
}
lines(apply(f1, 1, mean, na.rm=TRUE), col="black", lwd=1.5)
mtext("B", at=1, cex=1.5)

par(mar=c(4,4,2,1))
fa[fa==fa.int] <- NA
plot(fa[,1], type="l", col="lightgrey", ylim=c(0.3, 2.3),
    ylab=expression(italic(f)[italic(a)]), xlab="Time", axes=FALSE)
axis(1); axis(2)
for (s in 2:nsim){
  lines(fa[,s], col="lightgrey")
}
lines(apply(fa, 1, mean, na.rm=TRUE), col="black", lwd=1.5)
mtext("C", at=1, cex=1.5)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mean(N[1,T+1,] + N[2,T+1,])
# [1] 52.687

mean(f1[T,], na.rm=TRUE)
# [1] 1.237938
mean(fa[T,], na.rm=TRUE)
# [1] 1.768969

mean(N[1,T+1,] + N[2,T+1,] == 0)
# [1] 0.001
