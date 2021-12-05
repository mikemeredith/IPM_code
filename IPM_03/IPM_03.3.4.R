# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

# Run time approx. 3 mins

# 3.3 Classical analysis of a matrix population model
# ===================================================

# 3.3.4 Analysis of a matrix population model with demographic stochasticity
# --------------------------------------------------------------------------

rbinom(n=1, size=10, prob=0.55)
# [1] 4

rbinom(n=20, size=10, prob=0.55)
# [1] 4 6 7 6 7 4 4 7 9 3 5 5 6 6 1 6 2 6 1 6

size <- c(10, 25, 50, 75, 100, 250, 500, 750, 1000, 5000, 10000)
s <- matrix(NA, nrow=1000, ncol=length(size))
for (i in 1:length(size)){
  s[,i] <- rbinom(n=1000, size=size[i], prob=0.55) / size[i]
}
boxplot(s, ylab="Realized proportion of surviving birds",
    xlab="Population size (number of females)", names=size, outline=FALSE,
    col="dodgerblue", las=1, frame=FALSE)
abline(h=0.55, col="red", lwd=1)

# Define mean of the demographic parameters
mean.sj <- 0.3
mean.sa <- 0.55
mean.f1 <- 1.3
mean.fa <- 1.8

# Define the number of years over which we let the population evolve
T <- 200

# Define population matrix and initial stage-specific population sizes
N <- matrix(NA, nrow=2, ncol=T+1)
N[,1] <- c(10, 10)

# Project population
r <- numeric(T)
for (t in 1:T){
  N[1,t+1] <- rpois(1, mean.sj * mean.f1 * N[1,t] + mean.sj * mean.fa * N[2,t])
  N[2,t+1] <- rbinom(1, sum(N[,t]), mean.sa)
  if (sum(N[,t+1]) == 0) break                                # Stop calculation if pop. extinct
  r[t] <- log(sum(N[,t+1])) - log(sum(N[,t]))
}

mean(r)                                                       # Compute stochastic growth rate
# [1] 0.02552

# Plot graphs with estimated population growth rates and size
op <- par("mfrow", "cex")
layout(matrix(1:2, 1, 2, byrow=TRUE), widths=c(1, 1), heights=1, TRUE)
plot(r, type="l", lwd=1.5, ylab="Annual population growth rate", xlab="Time", axes=FALSE)
axis(1)
axis(2, las=1)
abline(h=mean(r))
mtext("A", at=1, cex=1.5)
plot(N[1,], type="l", lwd=1.5, ylab="Population size", xlab="Time", axes=FALSE,
    ylim=range(N, na.rm=TRUE))
axis(1)
axis(2, las=1)
mtext("B", at=3, cex=1.5)
lines(N[2,], lwd=1.5, col="red")
legend("topleft", lwd=c(1.5, 1.5), col=c("black", "red"), legend=c("1-year-old", ">1-year-old"), bty="n")
par(op)

# Define mean of the demographic parameters
mean.sj <- 0.3
mean.sa <- 0.55
mean.f1 <- 1.3
mean.fa <- 1.8

# Define the number of years with predictions and the Monte Carlo setting
T <- 200
nsim <- 100000

# Define population matrix and initial stage-specific population sizes
N <- array(NA, dim=c(2, T+1, nsim))
N[,1,] <- c(10, 10)
r <- matrix(NA, nrow=T, ncol=nsim)
alive <- matrix(NA, nrow=T, ncol=nsim)
mean.r <- numeric(nsim)

# Project population
for (s in 1:nsim){
  if(s %% 1000 == 0) {cat(paste("*** Simrep", s, "***\n")) } # Counter
  for (t in 1:T){
    N[1,t+1,s] <- rpois(1, mean.sj * mean.f1 * N[1,t,s] + mean.sj * mean.fa * N[2,t,s])
    N[2,t+1,s] <- rbinom(1, sum(N[,t,s]), mean.sa)
    if (sum(N[,t+1,s]) == 0) break
    r[t,s] <- log(sum(N[,t+1,s])) - log(sum(N[,t,s]))
    alive[t,s] <- t
  } #t
  mean.r[s] <- mean(r[min(alive[,s], na.rm=TRUE):max(alive[,s], na.rm=TRUE),s])
} #s

# Mean and SD of the population growth rate
mean(mean.r)
# [1] -0.009097617
sd(mean.r)
# [1] 0.05689681

# Mean and SD of the population growth rate for populations that did not go extinct
not.extinct <- which(!is.na(alive[T,]))
mean(mean.r[not.extinct])
# [1] 0.01966304
sd(mean.r[not.extinct])
# [1] 0.006095664

# Extinction probability (after T years)
sum(is.na(alive[T,])) / nsim
# [1] 0.28917

# Plot graph with population growth rates (Fig. 3.14 and 3.15)
op <- par(mar=c(4, 4, 2, 1), las=1, cex=1.1, "mfrow")
layout(matrix(c(1, 1, 2, 3), 2, 2, byrow=TRUE), widths=c(1, 1), heights=c(1, 1), TRUE)
# par(mar=c(4, 4, 2, 1), las=1, cex=1.1)
plot(r[,1], type="l", lwd=0.5, ylab="Annual population growth rate", xlab="Time",
    ylim=range(r[which(!is.na(alive))]), col="lightgrey", axes=FALSE)
axis(1); axis(2)
for (s in 2:nsim){
  lines(r[!is.na(alive[,s]),s], lwd=0.5, col="lightgrey")
}
lines(apply(r, 1, mean, na.rm=TRUE), lwd=1.5)
mtext("A", at=0, cex=1.5)
a <- hist(mean.r, nclass=100, col="dodgerblue", main="", xlab="Population growth rate",
    xlim=c(-0.5, 0.1), axes=FALSE)
axis(1)
axis(2, at=c(0, 10000, 20000, 30000, 40000), labels=c(0, 10, 20, 30, 40))
mtext("B", at=a$mids[26], cex=1.5)
a <- hist(mean.r[not.extinct], nclass=25, col="dodgerblue", main="",
    xlab="Population growth rate", axes=FALSE)
axis(1)
axis(2, at=c(0, 2000, 4000, 6000, 8000, 10000), labels=c(0, 2, 4, 6, 8, 10))
mtext("C", at=a$mids[1], cex=1.5)
par(op)
