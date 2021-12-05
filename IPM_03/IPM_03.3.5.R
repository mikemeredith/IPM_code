# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

# Run time approx. 3 mins

library(IPMbook)

# 3.3 Classical analysis of a matrix population model
# ===================================================

# 3.3.5 Analysis of a matrix population model with different sources of
#       stochasticity and parameter uncertainty
# ---------------------------------------------------------------------

# Define mean, measurement error and temporal variability of the demographic parameters
mean.sj <- 0.3          # Mean value of juv. survival
sd.sj.e <- 0.005        # Uncertainty of mean juv. survival as SD on natural scale
sd.sj.t <- 0.25         # Temporal variability of juv. survival as SD on logit scale
mean.sa <- 0.55         # Mean value of ad. survival
sd.sa.e <- 0.005        # Uncertainty of mean ad. survival as SD on natural scale
sd.sa.t <- 0.07         # Temporal variability of ad. survival as SD on logit scale
mean.f1 <- 1.3          # Mean value of productivity of 1y females
sd.f1.e <- 0.05         # Uncertainty of mean productivity as SD on natural scale
sd.f1.t <- 0.3          # Temporal variability of productivity as SD on natural scale
mean.fa <- 1.8          # Mean value of productivity of adult females
sd.fa.e <- 0.03         # Uncertainty of mean productivity as SD on natural scale
sd.fa.t <- 0.3          # Temporal variability of productivity as SD on natural scale

# Define the number of years with predictions and the Monte Carlo setting
T <- 200                # Number of years (projection time frame)
nsim <- 100000          # Number of replicate populations simulated

# Define population matrix and initial stage-specific population sizes
N <- array(NA, dim=c(2, T+1, nsim))
N[,1,] <- c(10, 10)
r <- matrix(NA, nrow=T, ncol=nsim)
alive <- matrix(NA, nrow=T, ncol=nsim)
mean.r <- numeric(nsim)

# Project population
for (s in 1:nsim){ # Loop over replicate populations
  if(s %% 100 == 0) {cat(paste("*** Simrep", s, "***\n")) }   # Counter
  # Generate a mean of the demographic rates (subject to measurement error)
  msj <- rbeta2(1, mean.sj, sd.sj.e)
  msa <- rbeta2(1, mean.sa, sd.sa.e)
  mf1 <- rnorm(1, mean.f1, sd.f1.e)
  mfa <- rnorm(1, mean.fa, sd.fa.e)

  # Generate annual demographic rates (subject to temporal variability)
  sj <- plogis(rnorm(T, qlogis(msj), sd.sj.t))
  sa <- plogis(rnorm(T, qlogis(msa), sd.sa.t))
  f1 <- pmax(0, rnorm(T, mf1, sd.f1.t))                       # Avoids negative values
  fa <- pmax(0, rnorm(T, mfa, sd.fa.t))                       # Dito

  # Project population (include demographic stochasticity)
  for (t in 1:T){                                             # Loop over years
    N[1,t+1,s] <- rpois(1, sj[t] * (f1[t] * N[1,t,s] + fa[t] * N[2,t,s]))
    N[2,t+1,s] <- rbinom(1, sum(N[,t,s]), sa[t])
    if (sum(N[,t+1,s]) == 0) break
    r[t,s] <- log(sum(N[,t+1,s])) - log(sum(N[,t,s]))
    alive[t,s] <- t
  } #t
  mean.r[s] <- mean(r[min(alive[,s], na.rm=TRUE):max(alive[,s], na.rm=TRUE),s])
} #s

mean(mean.r)
# [1] -0.01562592
sd(mean.r)
# [1] 0.06659048

not.extinct <- which(!is.na(alive[T,]))
mean(mean.r[not.extinct])
# [1] 0.02353411
sd(mean.r[not.extinct])
# [1] 0.01242122

# Extinction probability (after T years)
sum(is.na(alive[T,])) / nsim
# [1] 0.38023

# ~~~~ Plot graph with population growth rates (Fig. 3.15)  ~~~~
op <- par(mar=c(4, 4, 2, 1), las=1, cex=1.1, "mfrow")
layout(matrix(c(1, 1, 2, 3), 2, 2, byrow=TRUE), widths=c(1, 1), heights=c(1, 1), TRUE)
plot(r[,1], type="l", lwd=0.5, ylab="Annual population growth rate", xlab="Time",
    ylim=range(r[which(!is.na(alive))]), col="lightgrey", axes=FALSE)
axis(1); axis(2)
for (s in 2:nsim){
  lines(r[!is.na(alive[,s]),s], lwd=0.5, col="lightgrey")
}
lines(apply(r, 1, mean, na.rm=TRUE), lwd=1.5)
mtext("A", at=0, cex=1.5)
a <- hist(mean.r, nclass=100, col="dodgerblue", main="",
    xlab="Population growth rate", xlim=c(-0.5, 0.1), axes=FALSE)
axis(1)
axis(2, at=c(0, 5000, 10000, 15000, 20000, 25000, 30000), labels=c(0, 5, 10, 15, 20, 25, 30))
mtext("B", at=a$mids[51], cex=1.5)
a <- hist(mean.r[not.extinct], nclass=25, col="dodgerblue", main="",
    xlab="Population growth rate", axes=FALSE)
axis(1)
axis(2, at=c(0, 2000, 4000, 6000, 8000, 10000), labels=c(0, 2, 4, 6, 8, 10))
mtext("C", at=a$mids[1], cex=1.5)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
