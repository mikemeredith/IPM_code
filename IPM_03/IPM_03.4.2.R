# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

library(IPMbook) ; library(jagsUI)

# 3.4 Analysis of matrix population models with Markov Chain
#     Monte Carlo (MCMC) software
# ==========================================================

# 3.4.2 Analysis of a matrix population model with parameter uncertainty
# ----------------------------------------------------------------------

# Define mean and SD of the demographic parameters
mean.sj <- 0.3        # Point estimate of juv. survival
se.sj.e <- 0.03       # Uncertainty of juv. survival as SE on natural scale
mean.sa <- 0.55       # Point estimate of ad. survival
se.sa.e <- 0.015      # Uncertainty of ad. survival as SE on natural scale
mean.f1 <- 1.3        # Point estimate of productivity of 1y females
se.f1.e <- 0.3        # Uncertainty of productivity as SE on natural scale
mean.fa <- 1.8        # Point estimate of productivity of adult females
se.fa.e <- 0.1        # Uncertainty of productivity as SE on natural scale

# Bundle data
jags.data <- list(alpha.sj=getBeta2Par(mean.sj, se.sj.e)[1], beta.sj=getBeta2Par(mean.sj,
    se.sj.e)[2], alpha.sa=getBeta2Par(mean.sa, se.sa.e)[1], beta.sa=getBeta2Par(mean.sa, se.sa.e)[2],
    mean.f1=mean.f1, tau.f1=1/se.f1.e^2, mean.fa=mean.fa, tau.fa=1/se.fa.e^2, T=50)

# Write JAGS model file
cat(file="model3.txt", "
model {
  # Random number generators (RNGs)
  sj ~ dbeta(alpha.sj, beta.sj)           # These only *look* like priors
  sa ~ dbeta(alpha.sa, beta.sa)           # ... but they are not
  f1 ~ dnorm(mean.f1, tau.f1)             # ... as there is no estimation in this model
  fa ~ dnorm(mean.fa, tau.fa)

  # Initialize the population size nodes
  N[1,1] <- 1
  N[2,1] <- 1

  # Loop over time
  for (t in 1:T){
    # Population model
    N[1,t+1] <- sj * (f1 * N[1,t] + fa * N[2,t])
    N[2,t+1] <- sa * (N[1,t] + N[2,t])

    # Annual (realized) growth rate
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t])

    # Scaled stage distributions
    stage.distr[1,t] <- N[1,t+1] / (N[1,t+1] + N[2,t+1])
    stage.distr[2,t] <- N[2,t+1] / (N[1,t+1] + N[2,t+1])
  }
  lambda <- ann.growth.rate[T]
  stable.stage.distr <- stage.distr[,T]

  # Sensitivity and elasticity of lambda to changes in sj
  delta <- 0.001                          # size of perturbation
  N.star[1,1] <- 1
  N.star[2,1] <- 1
  for (t in 1:T){
    N.star[1,t+1] <- (sj + delta) * (f1 * N.star[1,t] + fa * N.star[2,t])
    N.star[2,t+1] <- sa * (N.star[1,t] + N.star[2,t])
    ann.growth.rate.star[t] <- (N.star[1,t+1] + N.star[2,t+1]) / (N.star[1,t] + N.star[2,t])
  }
  s.sj <- (ann.growth.rate.star[T] - ann.growth.rate[T]) / delta
  e.sj <- s.sj * sj / lambda

  # Calculation of net reproductive rate (R0)
  for (i in 1:100){
    u[i] <- pow(sa, i)
  }
  R0 <- sj * f1 + sj * fa * sum(u[])

  # Calculation of generation time (GT)
  GT <- log(R0) / log(lambda)
}
")

# Parameters monitored
parameters <- c("lambda", "stable.stage.distr", "s.sj", "e.sj", "R0", "GT")

# MCMC settings
ni <- 100000; nt <- 1; nb <- 0; nc <- 1; na <- 0

# Call JAGS (ART <1 min) and summarize results
out3 <- jags(jags.data, NULL, parameters, "model3.txt", n.adapt=na, n.chains=nc, n.thin=nt,
    n.iter=ni, n.burnin=nb, DIC=FALSE)
print(out3, 4)
#                         mean     sd   2.5%    50%  97.5% overlap0 f
# lambda                1.0223 0.0627 0.9102 1.0188 1.1553    FALSE 1
# stable.stage.distr[1] 0.4603 0.0323 0.3979 0.4600 0.5241    FALSE 1
# stable.stage.distr[2] 0.5397 0.0323 0.4759 0.5400 0.6021    FALSE 1
# s.sj                  1.4658 0.1935 1.1135 1.4567 1.8693    FALSE 1
# e.sj                  0.4273 0.0452 0.3438 0.4258 0.5203    FALSE 1
# R0                    1.0516 0.1490 0.7798 1.0453 1.3619    FALSE 1
# GT                    2.3842 0.1910 2.0426 2.3727 2.7911    FALSE 1

# ~~~~ Produce figure 3.17 ~~~~
op <- par(mar=c(4, 4, 3, 0), las=1, cex=1.1, "mfrow")
layout(matrix(1:4, 2, 2, byrow=TRUE), widths=c(1.1, 1), heights=c(1, 1), TRUE)
a <- hist(out3$sims.list$lambda, nclass=50, col="dodgerblue", main="",
    xlab=expression(paste("Asymptotic population growth rate (", lambda, ")")), prob=TRUE)
mtext("A", at=a$mids[2], cex=1.5)
par(mar=c(4, 2, 3, 2))
a <- hist(out3$sims.list$stable.stage[,1], nclass=50, col="dodgerblue", main="",
    xlab="Proportion of 1-year old individuals", prob=TRUE)
mtext("B", at=a$mids[2], cex=1.5)
par(mar=c(4, 4, 3, 0))
a <- hist(out3$sims.list$s.sj, nclass=50, col="dodgerblue", main="",
    xlab=expression('Sensitivity of '*lambda*' to variation of '*italic(s)[italic(j)]),
    prob=TRUE)
mtext("C", at=a$mids[2], cex=1.5)
par(mar=c(4, 2, 3, 2))
a <- hist(out3$sims.list$GT, nclass=40, col="dodgerblue", main="",
    xlab="Generation time", prob=TRUE)
mtext("D", at=a$mids[2], cex=1.5)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~ Produce figure 3.18 ~~~~
# Load output from section 3.3.2
load("IPM_03.3.2_output.RData")
op <- par(mar=c(4, 4, 3, 0), las=1, cex=1.1, "mfrow")
layout(matrix(1:4, 2, 2, byrow=TRUE), widths=c(1.1, 1), heights=c(1, 1), TRUE)
plot(density(lambda), main="",
    xlab=expression(paste("Asymptotic population growth rate (", lambda, ")")),
    col="blue", lwd=2, axes=FALSE)
lines(density(out3$sims.list$lambda), col="red", lwd=2, lty=2)
axis(1); axis(2)
# legend("topright", legend=c("classical MC", "Bayesian MCMC"),
#    lwd=c(2,2), col=c("blue","red"), bty="n")
mtext("A", at=0.8, cex=1.5)

par(mar=c(4, 2, 3, 2))
plot(density(stable.stage[,1]), main="", xlab="Proportion of 1-year old individuals",
    ylab="", col="blue", lwd=2, axes=FALSE)
axis(1); axis(2)
lines(density(out3$sims.list$stable.stage[,1]), col="red", lwd=2, lty=2)
mtext("B", at=0.35, cex=1.5)

par(mar=c(4, 4, 3, 0))
plot(density(sensitivity[,1]), main="",
    xlab=expression('Sensitivity of '*lambda*' to variation of '*italic(s)[italic(j)]),
    col="blue", lwd=2, axes=FALSE)
axis(1); axis(2)
lines(density(out3$sims.list$s.sj), col="red", lwd=2, lty=2)
mtext("C", at=0.83, cex=1.5)

par(mar=c(4, 2, 3, 2))
plot(density(GT), main="", xlab="Generation time", ylab="", col="blue",
    lwd=2, axes=FALSE)
axis(1); axis(2)
lines(density(out3$sims.list$GT), col="red", lwd=2, lty=2)
mtext("D", at=1.9, cex=1.5)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
