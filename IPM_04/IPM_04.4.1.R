# Schaub & Kéry (2022) Integrated Population Models
# Chapter 4 : Components of integrated population models
# ------------------------------------------------------

# Run time approx. 3 mins

library(IPMbook) ; library(jagsUI)

# 4.4 Models for productivity surveys
# ===================================

# 4.4.1 Poisson models for brood size data
# ----------------------------------------

# Choose constants in simulation
nbrood <- 1000                            # Number of broods with young counted
brood.mean <- 1.5                         # Average brood size
sd.brood <- 0.3                           # log-linear brood random effect

# Draw Poisson random numbers
set.seed(24)
expNyoung <- exp(log(brood.mean) + rnorm(nbrood, 0, sd.brood))
C <- rpois(nbrood, expNyoung)

# ~~~~ plot it ~~~~
plot(table(C), lwd=20, col='grey', lend='butt', frame=FALSE, xlab='Brood size',
    ylab='Number')
# ~~~~~~~~~~~~~~~~~

table(C)
# C
#   0   1   2   3   4   5   6   7   8
# 231 330 237 134  39  16   7   3   3

# Data bundle
jags.data <- list(C1=C, sumC1=sum(C), C1copy=C, nbrood=nbrood)
str(jags.data)
# List of 4
# $ C1    : int [1:1000] 0 0 1 1 4 1 0 2 0 1 ...
# $ sumC1 : int 1529
# $ C1copy: int [1:1000] 0 0 1 1 4 1 0 2 0 1 ...
# $ nbrood: num 1000

# Write JAGS model file
cat(file="model8.txt", "
model {
  # Priors and linear models
  rho1 ~ dunif(0, 5)                      # Mean brood size in model 1
  rho2 ~ dunif(0, 5)                      # Mean brood size in model 2
  rho3 ~ dunif(0, 5)                      # Mean brood size in model 3
  tau.rho3 <- pow(sd.rho3, -2)
  sd.rho3 ~ dunif(0, 3)                   # Brood-level overdispersion in model 3

  # Likelihoods for three separate models
  # Model 1: Poisson GLM for disaggregated data
  for (i in 1:nbrood){
    C1[i] ~ dpois(rho1)
  }

  # Model 2: Poisson GLM for aggregated data
  sumC1 ~ dpois(rho2 * nbrood)

  # Model 3: Poisson GLMM for aggregated data with brood-level overdispersion
  for (i in 1:nbrood){
    C1copy[i] ~ dpois(pois.mean[i])
    log(pois.mean[i]) <- logmean[i]
    logmean[i] ~ dnorm(log(rho3), tau.rho3)
  }
}
")

# Initial values
inits <- function(){list(rho1=runif(1, 0.5, 2.5), rho2=runif(1, 0.5, 2.5),
    rho3=runif(1, 0.5, 2.5))}

# Parameters monitored
parameters <- c("rho1", "rho2", "rho3", "sd.rho3")

# MCMC settings
ni <- 60000; nb <- 10000; nc <- 3; nt <- 10; na <- 1000

# Call JAGS from R (ART 3 min), check convergence and summarize posteriors
out11 <- jags(jags.data, inits, parameters, "model8.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out11)    # Not shown
print(out11, 3)

#              mean     sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# rho1        1.530  0.039    1.454    1.529    1.609    FALSE 1 1.000 15000
# rho2        1.530  0.039    1.453    1.530    1.606    FALSE 1 1.000 15000
# rho3        1.464  0.047    1.372    1.464    1.560    FALSE 1 1.002  1145
# sd.rho3     0.292  0.061    0.138    0.297    0.393    FALSE 1 1.031   156
