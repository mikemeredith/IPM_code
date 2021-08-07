# Schaub & Kery (2021) Integrated Population Models
# Chapter 4 : Components of integrated population models
# ------------------------------------------------------
# Code from proofs.

# Run time approx. 7 mins

library(IPMbook) ; library(jagsUI)

# 4.4 Models for productivity surveys
# ===================================

# 4.4.2 Zero-inflation in brood size data
# ---------------------------------------

# Choose constants in simulation
nbrood <- 1000                          # Number of broods with young counted
theta <- 0.7                            # Success probability in Bernoulli process
brood.mean <- 1.5                       # Average brood size in conditional Poisson process
sd.brood <- 0.3                         # Overdispersion in conditional Poisson process

# Simulate Bernoulli process dividing broods in failures and 'potential successes'
set.seed(46)
z <- rbinom(nbrood, 1, theta)           # z = 1 means 'potential success'

# Draw conditional Poisson random numbers with overdispersion
expNyoung <- z * exp(log(brood.mean) + rnorm(nbrood, 0, sd.brood))
Cx <- rpois(nbrood, expNyoung)

# ~~~~ plot the data ~~~~
plot(table(Cx), lwd=20, col='grey', lend='butt', frame=FALSE,
    xlab='Brood size', ylab='Number')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

table(Cx)
# Cx
#   0   1   2   3   4   5   6   7   9
# 466 234 167  76  31  14   9   1   2

# Data bundle
jags.data <- list(C=Cx, nbrood=nbrood)
str(jags.data)
# List of 2
# $ C : int [1:1000] 3 0 5 1 2 2 0 1 0 0 ...
# $ nbrood: num 1000

# Write JAGS model file
cat(file="model9.txt", "
model {
  # Priors and linear models
  theta ~ dunif(0, 1)                   # Success probabiliy in Bernoulli process
  rho ~ dunif(0, 5)                     # Mean brood size in conditional Poisson
  tau.rho <- pow(sd.rho, -2)
  sd.rho ~ dunif(0, 3)                  # Overdispersion in conditional Poisson

  # Likelihood
  # Note this is a zero-inflated Poisson log-normal GLMM .... cool !
  for (i in 1:nbrood){
    z[i] ~ dbern(theta) # Zero-inflation process
    C[i] ~ dpois(z[i] * cond.pois.mean[i])        # Conditional Poisson
    log(cond.pois.mean[i]) <- log.condmean[i]
    log.condmean[i] ~ dnorm(log(rho), tau.rho)    # Overdispersion
  }
}
")

# Initial values
inits <- function(){list(z = rep(1, nbrood), rho=runif(1, 0.5, 2.5))}

# Parameters monitored
parameters <- c("theta", "rho", "sd.rho")

# MCMC settings
ni <- 110000; nb <- 10000; nc <- 3; nt <- 100; na <- 1000

# Call JAGS from R (ART 6 min) and check convergence
out12 <- jags(jags.data, inits, parameters, "model9.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out12)                        # Not shown
print(out12, 3)
#              mean     sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# theta       0.705  0.033    0.648    0.704    0.774    FALSE 1 1.000  3000
# rho         1.456  0.107    1.244    1.456    1.656    FALSE 1 1.001  3000
# sd.rho      0.314  0.099    0.079    0.323    0.487    FALSE 1 1.006  3000
