# Schaub & Kéry (2022) Integrated Population Models
# Chapter 4 : Components of integrated population models
# ------------------------------------------------------

library(IPMbook) ; library(jagsUI)

# 4.3 Models for population size surveys
# ======================================

# 4.3.5 Transitioning from Gaussian to discrete-valued state-space
#       models for population counts
# ----------------------------------------------------------------

# Choose constants
nyear <- 25                             # Number of years
N1 <- 30                                # Initial abundance
mu.lam <- 1.02                          # Mean of the distribution of lambda
sig2.lam <- 0.02                        # Variance of the distribution of lambda

# Simulate true system state
N <- numeric(nyear)
N[1] <- N1                              # Set initial abundance
set.seed(1)                             # Initialize the RNGs
lambda <- rnorm(nyear-1, mu.lam, sqrt(sig2.lam)) # Draw random lambda
for (t in 1:(nyear-1)){
  N[t+1] <- rpois(1, lambda[t] * N[t])  # Propagate population size forwards
}

# Simulate observations
y <- rpois(nyear, N)                    # Observation error is now Poisson noise

# Data bundle
jags.data <- list(y=y, T=length(y))
str(jags.data)
# List of 2
# $ y: int [1:25] 26 26 40 38 24 21 16 33 25 31 ...
# $ T: int 25

# Write JAGS model file
cat(file = "model6.txt", "
model {
  # Priors and linear models
  mu.lam ~ dunif(0, 10)                 # Prior for mean growth rate
  sig.lam ~ dunif(0, 2)                 # Prior for sd of growth rate
  sig2.lam <- pow(sig.lam, 2)
  tau.lam <- pow(sig.lam, -2)

  # Likelihood
  # Model for the initial population size: uniform priors
  N[1] ~ dunif(0, 500)

  # Process model over time: our model of population dynamics
  for (t in 1:(T-1)){
    lambda[t] ~ dnorm(mu.lam, tau.lam)
    N[t+1] ~ dpois(N[t] * lambda[t])
  }

  # Observation process
  for (t in 1:T){
    y[t] ~ dpois(N[t])
  }
}
")

# Initial values
inits <- function(){list(sig.lam=runif(1, 0, 1), mu.lam=runif(1, 0.1, 2),
    N=round(runif(nyear, 20, 40)))}

# Parameters monitored
parameters <- c("mu.lam", "sig2.lam", "sig.lam", "N")

# MCMC settings
ni <- 20000; nb <- 10000; nc <- 3; nt <- 10; na <- 1000

# Call JAGS from R (ART <1 min), check convergence and summarize posteriors
out9 <- jags(jags.data, inits, parameters, "model6.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out9)     # Not shown
print(out9, 3)      # Not shown

# ~~~~ Plot of true and observed and estimated states (Fig. 4.8) ~~~~
ylim <- c(min(c(N, y)), max(c(N, y)))
plot(N, xlab='Year', ylab='Number', type='b', pch=16, ylim=ylim, col='red', axes=FALSE)
axis(1, at=1:nyear, tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
axis(2, las=1)
points(y, pch=16, col='black')
segments(1:nyear, N, 1:nyear, y, col='black')
points(out9$mean$N, pch=16, type='b', col='blue')
polygon(c(1:nyear, nyear:1), c(out9$q2.5$N, rev(out9$q97.5$N)),
    col=scales::alpha('blue', 0.15), border=NA)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
