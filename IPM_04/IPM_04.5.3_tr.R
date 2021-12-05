# Schaub & Kéry (2022) Integrated Population Models
# Chapter 4 : Components of integrated population models
# ------------------------------------------------------

# Run time for test script 70 secs, full run 11 mins

library(IPMbook) ; library(jagsUI)

# 4.5 Models for survival surveys
# ===============================

# 4.5.3 Dead-recovery data
# ------------------------

# 4.5.3.1 State-space formulation
# '''''''''''''''''''''''''''''''

# Choose constants in simulation
nmarked <- 100                          # Number of marked individuals at each occasion
nyears <- 11                            # Number of years
s <- 0.8                                # Survival probability
r <- 0.2                                # Recovery probability

# Determine occasion when an individual first captured and marked
f <- rep(1:(nyears-1), each=nmarked)
nind <- length(f)                       # Total number of marked individuals

# State or ecological process
# Simulate true system state
z <- array(NA, dim=c(nind, nyears))     # Empty alive/dead matrix

# Initial conditions: all individuals alive at f(i)
for (i in 1:nind){
  z[i,f[i]] <- 1
}

set.seed(3)                             # Initialize the RNGs in R
# Propagate alive/dead process forwards via transition rule
# Alive individuals survive with probability s
for (i in 1:nind){
  for (t in (f[i]+1):nyears){
    z[i,t] <- rbinom(1, 1, z[i,t-1] * s)
  } #t
} #i
z                                       # Look at the true state (not shown, but insightful)

# Observation process: simulate observations
y <- array(0, dim=c(nind, nyears))
for (i in 1:nind){
  y[i,f[i]] <- 1
  for(t in (f[i]+1):nyears){
    y[i,t] <- rbinom(1, 1, (z[i,t-1]-z[i,t]) * r)
  } #t
} #i

y                                       # Complete simulated data set(not shown)
# Compute the total number of dead recoveries
sum(rowSums(y, na.rm=TRUE)==2)
# [1] 131

# Data bundle
jags.data <- list(y=y, f=f, nind=nind, nyears=ncol(y))
str(jags.data)
# List of 4
# $ y     : num [1:1000, 1:11] 1 1 1 1 1 1 1 1 1 1 ...
# $ f     : int [1:1000] 1 1 1 1 1 1 1 1 1 1 ...
# $ nind  : int 1000
# $ nyears: num 11

# Write JAGS model file
cat(file="model20.txt", "
model {
  # Priors and linear models
  s.const ~ dunif(0, 1)                 # Vague prior for constant s
  r.const ~ dunif(0, 1)                 # Vague prior for constant r

  for (i in 1:nind){                    # Loop over individuals
    for (t in f[i]:(nyears-1)){         # Loop over time intervals/occasions
      s[i,t] <- s.const                 # Here model pattern in s ...
      r[i,t] <- r.const                 # ... and r
    } #t
  } #i

  # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):nyears){
      # State process
      z[i,t] ~ dbern(z[i,t-1] * s[i,t-1])
      # Observation process
      y[i,t] ~ dbern((z[i,t-1]- z[i,t]) * r[i,t-1])
    } #t
  } #i
}
")

# Initial values
inits <- function(){list(z=zInitDR(y))} # Function in IPMbook package

# Parameters monitored
parameters <- c("s.const", "r.const")

# MCMC settings
# ni <- 30000; nb <- 5000; nc <- 3; nt <- 10; na <- 5000
ni <- 3000; nb <- 500; nc <- 3; nt <- 1; na <- 500  # ~~~ for testing, 1 min

# Call JAGS from R (ART 13 min), check convergence and summarize posteriors
out23 <- jags(jags.data, inits, parameters, "model20.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out23) # Not shown
acf(out23$sims.list$s.const)            # Check out autocorrelation for s (not shown)
print(out23, 3)
#             mean     sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# s.const    0.777  0.037   0.709   0.775   0.852    FALSE 1 1.015  1992
# r.const    0.198  0.024   0.159   0.195   0.256    FALSE 1 1.015   811


# 4.5.3.2 Multinomial formulation
# '''''''''''''''''''''''''''''''

marr <- marrayDead(y)

# Data bundle
jags.data <- list(marr=marr, rel=rowSums(marr), nyears=ncol(y))
str(jags.data)
# List of 3
# $ marr  : num [1:10, 1:11] 4 0 0 0 0 0 0 0 0 0 ...
# $ rel   : num [1:10] 100 100 100 100 100 100 100 100 100 100
# $ nyears: int 11

# Write JAGS model file
cat(file="model21.txt", "
model {
  # Priors and linear models
  s.const ~ dunif(0, 1)                 # Vague prior for constant s
  r.const ~ dunif(0, 1)                 # Vague prior for constant r

  for (t in 1:(nyears-1)){              # Loop over time intervals/occasions
    s[t] <- s.const                     # Here model pattern in s ...
    r[t] <- r.const                     # ... and r
  } #t

  # Likelihood
  for (t in 1:(nyears-1)){
    marr[t,1:nyears] ~ dmulti(pi[t,], rel[t])
  }

  # Define the cell probabilities of the m-array
  for (t in 1:(nyears-1)){
    # Main diagonal
    pi[t,t] <- (1-s[t])*r[t]
    # Above main diagonal
    for (j in (t+1):(nyears-1)){
      pi[t,j] <- prod(s[t:(j-1)])*(1-s[j])*r[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pi[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recovery
  for (t in 1:(nyears-1)){
    pi[t,nyears] <- 1-sum(pi[t,1:(nyears-1)])
  } #t
}
")

# Initial values
inits <- function(){list(s.const=runif(1))}

# Parameters monitored
parameters <- c("s.const", "r.const")

# MCMC settings
ni <- 3000; nb <- 1000; nc <- 3; nt <- 1; na <- 3000

# Call JAGS from R (ART <1 min), check convergence and summarize posteriors
out24 <- jags(jags.data, inits, parameters, "model21.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out24) # Not shown
print(out24)
#             mean    sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# s.const    0.771 0.035   0.701   0.771   0.840    FALSE 1 1.001  3343
# r.const    0.195 0.022   0.157   0.193   0.244    FALSE 1 1.001  4078
