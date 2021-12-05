# Schaub & Kéry (2022) Integrated Population Models
# Chapter 4 : Components of integrated population models
# ------------------------------------------------------

# Run time testing 3 mins, full run 40 mins

library(IPMbook) ; library(jagsUI)

# 4.3 Models for population size surveys
# ======================================

# 4.3.6 The ‘demographic’ state-space model of Dail and Madsen
# ------------------------------------------------------------

# Choose constants in simulation
nsite <- 100                            # Number of sites
nvisit <- 2                             # Number of short-term replicate visits
nyear <- 12                             # Number of years
gamma <- 5                              # Initial expected abundance
phi <- 0.82                             # Apparent survival probability
rho <- 0.22                             # Recruitment rate (per-capita)
p <- 0.4                                # Detection probability (per individual)

# State or ecological process
# Simulate true system state
N <- array(NA, dim=c(nsite, nyear))     # Empty abundance matrix

# Set initial abundance (year 1)
set.seed(24)
N[,1] <- rpois(nsite, gamma)
head(N, 11)                             # Look at first 11 sites (not shown)

# Propagate abundance forwards via transition rule = pop.dyn.model
for (t in 1:(nyear-1)){
  N[,t+1] <- rbinom(nsite, N[,t], phi) + rpois(nsite, N[,t] * rho)
}
head(N, 11)                             # Look at first 11 sites (not shown)

# ~~~~ Plot the true abundance patterns (not shown) ~~~~
matplot(1:nyear, t(N), type='l', lty=1, lwd=3, xlab='Year of study',
    ylab='Abundance', frame=FALSE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Observation process: Simulate counts with imperfect detection
C <- array(NA, dim=c(nsite, nvisit, nyear))
for (i in 1:nsite){
  for (t in 1:nyear){
    C[i,,t] <- rbinom(nvisit, N[i,t], p)
  } #t
} #i

# Look at data in year 6 and compare with truth (for illustration; not shown)
cbind(Truth=N[,6], 'Count 1'=C[,1,6], 'Count 2'=C[,2,6])

trueN <- apply(N, 2, sum)                       # True number of pairs
obsN <- apply(apply(C, c(1, 3), max), 2, sum)   # Observed number of pairs

# ~~~~ Plot true and observed population size (not shown) ~~~~
matplot(1:nyear, cbind(trueN, obsN), type='b', lty=1, pch=16, col=c('red', 'black'),
    ylim=c(0, max(trueN)), lwd=3, xlab='Year of study', ylab='Number of pairs',
    frame=FALSE, cex=2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Data bundle
jags.data <- list(C=C, nsite=nsite, nvisit=nvisit, nyear=nyear)
str(jags.data)
# List of 4
# $ C     : int [1:100, 1:2, 1:12] 2 0 2 3 3 3 2 4 2 2 ...
# $ nsite : num 100
# $ nvisit: num 2
# $ nyear : num 12

# Write JAGS model file
cat(file = "model7.txt", "
model {
  # Priors and linear models
  gamma ~ dunif(0, 10)                  # Expected abundance in year 1
  phi ~ dunif(0, 1)                     # Apparent survival
  rho ~ dunif(0, 1)                     # Per-capita recruitment rate
  p ~ dunif(0, 1)                       # Per-individual detection probability

  # Likelihood
  for (i in 1:nsite){
    # Model for the initial population size
    N[i,1] ~ dpois(gamma)

    # Process model over time: our model of population dynamics
    for (t in 2:nyear){
      S[i,t] ~ dbinom(phi, N[i,t-1])    # Survivors
      R[i,t] ~ dpois(N[i,t-1] * rho)    # Recruits
      N[i,t] <- S[i,t] + R[i,t]
    } #t
  } #i

  # Observation process
  for (i in 1:nsite){
    for (j in 1:nvisit){
      for (t in 1:nyear){
        C[i,j,t] ~ dbinom(p, N[i,t])
      } #t
    } #j
  } #i

  # Derived quantities
  for (t in 1:nyear){
    ntot[t] <- sum(N[,t])               # Total abundance across all sites
  }
}
")

# Initial values
# Need to get inits for N that do not contradict model and data!
Nst <- apply(C, c(1,3), max) + 5
Nst[,2:nyear] <- NA
Rst <- apply(C, c(1,3), max) + 1        # Observed max. counts + 1 as inits
Rst[,1] <- NA
inits <- function(){list(N=Nst, R=Rst)}

# Parameters monitored
parameters <- c("gamma", "phi", "rho", "p", "ntot", "N")

# MCMC settings
# ni <- 200000; nb <- 100000; nc <- 3; nt <- 100; na <- 5000 # 45 min
ni <- 10000; nb <- 5000; nc <- 3; nt <- 5; na <- 500         # ~~~ for testing, 2 mins

# Call JAGS from R (ART 45 min), check convergence and summarize posteriors
out10 <- jags(jags.data, inits, parameters, "model7.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out10, c("gamma", "phi", "rho", "p", "ntot"))
print(out10)    # Not shown

# Compare truth and estimates
truth <- c('gamma'=gamma, 'phi'=phi, 'rho'=rho, 'p'=p, Ntot=apply(N, 2, sum))
print(cbind(truth, out10$summary[1:16,c(1:3, 7)]), 3)

#         truth    mean      sd    2.5%   97.5%
# gamma    5.00   4.919  0.3520   4.240   5.646
# phi      0.82   0.800  0.0264   0.744   0.847
# rho      0.22   0.240  0.0267   0.193   0.298
# p        0.40   0.396  0.0201   0.356   0.435
# Ntot1  476.00 490.508 26.9491 441.000 548.025
# Ntot2  509.00 516.968 27.4614 467.000 576.000
# Ntot3  537.00 546.495 28.5580 494.000 608.000
# Ntot4  562.00 559.367 29.8402 507.000 623.000
# Ntot5  585.00 576.243 31.1611 520.000 643.000
# Ntot6  623.00 601.033 32.5235 542.000 671.000
# Ntot7  646.00 636.714 33.3183 577.975 709.025
# Ntot8  655.00 659.010 34.3807 596.000 733.025
# Ntot9  675.00 681.968 35.7977 617.000 758.000
# Ntot10 691.00 694.363 37.7847 624.975 774.000
# Ntot11 736.00 728.986 39.6545 657.000 812.000
# Ntot12 773.00 758.156 42.9149 680.000 850.000
