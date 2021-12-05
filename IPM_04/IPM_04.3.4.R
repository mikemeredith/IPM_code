# Schaub & Kéry (2022) Integrated Population Models
# Chapter 4 : Components of integrated population models
# ------------------------------------------------------

# Run time approx. 90 secs

library(IPMbook) ; library(jagsUI)

# 4.3 Models for population size surveys
# ======================================

# 4.3.4 Correction of population count data for coverage
#       bias and detection bias
# ------------------------------------------------------

# Choose constants in simulation
nterritory <- 100                         # Number of territories (= sites)
nvisit <- 2                               # Number of short-term replicate visits
nyear <- 12                               # Number of years
psi <- 0.66                               # Initial occupancy probability
phi <- 0.9                                # Persistence probability
gamma <- 0.2                              # Colonization probability
p <- 0.4                                  # Detection probability (per visit)

# State or ecological process
z <- array(NA, dim=c(nterritory, nyear))  # Empty presence/absence matrix

# Set initial presence/absence (in year 1)
set.seed(24)
z[,1] <- rbinom(nterritory, 1, psi)

# Propagate presence/absence forwards via transition rule
for (i in 1:nterritory){
  for (t in 1:(nyear-1)){
    z[i,t+1] <- rbinom(1, 1, z[i,t] * phi + (1-z[i,t]) * gamma)
  } #t
} #i

# Plot the true presence/absence pattern (Fig. 4.6)
mapPalette <- colorRampPalette(c("white", "#C7B76BFF"))
image(x = 1:nyear, y=1:nterritory, z=t(z), col=mapPalette(10), axes=TRUE, xlab="Year",
    ylab="Territory")

# Observation process (1): Imperfect detection
y1 <- array(NA, dim=c(nterritory, nvisit, nyear))
for (i in 1:nterritory){
  for (t in 1:nyear){
    y1[i,,t] <- rbinom(nvisit, 1, z[i,t] * p)
  } #t
} #i

# Observation process (2): Incomplete coverage of territory
pvisit <- seq(0.3, 0.9, length.out=nyear) # Visitation probability
y2 <- y1                                  # Make a copy
for (i in 1:nterritory){
  for (t in 1:nyear){
    visit <- rbinom(1, 1, pvisit[t])      # 1 if visit takes place
    y2[i,,t] <- y1[i,,t] / visit
  } #t
} #i
y2[!is.finite(y2)] <- NA
y2[1:11,,1:3]                             # Look at data from first 11 sites and 3 years

# True population size (= number of occupied territories)
trueN <- apply(z, 2, sum)
# Observed population size if you had visited each territory twice a year
obsN1 <- apply(apply(y1, c(1, 3), max), 2, sum)
# Observed population size if you not had visited each territory twice a year
tmp <- apply(y2, c(1, 3), max, na.rm=TRUE)
tmp[tmp == -Inf] <- NA
obsN2 <- apply(tmp, 2, sum, na.rm=TRUE)

# ~~~~ Plot true and observed population size (Fig. 4.7 upper) ~~~~
op <- par(mfrow=c(2, 1), mar=c(4, 4, 4, 2))
library(scales)
co <- viridis_pal(option='E')(20)[c(4, 14)]
plot(trueN, type='b', lty=1, pch=16, col='red', ylim=c(0, 100),
    ylab='Number of territories', axes=FALSE)
points(obsN1, type='b', lty=1, pch=16, col=co[1])
points(obsN2, type='b', lty=1, pch=16, col=co[2])
axis(2, las=1)
axis(1)
legend(y=110, x=4, c("True population size", "Observed, imperfect detection (y1)",
    "Observed, additional coverage bias (y2)"), pch=rep(16, 3), col=c('red', co), bty='n')
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Fit the model to our first data set
# '''''''''''''''''''''''''''''''''''

# Data bundle
jags.data <- list(y=y1, nterritory=nterritory, nvisit=nvisit, nyear=nyear)
str(jags.data)
# List of 4
# $ y         : num [1:100, 1:2, 1:12] 0 0 0 0 0 0 0 0 0 0 ...
# $ nterritory: num 100
# $ nvisit    : num 2
# $ nyear     : num 12

# Write JAGS model file
cat(file = "model5.txt", "
model {
  # Priors and linear models
  psi ~ dunif(0, 1)
  phi ~ dunif(0, 1)
  gamma ~ dunif(0, 1)
  p ~ dunif(0, 1)

  # Likelihood
  for (i in 1:nterritory){
    # Model for initial occupancy state
    z[i,1] ~ dbern(psi)

    # Model for occupancy state dynamics
    for (t in 1:(nyear-1)){
      z[i, t+1] ~ dbern(z[i, t] * phi + (1 - z[i, t]) * gamma)
    } #t
  } #i

  # Observation process
  for (i in 1:nterritory){
    for (j in 1:nvisit) {
      for (t in 1:nyear){
        y[i,j,t] ~ dbern(z[i, t] * p)
      } #t
    } #j
  } #i

  # Derived quantities: population size
  # correcting for both coverage bias and nondetection bias
  for (t in 1:nyear){
    ntot[t] <- sum(z[,t])                 # Number of occupied territories = population size
  }
}
")

# Initial values
zst <- array(1, dim=c(nterritory, nyear))
inits <- function(){list(z=zst)}

# Parameters monitored
parameters <- c("psi", "phi", "gamma", "p", "ntot", "z")

# MCMC settings
ni <- 6000; nb <- 1000; nc <- 3; nt <- 5; na <- 1000

# Call JAGS from R (ART <1 min), check convergence and summarize posteriors
out7 <- jags(jags.data, inits, parameters, "model5.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out7)   # Not shown
print(out7, 3)

# Fit the same model to our second data set
# '''''''''''''''''''''''''''''''''''''''''

# Data bundle
jags.data <- list(y=y2, nterritory=nterritory, nvisit=nvisit, nyear=nyear)
str(jags.data) # Summary of data set (not shown)

# Call JAGS from R (ART <1 min), check convergence and summarize posteriors
out8 <- jags(jags.data, inits, parameters, "model5.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out8)     # Not shown
print(out8, 3)      # Not shown


# ~~~~ Plot true and observed population size (Fig. 4.7 upper) ~~~~
op <- par(mfrow=c(2, 1), mar=c(4, 4, 4, 2))
library(scales)
co <- viridis_pal(option='E')(20)[c(4, 14)]
plot(trueN, type='b', lty=1, pch=16, col='red', ylim=c(0, 100),
    ylab='Number of territories', axes=FALSE)
points(obsN1, type='b', lty=1, pch=16, col=co[1])
points(obsN2, type='b', lty=1, pch=16, col=co[2])
axis(2, las=1)
axis(1)
legend(y=110, x=4, c("True population size", "Observed, imperfect detection (y1)",
    "Observed, additional coverage bias (y2)"), pch=rep(16, 3), col=c('red', co), bty='n')

# Plot the estimates (lower panel of Fig. 4.7)
plot(trueN, type='b', lty=1, pch=16, col='red', ylim=c(0, 100),
    xlab='Year of study', ylab='Number of territories', axes=FALSE)
polygon(c(1:nyear, nyear:1), c(out7$q2.5$ntot, rev(out7$q97.5$ntot)),
    col=alpha(co[1], 0.25), border=NA)
polygon(c(1:nyear, nyear:1), c(out8$q2.5$ntot, rev(out8$q97.5$ntot)),
    col=alpha(co[2], 0.45), border=NA)
points(trueN, type='b', lty=1, pch=16, col='red')
points(out7$mean$ntot, type='b', lty=1, pch=15, col=co[1])
points(out8$mean$ntot, type='b', lty=1, pch=15, col=co[2])
axis(2, las=1)
axis(1)
legend(y=40, x=4,c('True population size', 'Population size estimate from y1',
    'Population size estimate from y2'), pch=c(16, rep(15, 2)), col=c('red', co), bty='n')
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
