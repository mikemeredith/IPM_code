# Schaub & Kéry (2022) Integrated Population Models
# Chapter 4 : Components of integrated population models
# ------------------------------------------------------

# Run time testing 6 mins, full run 12 mins

# ~~~ load results from earlier sections ~~~
load("IPM_04.3.1+2_output.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 4.3 Models for population size surveys
# ======================================

# 4.3.3 Use of estimates from another analysis in a
#       Gaussian state-space model or an IPM
# -------------------------------------------------

library(IPMbook); library(jagsUI)
?simMHB                                       # Check out help function

# Some examples of how the simMHB function works
# Explicit default values for all function arguments
str(dat <- simMHB(nsites=267, nsurveys=3, nyears=25, mean.lam=1, mean.beta=0.03,
    sd.lam=c(0.5, 0.05), mean.p=0.6, beta.p=0.1, show.plot=TRUE))
str(dat <- simMHB())                          # Same, implicit
str(dat <- simMHB(nsites=1000))               # More sites
str(dat <- simMHB(nsurveys=10))               # More surveys
str(dat <- simMHB(nyears=50))                 # More years
str(dat <- simMHB(mean.lam=5))                # Higher mean abundance
str(dat <- simMHB(mean.beta=-0.03))           # Population declines
str(dat <- simMHB(sd.lam=c(0, 0)))            # No site variability in lambda
str(dat <- simMHB(mean.p=1))                  # Perfect detection
str(dat <- simMHB(mean.p=0.6, beta.p=0))      # Constant p = 0.6
str(dat <- simMHB(mean.p=0.6, beta.p=-0.2))   # Declining p
str(dat <- simMHB(show.plot=FALSE))           # No plots (when used in simulations)

# Create one data set for our hypothetical surveys of jays
set.seed(1982)
str(MHBdata <- simMHB(nsites=267, nsurveys=3, nyears=25, mean.lam=1, mean.beta=0.03,
    sd.lam=c(0.5, 0.05), mean.p=0.6, beta.p=0, show.plot=TRUE))

# Bundle and summarize data set
jags.data <- with(MHBdata, (list(C=C, nsites=dim(C)[1], nsurveys=dim(C)[2], nyears=dim(C)[3])))
str(jags.data)
# List of 4
# $ C       : int [1:267, 1:3, 1:25] 0 0 0 2 0 0 0 1 0 3 ...
# $ nsites  : int 267
# $ nsurveys: int 3
# $ nyears  : int 25

# Write JAGS model file
cat(file = "model2.txt", "
model {
  # Priors
  for (t in 1:nyears){                              # Loop over years
    lambda[t] ~ dunif(0, 100)                       # Expected abundance
    p[t] ~ dunif(0, 1)                              # Detection probability
  }

  # Ecological model for true abundance
  for (i in 1:nsites){                              # Loop over 267 sites
    for (t in 1:nyears){                            # Loop over 25 years
      N[i,t] ~ dpois(lambda[t])
      # Observation model for replicated counts
      for (j in 1:nsurveys){                        # Loop over 3 occasions
        C[i,j,t] ~ dbin(p[t], N[i,t])
      } #j
    } #t
  } #i

  # Total abundance across all surveyed sites as a derived quantity
  for (t in 1:nyears){
    totalN[t] <- sum(N[,t])
  }
}
")

# Initial values
Nst <- apply(MHBdata$C, c(1, 3), max, na.rm=TRUE) + 1
Nst[Nst == '-Inf'] <- 1
nyears <- dim(MHBdata$C)[3]
inits <- function() list(N=Nst, lambda=runif(nyears))

# Parameters monitored
parameters <- c("lambda", "p", "totalN")

# MCMC settings
ni <- 5000; nb <- 1000; nc <- 3; nt <- 4; na <- 1000

# Call JAGS from R (ART 3 min), check convergence and summarize posteriors
out4 <- jags(jags.data, inits, parameters, "model2.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out4) # Not shown
print(out4, 2)

#                 mean      sd      2.5%       50%     97.5% overlap0 f  Rhat n.eff
# totalN[1]    258.159   6.900   246.000   258.000   273.000    FALSE 1 1.000  3000
# totalN[2]    231.374   8.227   217.000   231.000   250.000    FALSE 1 1.001  1133
# totalN[3]    238.130   6.574   227.000   238.000   252.000    FALSE 1 1.001  2609
# [... output truncated ...]
# totalN[23]   418.273  11.255   399.000   417.000   443.000    FALSE 1 1.000  3000
# totalN[24]   477.464  11.496   457.000   477.000   502.000    FALSE 1 1.001  2138
# totalN[25]   522.101  13.020   499.000   521.000   551.000    FALSE 1 1.001  1565

vcMat <- array(NA, dim=c(nyears, nyears))
for (i in 1:nyears){
  for (j in 1:nyears){
    vcMat[i,j] <- cov(out4$sims.list$totalN[,i], out4$sims.list$totalN[,j])
  } #j
} #i
print(round(vcMat,1), 1)                            # Look at that beast! (not shown)

# ~~~~ code to check ~~~~
diag(vcMat)                             # Check out (not shown)
apply(out4$sims.list$totalN, 2, var)    # Compare that variances in diagonal are right .... and it seems OK
all.equal(diag(vcMat), apply(out4$sims.list$totalN, 2, var))
# ~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code to plot the results ~~~~
# Compare estimates of total N and true values (not shown)
op <- par(cex = 1.25)
ylim <- range(c(out4$q2.5$totalN, out4$q97.5$totalN), 600)
plot(dat$totalN, xlab = 'Year', ylab = 'Number', type = 'b', pch = 16,
    ylim = ylim, col = 'red', axes = FALSE)
axis(1, at = 1:dat$nyear, tcl = -0.25, labels = NA)
axis(1, at = c(5, 10, 15, 20, 25), tcl = -0.5, labels = c('5', '10', '15', '20', '25'))
axis(2, las = 1)
#box()
points(out4$mean$totalN, pch = 16, col = 'black')
segments(1:dat$nyear, out4$q2.5$totalN, 1:dat$nyear, out4$q97.5$totalN, col = 'black')
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Plot posterior distributions of totalN (results not shown)
# par(mfrow=c(3, 2))
op <- par(mfrow=c(3, 2), ask=dev.interactive(orNone=TRUE))  # ~~~ better for testing
for (t in 1:25){
  hist(out4$sims.list$totalN[,t], freq=FALSE, col='grey', main=paste('totalN in year', t))
  lines(density(out4$sims.list$totalN[,t]), col='blue', lwd=3)
  curve(dnorm(x, mean=out4$mean$totalN[t], sd=out4$sd$totalN[t]),
  range(out4$sims.list$totalN[,t]), col='red', lwd=3, add=TRUE)
  # browser()
}
par(op)


# Data bundle
jags.data <- list(Nhat=out4$mean$totalN, var.Nhat=out4$sd$totalN^2, T=length(out4$mean$totalN))
str(jags.data)
# List of 3
# $ Nhat    : num [1:25(1d)] 258 231 238 264 237 ...
# $ var.Nhat: num [1:25(1d)] 47.6 67.7 43.2 79.3 38 ...
# $ T       : int 25

# Write JAGS model file
cat(file="model3.txt", "
model {
  # Priors and linear models
  mu.lam ~ dunif(0, 10)                             # Prior for mean growth rate
  sig.lam ~ dunif(0, 10)                            # Prior for sd of growth rate
  sig2.lam <- pow(sig.lam, 2)
  tau.lam <- pow(sig.lam, -2)
  sig.ystar ~ dunif(0, 10000)                       # Prior for sd of observation process
  sig2.ystar <- pow(sig.ystar, 2)
  tau.ystar <- pow(sig.ystar, -2)

  # Likelihood
  # Model for the initial population size: uniform priors
  N[1] ~ dunif(0, 500)

  # Process model over time: our model of population dynamics
  for (t in 1:(T-1)){
    lambda[t] ~ dnorm(mu.lam, tau.lam)
    N[t+1] <- N[t] * lambda[t]
  }

  # Observation process for the abundance estimates with their SEs
  for (t in 1:T){
    Nhat[t] ~ dnorm(ystar[t], tau.se[t])
    tau.se[t] <- 1/var.Nhat[t]                      # Assumed known and given by variance of Nhat
    ystar[t] ~ dnorm(N[t], tau.ystar)
  }
}
")

# Initial values
inits <- function(){list(sig.lam = runif(1, 0, 1))}

# Parameters monitored
parameters <- c("lambda", "mu.lam", "sig2.ystar", "sig2.lam", "sig.ystar", "sig.lam", "N", "ystar")

# MCMC settings
ni <- 50000; nb <- 25000; nc <- 3; nt <- 25; na <- 10000

# Call JAGS from R (ART <1 min), check convergence and summarize posteriors
out5 <- jags(jags.data, inits, parameters, "model3.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out5) # Not shown
print(out5, 3) # Not shown

# ~~~~ code for Fig. 4.5 ~~~~
library(scales)
co <- viridis_pal(option='E')(20)[c(6, 16)]
ylim <- range(out4$q2.5$totalN, out4$q97.5$totalN, 600)
plot(MHBdata$totalN, xlab='Year', ylab='Number', type='b', pch=16, ylim=ylim,
    col='red', axes=FALSE)
polygon(c(1:nyears, nyears:1), c(out5$q2.5$N, rev(out5$q97.5$N)),
    col=alpha(co[2], 0.35), border=NA)
axis(1, at=1:nyears, tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
axis(2, las=1)
points(out4$mean$totalN, pch=16, col=co[1])
segments(1:nyears, out4$q2.5$totalN, 1:nyears, out4$q97.5$totalN, col=co[1])
points(out5$mean$N, pch=16, col=co[2], type='b')
legend('topleft', c('True value of state',
    expression(paste(hat(bolditalic(N)), ': from N-mix model')),
    expression(paste(bolditalic(N), ': from Gaussian SSM'))),
    pch = 16, col=c('red', co), bty='n')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Data bundle
jags.data <- list(Nhat=out4$mean$totalN, vcMat.Nhat=vcMat, T=length(out4$mean$totalN))
str(jags.data)
# List of 3
# $ Nhat      : num [1:25(1d)] 258 231 238 264 237 ...
# $ vcMat.Nhat: num [1:25, 1:25] 47.606 0.943 0.71 -0.144 -0.673 ...
# $ T         : int 25

# Write JAGS model file
cat(file="model4.txt", "
model {
  # Priors and linear models
  mu.lam ~ dunif(0, 10)                             # Prior for mean growth rate
  sig.lam ~ dunif(0, 10)                            # Prior for sd of growth rate
  sig2.lam <- pow(sig.lam, 2)
  tau.lam <- pow(sig.lam, -2)
  sig.ystar ~ dunif(0, 10000)                       # Prior for sd of observation process
  sig2.ystar <- pow(sig.ystar, 2)
  tau.ystar <- pow(sig.ystar, -2)

  # Likelihood
  # Model for the initial population size: uniform priors
  N[1] ~ dunif(0, 500)

  # Process model over time: our model of population dynamics
  for (t in 1:(T-1)){
    lambda[t] ~ dnorm(mu.lam, tau.lam)
    N[t+1] <- N[t] * lambda[t]
  }

  # Observation process for the abundance estimates with their VC matrix
  Nhat[1:T] ~ dmnorm.vcov(ystar[1:T], vcMat.Nhat[1:T, 1:T])
  for (t in 1:T){
    ystar[t] ~ dnorm(N[t], tau.ystar)
  }
}
")

# Initial values
Nst <- c(out5$mean$N[1], rep(NA, 24))
lamst <- rnorm(24, out2$mean$lambda, 0.1)
inits <- function(){list(N=Nst, lambda=lamst, sig.lam=runif(1, 0, 1))} # Works fine

# Parameters monitored
parameters <- c("lambda", "mu.lam", "sig2.ystar", "sig2.lam", "sig.ystar", "sig.lam", "N", "ystar")

# MCMC settings
# ni <- 50000; nb <- 20000; nc <- 3; nt <- 10; na <- 5000   # 8 mins
ni <- 5000; nb <- 2000; nc <- 3; nt <- 1; na <- 500  # ~~~ for testing

# Call JAGS from R (ART 11 min), check convergence and summarize posteriors
out6 <- jags(jags.data, inits, parameters, "model4.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out6) # Not shown
print(out6, 3) # Not shown

# Compare posterior means and SDs under the two models
print(cbind('pm (model 2)'=out5$mean$N, 'psd (model 2)'=out5$sd$N, 'pm (model 3)'=out6$mean$N,
    'psd (model 3)'=out6$sd$N), 3)

#       pm (model 2) psd (model 2) pm (model 3) psd (model 3)
# [1,]           242         13.01          244         12.82
# [2,]           239          9.80          240          9.61
# [... output truncated ...]
# [24,]          472         13.06          473         12.64
# [25,]          506         17.85          507         17.67
