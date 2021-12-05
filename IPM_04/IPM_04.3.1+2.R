# Schaub & Kéry (2022) Integrated Population Models
# Chapter 4 : Components of integrated population models
# ------------------------------------------------------


# 4.3 Models for population size surveys
# ======================================

# 4.3.1 Gaussian state-space models
# ---------------------------------

# Choose constants
nyears <- 25                                        # Number of years
N1 <- 30                                            # Abundance at t = 1
mu.lam <- 1.02                                      # Mean of the distribution of lambda
sig2.lam <- 0.02                                    # Variance of the distribution of lambda
sig2.y <- 400                                       # Variance of observation error

# Simulate true system state
N <- numeric(nyears)
N[1] <- N1                                          # Set initial abundance
set.seed(1)                                         # Initialize the RNGs
lambda <- rnorm(nyears-1, mu.lam, sqrt(sig2.lam))   # Draw random lambdas
for (t in 1:(nyears-1)){
  N[t+1] <- lambda[t] * N[t]                        # Propagate population size forwards
}

# Simulate observations
eps <- rnorm(nyears, 0, sqrt(sig2.y))               # Draw random residuals
y <- N + eps                                        # Add residual (error) to value of true state

# ~~~~ code for Figure 4.1 ~~~~
op <- par(cex=1.25)
ylim <- c(min(c(N, y)), max(c(N, y)))
plot(N, xlab='Year', ylab='Number', type='b', pch=16, ylim=ylim, col='red', axes=FALSE)
points(y, pch=16, col='black')
segments(1:nyears, N, 1:nyears, y, col='black')
axis(1, at=1:nyears, tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
axis(2, las=1)
legend('topleft', c('True value of state', 'Observations'),
    pch=16, col=c('red', 'black'), bty='n', cex=1)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(IPMbook); library(jagsUI)

# Data bundle
jags.data <- list(y=y, T=length(y))
str(jags.data)
# List of 2
# $ y: num [1:25] 42.4 26.82 26.11 -3.06 23.27 ...
# $ T: int 25

# Write JAGS model file
cat(file="model1.txt", "
model {
  # Priors and linear models
  mu.lam ~ dunif(0, 10)                             # Prior for mean growth rate
  sig.lam ~ dunif(0, 1)                             # Prior for sd of growth rate
  sig2.lam <- pow(sig.lam, 2)
  tau.lam <- pow(sig.lam, -2)
  sig.y ~ dunif(0.1, 100)                           # Prior for sd of observation process
  sig2.y <- pow(sig.y, 2)
  tau.y <- pow(sig.y, -2)

  # Likelihood
  # Model for the initial population size: uniform priors
  N[1] ~ dunif(0, 500)

  # Process model over time: our model of population dynamics
  for (t in 1:(T-1)){
    lambda[t] ~ dnorm(mu.lam, tau.lam)
    N[t+1] <- N[t] * lambda[t]
  }

  # Observation process
  for (t in 1:T){
    y[t] ~ dnorm(N[t], tau.y)
  }
}
")

# Initial values
inits <- function(){list(sig.lam=runif(1, 0, 1))}

# Parameters monitored
parameters <- c("lambda", "mu.lam", "sig2.y", "sig2.lam", "sig.y", "sig.lam", "N")

# MCMC settings
ni <- 200000; nb <- 10000; nc <- 3; nt <- 100; na <- 5000

# Call JAGS from R (ART <1 min), check convergence and summarize posteriors
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out1) # Not shown
print(out1, 3)
#               mean     sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# lambda[1]    1.002  0.156   0.591   1.033   1.264    FALSE 1 1.001  2659
# lambda[2]    1.012  0.145   0.656   1.035   1.280    FALSE 1 1.002  3990
# [... output truncated ...]
# lambda[23]   1.071  0.108   0.871   1.057   1.331    FALSE 1 1.000  3275
# lambda[24]   0.952  0.151   0.577   1.004   1.149    FALSE 1 1.000  5700
# mu.lam       1.053  0.036   0.989   1.050   1.140    FALSE 1 1.001  5700
# sig2.y     246.009 96.637  98.641 229.619 477.841    FALSE 1 1.002  1145
# sig2.lam     0.024  0.038   0.000   0.009   0.132    FALSE 1 1.012  1073
# sig.y       15.403  2.958   9.932  15.153  21.860    FALSE 1 1.002  1157
# sig.lam      0.119  0.098   0.005   0.095   0.363    FALSE 1 1.003  1611
# N[1]        28.222  8.115  15.371  27.046  47.294    FALSE 1 1.000  5700
# N[2]        27.535  6.468  15.820  27.159  41.492    FALSE 1 1.001  5700
# [... output truncated ...]

# ~~~~ code for Figure 4.2 ~~~~
# Plot of true and observed and estimated states
library(scales)
op <- par(cex=1.25)
ylim <- c(min(c(N, y)), max(c(N, y)))
plot(N, xlab='Year', ylab='Number', type='n', pch=16, ylim=ylim, col='red', axes=FALSE)
polygon(c(1:nyears, nyears:1), c(out1$q2.5$N, rev(out1$q97.5$N)),
    col=alpha('blue', 0.15), border=NA)
axis(1, at=1:nyears, tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
axis(2, las=1)
points(N, pch=16, type='b', col='red')
points(y, pch=16, col='black')
segments(1:nyears, N, 1:nyears, y, col='black')
points(out1$mean$N, pch=16, type='b', col='blue')
legend('topleft', c('True value of state', 'Observations', 'Estimate of state'),
    pch=16, col = c('red', 'black', 'blue'), bty='n', lwd=rep(1, 3))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 4.3.2 Effects of ‘evil’ patterns in the measurement error
#       on a Gaussian state-space model
# ---------------------------------------------------------

# Simulate observations in Problem Case 1: p is constant at 0.8
set.seed(24)
p1 <- 0.8
y1 <- rbinom(length(N), round(N), p1)

# Simulate observations in Problem Case 2: p declines from 0.8 to 0.1
p2 <- seq(0.8, 0.1, length.out=length(N))
y2 <- rbinom(length(N), round(N), p2)

# ~~~~ code for Figure 4.3 ~~~~
# Plot both cases along with N
library(scales)
co <- viridis_pal(option='E')(20)[c(6, 16)]
op <- par(mar=c(4.5, 4.2, 1, 1), las=1)
ylim <- c(min(c(N, y)), max(c(N, y)))
plot(N, xlab='Year', ylab='Number', type='b', pch=16, ylim=ylim, col='red', axes=FALSE)
points(y1, pch=16, col=co[1], type='b')
points(y2, pch=16, col=co[2], type='b')
axis(1, at=1:nyears, tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
axis(2, las=1)
legend('topleft', c('True value of state', 'Observations (Problem case 1)',
    'Observations (Problem case 2)'), pch=c(16, 16, 16),
    col=c('red', co), bty='n', lwd=rep(1, 3))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Bundle the two new data sets
jags.data1 <- list(y=y1, T=length(y1))
jags.data2 <- list(y=y2, T=length(y2))

# MCMC settings
ni <- 100000; nb <- 50000; nc <- 3; nt <- 50; na <- 5000

# Fit the same Gaussian SSM to both new data sets
# Call JAGS from R (ART <1 min), check convergence and summarize posteriors
out2 <- jags(jags.data1, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
out3 <- jags(jags.data2, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out2) # Not shown
traceplot(out3) # Not shown
print(out2, 3) # Not shown
print(out3, 3) # Not shown

# ~~~~ code for Figure 4.4 ~~~~
# Plot of true and observed and estimated states
library(scales)
co <- viridis_pal(option='E')(20)[c(6, 16)]
ylim <- c(min(c(N, y)), max(c(N, y)))
plot(N, xlab='Year', ylab='Number', type='b', pch=16, ylim=ylim, col='red', axes=FALSE)
polygon(c(1:nyears, nyears:1), c(out2$q2.5$N, rev(out2$q97.5$N)),
    col=alpha(co[1], 0.3), border=NA)
polygon(c(1:nyears, nyears:1), c(out3$q2.5$N, rev(out3$q97.5$N)),
    col=alpha(co[2], 0.3), border=NA)
points(y1, pch=16, col=co[1], type='b')
points(y2, pch=16, col=co[2], type='b')
axis(1, at=1:nyears, tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20, 25), tcl=-0.5, labels=c(5, 10, 15, 20, 25))
axis(2, las=1)
legend('topleft', c('True value of state', 'Observations (Problem case 1)',
    'Observations (Problem case 2)'), pch=c(16, 16, 16),
    col=c('red', co), bty='n', lwd=rep(1, 3))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~ save output for use in subsequent sections ~~~
save(out1, out2, out3, file="IPM_04.3.1+2_output.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
