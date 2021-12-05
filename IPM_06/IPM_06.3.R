# Schaub & Kéry (2022) Integrated Population Models
# Chapter 6 : Benefits of integrated population modeling
# ------------------------------------------------------

# Run time approx. 2.5 mins

library(IPMbook) ; library(jagsUI)

# 6.3 Estimation of demographic parameters for which there is no explicit data
# ============================================================================

# Write JAGS model file
cat(file="model1.txt", "
model {
  # Priors and linear models
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }

  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- mean.f / 2 * mean.sj * (N[1,t] + N[2,t])
    N[2,t+1] <- mean.sa * (N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Capture-recapture data (CJS model with multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.occasions-1)){
    # Main diagonal
    q[t] <- 1 - p[t]                    # Probability of non-recapture
    pr.j[t,t] <- sj[t] * p[t]
    pr.a[t,t] <- sa[t] * p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- sj[t] * prod(sa[(t+1):j]) * prod(q[t:(j-1)]) * p[j]
      pr.a[t,j] <- prod(sa[t:j]) * prod(q[t:(j-1)]) * p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
  }

  # Derived parameters
  # Annual population growth rate
  for (t in 1:(n.occasions-1)){
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t])
  }
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
  }
}
")

library(IPMbook); library(jagsUI)
data(woodchat6)
str(woodchat6)
# List of 6
# $ ch   : num [1:947, 1:10] 1 1 1 1 1 0 1 1 1 0 ...
# $ age  : num [1:947] 2 2 2 2 2 2 2 2 2 2 ...
# $ count: num [1:10] 110 104 100 85 85 71 118 112 91 104
# $ J    : num [1:10] 147 144 132 131 178 178 235 169 177 186
# $ B    : num [1:10] 48 48 45 37 53 56 74 59 55 60
# $ f    : num [1:535] 4 7 5 5 1 5 4 1 3 0 ...

marr <- marrayAge(woodchat6$ch, woodchat6$age)

# Bundle data
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), C=woodchat6$count)

# Initial values
inits <- function(){list(mean.sj=runif(1, 0, 0.5))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "sigma", "N", "ann.growth.rate", "Ntot")

# MCMC settings
ni <- 40000; nb <- 10000; nc <- 3; nt <- 3; na <- 2000

# Call JAGS (ART <1 min), check convergence and summarize posteriors
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out1)
print(out1, 3)
#                       mean     sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# mean.sj              0.253  0.022   0.213   0.253   0.297    FALSE 1 1.000 11961
# mean.sa              0.567  0.022   0.523   0.567   0.610    FALSE 1 1.000 30000
# mean.p               0.604  0.031   0.543   0.604   0.664    FALSE 1 1.000 30000
# mean.f               3.413  0.384   2.717   3.394   4.238    FALSE 1 1.000 30000
# sigma               18.669  6.076  10.877  17.415  33.982    FALSE 1 1.000 14425
# N[1,1]              49.371 29.747   3.123  47.859 103.593    FALSE 1 1.003   752
# N[2,1]              50.932 29.696   3.375  51.263 103.410    FALSE 1 1.003   719
# [... output truncated ... ]
# N[1,10]             41.636  6.521  28.889  41.594  54.558    FALSE 1 1.000 30000
# N[2,10]             54.722  6.078  42.418  54.759  66.695    FALSE 1 1.000 30000
# ann.growth.rate[1]   0.996  0.023   0.948   0.996   1.040    FALSE 1 1.000 30000
# [... output truncated ... ]
# ann.growth.rate[9]   0.996  0.023   0.948   0.996   1.040    FALSE 1 1.000 30000
# Ntot[1]            100.304 11.905  78.043  99.764 125.839    FALSE 1 1.000 30000
# [... output truncated ... ]
# Ntot[10]            96.358 11.774  72.446  96.441 119.201    FALSE 1 1.000 30000


# Model for productivity and population count data only
# '''''''''''''''''''''''''''''''''''''''''''''''''''''

# Write JAGS model file
cat(file="model2.txt", "
model {
  # Priors and linear models
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)
  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- mean.f / 2 * mean.sj * (N[1,t] + N[2,t])
    N[2,t+1] <- mean.sa * (N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    count[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Productivity data (Poisson regression model)
  nJ ~ dpois(n.rep * mean.f)

  # Derived parameters
  # Annual population growth rate
  for (t in 1:(n.occasions-1)){
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t])
  }
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
  }
}
")

# Bundle data
jags.data <- list(n.occasions=length(woodchat6$count), nJ=sum(woodchat6$J),
    n.rep=sum(woodchat6$B), count=woodchat6$count)

# Initial values
inits <- function(){list(mean.sj=runif(1, 0, 0.5))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.f", "sigma", "N", "ann.growth.rate", "Ntot")

# MCMC settings
ni <- 1000000; nb <- 50000; nc <- 3; nt <- 500; na <- 2000

# Call JAGS (ART 1 min), check convergence and summarize posteriors
out2 <- jags(jags.data, inits, parameters, "model2.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out2)
print(out2, 3)
#                       mean     sd   2.5%    50%   97.5% overlap0 f  Rhat n.eff
# mean.sj              0.317  0.183  0.018  0.318   0.619    FALSE 1 1.001  1201
# mean.sa              0.500  0.286  0.029  0.500   0.966    FALSE 1 1.001  1205
# mean.f               3.136  0.076  2.991  3.134   3.289    FALSE 1 1.000  5700
# sigma               18.668  5.912 10.963 17.410  33.405    FALSE 1 1.000  4729
# N[1,1]              49.994 29.343  3.387 48.571 103.446    FALSE 1 1.000  5265
# N[2,1]              50.019 29.069  2.962 50.145 101.707    FALSE 1 1.000  4450
# [ ... output truncated ... ]
# N[1,10]             48.257 28.429  2.759 47.690  99.152    FALSE 1 1.001  1426
# N[2,10]             48.477 28.350  2.683 48.059  99.484    FALSE 1 1.001  1260
# ann.growth.rate[1]   0.997  0.022  0.950  0.997   1.041    FALSE 1 1.000  5700
# [ ... output truncated ... ]
# ann.growth.rate[9]   0.997  0.022  0.950  0.997   1.041    FALSE 1 1.000  5700
# Ntot[1]            100.014 11.773 77.965 99.496 124.997    FALSE 1 1.001  5700
# [ ... output truncated ... ]
# Ntot[10]            96.734 11.539 73.067 96.764 119.929    FALSE 1 1.000  5700


# Model with population count data alone
# ''''''''''''''''''''''''''''''''''''''

# Write JAGS model file
cat(file="model3.txt", "
model {
  # Priors and linear models
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)
  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- mean.f / 2 * mean.sj * (N[1,t] + N[2,t])
    N[2,t+1] <- mean.sa * (N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    count[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Derived parameters
  # Annual population growth rate
  for (t in 1:(n.occasions-1)){
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t])
  }
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
  }
}
")

# Bundle data
jags.data <- list(n.occasions=length(woodchat6$count), count=woodchat6$count)

# Initial values
inits <- function(){list(mean.sj=runif(1, 0, 0.5))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.f", "sigma", "N", "ann.growth.rate", "Ntot")

# MCMC settings
ni <- 1000000; nb <- 50000; nc <- 3; nt <- 500; na <- 2000

# Call JAGS (ART 1 min), check convergence and summarize posteriors
out3 <- jags(jags.data, inits, parameters, "model3.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out3)
print(out3, 3)
#                       mean     sd   2.5%    50%   97.5% overlap0 f  Rhat n.eff
# mean.sj              0.352  0.267  0.018  0.277   0.937    FALSE 1 1.000  5700
# mean.sa              0.586  0.285  0.039  0.620   0.982    FALSE 1 1.002   951
# mean.f               3.447  2.623  0.166  2.737   9.228    FALSE 1 1.003   729
# sigma               18.860  6.173 11.089 17.473  34.453    FALSE 1 1.000  5700
# N[1,1]              49.962 29.036  3.965 48.763 102.803    FALSE 1 1.000  4570
# N[2,1]              50.352 29.215  3.399 50.355 102.760    FALSE 1 1.000  4231
# [ ... output truncated ... ]
# N[1,10]             39.787 28.355  1.320 35.676  96.392    FALSE 1 1.002   970
# N[2,10]             56.734 28.456  3.851 59.533 102.348    FALSE 1 1.002  1048
# ann.growth.rate[1]   0.996  0.023  0.948  0.997   1.039    FALSE 1 1.001  5700
# [ ... output truncated ... ]
# ann.growth.rate[9]   0.996  0.023  0.948  0.997   1.039    FALSE 1 1.001  5700
# Ntot[1]            100.314 12.037 78.599 99.593 126.060    FALSE 1 1.001  5700
# [ ... output truncated ... ]
# Ntot[10]            96.521 11.738 73.273 96.525 119.578    FALSE 1 1.000  5700


# ~~~~ extra code for figure 6.3 ~~~~
cat(file="model14.txt", "
model {
  # Priors and linear models
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }

  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)

  # Process model over time
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- mean.f/2 * mean.sj * (N[1,t] + N[2,t])
    N[2,t+1] <- mean.sa * (N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Capture-recapture data (CJS model with multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.occasions-1)){
    # Main diagonal
    q[t] <- 1-p[t]   # Probability of non-recapture
    pr.j[t,t] <- sj[t]*p[t]
    pr.a[t,t] <- sa[t]*p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- sj[t]*prod(sa[(t+1):j])*prod(q[t:(j-1)])*p[j]
      pr.a[t,j] <- prod(sa[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
  }

  # Productivity data (Poisson regression model)
  nJ ~ dpois(n.rep * mean.f)

  # Derived parameters
  # Annual population growth rate
  for (t in 1:(n.occasions-1)){
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t])
  }
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
  }
}
")

# Bundle data
jags.data <- list(marr.j=marr[,,1], marr.a=marr [,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr [,,2]), nJ=sum(woodchat6$J),
    n.rep=sum(woodchat6$B), C= woodchat6$count)

# Initial values
inits <- function(){list(mean.sj=runif(1, 0, 0.5))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "sigma", "N",
    "ann.growth.rate", "Ntot")

# MCMC settings
ni <- 40000; nb <- 10000; nc <- 3; nt <- 3; na <- 2000

# Call JAGS (ART <1 min)
out14 <- jags(jags.data, inits, parameters, "model14.txt",
    n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=TRUE)

# Produce figure 6.3
mag <- 1
cex.tif <- mag * 1
lwd.tif <- mag
lwd.tif2 <- mag*1.25
op <- par(las=1, mar=c(4,4,1,1), cex=cex.tif, "mfrow")
layout(matrix(1:6, 2, 3, byrow=TRUE), widths=c(1.05, 1, 1), heights=c(1, 1), TRUE)
co <- c("red", "dodgerblue", "darkolivegreen", "orange")
plot(density(out14$sims.list$mean.sj), xlim=c(0, 1), main="",
    xlab=expression('Juvenile survival ('*italic(s)[italic(j)]*')'),
    lwd=lwd.tif2, col=co[1], axes=FALSE)
lines(density(out1$sims.list$mean.sj), col=co[2], lwd=lwd.tif2)
lines(density(out2$sims.list$mean.sj), col=co[3], lwd=lwd.tif2)
lines(density(out3$sims.list$mean.sj), col=co[4], lwd=lwd.tif2)
axis(1, lwd=lwd.tif)
axis(2, lwd=lwd.tif)
legend("topright", legend=c("CMR, Prod., Count", "CMR, Count", "Prod., Count", "Count"),
    col=co, lwd=rep(lwd.tif2,4), bty="n")

plot(density(out14$sims.list$mean.sa), xlim=c(0, 1), main="",
    xlab=expression('Adult survival ('*italic(s)[italic(a)]*')'),
    lwd=lwd.tif2, col=co[1], ylab=NA, axes=FALSE)
lines(density(out1$sims.list$mean.sa), col=co[2], lwd=lwd.tif2)
lines(density(out2$sims.list$mean.sa), col=co[3], lwd=lwd.tif2)
lines(density(out3$sims.list$mean.sa), col=co[4], lwd=lwd.tif2)
axis(1, lwd=lwd.tif)
axis(2, lwd=lwd.tif)

plot(density(out14$sims.list$mean.f), xlim=c(0, 10), main="",
    xlab=expression('Productivity ('*italic(f)*')'),
    lwd=lwd.tif2, col=co[1], ylab=NA, axes=FALSE)
lines(density(out1$sims.list$mean.f), col=co[2], lwd=lwd.tif2)
lines(density(out2$sims.list$mean.f), col=co[3], lwd=lwd.tif2)
lines(density(out3$sims.list$mean.f), col=co[4], lwd=lwd.tif2)
axis(1, lwd=lwd.tif)
axis(2, lwd=lwd.tif)

plot(density(out14$sims.list$N[,1,10]), xlim=c(0, 130), main="",
    xlab=expression('Number of first year indviduals ('*italic(t)*' = 10)'),
    lwd=lwd.tif2, col=co[1], axes=FALSE)
lines(density(out1$sims.list$N[,1,10]), col=co[2], lwd=lwd.tif2)
lines(density(out2$sims.list$N[,1,10]), col=co[3], lwd=lwd.tif2)
lines(density(out3$sims.list$N[,1,10]), col=co[4], lwd=lwd.tif2)
axis(1, lwd=lwd.tif)
axis(2, lwd=lwd.tif)

plot(density(out14$sims.list$N[,2,10]), xlim=c(0, 130), main="",
    xlab=expression('Number of adults ('*italic(t)*' = 10)'),
    lwd=lwd.tif2, col=co[1], ylab=NA, axes=FALSE)
lines(density(out1$sims.list$N[,2,10]), col=co[2], lwd=lwd.tif2)
lines(density(out2$sims.list$N[,2,10]), col=co[3], lwd=lwd.tif2)
lines(density(out3$sims.list$N[,2,10]), col=co[4], lwd=lwd.tif2)
axis(1, lwd=lwd.tif)
axis(2, lwd=lwd.tif)

plot(density(out14$sims.list$Ntot[,10]), xlim=c(60, 140), main="",
    xlab=expression('Total population size ('*italic(t)*' = 10)'),
    lwd=lwd.tif2, col=co[1], ylab=NA, axes=FALSE)
lines(density(out1$sims.list$Ntot[,10]), col=co[2], lwd=lwd.tif2)
lines(density(out2$sims.list$Ntot[,10]), col=co[3], lwd=lwd.tif2)
lines(density(out3$sims.list$Ntot[,10]), col=co[4], lwd=lwd.tif2)
axis(1, lwd=lwd.tif)
axis(2, lwd=lwd.tif)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code for the simulations for figure 6.4 is in the script "IPM_06.3_sims.R"

# ~~~~ extra code for figure 6.5 ~~~~
data(woodchat6)
marr <- marrayAge(woodchat6$ch, woodchat6$age)

# The productivity data are not aggregated anymore, because it is easier to thin them if the outcome of each brood is used as data. Consequently, the Poisson regression model to analyze these data has been reformulated.

# Bundle data
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), C=woodchat6$count,
    f=woodchat6$f, n.rep=length(woodchat6$f))

# Write JAGS model file
cat(file="model15.txt", "
model {
  # Priors and linear models
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }

  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)

  # Process model over time: our model for population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- mean.f / 2 * mean.sj * (N[1,t] + N[2,t])
    N[2,t+1] <- mean.sa * (N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Capture-recapture data (CJS model with multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.occasions-1)){
    # Main diagonal
    q[t] <- 1 - p[t]   # Probability of non-recapture
    pr.j[t,t] <- sj[t] * p[t]
    pr.a[t,t] <- sa[t] * p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- sj[t] * prod(sa[(t+1):j]) * prod(q[t:(j-1)]) * p[j]
      pr.a[t,j] <- prod(sa[t:j]) * prod(q[t:(j-1)]) * p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
  }

  # Productivity data (Poisson regression model)
  for (i in 1:n.rep){
    f[i] ~ dpois(mean.f)
  }

  # Derived parameters
  # Annual population growth rate
  for (t in 1:(n.occasions-1)){
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t])
  }
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
  }
}
")

# Initial values
inits <- function(){list(mean.sj=runif(1, 0, 0.5))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "sigma", "N", "ann.growth.rate")

# MCMC settings
ni <- 20000; nb <- 10000; nc <- 3; nt <- 2; na <- 2000

# Call JAGS (ART <1 min)
m1 <- jags(jags.data, inits, parameters, "model15.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

# only 90% of the productivity data are available
n.fl <- round(length(woodchat6$f) * 0.9)
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), C=woodchat6$count,
    f=sample(woodchat6$f, n.fl, replace=TRUE), n.rep=n.fl)

m2 <- jags(jags.data, inits, parameters, "model15.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

# only 80% of the productivity data are available
n.fl <- round(length(woodchat6$f) * 0.8)
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), C=woodchat6$count,
    f=sample(woodchat6$f, n.fl, replace=TRUE), n.rep=n.fl)
m3 <- jags(jags.data, inits, parameters, "model15.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

# only 70% of the productivity data are available
n.fl <- round(length(woodchat6$f) * 0.7)
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), C=woodchat6$count,
    f=sample(woodchat6$f, n.fl, replace=TRUE), n.rep=n.fl)
m4 <- jags(jags.data, inits, parameters, "model15.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

# only 60% of the productivity data are available
n.fl <- round(length(woodchat6$f) * 0.6)
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), C=woodchat6$count,
    f=sample(woodchat6$f, n.fl, replace=TRUE), n.rep=n.fl)
m5 <- jags(jags.data, inits, parameters, "model15.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

# only 50% of the productivity data are available
n.fl <- round(length(woodchat6$f) * 0.5)
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), C=woodchat6$count,
    f=sample(woodchat6$f, n.fl, replace=TRUE), n.rep=n.fl)
m6 <- jags(jags.data, inits, parameters, "model15.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

# only 40% of the productivity data are available
n.fl <- round(length(woodchat6$f) * 0.4)
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), C=woodchat6$count,
    f=sample(woodchat6$f, n.fl, replace=TRUE), n.rep=n.fl)
m7 <- jags(jags.data, inits, parameters, "model15.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

# only 30% of the productivity data are available
n.fl <- round(length(woodchat6$f) * 0.3)
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), C=woodchat6$count,
    f=sample(woodchat6$f, n.fl, replace=TRUE), n.rep=n.fl)
m8 <- jags(jags.data, inits, parameters, "model15.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

# only 20% of the productivity data are available
n.fl <- round(length(woodchat6$f) * 0.2)
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), C=woodchat6$count,
    f=sample(woodchat6$f, n.fl, replace=TRUE), n.rep=n.fl)
m9 <- jags(jags.data, inits, parameters, "model15.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

# only 10% of the productivity data are available
n.fl <- round(length(woodchat6$f) * 0.1)
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), C=woodchat6$count,
    f=sample(woodchat6$f, n.fl, replace=TRUE), n.rep=n.fl)
m10 <- jags(jags.data, inits, parameters, "model15.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

# only 5% of the productivity data are available
n.fl <- round(length(woodchat6$f) * 0.05)
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), C=woodchat6$count,
    f=sample(woodchat6$f, n.fl, replace=TRUE), n.rep=n.fl)
m11 <- jags(jags.data, inits, parameters, "model15.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

# no productivity data
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), C=woodchat6$count)
m12 <- jags(jags.data, inits, parameters, "model1.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

sd.prod <- c(m1$sd$mean.f, m2$sd$mean.f, m3$sd$mean.f, m4$sd$mean.f, m5$sd$mean.f,
    m6$sd$mean.f, m7$sd$mean.f, m8$sd$mean.f, m9$sd$mean.f, m10$sd$mean.f,
    m11$sd$mean.f, m12$sd$mean.f)

plot(c(100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0), sd.prod, type="b",
    pch=16, axes=FALSE, ylab="SD (productivity)",
    xlab="Percentage of productivity data (%)")
axis(2, las=1)
axis(1)

save(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, file="Data Fig 6.5.Rdata")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
