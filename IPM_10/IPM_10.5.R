# Schaub & Kéry (2022) Integrated Population Models
# Chapter 10 : Population viability analysis
# ------------------------------------------

# Run time approx 9 mins


# 10.5 A PVA for simulated woodchat shrike data
# =============================================

# 10.5.1 Estimation of extinction probability and related quantities
# ------------------------------------------------------------------

# Load woodchat shrike data and produce data overview
library(IPMbook); library(jagsUI)
data(woodchat10)
str(woodchat10)
# List of 5
# $ marr.a: num [1:19, 1:20] 8 0 0 0 0 0 0 0 0 0 ...
# $ marr.j: num [1:19, 1:20] 1 0 0 0 0 0 0 0 0 0 ...
# $ J     : num [1:20] 11 12 17 9 12 30 17 11 10 20 ...
# $ B     : num [1:20] 9 10 9 4 10 13 11 11 8 10 ...
# $ count : num [1:20] 10 12 9 4 11 14 18 11 8 13 ...

# Bundle data
K <- 15                                           # Number of years with predictions
jags.data <- list(marr.j=woodchat10$marr.j, marr.a=woodchat10$marr.a,
    n.occasions=ncol(woodchat10$marr.j), rel.j=rowSums(woodchat10$marr.j),
    rel.a=rowSums(woodchat10$marr.a), J=woodchat10$J, B=woodchat10$B, count=woodchat10$count,
    pNinit=dUnif(1, 50), K=K)
str(jags.data)
# List of 10
# $ marr.j     : num [1:19, 1:20] 1 0 0 0 0 0 0 0 0 0 ...
# $ marr.a     : num [1:19, 1:20] 8 0 0 0 0 0 0 0 0 0 ...
# $ n.occasions: int 20
# $ rel.j      : num [1:19] 30 41 39 42 44 26 39 35 31 47 ...
# $ rel.a      : num [1:19] 18 16 14 20 21 21 20 17 19 18 ...
# $ J          : num [1:20] 11 12 17 9 12 30 17 11 10 20 ...
# $ B          : num [1:20] 9 10 9 4 10 13 11 11 8 10 ...
# $ count      : num [1:20] 10 12 9 4 11 14 18 11 8 13 ...
# $ pNinit     : num [1:50] 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 ...
# $ K          : num 15

# Write JAGS model file
cat(file="model1.txt", "
model {
  # Priors and linear models
  mean.logit.sj <- logit(mean.sj)
  mean.sj ~ dunif(0, 1)
  mean.logit.sa <- logit(mean.sa)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.log.f <- log(mean.f)
  mean.f ~ dunif(0, 10)

  for (t in 1:(n.occasions-1)){
    p[t] <- mean.p
  }

  for (t in 1:(n.occasions-1+K)){                 # Here we extend the loop to K more years
    logit.sj[t] <- mean.logit.sj + eps.sj[t]
    eps.sj[t] ~ dnorm(0, tau.sj)
    sj[t] <- ilogit(logit.sj[t])
    logit.sa[t] <- mean.logit.sa + eps.sa[t]
    eps.sa[t] ~ dnorm(0, tau.sa)
    sa[t] <- ilogit(logit.sa[t])
  }

  for (t in 1:(n.occasions+K)){                   # Extended loop also here
    log.f[t] <- mean.log.f + eps.f[t]
    eps.f[t] ~ dnorm(0, tau.f)
    f[t] <- exp(log.f[t])
  }

  sigma.sj ~ dunif(0, 10)
  tau.sj <- pow(sigma.sj, -2)
  sigma.sa ~ dunif(0, 10)
  tau.sa <- pow(sigma.sa, -2)
  sigma.f ~ dunif(0, 10)
  tau.f <- pow(sigma.f, -2)

  sigma ~ dunif(0.5, 50)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dcat(pNinit)
  N[2,1] ~ dcat(pNinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1+K)){                 # Note extended loop
    N[1,t+1] ~ dpois(sj[t] * f[t] * (N[1,t] + N[2,t]))
    N[2,t+1] ~ dbin(sa[t], (N[1,t] + N[2,t]))
  }

  # Observation model
  for (t in 1:n.occasions){
    count[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Productivity data (Poisson regression model)
  for (t in 1:n.occasions){
    J[t] ~ dpois(f[t] * B[t])
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
    q[t] <- 1-p[t]                                # Probability of non-recapture
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
  # Total population size
  for (t in 1:(n.occasions+K)){
    Ntot[t] <- N[1,t] + N[2,t]
  }
  # Check whether the population is extinct in the future
  for (t in 1:K){
    extinct[t] <- equals(Ntot[n.occasions+t], 0)
  }
}
")

# Initial values
inits <- function(){list(mean.sj=runif(1, 0, 0.5), mean.sa=runif(1, 0.4, 0.6),
    mean.f=runif(1, 1.3, 2))}

# Parameters monitored
parameters <- c("mean.sj", "sigma.sj", "mean.sa", "sigma.sa", "mean.p", "mean.f", "sigma.f", "sigma",
    "sj", "sa", "f", "N", "Ntot", "extinct")

# MCMC settings
ni <- 20000; nb <- 5000; nc <- 3; nt <- 3; na <- 1000

# Call JAGS (ART 1 min), check convergence and summarize posteriors
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out1)
print(out1, 3)
#                mean     sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# mean.sj       0.284  0.029   0.230   0.284   0.343    FALSE 1 1.000 15000
# sigma.sj      0.413  0.160   0.122   0.402   0.764    FALSE 1 1.002  1046
# mean.sa       0.542  0.027   0.490   0.542   0.594    FALSE 1 1.000 11801
# sigma.sa      0.168  0.136   0.006   0.136   0.501    FALSE 1 1.006   566
# mean.p        0.633  0.033   0.568   0.633   0.696    FALSE 1 1.000 15000
# mean.f        1.529  0.093   1.347   1.529   1.714    FALSE 1 1.001  3187
# sigma.f       0.119  0.077   0.009   0.110   0.291    FALSE 1 1.005   407
# sigma         2.342  0.847   0.928   2.253   4.262    FALSE 1 1.001 14429
# sj[1]         0.231  0.060   0.114   0.230   0.349    FALSE 1 1.000  8420
# [ ...  output truncated... ]
# sj[34]        0.292  0.092   0.131   0.284   0.510    FALSE 1 1.000 15000
# sa[1]         0.555  0.052   0.464   0.550   0.678    FALSE 1 1.001 14315
# [ ...  output truncated... ]
# sa[34]        0.542  0.058   0.418   0.543   0.662    FALSE 1 1.001 15000
# f[1]          1.491  0.186   1.101   1.497   1.869    FALSE 1 1.001 12053
# [ ...  output truncated... ]
# f[35]         1.544  0.250   1.083   1.529   2.104    FALSE 1 1.005  3272
# N[1,1]        5.624  3.198   1.000   5.000  12.000    FALSE 1 1.000 15000
# N[2,1]        5.604  3.188   1.000   5.000  12.000    FALSE 1 1.000 15000
# [ ...  output truncated... ]
# N[1,35]       7.003 16.365   0.000   3.000  38.000     TRUE 1 1.012  9918
# N[2,35]       7.872 15.044   0.000   4.000  41.000     TRUE 1 1.003  8881
# Ntot[1]      11.229  2.046   8.000  11.000  16.000    FALSE 1 1.000 15000
# [ ...  output truncated... ]
# Ntot[35]     14.875 30.282   0.000   7.000  78.000     TRUE 1 1.006  8747
# extinct[1]    0.000  0.008   0.000   0.000   0.000    FALSE 1 1.291 15000
# [ ...  output truncated... ]
# extinct[15]   0.229  0.421   0.000   0.000   1.000     TRUE 1 1.000 15000

# ~~~ Save output for use later ~~~
save(out1, file="DataFig10.2.Rdata")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ plots of the predicted population sizes and demographic rates ~~~~
# Figure 10.2

n.years <- length(out1$mean$Ntot)

op <- par(mfrow=c(2, 2), mar=c(2, 4, 2.5, 1), cex=1.1)
plot(0, 0, ylim=range(c(out1$q2.5$Ntot, out1$q97.5$Ntot)), xlim=c(0.5, n.years),
    ylab=expression('Total population size ('*italic(N)*')'), xlab=NA, las=1,
    col="black", type="l", axes=FALSE)
axis(2, las=1)
axis(1, at=seq(5, n.years, 5), labels=NA)
axis(1, at=1:n.years, labels=NA, tcl=-0.25)
polygon(x=c(1:n.years, n.years:1), y=c(out1$q2.5$Ntot, rev(out1$q97.5$Ntot)),
    col="gray90", border="gray90")
lines(out1$q50$Ntot, col="blue")
abline(v=20, lty=2)

plot(0, 0, ylim=range(c(out1$q2.5$f, out1$q97.5$f)), xlim=c(0.5, n.years),
    ylab=expression('Productivity ('*italic(f)*')'), xlab=NA, las=1,
    col="black", type="l", axes=FALSE)
axis(2, las=1)
axis(1, at=seq(5, n.years, 5), labels=NA)
axis(1, at=1:n.years, labels=NA, tcl=-0.25)
polygon(x=c(1:n.years, n.years:1), y=c(out1$q2.5$f, rev(out1$q97.5$f)),
    col="gray90", border="gray90")
lines(out1$mean$f, col="blue")
abline(v=20, lty=2)

par(mar=c(4.5, 4, 0, 1))
plot(0, 0, ylim=range(c(out1$q2.5$sj, out1$q97.5$sj)), xlim=c(0.5, n.years),
    ylab=expression('Juvenile survival ('*italic(s)[italic(j)]*')'), xlab="Year",
    las=1, col="black", type="n", axes=FALSE)
axis(2, las=1)
axis(1, at=seq(5, n.years, 5), labels=c(5, 10, 15, 20, 25, 30, NA))
axis(1, at=1:n.years, labels=NA, tcl=-0.25)
polygon(x=c(1:(n.years-1)+0.5, (n.years-1):1+0.5),
    y=c(out1$q2.5$sj, rev(out1$q97.5$sj)), col="gray90", border="gray90")
lines(x=(1:(n.years-1))+0.5, y=out1$mean$sj, col="blue")
abline(v=20, lty=2)

plot(0, 0, ylim=range(c(out1$q2.5$sa, out1$q97.5$sa)), xlim=c(0.5, n.years),
    ylab=expression('Adult survival ('*italic(s)[italic(a)]*')'), xlab="Year",
    las=1, col="black", type="l", axes=FALSE)
axis(2, las=1)
axis(1, at=seq(5, n.years, 5), labels=c(5, 10, 15, 20, 25, 30, NA))
axis(1, at=1:n.years, labels=NA, tcl=-0.25)
polygon(x=c(1:(n.years-1)+0.5, (n.years-1):1+0.5),
    y=c(out1$q2.5$sa, rev(out1$q97.5$sa)), col="gray90", border="gray90")
lines(x=(1:(n.years-1))+0.5, y=out1$mean$sa, col="blue")
abline(v=20, lty=2)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code for figure 10.3 ~~~~
plot(out1$mean$extinct, type="l", ylab="Extinction probability", lwd=3,
    xlab="Year of forecast", frame=FALSE, axes=FALSE)
axis(1, at=1:K, tck=-0.0125, labels=FALSE)
axis(1, at=c(1, 3, 5, 7, 9, 11, 13, 15), labels=c(1, 3, 5, 7, 9, 11, 13, 15), tck=-0.025)
axis(2, las=1)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

D <- c(1, 3, 5, 7)
T <- length(woodchat10$count)
library(scales)
color <- viridis_pal(option='E')(20)[c(18,13,7,1)]
plot(y=rep(0,K), x=1:K, type="n", ylim=c(0, 0.5), ylab="Quasi-extinction probability",
    xlab="Year of forecast", axes=FALSE)
axis(2, las=1)
axis(1, at=1:K, tck=-0.0125, labels=FALSE)
axis(1, at=c(1, 3, 5, 7, 9, 11, 13, 15), labels=c(1, 3, 5, 7, 9, 11, 13, 15), tck=-0.025)
for (i in 1:length(D)){
  qextinct <- out1$sims.list$Ntot[,(T+1):(T+K)] <= D[i]
  lines(apply(qextinct, 2, mean), lwd=2, col=color[i])
}
legend("topleft", legend=c("D = 7", "D = 5", "D = 3", "D = 1"), col=rev(color), lwd=2, bty="n")

# Probability of smaller future population size
round(mean(out1$sims.list$Ntot[,T] > out1$sims.list$Ntot[,T+K]), 2)
# [1] 0.65

# Probability of the same population size in the future
round(mean(out1$sims.list$Ntot[,T] == out1$sims.list$Ntot[,T+K]), 2)
# [1] 0.02

# Probability of a larger future population size
round(mean(out1$sims.list$Ntot[,T] < out1$sims.list$Ntot[,T+K]), 2)
# [1] 0.33

D <- 3
T <- length(woodchat10$count)
qextinct <- out1$sims.list$Ntot[,(T+1):(T+K)] <= D
ext <- qextinct[,K]                               # extinct by year K, TRUE/FALSE
time.to.extinction <- K + 1 - apply(qextinct[ext,], 1, sum)
ns <- sum(ext)

# Fig 10.5
a <- barplot(table(time.to.extinction) / ns, col="grey", ylab="Relative frequency",
    xlab="Conditional time to extinction (years)", axes=FALSE, names.arg=c(1, NA, 3, NA, 5, NA, 7,
        NA, 9, NA, 11, NA, 13, NA, 15), border=NA)
axis(1, at=a, labels=NA)
axis(2, las=1)


# 10.5.2 Comparison of different management options
# -------------------------------------------------

# Write JAGS model file
cat(file="model2.txt", "
model {
  # Priors and linear models
  mean.logit.sj <- logit(mean.sj)
  mean.sj ~ dunif(0, 1)
  mean.logit.sa <- logit(mean.sa)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.log.f <- log(mean.f)
  mean.f ~ dunif(0, 10)

  for (t in 1:(n.occasions-1)){
    p[t] <- mean.p
  }

  sigma.sj ~ dunif(0, 10)
  tau.sj <- pow(sigma.sj, -2)
  sigma.sa ~ dunif(0, 10)
  tau.sa <- pow(sigma.sa, -2)
  sigma.f ~ dunif(0, 10)
  tau.f <- pow(sigma.f, -2)

  sigma ~ dunif(0.5, 50)
  tau <- pow(sigma, -2)

  # Models for demographic rates
  # Control, no change (option 4)
  for (t in 1:(n.occasions-1+K)){
    logit.sj[t] <- mean.logit.sj + eps.sj[t]
    eps.sj[t] ~ dnorm(0, tau.sj)
    sj[t] <- ilogit(logit.sj[t])
    logit.sa[t,1] <- mean.logit.sa + eps.sa[t,1]
    eps.sa[t,1] ~ dnorm(0, tau.sa)
    sa[t,1] <- ilogit(logit.sa[t,1])
  }
  for (t in 1:(n.occasions+K)){
    log.f[t,1] <- mean.log.f + eps.f[t,1]
    eps.f[t,1] ~ dnorm(0, tau.f)
    f[t,1] <- exp(log.f[t,1])
  }

  # Option 1: increase of productivity
  # Past: identical to control
  for (t in 1:n.occasions){
    log.f[t,2] <- log.f[t,1]
    eps.f[t,2] <- eps.f[t,1]
    f[t,2] <- f[t,1]
  }
  # Future: increase mean productivity by 20%
  for (t in (n.occasions+1):(n.occasions+K)){
    log.f[t,2] <- mean.log.f + log(1.2) + eps.f[t,2]
    eps.f[t,2] ~ dnorm(0, tau.f)
    f[t,2] <- exp(log.f[t,2])
  }

  # Option 2: reduction of temporal variability in adult survival
  # Past: identical to control
  for (t in 1:(n.occasions-1)){
    logit.sa[t,2] <- logit.sa[t,1]
    eps.sa[t,2] <- eps.sa[t,1]
    sa[t,2] <- sa[t,1]
  }
  # Future: reduction of temporal variability by half
  for (t in n.occasions:(n.occasions-1+K)){
    logit.sa[t,2] <- mean.logit.sa + eps.sa[t,2]
    eps.sa[t,2] ~ dnorm(0, tau.sa*2)              # Temporal precision increased
    sa[t,2] <- ilogit(logit.sa[t,2])
  }

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  N[1,1,1] ~ dcat(pNinit)
  N[2,1,1] ~ dcat(pNinit)

  # Process model over time: our model of population dynamics
  # Control, no change (option 4)
  for (t in 1:(n.occasions-1+K)){
    N[1,t+1,1] ~ dpois(f[t,1] * sj[t] * (N[1,t,1] + N[2,t,1]))
    N[2,t+1,1] ~ dbin(sa[t,1], (N[1,t,1] + N[2,t,1]))
  }

  # Option 1: increase of productivity
  # Past
  for (t in 1:n.occasions){
    N[1,t,2] <- N[1,t,1]
    N[2,t,2] <- N[2,t,1]
  }
  # Future
  for (t in n.occasions:(n.occasions-1+K)){
    N[1,t+1,2] ~ dpois(f[t,2] * sj[t] * (N[1,t,2] + N[2,t,2]))
    N[2,t+1,2] ~ dbin(sa[t,1], (N[1,t,2] + N[2,t,2]))
  }

  # Option 2: decrease of temporal variability in adult survival
  # Past
  for (t in 1:n.occasions){
    N[1,t,3] <- N[1,t,1]
    N[2,t,3] <- N[2,t,1]
  }
  # Future
  for (t in n.occasions:(n.occasions-1+K)){
    N[1,t+1,3] ~ dpois(f[t,1] * sj[t] * (N[1,t,3] + N[2,t,3]))
    N[2,t+1,3] ~ dbin(sa[t,2], (N[1,t,3] + N[2,t,3]))
  }

  # Option 3: release of 3 females annually during 5 years
  # Past
  for (t in 1:n.occasions){
    N[1,t,4] <- N[1,t,1]
    N[2,t,4] <- N[2,t,1]
  }
  # Future, phase with releases
  for (t in n.occasions:(n.occasions+5)){
    N[1,t+1,4] ~ dpois(f[t,1] * sj[t] * (N[1,t,4] + N[2,t,4] + 3))
    N[2,t+1,4] ~ dbin(sa[t,1], (N[1,t,4] + N[2,t,4] + 3))
  }
  # Future, after the phase with releases
  for (t in (n.occasions+6):(n.occasions-1+K)){
    N[1,t+1,4] ~ dpois(f[t,1] * sj[t] * (N[1,t,4] + N[2,t,4]))
    N[2,t+1,4] ~ dbin(sa[t,1], (N[1,t,4] + N[2,t,4]))
  }

  # Observation model
  for (t in 1:n.occasions){
    count[t] ~ dnorm(N[1,t,1] + N[2,t,1], tau)
  }

  # Productivity data (Poisson regression model)
  for (t in 1:n.occasions){
    J[t] ~ dpois(f[t,1]*B[t])
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
    q[t] <- 1-p[t]                                # Probability of non-recapture
    pr.j[t,t] <- sj[t] * p[t]
    pr.a[t,t] <- sa[t,1] * p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- sj[t] * prod(sa[(t+1):j,1]) * prod(q[t:(j-1)]) * p[j]
      pr.a[t,j] <- prod(sa[t:j,1]) * prod(q[t:(j-1)]) * p[j]
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
  } #t

  # Derived parameters
  for (t in 1:(n.occasions+K)){
    Ntot[t,1] <- N[1,t,1] + N[2,t,1] # Total population sizes control
    Ntot[t,2] <- N[1,t,2] + N[2,t,2] # Total population sizes option 1
    Ntot[t,3] <- N[1,t,3] + N[2,t,3] # Total population sizes option 2
    Ntot[t,4] <- N[1,t,4] + N[2,t,4] # Total population sizes option 3
  }
}
")

# Initial values
inits <- function(){list(mean.sj=runif(1, 0, 0.5), mean.sa=runif(1, 0.4, 0.6),
    mean.f=runif(1, 1.3, 2))}

# Parameters monitored
parameters <- c("mean.sj", "sigma.sj", "mean.sa", "sigma.sa", "mean.p", "mean.f", "sigma.f", "sigma",
    "sj", "sa", "f", "N", "Ntot")

# MCMC settings
ni <- 20000; nb <- 5000; nc <- 3; nt <- 3; na <- 1000

# Call JAGS (ART 1 min), check convergence and summarize posteriors
out2 <- jags(jags.data, inits, parameters, "model2.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out2)
print(out2, 3)

# Fig. 10.6
D <- 0                                            # Define extinction threshold (here for absolute extinction)
T <- length(woodchat10$count)
extinct <- (out2$sims.list$Ntot[, (T+1):(T+K), ] <= D)
ext.prob <- apply(extinct, 2:3, mean)
color <- c("black", "skyblue", "seagreen", "orange")
op <- par(las=1)
matplot(ext.prob, type="l", ylab="Extinction probability", lwd=2, lty=1, xlab="Year of forecast",
    axes=FALSE, col=color)
axis(1, at=1:K, tck=-0.0125, labels=FALSE)
axis(1, at=c(1, 3, 5, 7, 9, 11, 13, 15), labels=c(1, 3, 5, 7, 9, 11, 13, 15), tck=-0.025)
axis(2, las=1)
legend("topleft", lty=rep(1, 4), lwd=rep(2, 4), col=color, legend=c("Do nothing (control)",
    "Increase productivity", "Reduce variability", "Translocation"), bty="n")
par(op)

last.year <- jags.data$n.occasions + jags.data$K

# Probability that option 1 (increased productivity) is better than control
round(mean(out2$sims.list$Ntot[,last.year,2] >
    out2$sims.list$Ntot[,last.year,1]), 2)
# [1] 0.77

# Probability that option 2 (reduced variability) is better than control
round(mean(out2$sims.list$Ntot[,last.year,3] >
    out2$sims.list$Ntot[,last.year,1]), 2)
# [1] 0.45

# Probability that option 3 (translocation) is better than control
round(mean(out2$sims.list$Ntot[,last.year,4] >
    out2$sims.list$Ntot[,last.year,1]), 2)
# [1] 0.80

# Probability that translocation is better than increased productivity
round(mean(out2$sims.list$Ntot[,last.year,4] >
    out2$sims.list$Ntot[,last.year,2]), 2)
# [1] 0.43


# Write JAGS model file
cat(file="model3.txt", "
model {
  # Priors and linear models
  mean.logit.sj <- logit(mean.sj)
  mean.sj ~ dunif(0, 1)
  mean.logit.sa <- logit(mean.sa)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.log.f <- log(mean.f)
  mean.f ~ dunif(0, 10)

  for (t in 1:(n.occasions-1)){
    p[t] <- mean.p
  }

  for (t in 1:(n.occasions-1+K)){                 # Extend the loop to K more years
    logit.sj[t] <- mean.logit.sj + eps.sj[t]
    eps.sj[t] ~ dnorm(0, tau.sj)
    sj[t] <- ilogit(logit.sj[t])
    logit.sa[t] <- mean.logit.sa + eps.sa[t]
    eps.sa[t] ~ dnorm(0, tau.sa)
    sa[t] <- ilogit(logit.sa[t])
  }

  # Productivity in the past
  for (t in 1:(n.occasions-1)){
    log.f[t] <- mean.log.f + eps.f[t]
    eps.f[t] ~ dnorm(0, tau.f)
    f[t] <- exp(log.f[t])
  }
  # Productivity in the future: increased by rep.inc
  for (t in n.occasions:(n.occasions+K)){
    log.f[t] <- mean.log.f + eps.f[t] + log(rep.inc)
    eps.f[t] ~ dnorm(0, tau.f)
    f[t] <- exp(log.f[t])
  }

  sigma.sj ~ dunif(0, 10)
  tau.sj <- pow(sigma.sj, -2)
  sigma.sa ~ dunif(0, 10)
  tau.sa <- pow(sigma.sa, -2)
  sigma.f ~ dunif(0, 10)
  tau.f <- pow(sigma.f, -2)

  sigma ~ dunif(0.5, 50)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dcat(pNinit)
  N[2,1] ~ dcat(pNinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1+K)){ # Extend the loop to K more years
    N[1,t+1] ~ dpois(sj[t] * f[t] * (N[1,t] + N[2,t]))
    N[2,t+1] ~ dbin(sa[t], (N[1,t] + N[2,t]))
  }

  # Observation model
  for (t in 1:n.occasions){
    count[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Productivity data (Poisson regression model)
  for (t in 1:n.occasions){
    J[t] ~ dpois(f[t] * B[t])
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
    q[t] <- 1-p[t]                                # Probability of non-recapture
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
  } #t

  # Derived parameters
  # Total population size
  for (t in 1:(n.occasions+K)){
    Ntot[t] <- N[1,t] + N[2,t]
  }
}
")

# Parameters monitored
parameters <- c("mean.sj", "sigma.sj", "mean.sa", "sigma.sa", "mean.p", "mean.f", "sigma.f", "sigma",
    "sj", "sa", "f", "N", "Ntot")

# MCMC settings
ni <- 20000; nb <- 5000; nc <- 3; nt <- 3; na <- 1000

# Define vector with the levels of increased productivity
rep.inc <- seq(1.4, 1.6, 0.025)

# Define matrix to store results
postN <- matrix(NA, nrow=(ni - nb) / nt * nc, ncol=length(rep.inc))

# Fit IPM for each level of increased productivity
for (i in 1:length(rep.inc)){                     # Takes about 10 min total

  cat(paste('*** Running model for', (rep.inc[i]-1)*100, '% productivity increase ***\n'))

  # Bundle data
  K <- 15                                         # Number of years with predictions
  jags.data <- list(marr.j=woodchat10$marr.j, marr.a=woodchat10$marr.a,
      n.occasions=ncol(woodchat10$marr.j), rel.j=rowSums(woodchat10$marr.j),
      rel.a=rowSums(woodchat10$marr.a), J=woodchat10$J, B=woodchat10$B, count=woodchat10$count,
      pNinit=dUnif(1, 50), K=K, rep.inc=rep.inc[i])

  # Initial values
  inits <- function(){list(mean.sj=runif(1, 0, 0.5), mean.sa=runif(1, 0.4, 0.6),
      mean.f=runif(1, 1.3, 2))}

  # Call JAGS (ART 1 min for each model)
  out3 <- jags(jags.data, inits, parameters, "model3.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
      n.thin=nt, parallel=TRUE)

  # Save posterior distribution of pop. size in 15 years
  postN[,i] <- out3$sims.list$Ntot[,20+K]
}

# Calculate probability that population size is >20
prob <- colMeans(postN > 20)

# Produce figure 10.7
plot(x=rep.inc, y=prob, type="b", xlab="Increase in productivity (%)",
    ylab=expression(Prob (N[35] > 20)), axes=FALSE, pch=16)
axis(2, las=1)
axis(1, at=rep.inc[c(1,3,5,7,9)], labels=c(40, 45, 50, 55, 60))
axis(1, at=rep.inc, labels=NA, tcl=-0.25)
segments(0, 0.9, 10, 0.9, col="red")
