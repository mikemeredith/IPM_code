# Schaub & Kéry (2022) Integrated Population Models
# Chapter 12 : Peregrine falcon
# -----------------------------

# Run time for test script 6 mins, full run 1 hr

# 12.4 Component data likelihoods
# =============================================

library(IPMbook); library(jagsUI)
data(peregrine)
str(peregrine)
# List of 3
# $ count       : num [1:43, 1:2] 1965 1966 1967 1968 1969 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : NULL
# .. ..$ : chr [1:2] "Year" "Breeding_pairs"
# $ productivity: num [1:43, 1:3] 1965 1966 1967 1968 1969 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : NULL
# .. ..$ : chr [1:3] "Year" "No_surveyed_brood" "No_fledglings"
# $ recoveries  : int [1:1810, 1:43] 0 0 0 0 0 0 0 0 0 0 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : NULL
# .. ..$ : chr [1:43] "1965" "1966" "1967" "1968" ...

# ~~~~ code for Fig. 12.2 ~~~~
m.dead <- marrayDead(peregrine$recoveries)
rel <- c(rowSums(m.dead),75)
library(scales)
co <- viridis_pal(option='E')(20)[2]
op <- par(mfrow=c(1, 3), las=1, mar=c(3,4,2,1))
plot(peregrine$count[,2], xlab=NA, ylab='Number', axes=FALSE, type='l',
    lty=1, lwd=3, col=co, ylim=c(0, 250))
mtext("Count of pairs", side=3, line=0.8)
axis(2)
axis(1, at=1:43, tcl=0, labels=NA)
u <- seq(1,43, by=10)
axis(1, at=u, tcl=-0.5, labels=(1965:2007)[u])
u <- seq(1,43, by=5)
axis(1, at=u, tcl=-0.25, labels=NA)

plot(peregrine$productivity[,3] / peregrine$productivity[,2], xlab=NA,
    ylab='Mean brood size', axes=FALSE, type='l', lty=1, lwd=3, col=co)
mtext("Productivity data", side=3, line=0.8)
axis(2)
axis(1, at=1:43, tcl=0, labels=NA)
u <- seq(1,43, by=10)
axis(1, at=u, tcl=-0.5, labels=(1965:2007)[u])
u <- seq(1,43, by=5)
axis(1, at=u, tcl=-0.25, labels=NA)

a <- barplot(rel, axes=FALSE, col=co, border=NA, ylab='Number',)
mtext("Nestlings ringed", side=3, line=0.8)
axis(2)
axis(1, at=a, tcl=0, labels=NA)
u <- seq(1,43, by=10)
axis(1, at=a[u], tcl=-0.5, labels=(1965:2007)[u])
u <- seq(1,43, by=5)
axis(1, at=a[u], tcl=-0.25, labels=NA)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 12.4.1 Population count data (no code)
# 12.4.2 Productivity data (no code)

# 12.4.3 Dead-recovery data
# -------------------------

# Create the dead-recovery m-array
m.dead <- marrayDead(peregrine$recoveries)

# Compute the annual number of marked individuals
rel <- rowSums(m.dead)


# 12.5 The integrated population model
# ====================================

# Bundle data and produce data overview
jags.data <- with(peregrine, list(nyears=ncol(peregrine$recoveries), marr=m.dead, rel=rel,
    y=count[,2], B=productivity[,2], J=productivity[,3], pNinit=dUnif(1, 70)))
str(jags.data)
# List of 7
# $ nyears: int 43
# $ marr  : num [1:42, 1:43] 0 0 0 0 0 0 0 0 0 0 ...
# $ rel   : num [1:42] 2 2 3 3 3 1 0 0 0 1 ...
# $ y     : num [1:43] 50 43 39 24 21 23 20 18 20 23 ...
# $ B     : num [1:43] 29 36 24 14 13 14 14 16 18 20 ...
# $ J     : num [1:43] 32 38 41 24 12 10 17 15 23 35 ...
# $ pNinit: num [1:70] 0.0143 0.0143 0.0143 0.0143 0.0143 ...

# IPM1: temporal random effects on the demographic rates
# ''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Write JAGS model file
cat(file = "model1.txt", "
model {
  # Priors and linear models
  for (t in 1:(nyears-1)){
    logit.s[1,t] ~ dnorm(l.mean.s[1], tau.s[1])
    s[1,t] <- ilogit(logit.s[1,t])
    logit.s[2,t] ~ dnorm(l.mean.s[2], tau.s[2])
    s[2,t] <- ilogit(logit.s[2,t])
    alpha[t] <- mean.alpha
    logit(r[t]) <- beta[1] + beta[2] * t          # Linear trend in recovery
  }
  for (u in 1:2){
    mean.s[u] ~ dunif(0, 1)                       # Priors for mean age-dep. survival
    l.mean.s[u] <- logit(mean.s[u])
    sigma.s[u] ~ dunif(0, 5)
    tau.s[u] <- pow(sigma.s[u], -2)
  }
  mean.alpha ~ dunif(0, 1)                        # Prior for prob. start reproduction at age 2y
  for (i in 1:2){                                 # Priors for betas
    beta[i] ~ dnorm(0, 0.001)
  }
  for (t in 1:nyears){
    log.rho[t] ~ dnorm(l.mean.rho, tau.rho)
    rho[t] <- exp(log.rho[t])
  }
  mean.rho ~ dunif(0, 5)                          # Prior for mean productivity
  l.mean.rho <- log(mean.rho)
  sigma.rho ~ dunif(0, 5)
  tau.rho <- pow(sigma.rho, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  for (a in 1:4){
    N[a,1] ~ dcat(pNinit)
  }

  # Process model over time: our model of population dynamics
  for (t in 1:(nyears-1)){
    N[1,t+1] ~ dpois(rho[t] / 2 * s[1,t] * (N[3,t] + N[4,t]))
    N[2,t+1] ~ dbin(s[2,t] * (1-alpha[t]), N[1,t])
    N[3,t+1] ~ dbin(s[2,t] * alpha[t], N[1,t])
    N[4,t+1] ~ dbin(s[2,t], (N[2,t] + N[3,t] + N[4,t]))
  }

  # Observation model
  for (t in 1:nyears){
    NB[t] <- N[3,t] + N[4,t]
    y[t] ~ dpois(NB[t])

    # GOF for population count data: mean absolute percentage error
    y.pred[t] ~ dpois(NB[t])
    disc.y[t] <- pow(((y[t] - NB[t]) / y[t]) * ((y[t] - NB[t]) / (y[t] +
        0.001)), 0.5)                             # Add a small number to avoid potential division by 0
    discN.y[t] <- pow(((y.pred[t] - NB[t]) / (y.pred[t] + 0.001)) *
        ((y.pred[t] - NB[t]) / (y.pred[t] + 0.001)), 0.5)
  }
  fit.y <- 100 / nyears * sum(disc.y)
  fitN.y <- 100 / nyears * sum(discN.y)

  # Dead-recovery data (multinomial model)
  # Define the multinomial likelihood
  for (t in 1:(nyears-1)){
    marr[t,1:nyears] ~ dmulti(pr[t,], rel[t])
  }
  # Define the cell probabilities of the m-array
  for (t in 1:(nyears-1)){
    # Main diagonal
    pr[t,t] <- (1-s[1,t]) * r[t]
    # Further than one above main diagonal
    for (j in (t+2):(nyears-1)){
      pr[t,j] <- s[1,t] * prod(s[2,(t+1):(j-1)]) * (1-s[2,j]) * r[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr[t,j] <- 0
    } #j
  } #t
  for (t in 1:(nyears-2)){
  # One above main diagonal
    pr[t,t+1] <- s[1,t] * (1-s[2,t+1]) * r[t+1]
  } #t
  # Last column: probability of non-recovery
  for (t in 1:(nyears-1)){
    pr[t,nyears] <- 1-sum(pr[t,1:(nyears-1)])
  } #t

  # GOF for dead-recovery data: Freeman-Tukey test statistics
  for (t in 1:(nyears-1)){
    # Simulated m-arrays
    marr.pred[t,1:nyears] ~ dmulti(pr[t,], rel[t])
    # Expected values and test statistics
    for (j in 1:nyears){
      marr.E[t,j] <- pr[t,j] * rel[t]
      E.org[t,j] <- pow((pow(marr[t,j], 0.5) - pow(marr.E[t,j], 0.5)), 2)
      E.new[t,j] <- pow((pow(marr.pred[t,j], 0.5) - pow(marr.E[t,j], 0.5)), 2)
    } #j
  } #t
  fit.DR <- sum(E.org)
  fitN.DR <- sum(E.new)

  # Productivity data (Poisson regression)
  for (t in 1:nyears){
    J[t] ~ dpois(B[t] * rho[t])

    # GOF for productivity data: deviance
    J.pred[t] ~ dpois(B[t] * rho[t])
    J.exp[t] <- B[t] * rho[t]
    dev[t] <- J[t] * log(J[t] / J.exp[t]) - (J[t] - J.exp[t])
    devN[t] <- J.pred[t] * log(J.pred[t] / J.exp[t]) - (J.pred[t] - J.exp[t])
  }
  fit.J <- sum(dev)
  fitN.J <- sum(devN)
}
")

# Initial values
inits <- function(){list(mean.s=runif(2, 0.6, 0.8))}

# Parameters monitored
parameters <- c("mean.s", "mean.alpha", "mean.rho", "beta", "sigma.rho", "sigma.s", "s", "rho", "r",
    "N", "NB", "fit.y", "fitN.y", "fit.DR", "fitN.DR", "fit.J", "fitN.J")

# MCMC settings
# ni <- 100000; nb <- 20000; nc <- 3; nt <- 80; na <- 5000
ni <- 10000; nb <- 2000; nc <- 3; nt <- 8; na <- 500  # ~~~ for testing, 2 mins

# Call JAGS from R (ART 56 min) and check convergence
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out1)


# IPM2: random-walk smoothers for all demographic rates
# '''''''''''''''''''''''''''''''''''''''''''''''''''''

# Write JAGS model file
cat(file = "model2.txt", "
model {
  # Priors and linear models
  s[1,1] ~ dunif(0, 1)                            # Prior of first-year survival in year 1
  logit.s[1,1] <- logit(s[1,1])
  s[2,1] ~ dunif(0, 1)                            # Prior of adult survival in year 1
  logit.s[2,1] <- logit(s[2,1])
  rho[1] ~ dunif(0, 5)                            # Prior of productivity in year 1
  log.rho[1] <- log(rho[1])

  for (t in 2:(nyears-1)){                        # Autoregressive models for s
    logit.s[1,t] ~ dnorm(logit.s[1,t-1], tau.s[1])
    s[1,t] <- ilogit(logit.s[1,t])
    logit.s[2,t] ~ dnorm(logit.s[2,t-1], tau.s[2])
    s[2,t] <- ilogit(logit.s[2,t])
  }
  for (u in 1:2){
    sigma.s[u] ~ dunif(0, 5)
    tau.s[u] <- pow(sigma.s[u], -2)
  }

  for (t in 2:nyears){                            # Autoregressive models for rho
    log.rho[t] ~ dnorm(log.rho[t-1], tau.rho)
    rho[t] <- exp(log.rho[t])
  }
  sigma.rho ~ dunif(0, 5)
  tau.rho <- pow(sigma.rho, -2)

  for (t in 1:(nyears-1)){                        # Models for alpha and r
    alpha[t] <- mean.alpha
    logit(r[t]) <- beta[1] + beta[2] * t
  }
  mean.alpha ~ dunif(0, 1)                        # Prior for prob. start reproduction at age 2y
  for (i in 1:2){
    beta[i] ~ dnorm(0, 0.001)
  }

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  for (a in 1:4){
    N[a,1] ~ dcat(pNinit)
  }

  # Process model over time: our model of population dynamics
  for (t in 1:(nyears-1)){
    N[1,t+1] ~ dpois(rho[t] / 2 * s[1,t] * (N[3,t] + N[4,t]))
    N[2,t+1] ~ dbin(s[2,t] * (1-alpha[t]), N[1,t])
    N[3,t+1] ~ dbin(s[2,t] * alpha[t], N[1,t])
    N[4,t+1] ~ dbin(s[2,t], (N[2,t] + N[3,t] + N[4,t]))
  }

  # Observation model
  for (t in 1:nyears){
    NB[t] <- N[3,t] + N[4,t]
    y[t] ~ dpois(NB[t])

    # GOF for population count data: mean absolute percentage error
    y.pred[t] ~ dpois(NB[t])
    disc.y[t] <- pow(((y[t] - NB[t]) / y[t]) * ((y[t] - NB[t]) / (y[t] +
        0.001)), 0.5)                             # Add a small number to avoid potential division by 0
    discN.y[t] <- pow(((y.pred[t] - NB[t]) / (y.pred[t] + 0.001)) *
        ((y.pred[t] - NB[t]) / (y.pred[t] + 0.001)), 0.5)
  }
  fit.y <- 100 / nyears * sum(disc.y)
  fitN.y <- 100 / nyears * sum(discN.y)

  # Dead-recovery data (multinomial model)
  # Define the multinomial likelihood
  for (t in 1:(nyears-1)){
    marr[t,1:nyears] ~ dmulti(pr[t,], rel[t])
  }
  # Define the cell probabilities of the m-array
  for (t in 1:(nyears-1)){
    # Main diagonal
    pr[t,t] <- (1-s[1,t]) * r[t]
    # Further than one above main diagonal
    for (j in (t+2):(nyears-1)){
      pr[t,j] <- s[1,t] * prod(s[2,(t+1):(j-1)]) * (1-s[2,j]) * r[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr[t,j] <- 0
    } #j
  } #t
  for (t in 1:(nyears-2)){
  # One above main diagonal
    pr[t,t+1] <- s[1,t] * (1-s[2,t+1]) * r[t+1]
  } #t
  # Last column: probability of non-recovery
  for (t in 1:(nyears-1)){
    pr[t,nyears] <- 1-sum(pr[t,1:(nyears-1)])
  } #t

  # GOF for dead-recovery data: Freeman-Tukey test statistics
  for (t in 1:(nyears-1)){
    # Simulated m-arrays
    marr.pred[t,1:nyears] ~ dmulti(pr[t,], rel[t])

    # Expected values and test statistics
    for (j in 1:nyears){
      marr.E[t,j] <- pr[t,j] * rel[t]
      E.org[t,j] <- pow((pow(marr[t,j], 0.5) - pow(marr.E[t,j], 0.5)), 2)
      E.new[t,j] <- pow((pow(marr.pred[t,j], 0.5) - pow(marr.E[t,j], 0.5)), 2)
    } #j
  } #t
  fit.DR <- sum(E.org)
  fitN.DR <- sum(E.new)

  # Productivity data (Poisson regression)
  for (t in 1:nyears){
    J[t] ~ dpois(B[t] * rho[t])

    # GOF for productivity data: deviance
    J.pred[t] ~ dpois(B[t] * rho[t])
    J.exp[t] <- B[t] * rho[t]
    dev[t] <- J[t] * log(J[t] / J.exp[t]) - (J[t] - J.exp[t])
    devN[t] <- J.pred[t] * log(J.pred[t] / J.exp[t]) - (J.pred[t] - J.exp[t])
  }
  fit.J <- sum(dev)
  fitN.J <- sum(devN)
}
")

# Initial values
s <- matrix(NA, nrow=2, ncol=jags.data$nyears-1)
s[1,1] <- runif(1, 0.5, 0.6)
s[2,1] <- runif(1, 0.7, 0.9)
inits <- function(){list(s=s)}

# Parameters monitored
parameters <-c("mean.alpha", "beta", "s", "rho", "r", "sigma.s", "sigma.rho", "N", "NB", "fit.y",
    "fitN.y", "fit.DR", "fitN.DR", "fit.J", "fitN.J")

# MCMC settings
# ni <- 200000; nb <- 50000; nc <- 3; nt <- 150; na <- 5000
ni <- 20000; nb <- 5000; nc <- 3; nt <- 15; na <- 500  # ~~~ for testing, 3 mins

# Call JAGS from R (ART 103 min) and check convergence
out2 <- jags(jags.data, inits, parameters, "model2.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out2)


# 12.6 Results
# ============

print(out1, 3)
#                mean     sd     2.5%      50%    97.5% overlap0     f  Rhat n.eff
# mean.s[1]     0.595  0.049    0.501    0.594    0.692    FALSE 1.000 1.001  1801
# mean.s[2]     0.780  0.016    0.748    0.780    0.810    FALSE 1.000 1.009   229
# mean.alpha    0.476  0.184    0.144    0.468    0.861    FALSE 1.000 1.056    41
# mean.rho      1.424  0.048    1.330    1.424    1.521    FALSE 1.000 1.000  3000
# beta[1]      -1.713  0.348   -2.413   -1.707   -1.049    FALSE 1.000 1.001  1581
# beta[2]      -0.019  0.011   -0.040   -0.019    0.004     TRUE 0.952 1.001  1690
# sigma.rho     0.191  0.028    0.143    0.188    0.250    FALSE 1.000 1.000  3000
# s[1,1]        0.479  0.173    0.124    0.496    0.788    FALSE 1.000 1.000  3000
# s[2,1]        0.762  0.048    0.637    0.770    0.834    FALSE 1.000 1.004   763
# [... output truncated ...]
# s[1,42]       0.672  0.141    0.385    0.674    0.929    FALSE 1.000 1.000  3000
# s[2,42]       0.774  0.031    0.705    0.775    0.832    FALSE 1.000 1.001  1505
# rho[1]        1.200  0.151    0.928    1.199    1.514    FALSE 1.000 1.000  3000
# rho[2]        1.157  0.141    0.899    1.151    1.447    FALSE 1.000 1.002   769
# [... output truncated ...]
# rho[42]       1.142  0.074    1.002    1.141    1.292    FALSE 1.000 1.001  2347
# rho[43]       1.308  0.076    1.165    1.307    1.462    FALSE 1.000 1.000  3000
# r[1]          0.155  0.044    0.082    0.151    0.252    FALSE 1.000 1.001  1842
# r[2]          0.153  0.042    0.083    0.149    0.245    FALSE 1.000 1.001  1832
# [... output truncated ...]
# r[41]         0.078  0.011    0.058    0.077    0.101    FALSE 1.000 1.002  2152
# r[42]         0.076  0.011    0.056    0.076    0.100    FALSE 1.000 1.002  2065
# N[1,1]        8.257  7.172    1.000    6.000   26.025    FALSE 1.000 1.001  3000
# N[2,1]        7.489  6.230    1.000    6.000   23.000    FALSE 1.000 1.004  1132
# [... output truncated ...]
# N[3,43]      24.739 12.885    4.000   23.000   53.000    FALSE 1.000 1.033    67
# N[4,43]     216.995 16.523  184.000  217.000  249.000    FALSE 1.000 1.024    88
# NB[1]        44.866  5.650   34.975   45.000   57.000    FALSE 1.000 1.000  3000
# NB[2]        40.791  4.452   33.000   41.000   50.000    FALSE 1.000 1.001  1331
# [... output truncated ...]
# NB[42]      242.664 10.881  222.000  243.000  263.025    FALSE 1.000 1.000  3000
# NB[43]      241.734 12.494  218.000  241.000  267.000    FALSE 1.000 1.000  3000
# fit.y         7.623  1.272    5.411    7.542   10.459    FALSE 1.000 1.001  1935
# fitN.y        9.571  1.475    7.031    9.419   12.804    FALSE 1.000 1.000  3000
# fit.DR       83.500  4.303   75.718   83.238   92.400    FALSE 1.000 1.000  3000
# fitN.DR      79.904  6.026   68.750   79.715   92.583    FALSE 1.000 1.002  1186
# fit.J        21.533  4.343   13.864   21.347   30.728    FALSE 1.000 1.000  2079
# fitN.J       21.365  4.645   13.245   21.051   31.076    FALSE 1.000 1.002  1062


print(out2, 3)
#                mean     sd     2.5%      50%    97.5% overlap0     f  Rhat n.eff
# mean.alpha    0.396  0.232    0.023    0.379    0.886    FALSE 1.000 1.003   568
# beta[1]      -1.683  0.354   -2.412   -1.675   -0.995    FALSE 1.000 1.000  3000
# beta[2]      -0.020  0.011   -0.041   -0.020    0.003     TRUE 0.956 1.000  3000
# s[1,1]        0.579  0.136    0.253    0.592    0.826    FALSE 1.000 1.005   401
# s[2,1]        0.615  0.091    0.436    0.616    0.784    FALSE 1.000 1.004   648
# [... output truncated ...]
# s[1,42]       0.587  0.087    0.416    0.585    0.783    FALSE 1.000 1.000  3000
# s[2,42]       0.782  0.033    0.714    0.782    0.846    FALSE 1.000 1.000  3000
# rho[1]        1.139  0.154    0.866    1.130    1.465    FALSE 1.000 1.000  3000
# rho[2]        1.167  0.130    0.929    1.161    1.430    FALSE 1.000 1.000  3000
# [... output truncated ...]
# rho[42]       1.122  0.067    0.994    1.122    1.251    FALSE 1.000 1.000  3000
# rho[43]       1.271  0.076    1.127    1.269    1.425    FALSE 1.000 1.000  2300
# r[1]          0.159  0.045    0.083    0.155    0.261    FALSE 1.000 1.000  3000
# r[2]          0.156  0.043    0.083    0.152    0.253    FALSE 1.000 1.000  3000
# [... output truncated ...]

# r[41]         0.077  0.011    0.058    0.077    0.099    FALSE 1.000 1.000  3000
# r[42]         0.076  0.011    0.056    0.075    0.099    FALSE 1.000 1.000  3000
# sigma.s[1]    0.164  0.154    0.006    0.117    0.581    FALSE 1.000 1.002  1579
# sigma.s[2]    0.177  0.080    0.025    0.173    0.354    FALSE 1.000 1.001  3000
# sigma.rho     0.174  0.031    0.120    0.172    0.238    FALSE 1.000 1.000  3000
# N[1,1]       23.868 18.416    1.000   19.000   65.000    FALSE 1.000 1.001  1499
# N[2,1]       17.664 14.066    1.000   14.000   53.000    FALSE 1.000 1.000  3000
# [... output truncated ...]

# N[3,43]      22.841 14.379    1.000   22.000   54.000    FALSE 1.000 1.002   727
# N[4,43]     218.163 18.226  182.000  218.000  254.000    FALSE 1.000 1.003   703
# NB[1]        49.149  6.600   37.000   49.000   63.000    FALSE 1.000 1.000  3000
# NB[2]        43.923  5.349   34.000   44.000   55.000    FALSE 1.000 1.001  1598
# [... output truncated ...]
# NB[42]      240.518 10.260  221.000  240.000  261.000    FALSE 1.000 1.001  2142
# NB[43]      241.003 12.953  216.000  241.000  267.000    FALSE 1.000 1.000  3000
# fit.y         6.866  1.041    5.059    6.757    9.172    FALSE 1.000 1.000  3000
# fitN.y        9.696  1.463    7.160    9.615   12.894    FALSE 1.000 1.000  3000
# fit.DR       88.268  3.838   81.488   88.097   96.127    FALSE 1.000 1.001  2423
# fitN.DR      81.585  6.074   70.321   81.355   93.988    FALSE 1.000 1.000  3000
# fit.J        22.889  4.837   14.527   22.505   33.610    FALSE 1.000 1.001  3000
# fitN.J       21.520  4.667   13.436   21.183   31.861    FALSE 1.000 1.002  1181


# ~~~~ Code to calculate the Bayesian p-values ~~~~
# IPM1
mean(out1$sims.list$fit.y < out1$sims.list$fitN.y)
mean(out1$sims.list$fit.DR < out1$sims.list$fitN.DR)
mean(out1$sims.list$fit.J < out1$sims.list$fitN.J)

# IPM2
mean(out2$sims.list$fit.y < out2$sims.list$fitN.y)
mean(out2$sims.list$fit.DR < out2$sims.list$fitN.DR)
mean(out2$sims.list$fit.J < out2$sims.list$fitN.J)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Fig. 12.4 ~~~~

library(scales)
cl <- c('forestgreen', 'black')
op <- par(mfrow=c(2, 2), mar=c(3, 4.5, 2, 0.5), las=1)

time <- seq(1965, 2007, by=5)
year <- 1965:2007
nyear <- length(year)

plot(x=(1:(nyear-1))+0.5, y=out1$mean$s[1,], type="n", pch=16, ylim=c(0, 1),
    ylab=expression(paste("Juvenile survival ( ", italic(s)[1], ")")),
    xlab=NA, axes=F)
polygon(x=c((1:(nyear-1))+0.5, rev((1:(nyear-1))+0.5)),
    y=c(out1$q2.5$s[1,], rev(out1$q97.5$s[1,])), col=alpha(cl[1], 0.3),
    border=NA)
lines(x=(1:(nyear-1))+0.5, y=out1$mean$s[1,], col='white')
axis(1, at=1:nyear, tcl=-0.25, label=NA)
axis(1, at=seq(1, 41, by=5), label=time)
axis(2)
segments((1:(nyear-1))+0.5, out2$q2.5$s[1,], (1:(nyear-1))+0.5,
    out2$q97.5$s[1,], col=cl[2])
points(x=(1:(nyear-1))+0.5, y=out2$mean$s[1,], type="b", pch=16, col=cl[2])

plot(x=(1:(nyear-1))+0.5, y=out1$mean$s[2,], type="n", pch=16, ylim=c(0, 1),
    ylab=expression(paste("Adult survival ( ", italic(s)[2], ")")),
    xlab=NA, axes=F)
polygon(x=c((1:(nyear-1))+0.5, rev((1:(nyear-1))+0.5)),
    y=c(out1$q2.5$s[2,], rev(out1$q97.5$s[2,])), col=alpha(cl[1], 0.3),
    border=NA)
lines(x=(1:(nyear-1))+0.5, y=out1$mean$s[2,], col='white')
axis(1, at=1:nyear, tcl=-0.25, label=NA)
axis(1, at=seq(1, 41, by=5), label=time)
axis(2)
segments((1:(nyear-1))+0.5, out2$q2.5$s[2,], (1:(nyear-1))+0.5,
    out2$q97.5$s[2,], col=cl[2])
points(x=(1:(nyear-1))+0.5, y=out2$mean$s[2,], type="b", pch=16, col =cl[2])
legend('bottomright', pch=c(NA, 16), lwd=c(1, NA), col=cl, bty='n',
    legend=c(expression(IPM[1]: 'Unstructured random'),
    expression(IPM[2]: 'Random walk')))

plot(x=(1:nyear)+0.5, y=out1$mean$rho, type="n", pch=16, ylim=c(0.7, 2.4),
    ylab=expression(paste("Productivity ( ", rho, ")")), xlab=NA, axes=F)
polygon(x=c((1:nyear)+0.5, rev((1:nyear)+0.5)),
    y=c(out1$q2.5$rho, rev(out1$q97.5$rho)), col=alpha(cl[1], 0.3),
    border=NA)
lines(x=(1:nyear)+0.5, y=out1$mean$rho, col='white')
axis(1, at=1:nyear, tcl=-0.25, label=NA)
axis(1, at=seq(1, 41, by=5), label=time)
axis(2)
axis(2, at=c(0.75, 1.25, 1.75, 2.25), tcl=-0.25, labels=NA)
segments((1:nyear)+0.5, out2$q2.5$rho, (1:nyear)+0.5,
    out2$q97.5$rho, col=cl[2])
points(x=(1:nyear)+0.5, y=out2$mean$rho, type="b", pch=16, col =cl[2])

plot(x=(1:(nyear-1))+0.5, y=out1$mean$r, type="n", pch=16, ylim=c(0, 0.4),
    ylab=expression(paste('Recovery probability ( ', italic(r),')')),
    xlab=NA, axes=F)
polygon(x=c((1:(nyear-1))+0.5, rev((1:(nyear-1))+0.5)),
    y=c(out1$q2.5$r, rev(out1$q97.5$r)), col=alpha(cl[1], 0.3), border=NA)
lines(x=(1:(nyear-1))+0.5, y=out1$mean$r, col='white')
axis(1, at=1:nyear, tcl=-0.25, label=NA)
axis(1, at=seq(1, 41, by=5), label=time)
axis(2)
segments((1:(nyear-1))+0.5, out2$q2.5$r, (1:(nyear-1))+0.5, out2$q97.5$r, col=cl[2])
points(x=(1:(nyear-1))+0.5, y=out2$mean$r, type="b", pch=16, col =cl[2], cex=0.9)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Fig. 12.5 ~~~~
library(scales)
cl <- c('forestgreen', 'black')
op <- par(las=1, mar=c(3,4,1,1))
year <- 1965:2007
nyear <- length(year)
time <- seq(1965, 2007, by=5)

plot(x=1:nyear, y=out1$mean$NB, type="n", pch=16, ylim=c(0, 270),
    ylab='Number', xlab=NA, axes=F)
polygon(x=c(1:nyear, rev(1:nyear)), y=c(out1$q2.5$NB, rev(out1$q97.5$NB)),
    col=alpha(cl[1], 0.3), border=NA)
lines(x=1:nyear, y=out1$mean$NB, col='white')
axis(1, at=1:nyear, tcl=-0.25, label=NA)
axis(1, at=seq(1, 43, by=5), label=time)
axis(2)
segments(1:nyear, out2$q2.5$NB, 1:nyear, out2$q97.5$NB, col=cl[2])
points(x=1:nyear, y=peregrine$count[,2], type='p', pch=16, col='red')
points(x=1:nyear, y=out2$mean$NB, type="p", pch=16, col =cl[2], cex=1.25)
legend(x=1.5, y=250, pch=c(16, NA, 16), lwd=c(NA, 1, NA),
    legend=c('Breeding counts', expression(IPM[1]: 'Unstructured random'),
    expression(IPM[2]: 'Random walk')), col=c('red',cl),
    bty='n', pt.cex=c(1, 1, 1.25))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Estimate age/stage composition of population averaged over both models
str(out1$sims.list$N)                             # Draws of N from model 1
str(out2$sims.list$N)                             # Draws of N from model 2
Ndraws <- array(NA, dim=c(6000, 4, 43))           # Create empty array
Ndraws[1:3000,,] <- out1$sims.list$N              # Fill in draws from IPM1
Ndraws[3001:6000,,] <- out2$sims.list$N           #... and from IPM2
stage.comp1 <- apply(Ndraws, c(2,3), mean)        # Compute posterior mean
dimnames(stage.comp1)<- list(c("N1 (Nonbr.)", "N2 (Nonbr.)", "N3 (Br.)", "N4 (Br)"), 1965:2007)

floaters <- apply(apply(Ndraws[,1:2,], c(1,3), sum), 2, mean)
breeders <- apply(apply(Ndraws[,3:4,], c(1,3), sum), 2, mean)
stage.comp2 <- rbind(floaters, breeders)
dimnames(stage.comp2)<- list(c("Floaters", "Breeders"), 1965:2007)
# Total population size
all <- apply(apply(Ndraws, c(1,3), sum), 2, mean)
# Number of floaters
prop.floaters <- (Ndraws[,1,] + Ndraws[,2,]) / apply(Ndraws, c(1,3), sum)
round(range(all))
# [1] 33 374

# ~~~~ Fig. 12.6 ~~~~
co1 <- viridis_pal(option='E')(20)[c(2,9,15,19)]
co2 <- viridis_pal(option='E')(20)[c(2,19)]

qu <- function(x) quantile(x, c(0.025, 0.975))
year <- 1965:2007
nyear <- length(year)
time <- seq(1965, 2007, by = 5)

op <- par(las=1, mfrow=c(3, 1), mar=c(4,4,2,1))

a <- barplot(stage.comp1[4:1,], axes=FALSE, col=co1, border=NA, ylab='Number')
mtext("'Raw' stage composition of population", side=3, line=0.2)
axis(2)
axis(1, at=a, tcl=-0.25, labels=NA)
legend(x=0.8, y=350, pch=rep(15,4), col=co1, bty='n',
    legend=c(expression(N[4]), expression(N[3]), expression(N[2]),
    expression(N[1])), pt.cex=1.5)

a <- barplot(stage.comp2[c(2,1),], axes=FALSE, col=co2, border=NA, ylab='Number')
mtext("Population composition in terms of breeders and floaters", side=3, line=0.2)
axis(2)
axis(1, at=a, tcl=-0.25, labels=NA)
legend(x=0.8, y=350, pch=rep(15,2), col=co2, bty='n',
    legend=c('Breeders', 'Floaters'), pt.cex=1.5)

plot(apply(prop.floaters, 2, mean), type='p', pch=16, axes=FALSE,
    ylab='Proportion of floaters', xlab=NA, ylim=c(0,0.7), cex=1.4)
segments(1:43, apply(prop.floaters, 2, qu)[1,], 1:43, apply(prop.floaters, 2, qu)[2,])
mtext("Proportion of floaters", side = 3, line=0.2)
axis(2)
axis(1, at=1:43, tcl=-0.25, labels=NA)
u <- seq(1,43, by=10)
axis(1, at=u, tcl=-0.5, labels=(1965:2007)[u])
u <- seq(1,43, by=5)
axis(1, at=u, tcl=-0.5, labels=NA)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
