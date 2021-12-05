# Schaub & Kéry (2022) Integrated Population Models
# Chapter 14 : Hoopoe
# -------------------

# Run time for test script 2 mins, full run 15 mins

# 14.4 Component data likelihoods
# =============================================

# Load hoopoe data and produce data overview
library(IPMbook); library(jagsUI)
data(hoopoe)
str(hoopoe)
# List of 5
# $ ch      : int [1:3844, 1:16] 0 0 0 0 0 0 0 0 0 0 ...
# $ age     : int [1:3844] 1 1 1 1 1 1 1 1 1 1 ...
# $ count   : num [1:16] 34 46 68 93 88 87 85 78 82 84 ...
# $ reproAgg:List of 4
# ..$ J1: num [1:16] 73 188 243 320 261 222 206 154 278 220 ...
# ..$ J2: num [1:16] 51 83 101 182 256 226 206 278 237 244 ...
# ..$ B1: num [1:16] 15 23 36 50 44 37 39 30 48 48 ...
# ..$ B2: num [1:16] 6 13 12 23 38 38 32 41 39 41 ...
# $ reproInd:List of 3
# ..$ f   : int [1:1092] 6 11 11 3 15 6 11 12 7 6 ...
# ..$ id  : num [1:1092] 483 486 530 531 532 533 534 535 536 538 ...
# ..$ year: int [1:1092] 2002 2002 2002 2002 2002 2002 2002 2002 2002 ...


# 14.4.1 Population count data (no code)

# 14.4.2 Capture-recapture data
# -----------------------------


# Produce m-array
marr <- marrayAge(hoopoe$ch, hoopoe$age, 2)
marr.j <- marr[,,1]
marr.a <- marr[,,2]


# 14.4.3 Productivity data
# ------------------------

# ~~~~ investigate label-switching ~~~~
# Bundle the data
jags.data <- with(hoopoe$reproInd, list(marr.j=marr.j, marr.a=marr.a,
    n.years=ncol(marr.j), rel.j=rowSums(marr.j), rel.a=rowSums(marr.a),
    nind=length(unique(id)), id=id, nrep=length(id), f=f, year=year-2001,
    C=hoopoe$count, pNinit=dUnif(1, 50)))

# This model has no identifiability constraints built in
# Write JAGS model file
cat(file="model2.txt", "
model {
  # Priors and linear models
  # For mean and variance parameters
  mean.logit.phij <- logit(mean.phij)
  mean.phij ~ dunif(0, 1)
  mean.logit.phia <- logit(mean.phia)
  mean.phia ~ dunif(0, 1)
  mean.log.omega <- log(mean.omega)
  mean.omega ~ dunif(0, 50)
  mean.pa ~ dunif(0, 1)
  mean.pj ~ dunif(0, 1)

  sigma.phij ~ dunif(0, 3)
  tau.phij <- pow(sigma.phij, -2)
  sigma.phia ~ dunif(0, 3)
  tau.phia <- pow(sigma.phia, -2)
  sigma.omega ~ dunif(0.001, 5)
  tau.omega <- pow(sigma.omega, -2)

  sigma.rep ~ dunif(0.001, 5)
  tau.rep <- pow(sigma.rep, -2)

  for (g in 1:2){
    alpha[g] ~ dnorm(0, 0.01)     # Not too low precision, for convergence
    sigma.rep.t[g] ~ dunif(0.001, 2)
    tau.rep.t[g] <- pow(sigma.rep.t[g], -2)
  }

  # Mixture parameter
  gamma ~ dunif(0, 1)

  # Residual (observation) error
  sigma ~ dunif(0.4, 20)          # lower bound >0
  tau <- pow(sigma, -2)

  # Linear models for demographic parameters
  # Survival and recapture
  for (t in 1:(n.years-1)){
    logit.phij[t] ~ dnorm(mean.logit.phij, tau.phij)
    phij[t] <- ilogit(logit.phij[t])
    logit.phia[t] ~ dnorm(mean.logit.phia, tau.phia)
    phia[t] <- ilogit(logit.phia[t])
    pa[t] <- mean.pa
    pj[t] <- mean.pj
  }

  # Reproduction parameters
  for (g in 1:2){                 # Loop over two groups
    for (t in 1:n.years){
      log(rep[g,t]) <- alpha[g] + eps[g,t]
      eps[g,t] ~ dnorm(0, tau.rep.t[g])
    } #t
  } #g

  # Immigration parameters
  for (t in 1:n.years){
    log.omega[t] ~ dnorm(mean.log.omega, tau.omega)
    omega[t] <- exp(log.omega[t])
  }

  # Population count data (state-space model)
  # Model for initial stage-spec. population sizes
  ini[1] ~ dcat(pNinit)               # Local recruits
  ini[2] ~ dcat(pNinit)               # Surviving adults
  ini[3] ~ dpois(omega[1])            # Immigrants
  R[1,1] ~ dbin(1-gamma, ini[1])
  R[2,1] <- ini[1]-R[1,1]
  S[1,1] ~ dbin(1-gamma, ini[2])
  S[2,1] <- ini[2]-S[1,1]
  I[1,1] ~ dbin(1-gamma, ini[3])
  I[2,1] <- ini[3]-I[1,1]

  # Process model over time: our model of population dynamics
  for (t in 2:n.years){
    J[1,t] ~ dpois(rep[1,t-1] / 2 * phij[t-1] * N[1,t-1]) # Local recruits produced by females of group 1 (L)
    J[2,t] ~ dpois(rep[2,t-1] / 2 * phij[t-1] * N[2,t-1]) # Local recruits produced by females of group 2 (H)
    R[1,t] ~ dbin(1-gamma, J[1,t] + J[2,t]) # Allocate recruits to group 1
    R[2,t] <- J[1,t] + J[2,t] - R[1,t]    # Allocate recruits to group 2 (H)
    S[1,t] ~ dbin(phia[t-1], N[1,t-1])    # Surviving adults, group 1 (L)
    S[2,t] ~ dbin(phia[t-1], N[2,t-1])    # Surviving adults, group 2 (H)
    Itot[t] ~ dpois(omega[t])             # Total number of immigrants
    I[1,t] ~ dbin(1-gamma, Itot[t])       # Immigrants of group 1 (L)
    I[2,t] <- Itot[t]-I[1,t]              # Immigrants of group 2 (H)
  }

  # Observation model
  for (t in 1:n.years){
    N[1,t] <- S[1,t] + R[1,t] + I[1,t]    # Total no of females in group 1
    N[2,t] <- S[2,t] + R[2,t] + I[2,t]    # Total no of females in group 2
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Capture-recapture data (CJS model with multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.years-1)){
    marr.j[t,1:n.years] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.years] ~ dmulti(pr.a[t,], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.years-1)){
    # Main diagonal
    qj[t] <- 1-pj[t]
    qa[t] <- 1-pa[t]
    pr.j[t,t] <- phij[t] * pj[t]
    pr.a[t,t] <- phia[t] * pa[t]
    # Above main diagonal
    for (j in (t+1):(n.years-1)){
      pr.j[t,j] <- phij[t] * prod(phia[(t+1):j]) * qj[t] * prod(qa[t:(j-1)]) * pa[j] / qa[t]
      pr.a[t,j] <- prod(phia[t:j]) * prod(qa[t:(j-1)]) * pa[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(n.years-1)){
    pr.j[t,n.years] <- 1-sum(pr.j[t,1:(n.years-1)])
    pr.a[t,n.years] <- 1-sum(pr.a[t,1:(n.years-1)])
  }

  # Productivity data (finite-mixture model with Gaussian error)
  # Mixture
  for (i in 1:nind){              # Loop over all ID of mothers with known annual prod.
    k[i] ~ dbern(gamma)
    h[i] <- k[i] + 1
  }

  # Productivity
  for (i in 1:nrep){              # Loop over all known-ID annual productivity data
    f[i] ~ dnorm(rep[h[id[i]], year[i]], tau.rep)
  }
}
")

# Initial values
inits <- function(){list(mean.omega=runif(1, 8, 12))}

# Parameters to be monitored
parameters <- c("mean.phij", "mean.phia", "alpha", "gamma", "mean.omega",
    "mean.pj", "mean.pa", "sigma.rep", "sigma.phij", "sigma.phia",
    "sigma.omega", "sigma.rep.t", "sigma", "phij", "phia", "rep",
    "R", "S", "I", "N")

# MCMC settings
# ni <- 40000; nb <- 10000; nc <- 3; nt <- 10; na <- 5000
ni <- 4000; nb <- 1000; nc <- 3; nt <- 1; na <- 5000  # ~~~ for testing

# Call JAGS (ART 8 min), check convergence and summarize posteriors
set.seed(2)                                       # ~~~ to get the example of switching
out2 <- jags(jags.data, inits, parameters, "model2.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

traceplot(out2, c('alpha', 'gamma'), layout=c(1,3))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 14.5 The integrated population model
# ====================================

# Bundle the data
jags.data <- with(hoopoe$reproInd, list(marr.j=marr.j, marr.a=marr.a, n.years=ncol(marr.j),
    rel.j=rowSums(marr.j), rel.a=rowSums(marr.a), nind=length(unique(id)), id=id, nrep=length(id),
    f=f, year=year-2001, C=hoopoe$count, pNinit=dUnif(1, 50)))
str(jags.data)                                    # Remind ourselves of how the data look like
# List of 12
# $ marr.j : num [1:15, 1:16] 6 0 0 0 0 0 0 0 0 0 ...
# $ marr.a : num [1:15, 1:16] 6 0 0 0 0 0 0 0 0 0 ...
# $ n.years: int 16
# $ rel.j  : num [1:15] 104 165 256 285 265 264 237 222 260 242 ...
# $ rel.a  : num [1:15] 29 36 57 77 82 83 78 73 86 90 ...
# $ nind   : int 767
# $ id     : num [1:1092] 483 486 530 531 532 533 534 535 536 538 ...
# $ nrep   : int 1092
# $ f      : int [1:1092] 6 11 11 3 15 6 11 12 7 6 ...
# $ year   : num [1:1092] 1 1 1 1 1 1 1 1 1 1 ...
# $ C      : num [1:16] 34 46 68 93 88 87 85 78 82 84 ...
# $ pNinit : num [1:50] 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 ...

# Write JAGS model file
cat(file="model1.txt", "
model {
  # Priors and linear models
  # For mean and variance parameters
  mean.logit.phij <- logit(mean.phij)
  mean.phij ~ dunif(0, 1)
  mean.logit.phia <- logit(mean.phia)
  mean.phia ~ dunif(0, 1)
  mean.log.omega <- log(mean.omega)
  mean.omega ~ dunif(0, 50)
  mean.pa ~ dunif(0, 1)
  mean.pj ~ dunif(0, 1)

  sigma.phij ~ dunif(0, 3)
  tau.phij <- pow(sigma.phij, -2)
  sigma.phia ~ dunif(0, 3)
  tau.phia <- pow(sigma.phia, -2)
  sigma.omega ~ dunif(0.001, 5)
  tau.omega <- pow(sigma.omega, -2)

  sigma.rep ~ dunif(0.001, 5)
  tau.rep <- pow(sigma.rep, -2)

  for (g in 1:2){
    sigma.rep.t[g] ~ dunif(0.001, 2)
    tau.rep.t[g] <- pow(sigma.rep.t[g], -2)
  }

  # Enforce order on the alphas: alpha[1] < alpha[2]
  alpha[1] ~ dnorm(0, 0.01)                       # Not too low precision, for convergence
  alpha[2] <- alpha[1] + diff.alpha
  diff.alpha ~ dnorm(0, 0.01)T(0,)                # Diff. must be positive; alpha[2] is
                                                  # defined to have higher values

  # Mixture parameter
  gamma ~ dunif(0, 0.5)

  # Residual (observation) error
  sigma ~ dunif(0.4, 20)                          # Lower bound >0
  tau <- pow(sigma, -2)

  # Linear models for demographic parameters
  # Survival and recapture
  for (t in 1:(n.years-1)){
    logit.phij[t] ~ dnorm(mean.logit.phij, tau.phij)
    phij[t] <- ilogit(logit.phij[t])
    logit.phia[t] ~ dnorm(mean.logit.phia, tau.phia)
    phia[t] <- ilogit(logit.phia[t])
    pa[t] <- mean.pa
    pj[t] <- mean.pj
  }

  # Reproduction parameters
  for (g in 1:2){ # Loop over two groups
    for (t in 1:n.years){
      log(rep[g,t]) <- alpha[g] + eps[g,t]
      eps[g,t] ~ dnorm(0, tau.rep.t[g])
    } #t
  } #g
  # Immigration parameters
  for (t in 1:n.years){
    log.omega[t] ~ dnorm(mean.log.omega, tau.omega)
    omega[t] <- exp(log.omega[t])
  }

  # Population count data (state-space model)
  # Model for initial stage-spec. population sizes
  ini[1] ~ dcat(pNinit)                           # Local recruits
  ini[2] ~ dcat(pNinit)                           # Surviving adults
  ini[3] ~ dpois(omega[1])                        # Immigrants
  R[1,1] ~ dbin(1-gamma, ini[1])
  R[2,1] <- ini[1]-R[1,1]
  S[1,1] ~ dbin(1-gamma, ini[2])
  S[2,1] <- ini[2]-S[1,1]
  I[1,1] ~ dbin(1-gamma, ini[3])
  I[2,1] <- ini[3]-I[1,1]

  # Process model over time: our model of population dynamics
  for (t in 2:n.years){
    J[1,t] ~ dpois(rep[1,t-1] / 2 * phij[t-1] * N[1,t-1]) # Local recruits
    # produced by females of group 1 (L)
    J[2,t] ~ dpois(rep[2,t-1] / 2 * phij[t-1] * N[2,t-1]) # Local recruits
    # produced by females of group 2 (H)
    R[1,t] ~ dbin(1-gamma, J[1,t] + J[2,t])               # Allocate recruits to group 1 (L)
    R[2,t] <- J[1,t] + J[2,t] - R[1,t]                    # Allocate recruits to group 2 (H)
    S[1,t] ~ dbin(phia[t-1], N[1,t-1])                    # Surviving adults, group 1 (L)
    S[2,t] ~ dbin(phia[t-1], N[2,t-1])                    # Surviving adults, group 2 (H)
    Itot[t] ~ dpois(omega[t])                             # Total number of immigrants
    I[1,t] ~ dbin(1-gamma, Itot[t])                       # Immigrants of group 1 (L)
    I[2,t] <- Itot[t]-I[1,t]                              # Immigrants of group 2 (H)
  }

  # Observation model
  for (t in 1:n.years){
    N[1,t] <- S[1,t] + R[1,t] + I[1,t]                    # Total no of females in group 1
    N[2,t] <- S[2,t] + R[2,t] + I[2,t]                    # Total no of females in group 2
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Capture-recapture data (CJS model with multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.years-1)){
    marr.j[t,1:n.years] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.years] ~ dmulti(pr.a[t,], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.years-1)){
    # Main diagonal
    qj[t] <- 1-pj[t]
    qa[t] <- 1-pa[t]
    pr.j[t,t] <- phij[t] * pj[t]
    pr.a[t,t] <- phia[t] * pa[t]
    # Above main diagonal
    for (j in (t+1):(n.years-1)){
      pr.j[t,j] <- phij[t] * prod(phia[(t+1):j]) * qj[t] * prod(qa[t:(j-1)]) * pa[j] / qa[t]
      pr.a[t,j] <- prod(phia[t:j]) * prod(qa[t:(j-1)]) * pa[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(n.years-1)){
    pr.j[t,n.years] <- 1-sum(pr.j[t,1:(n.years-1)])
    pr.a[t,n.years] <- 1-sum(pr.a[t,1:(n.years-1)])
  }

  # Productivity data (finite-mixture model with Gaussian error)
  # Mixture
  for (i in 1:nind){                              # Loop over all ID of mothers with known annual productivity
    k[i] ~ dbern(gamma)
    h[i] <- k[i] + 1
  }

  # Productivity
  for (i in 1:nrep){                              # Loop over all known-ID annual productivity data
    f[i] ~ dnorm(rep[h[id[i]], year[i]], tau.rep)
  }
}
")

# Initial values
inits <- function(){list(mean.omega=runif(1, 8, 12))}

# Parameters to be monitored
parameters <- c("mean.phij", "mean.phia", "alpha", "gamma", "diff.alpha", "mean.omega", "mean.pj",
    "mean.pa", "sigma.rep", "sigma.phij", "sigma.phia", "sigma.omega", "sigma.rep.t", "sigma",
    "phij", "phia", "rep", "R", "S", "I", "N")

# MCMC settings
# ni <- 40000; nb <- 10000; nc <- 10; nt <- 10; na <- 5000
ni <- 4000; nb <- 1000; nc <- 3; nt <- 1; na <- 500  # ~~~ for testing, < 1 min

# Call JAGS (ART 15 min), check convergence and summarize posteriors
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out1)

# ~~~~ Fig. 14.3 ~~~~
traceplot(out1, c('alpha', 'gamma'), layout=c(1,3))
# ~~~~~~~~~~~~~~~~~~~

# ~~~~ save the results ~~~~
save(out1, file ="Hoopoe.Results.Rdata")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

# 14.6 Results
# ============

print(out1, 3)
                   # mean     sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# mean.phij         0.122  0.012    0.099    0.122    0.147    FALSE 1 1.001  5725
# mean.phia         0.363  0.020    0.324    0.363    0.403    FALSE 1 1.001  7546
# alpha[1]          1.673  0.052    1.567    1.674    1.771    FALSE 1 1.001  5950
# alpha[2]          2.132  0.117    1.878    2.139    2.346    FALSE 1 1.002  3125
# gamma             0.154  0.055    0.051    0.152    0.269    FALSE 1 1.005  1315
# diff.alpha        0.459  0.125    0.179    0.470    0.679    FALSE 1 1.001  5717
# mean.omega       19.886  3.032   13.968   19.885   25.799    FALSE 1 1.002  2672
# mean.pj           0.603  0.039    0.526    0.604    0.680    FALSE 1 1.001  6456
# mean.pa           0.795  0.033    0.727    0.796    0.855    FALSE 1 1.000 14285
# sigma.rep         2.931  0.092    2.759    2.927    3.122    FALSE 1 1.001  6303
# sigma.phij        0.278  0.101    0.105    0.268    0.505    FALSE 1 1.000  9746
# sigma.phia        0.164  0.103    0.008    0.154    0.395    FALSE 1 1.002  3011
# sigma.omega       0.178  0.143    0.008    0.145    0.537    FALSE 1 1.003  1908
# sigma.rep.t[1]    0.149  0.045    0.072    0.144    0.251    FALSE 1 1.001  4207
# sigma.rep.t[2]    0.281  0.104    0.131    0.265    0.529    FALSE 1 1.002  7348
# sigma             2.807  2.163    0.430    2.271    8.149    FALSE 1 1.002  3828
# phij[1]           0.114  0.023    0.071    0.114    0.162    FALSE 1 1.001  6568
# [... output truncated ...]
# phij[15]          0.104  0.022    0.062    0.104    0.149    FALSE 1 1.001  7592
# phia[1]           0.354  0.041    0.267    0.356    0.435    FALSE 1 1.001  7006
# [... output truncated ...]
# phia[15]          0.394  0.045    0.325    0.385    0.503    FALSE 1 1.001  7654
# rep[1,1]          5.633  0.566    4.559    5.617    6.789    FALSE 1 1.001  4721
# rep[2,1]         10.137  1.531    7.092   10.145   13.176    FALSE 1 1.000 15369
# [... output truncated ...]
# rep[1,16]         5.613  0.487    4.688    5.604    6.591    FALSE 1 1.000 16855
# rep[2,16]         8.915  1.999    4.904    9.020   12.572    FALSE 1 1.001  6217
# R[1,1]            7.361  5.144    1.000    6.000   19.000    FALSE 1 1.003  1731
# R[2,1]            1.375  1.534    0.000    1.000    5.000     TRUE 1 1.001  5410
# [... output truncated ...]
# R[1,16]          10.498  3.704    4.000   10.000   18.000    FALSE 1 1.001  4311
# R[2,16]           1.888  1.561    0.000    2.000    6.000     TRUE 1 1.001  7532
# S[1,1]            7.452  5.259    1.000    6.000   20.000    FALSE 1 1.003  2333
# S[2,1]            1.390  1.551    0.000    1.000    5.000     TRUE 1 1.002  3545
# [... output truncated ...]
# S[2,15]           3.085  2.061    0.000    3.000    8.000     TRUE 1 1.002  3401
# S[1,16]          15.714  3.678    9.000   16.000   23.000    FALSE 1 1.001  4489
# S[2,16]           2.835  1.950    0.000    3.000    7.000     TRUE 1 1.001  4372
# I[1,1]           14.758  4.874    6.000   15.000   25.000    FALSE 1 1.001  4125
# I[2,1]            2.770  2.011    0.000    2.000    7.000     TRUE 1 1.001  4911
# [... output truncated ...]
# I[1,16]          14.995  4.250    7.000   15.000   24.000    FALSE 1 1.001  4080
# I[2,16]           2.738  1.948    0.000    2.000    7.000     TRUE 1 1.001  5716
# N[1,1]           29.571  4.194   22.000   29.000   39.000    FALSE 1 1.003  3959
# N[2,1]            5.535  2.950    1.000    5.000   12.000    FALSE 1 1.002  3206
# [... output truncated ...]
# N[1,16]          41.206  4.570   32.000   41.000   50.000    FALSE 1 1.002  3316
# N[2,16]           7.462  3.700    1.000    7.000   16.000    FALSE 1 1.002  3310

# ~~~~ code for Fig. 14.4 ~~~~
library(scales)
cl <- viridis_pal(option='E')(20)[c(18,2)]
n.years <- ncol(out1$mean$rep)
qu <- function(x) quantile(x, c(0.025, 0.975))
op <- par(mfrow=c(2,1), mar=c(1.5, 4, 2, 0), las=1)
d <- 0.1
plot(x=(1:n.years)+d, y=out1$mean$rep[1,], type="b", pch=16, ylim=c(3, 15),
    xlab="", ylab="Productivity", axes=FALSE, col=cl[1])
segments((1:n.years)+d, out1$q2.5$rep[1,], (1:n.years)+d, out1$q97.5$rep[1,], col=cl[1])
axis(2, las=1)
axis(1, at=1:n.years, tcl=-0.25, labels=NA)
axis(1, at=c(1, 4, 7, 10, 13, 16), labels=NA)
points(y=out1$mean$rep[2,], x=(1:n.years)-d, type="b", pch=17, col=cl[2])
segments((1:n.years)-d, out1$q2.5$rep[2,], (1:n.years)-d, out1$q97.5$rep[2,], col=cl[2])
legend("topright", pch=c(17, 16), col=rev(cl),
    legend=c("High quality group", "Low quality group"), bty="n")

par(mar=c(4, 4, 1.5, 0))
plot(x=(1:n.years)+d, y=out1$mean$N[1,], pch=16,
    ylim=range(c(out1$q2.5$N, out1$q97.5$N)), ylab="Population size",
    xlab="", axes=FALSE, col=cl[1], type="b")
segments((1:n.years)+d, out1$q2.5$N[1,], (1:n.years)+d, out1$q97.5$N[1,], col=cl[1])
points(x=(1:n.years)-d, y=out1$mean$N[2,], pch=17, col=cl[2], type="b")
segments((1:n.years)-d, out1$q2.5$N[2,], (1:n.years)-d, out1$q97.5$N[2,], col=cl[2])
axis(2, las=1)
axis(1, at=1:n.years, tcl=-0.25, labels=NA)
axis(1, at=c(1, 4, 7, 10, 13, 16),
    labels=c("2002", "2005", "2008", "2011", "2014", "2017"))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
