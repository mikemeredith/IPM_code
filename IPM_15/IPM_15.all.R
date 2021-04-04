# Schaub & Kery (2021) Integrated Population Modeling
# Chapter 15 : Black grouse
# -------------------------
# Code from MS submitted to publisher.

# Run time 3 mins

library(IPMbook) ; library(jagsUI)

# 15.4 Component data likelihoods
# =============================================

data(grouse)
str(grouse)

# List of 10
 # $ ch       : int [1:96, 1:128] 1 1 NA NA NA NA NA NA NA NA ...
  # ..- attr(*, "dimnames")=List of 2
  # .. ..$ : NULL
  # .. ..$ : chr [1:128] "Apr98" "May98" "June98" "July98" ...
 # $ age      : num [1:96, 1:127] 2 2 NA NA NA NA NA NA NA NA ...
  # ..- attr(*, "dimnames")=List of 2
  # .. ..$ : NULL
  # .. ..$ : chr [1:127] "Apr98" "May98" "June98" "July98" ...
 # $ sex      : int [1:96] 2 2 2 2 2 2 2 2 2 2 ...
 # $ season   : int [1:128] 1 1 2 2 2 3 3 4 4 4 ...
 # $ count.sp : int [1:20] 36 45 52 64 46 49 50 35 36 39 ...
 # $ count.lsM: int [1:20] 15 11 14 14 14 5 16 19 19 16 ...
 # $ count.lsF: int [1:20] 23 32 29 20 19 18 28 21 31 39 ...
 # $ count.lsC: int [1:20] 44 86 76 28 18 24 99 53 101 82 ...
 # $ v        : int [1:20] 30 39 21 10 11 6 67 32 65 55 ...
 # $ u        : int [1:20] 18 22 11 5 4 2 25 15 29 26 ...

# 15.4.1 Population count data (no code)
# 15.4.2 Radio tracking data (no code)
# 15.4.3 Modeling productivity and chick sex ratio (no code)


# 15.5 The integrated population model
# ====================================

# Bundle data
getLast <- function(x) max(which(!is.na(x))) # Last occasion function

jags.data <- with(grouse, list(ch=ch, age=age, sex=sex, season=season,
    f=getFirst(ch), k=apply(ch, 1, getLast), C.sp=count.sp, C.lsM=count.lsM,
    C.lsF=count.lsF, C.lsC=count.lsC, u=u, v=v, nind=nrow(ch),
    ny=length(count.sp), pinit=dUnif(1, 100)))

str(jags.data)

# List of 15
 # $ ch    : int [1:96, 1:128] 1 1 NA NA NA NA NA NA NA NA ...
  # ..- attr(*, "dimnames")=List of 2
  # .. ..$ : NULL
  # .. ..$ : chr [1:128] "Apr98" "May98" "June98" "July98" ...
 # $ age   : num [1:96, 1:127] 2 2 NA NA NA NA NA NA NA NA ...
  # ..- attr(*, "dimnames")=List of 2
  # .. ..$ : NULL
  # .. ..$ : chr [1:127] "Apr98" "May98" "June98" "July98" ...
 # $ sex   : int [1:96] 2 2 2 2 2 2 2 2 2 2 ...
 # $ season: int [1:128] 1 1 2 2 2 3 3 4 4 4 ...
 # $ f     : int [1:96] 1 1 2 2 2 2 14 14 14 14 ...
 # $ k     : int [1:96] 8 15 9 4 6 27 25 18 23 25 ...
 # $ C.sp  : int [1:20] 36 45 52 64 46 49 50 35 36 39 ...
 # $ C.lsM : int [1:20] 15 11 14 14 14 5 16 19 19 16 ...
 # $ C.lsF : int [1:20] 23 32 29 20 19 18 28 21 31 39 ...
 # $ C.lsC : int [1:20] 44 86 76 28 18 24 99 53 101 82 ...
 # $ u     : int [1:20] 18 22 11 5 4 2 25 15 29 26 ...
 # $ v     : int [1:20] 30 39 21 10 11 6 67 32 65 55 ...
 # $ nind  : int 96
 # $ ny    : int 20
 # $ pinit : num [1:100] 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 ...


# Write JAGS model file
cat(file="model1.txt", "
model {
  # Priors and linear models
  # Survival
  for (a in 1:2){                     # age
    for (j in 1:2){                   # sex
      for (m in 1:4){                 # season
        s[a,j,m] ~ dbeta(1, 1)
      } #m
    } #j
  } #a

  # Productivity
  for (t in 1:ny){
    log.rho[t] ~ dnorm(log.mean.rho, tau.rho)
    rho[t] <- exp(log.rho[t])
  }
  mean.rho ~ dunif(0, 5)
  log.mean.rho <- log(mean.rho)
  sigma.rho ~ dunif(0, 1)
  tau.rho <- pow(sigma.rho, -2)

  # Sex ratio of the 5 weeks old chicks
  gamma ~ dbeta(1, 1)

  # Resighting probabilities for the late summer counts
  for (t in 1:ny){
    logit.p[t] ~ dnorm(lmean.p, tau.p)
    p[t] <- ilogit(logit.p[t])
  }
  mean.p ~ dbeta(1, 1)
  lmean.p <- logit(mean.p)
  sigma.p ~ dunif(0, 5)
  tau.p <- pow(sigma.p, -2)

  # Population count data (state-space model)
  # Model for the initial population size in spring: uniform priors
  Nsp[1,1,1] ~ dcat(pinit[])
  Nsp[2,1,1] ~ dcat(pinit[])
  Nsp[1,2,1] ~ dcat(pinit[])
  Nsp[2,2,1] ~ dcat(pinit[])

  # Process model: population sizes in spring
  for (t in 1:(ny-1)){
    # Females
    Nsp[1,1,t+1] ~ dbin(s[1,1,2]^0.5 * s[1,1,3]^2 * s[1,1,4]^5 * s[2,1,1], Nls[1,1,t])  # 1y
    Nsp[2,1,t+1] ~ dbin(s[2,1,2]^0.5 * s[2,1,3]^2 * s[2,1,4]^5 * s[2,1,1], Nls[2,1,t])  # >1y

    # Males
    Nsp[1,2,t+1] ~ dbin(s[1,2,2]^0.5 * s[1,2,3]^2 * s[1,2,4]^5 * s[2,2,1], Nls[1,2,t])  # 1y
    Nsp[2,2,t+1] ~ dbin(s[2,2,2]^0.5 * s[2,2,3]^2 * s[2,2,4]^5 * s[2,2,1], Nls[2,2,t])  # >1y
  }

  # Process model: population sizes in late summer
  for (t in 1:ny){
    # Total number of chicks
    F[t] ~ dpois((Nsp[1,2,t] + Nsp[2,2,t]) * s[2,2,1] * rho[t])
    # Allocate chicks to a sex
    Nls[1,1,t] ~ dbin(gamma, F[t])            # Female chicks
    Nls[1,2,t] <- F[t] - Nls[1,1,t]           # Male chicks
    # Survival
    Nls[2,1,t] ~ dbin(s[2,1,1] * s[2,1,2]^2.5, (Nsp[1,1,t] + Nsp[2,1,t]))    # ≥1y females
    Nls[2,2,t] ~ dbin(s[2,2,1] * s[2,2,2]^2.5, (Nsp[1,2,t] + Nsp[2,2,t]))    # ≥1y males
  }

  # Observation models
  for (t in 1:ny){
    # Spring counts of 1y and >1y males
    C.sp[t] ~ dpois((Nsp[1,2,t] + Nsp[2,2,t]))
    # Late summer counts of ≥1y females and of ≥1y males
    C.lsF[t] ~ dbin(p[t], Nls[2,1,t])
    C.lsM[t] ~ dbin(p[t], Nls[2,2,t])
  }

  # Radio tracking data (Bernoulli model)
  for (i in 1:nind){
    for (t in (f[i]+1):k[i]){
      ch[i,t] ~ dbern(s[age[i,t-1], sex[i], season[t-1]])
    } #t
  } #i

  # Productivity data (Poisson model)
  for (t in 1:ny){
    C.lsC[t] ~ dpois(C.lsF[t] * rho[t])
  }

  # Chick sex ratio data (binomial model)
  for (t in 1:ny){
    u[t] ~ dbin(gamma, v[t])
  }

  # Derived quantities: annual survival
  ann.s[1,1] <- s[1,1,2] * s[1,1,3]^2 * s[1,1,4]^5 * s[2,1,1]      # juv females
  ann.s[2,1] <- s[2,1,1]^2 * s[2,1,2]^3 * s[2,1,3]^2 * s[2,1,4]^5  # ad females
  ann.s[1,2] <- s[1,2,2] * s[1,2,3]^2 * s[1,2,4]^5 * s[2,2,1]      # juv males
  ann.s[2,2] <- s[2,2,1]^2 * s[2,2,2]^3 * s[2,2,3]^2 * s[2,2,4]^5  # ad males
}
")

# Initial values
inits <- function(){
  Nsp <- array(NA, dim=c(2, 2, 20))
  Nsp[,,1] <- round(runif(4, 30, 40))
  s <- array(runif(2*2*4, 0.90, 0.97), dim=c(2, 2, 4))
  list(Nsp=Nsp, s=s)
}

# Parameters monitored
parameters <- c("s", "ann.s", "mean.rho", "gamma", "mean.p", "sigma.rho",
    "sigma.p", "rho", "p", "Nsp", "Nls")

# MCMC settings
ni <- 30000; nb <- 10000; nc <- 3; nt <- 5; na <- 5000

# Call JAGS (ART 5 min), check convergence and summarize posteriors
out1 <- jags(jags.data, inits, parameters, "model1.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)
op <- par(mfrow=c(3, 3)); traceplot(out1)
par(op)

# ~~~~ model 1 plus posterior predictive checks ~~~~
# Write JAGS model file
cat(file="model2.txt", "
model {
  # Priors and linear models
  # Survival
  for (a in 1:2){                          # age
    for (j in 1:2){                       # sex
      for (m in 1:4){                    # season
        s[a,j,m] ~ dbeta(1, 1)
      } #m
    } #j
  } #a

  # Productivity
  for (t in 1:ny){
    log.rho[t] ~ dnorm(log.mean.rho, tau.rho)
    rho[t] <- exp(log.rho[t])
  }
  mean.rho ~ dunif(0, 5)
  log.mean.rho <- log(mean.rho)
  sigma.rho ~ dunif(0, 1)
  tau.rho <- pow(sigma.rho, -2)

  # Sex ratio of the 5 weeks old chicks
  gamma ~ dbeta(1, 1)

  # Resighting probabilities for the late summer counts
  for (t in 1:ny){
    logit.p[t] ~ dnorm(lmean.p, tau.p)
    p[t] <- ilogit(logit.p[t])
  }
  mean.p ~ dbeta(1, 1)
  lmean.p <- logit(mean.p)
  sigma.p ~ dunif(0, 5)
  tau.p <- pow(sigma.p, -2)

  # Population count data (state-space model)
  # Model for the initial population size in spring: uniform priors
  Nsp[1,1,1] ~ dcat(pinit[])
  Nsp[2,1,1] ~ dcat(pinit[])
  Nsp[1,2,1] ~ dcat(pinit[])
  Nsp[2,2,1] ~ dcat(pinit[])

  # Process model: population sizes in spring
  for (t in 1:(ny-1)){
    # Females
    Nsp[1,1,t+1] ~ dbin(s[1,1,2]^0.5 * s[1,1,3]^2 * s[1,1,4]^5 * s[2,1,1], Nls[1,1,t])  # 1y
    Nsp[2,1,t+1] ~ dbin(s[2,1,2]^0.5 * s[2,1,3]^2 * s[2,1,4]^5 * s[2,1,1], Nls[2,1,t])  # >1y

    # Males
    Nsp[1,2,t+1] ~ dbin(s[1,2,2]^0.5 * s[1,2,3]^2 * s[1,2,4]^5 * s[2,2,1], Nls[1,2,t])  # 1y
    Nsp[2,2,t+1] ~ dbin(s[2,2,2]^0.5 * s[2,2,3]^2 * s[2,2,4]^5 * s[2,2,1], Nls[2,2,t])  # >1y
  }

  # Process model: population sizes in late summer
  for (t in 1:ny){
    # Total number of chicks
    F[t] ~ dpois((Nsp[1,2,t] + Nsp[2,2,t]) * s[2,2,1] * rho[t])
    # Allocate chicks to a sex
    Nls[1,1,t] ~ dbin(gamma, F[t])            # Female chicks
    Nls[1,2,t] <- F[t] - Nls[1,1,t]           # Male chicks
    # Survival
    Nls[2,1,t] ~ dbin(s[2,1,1] * s[2,1,2]^2.5, (Nsp[1,1,t] + Nsp[2,1,t]))    # ≥1y females
    Nls[2,2,t] ~ dbin(s[2,2,1] * s[2,2,2]^2.5, (Nsp[1,2,t] + Nsp[2,2,t]))    # ≥1y males
  }

  # Observation models
  for (t in 1:ny){
    # Spring counts of 1y and >1y males
    C.sp[t] ~ dpois((Nsp[1,2,t] + Nsp[2,2,t]))
    # Late summer counts of ≥1y females and of ≥1y males
    C.lsF[t] ~ dbin(p[t], Nls[2,1,t])
    C.lsM[t] ~ dbin(p[t], Nls[2,2,t])
    # Elements for GOF
    C.sp.pred[t] ~ dpois((Nsp[1,2,t] + Nsp[2,2,t]))
    C.ls.pred[1,t] ~ dbin(p[t], Nls[2,1,t])
    C.ls.pred[2,t] ~ dbin(p[t], Nls[2,2,t])
    # Discrepancy measures: absolute percentage error
    d.Csp.org[t] <- (C.sp[t] - (Nsp[1,2,t] + Nsp[2,2,t])) / (C.sp[t] + 0.001)
    d.Cls.org[1,t] <- (C.lsF[t] - p[t] * Nls[2,1,t]) / (C.lsF[t] + 0.001)
    d.Cls.org[2,t] <- (C.lsM[t] - p[t] * Nls[2,2,t]) / (C.lsM[t] + 0.001)
    d.Csp.new[t] <- (C.sp.pred[t] - (Nsp[1,2,t] + Nsp[2,2,t])) / (C.sp.pred[t] + 0.001)
    d.Cls.new[1,t] <- (C.ls.pred[1,t] - p[t] * Nls[2,1,t]) / (C.ls.pred[1,t] + 0.001)
    d.Cls.new[2,t] <- (C.ls.pred[2,t] - p[t] * Nls[2,2,t]) / (C.ls.pred[2,t] + 0.001)
  }
  fit1[1] <- 100 / ny * sum(d.Csp.org)
  fit1[2] <- 100 / ny * sum(d.Csp.new)
  fit2[1] <- 100 / ny * sum(d.Cls.org[1,])
  fit2[2] <- 100 / ny * sum(d.Cls.new[1,])
  fit3[1] <- 100 / ny * sum(d.Cls.org[2,])
  fit3[2] <- 100 / ny * sum(d.Cls.new[2,])

  # Radio tracking data (Bernoulli model)
  for (i in 1:nind){
    for (t in (f[i]+1):k[i]){
      ch[i,t] ~ dbern(s[age[i,t-1], sex[i], season[t-1]])
    } #t
  } #i

  # Productivity data (Poisson model)
  for (t in 1:ny){
    C.lsC[t] ~ dpois(C.lsF[t] * rho[t])

    # Elements for GOF
    C.lsC.pred[t] ~ dpois(C.lsF[t] * rho[t])
    C.lsC.exp[t] <- C.lsF[t] * rho[t]
    # Discrepancy measure: deviance
    dev.org[t] <- C.lsC[t] * log(C.lsC[t] / C.lsC.exp[t])
    dev.new[t] <- C.lsC.pred[t] * log(C.lsC.pred[t] / C.lsC.exp[t])
  }
  fit4[1] <- 2*sum(dev.org)
  fit4[2] <- 2*sum(dev.new)

  # Chick sex ratio data (binomial model)
  for (t in 1:ny){
    u[t] ~ dbin(gamma, v[t])

    # Elements for GOF
    u.pred[t] ~ dbin(gamma, v[t])
    u.exp[t] <- gamma * v[t]
    # Discrepancy measure: Freeman-Tukey statistics
    Chi2.org[t] <- pow((pow(u[t], 0.5) - pow(u.exp[t], 0.5)), 2)
    Chi2.new[t] <- pow((pow(u.pred[t], 0.5) - pow(u.exp[t], 0.5)), 2)
  }
  fit5[1] <- sum(Chi2.org)
  fit5[2] <- sum(Chi2.new)

  # Derived quantities: annual survival
  ann.s[1,1] <- s[1,1,2] * s[1,1,3]^2 * s[1,1,4]^5 * s[2,1,1]      # juv females
  ann.s[2,1] <- s[2,1,1]^2 * s[2,1,2]^3 * s[2,1,3]^2 * s[2,1,4]^5  # ad females
  ann.s[1,2] <- s[1,2,2] * s[1,2,3]^2 * s[1,2,4]^5 * s[2,2,1]      # juv males
  ann.s[2,2] <- s[2,2,1]^2 * s[2,2,2]^3 * s[2,2,3]^2 * s[2,2,4]^5  # ad males
}
")


# Initial values
inits <- function(){
  Nsp <- array(NA, dim=c(2, 2, 20))
  Nsp[,,1] <- round(runif(4, 30, 40))
  s <- array(runif(2*2*4, 0.90, 0.97), dim=c(2, 2, 4))
  list(Nsp=Nsp, s=s)
}

# Parameters monitored
parameters <- c("s", "ann.s", "mean.rho", "gamma", "mean.p", "sigma.rho",
    "sigma.p", "rho", "p", "Nsp", "Nls", "fit1", "fit2", "fit3", "fit4", "fit5")

# MCMC settings
ni <- 30000; nb <- 10000; nc <- 3; nt <- 5; na <- 5000

# Call JAGS (ART 5 min), check convergence and summarize posteriors
out2 <- jags(jags.data, inits, parameters, "model2.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

save(out1, out2, file="Grouse.Results.Rdata")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 15.6 Results
# ============

print(out1, 3)
                # mean     sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# s[1,1,1]       0.499  0.288    0.026    0.497    0.975    FALSE 1 1.000  9165
# s[2,1,1]       0.949  0.016    0.914    0.950    0.975    FALSE 1 1.001  2959
# s[1,2,1]       0.498  0.288    0.022    0.499    0.974    FALSE 1 1.000  7815
# s[2,2,1]       0.898  0.029    0.837    0.900    0.951    FALSE 1 1.005   506
# s[1,1,2]       0.940  0.055    0.792    0.955    0.998    FALSE 1 1.000  8495
# s[2,1,2]       0.954  0.012    0.929    0.954    0.974    FALSE 1 1.001  4302
# s[1,2,2]       0.860  0.087    0.653    0.875    0.981    FALSE 1 1.001  2753
# s[2,2,2]       0.969  0.017    0.928    0.972    0.994    FALSE 1 1.007   315
# s[1,1,3]       0.908  0.044    0.811    0.914    0.978    FALSE 1 1.000  7056
# s[2,1,3]       0.932  0.020    0.889    0.934    0.966    FALSE 1 1.000  8804
# s[1,2,3]       0.893  0.051    0.781    0.899    0.974    FALSE 1 1.001  2499
# s[2,2,3]       0.923  0.033    0.849    0.927    0.974    FALSE 1 1.002  1521
# s[1,1,4]       0.955  0.022    0.906    0.957    0.989    FALSE 1 1.002  1374
# s[2,1,4]       0.970  0.010    0.948    0.971    0.986    FALSE 1 1.003   937
# s[1,2,4]       0.947  0.026    0.889    0.951    0.987    FALSE 1 1.009   268
# s[2,2,4]       0.962  0.016    0.924    0.964    0.987    FALSE 1 1.003  1418
# ann.s[1,1]     0.585  0.077    0.433    0.583    0.732    FALSE 1 1.001  3071
# ann.s[2,1]     0.582  0.041    0.504    0.582    0.661    FALSE 1 1.005   569
# ann.s[1,2]     0.471  0.074    0.325    0.471    0.612    FALSE 1 1.009   231
# ann.s[2,2]     0.517  0.063    0.397    0.517    0.643    FALSE 1 1.013   160
# mean.rho       2.188  0.164    1.882    2.182    2.530    FALSE 1 1.000 12000
# gamma          0.494  0.021    0.454    0.494    0.535    FALSE 1 1.002  1339
# mean.p         0.426  0.039    0.357    0.423    0.509    FALSE 1 1.017   128
# sigma.rho      0.295  0.063    0.194    0.289    0.439    FALSE 1 1.001  2281
# sigma.p        0.346  0.112    0.155    0.335    0.599    FALSE 1 1.003   771
# rho[1]         2.053  0.261    1.574    2.041    2.596    FALSE 1 1.000 12000
# [... output truncated ...]
# rho[20]        2.090  0.249    1.630    2.080    2.607    FALSE 1 1.001  2093
# p[1]           0.373  0.061    0.261    0.371    0.502    FALSE 1 1.002   935
# [... output truncated ...]
# p[20]          0.390  0.052    0.298    0.387    0.504    FALSE 1 1.008   261
# Nsp[1,1,1]    44.709 25.655    3.000   44.000   93.000    FALSE 1 1.002   985
# Nsp[2,1,1]    42.936 25.607    3.000   41.000   93.000    FALSE 1 1.003   770
# Nsp[1,2,1]    21.090 12.434    1.000   20.000   44.000    FALSE 1 1.007   285
# Nsp[2,2,1]    21.632 12.376    1.000   22.000   43.000    FALSE 1 1.008   260
# [... output truncated ...]
# Nsp[1,1,20]   47.030  9.908   29.000   47.000   68.000    FALSE 1 1.002  1025
# Nsp[2,1,20]   43.315  7.283   30.000   43.000   58.025    FALSE 1 1.011   194
# Nsp[1,2,20]   42.234  7.408   28.000   42.000   57.025    FALSE 1 1.003   629
# Nsp[2,2,20]   33.473  5.893   22.000   33.000   45.000    FALSE 1 1.010   194
# Nls[1,1,1]    41.313  8.804   25.000   41.000   60.000    FALSE 1 1.001  4464
# Nls[2,1,1]    73.558 14.826   48.000   73.000  106.000    FALSE 1 1.001  2510
# Nls[1,2,1]    39.570  7.463   26.000   39.000   55.000    FALSE 1 1.000  8933
# Nls[2,2,1]    35.942  4.451   28.000   36.000   45.000    FALSE 1 1.002   857
# [... output truncated ...]
# Nls[1,1,20]   70.303 13.962   45.000   70.000  100.000    FALSE 1 1.002   939
# Nls[2,1,20]   75.709 10.508   56.000   75.000   98.000    FALSE 1 1.011   191
# Nls[1,2,20]   71.980 14.120   47.000   71.000  102.000    FALSE 1 1.001  1828
# Nls[2,2,20]   62.821  6.975   50.000   63.000   77.000    FALSE 1 1.006   316


# ~~~~ posterior predictive checks ~~~~
# (only works if you ran model 2)
# Fig. 15.3
op <- par(mar=c(4, 4.5, 4, 0.5), las=1, "mfrow")
layout(matrix(1:6, 2, 3, byrow=TRUE), widths=c(1, 1, 1), heights=c(1, 1), TRUE)

z <- range(c(out2$sims.list$fit1[,1], out2$sims.list$fit1[,2]))
plot(x=out2$sims.list$fit1[,1], y=out2$sims.list$fit1[,2], pch=16,
    ylim=z, xlim=z, ylab="Discrepancy replicate data",
    xlab=" Discrepancy actual data", main="Spring counts",
    col=rgb(0, 0, 0, 0.3), axes=FALSE, cex=0.8)
segments(z[1], z[1], z[2], z[2])
pb <- round(mean(out2$sims.list$fit1[,1] < out2$sims.list$fit1[,2]), 2)
legend('bottomright', legend = bquote(p[B]==.(pb)), bty = "n")
axis(1); axis(2)

z <- range(c(out2$sims.list$fit2[,1], out2$sims.list$fit2[,2]))
z <- c(-30, z[2])
plot(x=out2$sims.list$fit2[,1], y=out2$sims.list$fit2[,2], pch=16,
    ylim=z, xlim=z, ylab="Discrepancy replicate data",
    xlab=" Discrepancy actual data", main="Summer counts females",
    col=rgb(0, 0, 0, 0.3), axes=FALSE, cex=0.8)
segments(z[1], z[1], z[2], z[2])
pb <- round(mean(out2$sims.list$fit2[,1] < out2$sims.list$fit2[,2]), 2)
legend('bottomright', legend = bquote(p[B]==.(pb)), bty = "n")
axis(1); axis(2)

z <- range(c(out2$sims.list$fit3[,1], out2$sims.list$fit3[,2]))
z <- c(-40, z[2])
plot(x=out2$sims.list$fit3[,1], y=out2$sims.list$fit3[,2], pch=16, ylim=z,
    xlim=z, ylab="Discrepancy replicate data", xlab=" Discrepancy actual data",
    main="Summer counts males", col=rgb(0, 0, 0, 0.3), axes=FALSE, cex=0.8)
segments(z[1], z[1], z[2], z[2])
pb <- round(mean(out2$sims.list$fit3[,1] < out2$sims.list$fit3[,2]), 2)
legend('bottomright', legend = bquote(p[B]==.(pb)), bty = "n")
axis(1); axis(2)

z <- range(c(out2$sims.list$fit4[,1], out2$sims.list$fit4[,2]))
plot(x=out2$sims.list$fit4[,1], y=out2$sims.list$fit4[,2], pch=16, ylim=z,
    xlim=z, ylab="Discrepancy replicate data",
    xlab=" Discrepancy actual data", main="Productivity data",
    col=rgb(0, 0, 0, 0.3), axes=FALSE, cex=0.8)
segments(z[1], z[1], z[2], z[2])
pb <- round(mean(out2$sims.list$fit4[,1] < out2$sims.list$fit4[,2]), 2)
legend('bottomright', legend = bquote(p[B]==.(pb)), bty = "n")
axis(1); axis(2)

z <- range(c(out2$sims.list$fit5[,1], out2$sims.list$fit5[,2]))
plot(x=out2$sims.list$fit5[,1], y=out2$sims.list$fit5[,2], pch=16, ylim=z,
    xlim=z, ylab="Discrepancy replicate data",
    xlab=" Discrepancy actual data", main="Chick sex ratio data",
    col=rgb(0, 0, 0, 0.3), axes=FALSE, cex=0.8)
segments(z[1], z[1], z[2], z[2])
pb <- round(mean(out2$sims.list$fit5[,1] < out2$sims.list$fit5[,2]), 2)
legend('bottomright', legend = bquote(p[B]==.(pb)), bty = "n")
axis(1); axis(2)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~ code for Fig. 15.4 ~~~~
library(scales)
library(RColorBrewer)
co <- c(brewer.pal(9, "Blues")[8], brewer.pal(9, "Reds")[7])

# Calculate population size in spring
qu <- function(x) quantile(x, p=c(0.025, 0.975))
NspF <- out1$sims.list$Nsp[,1,1,] + out1$sims.list$Nsp[,2,1,]
NspM <- out1$sims.list$Nsp[,1,2,] + out1$sims.list$Nsp[,2,2,]
NlsF <- out1$sims.list$Nls[,1,1,] + out1$sims.list$Nls[,2,1,]
NlsM <- out1$sims.list$Nls[,1,2,] + out1$sims.list$Nls[,2,2,]
year <- c(1997, 2000, 2003, 2006, 2009, 2012, 2015)

op <- par(mfrow=c(1,2), las=1, mar=c(4.5,4,2,1))

d <- 0.1
plot(y=apply(NspM, 2, mean), x=(1:20)-d, type="b", pch=16,
    ylab="Population size", xlab=NA, ylim=c(0, 180), col="blue",
    axes=FALSE, main="Spring")
polygon(c((1:20)-d, (20:1)-d), c(apply(NspM, 2, qu)[1,], rev(apply(NspM, 2, qu)[2,])),
    col=alpha(co[1], alpha=0.25), border=NA)
polygon(c((1:20)+d, (20:1)+d), c(apply(NspF, 2, qu)[1,], rev(apply(NspF, 2, qu)[2,])),
    col=alpha(co[2], alpha=0.25), border=NA)
points(y=apply(NspM, 2, mean), x=(1:20)-d, type="b", pch=16, col=co[1])
points(y=apply(NspF, 2, mean), x=(1:20)+d, type="b", pch=18, col=co[2])
points(grouse$count.sp, col="black", type="b", pch=1)
axis(2)
axis(1, at=1:20, labels=NA, tcl=-0.25)
axis(1, at=c(1, 4, 7, 10, 13, 16, 19), labels=year)
legend("topright", pch=c(18,16,1), col=c(rev(co), "black"),
    legend=c("Females", "Males", "Males count"), bty="n")

plot(y=apply(NlsM, 2, mean), x=(1:20)-d, type="b", pch=16, ylab=NA, xlab=NA,
    ylim=c(0, 180), col=co[1], axes=FALSE, main="Late summer")
polygon(c((1:20)-d, (20:1)-d), c(apply(NlsM, 2, qu)[1,], rev(apply(NlsM, 2, qu)[2,])),
    col=alpha(co[1], alpha=0.25), border=NA)
polygon(c((1:20)+d, (20:1)+d), c(apply(NlsF, 2, qu)[1,], rev(apply(NlsF, 2, qu)[2,])),
    col=alpha(co[2], alpha=0.25), border=NA)
points(y=apply(NlsM, 2, mean), x=(1:20)-d, type="b", pch=16, col=co[1])
points(y=apply(NlsF, 2, mean), x=(1:20)+d, type="b", pch=18, col=co[2])
axis(2, labels=NA)
axis(1, at=1:20, labels=NA, tcl=-0.25)
axis(1, at=c(1, 4, 7, 10, 13, 16, 19), labels=year)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Fig. 15.5 ~~~~
library(RColorBrewer)
co <- brewer.pal(9, "Reds")[7]

quant <- function(x) quantile(x, p=c(0.025, 0.975))
NF <- cbind(NspF[,1], NlsF[,1], NspF[,2], NlsF[,2], NspF[,3], NlsF[,3],  #### FIXME
    NspF[,4], NlsF[,4], NspF[,5], NlsF[,5], NspF[,6], NlsF[,6],
    NspF[,7], NlsF[,7], NspF[,8], NlsF[,8], NspF[,9], NlsF[,9],
    NspF[,10], NlsF[,10], NspF[,11], NlsF[,11], NspF[,12], NlsF[,12],
    NspF[,13], NlsF[,13], NspF[,14], NlsF[,14], NspF[,15], NlsF[,15],
    NspF[,16], NlsF[,16], NspF[,17], NlsF[,17], NspF[,18], NlsF[,18],
    NspF[,19], NlsF[,19], NspF[,20], NlsF[,20])
year <- c(1997, 2000, 2003, 2006, 2009, 2012, 2015)

op <- par(mar=c(3, 4.2, 1, 1), las=1)
plot(apply(NF, 2, mean), type="b", pch=rep(c(18,16),20),
    ylab="Female population size", xlab=NA, ylim=c(40, 180), col=co, axes=FALSE)
polygon(x=c(1:40, 40:1), y=c(apply(NF, 2, quant)[1,], rev(apply(NF, 2, quant)[2,])),
    col=alpha(co, alpha=0.25), border=NA)
points(apply(NF, 2, mean), type="b", pch=rep(c(18,16),20), col=co)
axis(2)
axis(1, at=seq(1, 40, 2), labels=NA, tcl=-0.25)
axis(1, at=c(1, 7, 13, 19, 25, 31, 37), labels=year)
legend("topleft", pch=c(18,16), col=c(co, co),
    legend=c("Spring", "Late summer"), bty="n")
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Fig. 15.6 ~~~~
library(RColorBrewer)
co <- brewer.pal(9, "Reds")[7]
year <- c(1997, 2000, 2003, 2006, 2009, 2012, 2015)

op <- par(las=1, mar=c(3, 4.2, 1, 1))
plot(y=out1$mean$rho, x=1:20, type="b", pch=16, ylab="Productivity",
    ylim=c(0.5,4), axes=FALSE, xlab=NA)
segments(1:20, out1$q2.5$rho, 1:20, out1$q97.5$rho)
lines(x=c(0.5, 20.5), y=rep(out1$mean$mean.rho, 2), col=co)
polygon(x=c(0.5,20.5,20.5,0.5),
    y=c(rep(out1$q2.5$mean.rho, 2), rep(out1$q97.5$mean.rho, 2)),
    col=alpha(co, alpha=0.25), border=NA)
axis(2)
axis(1, at=1:20, labels=NA, tcl=-0.25)
axis(1, at=c(1, 4, 7, 10, 13, 16, 19), labels=year)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Fig. 15.7 ~~~~
library(denstrip)
op <- par(mar=c(4.5, 6, 1, 1))

plot(0, ylim=c(0.8, 4.2), xlim=c(0,1), axes=FALSE, pch=NA,
    xlab="Annual survival", ylab=NA)
denstrip(out1$sims.list$ann.s[,1,1], at=4,
    ticks=c(out1$mean$ann.s[1,1], out1$q2.5$ann.s[1,1], out1$q97.5$ann.s[1,1]),
    twd=c(7,2.5,2.5), tlen=c(2,2,2), width=1/5)
denstrip(out1$sims.list$ann.s[,2,1], at=3,
    ticks=c(out1$mean$ann.s[2,1], out1$q2.5$ann.s[2,1], out1$q97.5$ann.s[2,1]),
    twd=c(7,2.5,2.5), tlen=c(2,2,2), width=1/5)
denstrip(out1$sims.list$ann.s[,1,2], at=2,
    ticks=c(out1$mean$ann.s[1,2], out1$q2.5$ann.s[1,2], out1$q97.5$ann.s[1,2]),
    twd=c(7,2.5,2.5), tlen=c(2,2,2), width=1/5)
denstrip(out1$sims.list$ann.s[,2,2], at=1,
    ticks=c(out1$mean$ann.s[2,2], out1$q2.5$ann.s[2,2], out1$q97.5$ann.s[2,2]),
    twd=c(7,2.5,2.5), tlen=c(2,2,2), width=1/5)
axis(1)
axis(2, las=1, at=c(4,3,2,1),
    labels=c("Juvenile\n females", "Adult\n females", "Juvenile\n males",
        "Adult\n males"))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
