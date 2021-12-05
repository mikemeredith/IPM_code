# Schaub & Kéry (2022) Integrated Population Models
# Chapter 15 : Black grouse
# -------------------------

# Run time 8 mins

# 15.4 Component data likelihoods
# =============================================

library(IPMbook); library(jagsUI)
data(grouse)
str(grouse)
# List of 10
# $ ch       : int [1:96, 1:128] 1 1 NA NA NA NA NA NA NA NA ...
# $ age      : num [1:96, 1:127] 2 2 NA NA NA NA NA NA NA NA ...
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
getLast <- function(x) max(which(!is.na(x)))              # Last occasion function
jags.data <- with(grouse, list(ch=ch, age=age, sex=sex, season=season, f=getFirst(ch),
    k=apply(ch, 1, getLast), C.sp=count.sp, C.lsM=count.lsM, C.lsF=count.lsF, C.lsC=count.lsC,
    u=u, v=v, nind=nrow(ch), ny=length(count.sp), pinit=dUnif(1, 100)))
str(jags.data)
# List of 15
# $ ch    : int [1:96, 1:128] 1 1 NA NA NA NA NA NA NA NA ...
# $ age   : num [1:96, 1:127] 2 2 NA NA NA NA NA NA NA NA ...
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
  for (a in 1:2){                                         # age
    for (j in 1:2){                                       # sex
      for (m in 1:4){                                     # season
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
    Nsp[1,1,t+1] ~ dbin(s[1,1,2]^0.5 * s[1,1,3]^2 * s[1,1,4]^5 * s[2,1,1], Nls[1,1,t]) # 1y
    Nsp[2,1,t+1] ~ dbin(s[2,1,2]^0.5 * s[2,1,3]^2 * s[2,1,4]^5 * s[2,1,1], Nls[2,1,t]) # >1y

    # Males
    Nsp[1,2,t+1] ~ dbin(s[1,2,2]^0.5 * s[1,2,3]^2 * s[1,2,4]^5 * s[2,2,1], Nls[1,2,t]) # 1y
    Nsp[2,2,t+1] ~ dbin(s[2,2,2]^0.5 * s[2,2,3]^2 * s[2,2,4]^5 * s[2,2,1], Nls[2,2,t]) # >1y
  }

  # Process model: population sizes in late summer
  for (t in 1:ny){
    # Total number of chicks
    F[t] ~ dpois(Nls[2,1,t] * rho[t])
    # Allocate chicks to a sex
    Nls[1,1,t] ~ dbin(gamma, F[t])                        # Female chicks
    Nls[1,2,t] <- F[t] - Nls[1,1,t]                       # Male chicks
    # Survival
    Nls[2,1,t] ~ dbin(s[2,1,1] * s[2,1,2]^2.5, (Nsp[1,1,t] + Nsp[2,1,t])) # ≥1y females
    Nls[2,2,t] ~ dbin(s[2,2,1] * s[2,2,2]^2.5, (Nsp[1,2,t] + Nsp[2,2,t])) # ≥1y males
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
  ann.s[1,1] <- s[1,1,2] * s[1,1,3]^2 * s[1,1,4]^5 * s[2,1,1]       # juv females
  ann.s[2,1] <- s[2,1,1]^2 * s[2,1,2]^3 * s[2,1,3]^2 * s[2,1,4]^5   # ad females
  ann.s[1,2] <- s[1,2,2] * s[1,2,3]^2 * s[1,2,4]^5 * s[2,2,1]       # juv males
  ann.s[2,2] <- s[2,2,1]^2 * s[2,2,2]^3 * s[2,2,3]^2 * s[2,2,4]^5   # ad males
}
")

# Initial values
inits <- function(){
  Nsp <- array(NA, dim=c(2, 2, 20))
  Nsp[,,1] <- round(runif(4, 30, 40))
  s <- array(runif(2*2*4, 0.94, 0.97), dim=c(2, 2, 4))
  list(Nsp=Nsp, s=s)
}

# Parameters monitored
parameters <- c("s", "ann.s", "mean.rho", "gamma", "mean.p", "sigma.rho", "sigma.p", "rho", "p",
    "Nsp", "Nls")

# MCMC settings
# ni <- 120000; nb <- 20000; nc <- 3; nt <- 20; na <- 10000
ni <- 30000; nb <- 10000; nc <- 3; nt <- 5; na <- 5000  # ~~~ for testing

# Call JAGS (ART 8 min), check convergence and summarize posteriors
set.seed(101)                                             # Reproducible results without crashing
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out1)


# ~~~~ model 1 plus posterior predictive checks ~~~~
# Write JAGS model file
cat(file="model2.txt", "
model {
  # Priors and linear models
  # Survival
  for (a in 1:2){                          # age
    for (j in 1:2){                        # sex
      for (m in 1:4){                      # season
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
    # F[t] ~ dpois((Nsp[1,2,t] + Nsp[2,2,t]) * s[2,2,1] * rho[t]) # Correction 2021-07-28
    F[t] ~ dpois(Nls[2,1,t] * rho[t])
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
  # s <- array(runif(2*2*4, 0.90, 0.97), dim=c(2, 2, 4)) # Correction 2021-07-28
  s <- array(runif(2*2*4, 0.94, 0.97), dim=c(2, 2, 4))
  list(Nsp=Nsp, s=s)
}

# Parameters monitored
parameters <- c("s", "ann.s", "mean.rho", "gamma", "mean.p", "sigma.rho",
    "sigma.p", "rho", "p", "Nsp", "Nls", "fit1", "fit2", "fit3", "fit4", "fit5")

# MCMC settings
ni <- 120000; nb <- 20000; nc <- 3; nt <- 20; na <- 10000

# Call JAGS (ART 5 min), check convergence and summarize posteriors
set.seed(101) # Reproducible results without crashing
out2 <- jags(jags.data, inits, parameters, "model2.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

save(out1, out2, file="Grouse.Results.Rdata")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 15.6 Results
# ============

print(out1, 3)
#                 mean     sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# s[1,1,1]       0.501  0.290    0.024    0.503    0.976    FALSE 1 1.000 15000
# s[2,1,1]       0.943  0.018    0.904    0.944    0.972    FALSE 1 1.000 15000
# s[1,2,1]       0.501  0.288    0.025    0.506    0.974    FALSE 1 1.000 15000
# s[2,2,1]       0.880  0.035    0.807    0.883    0.942    FALSE 1 1.000 15000
# s[1,1,2]       0.930  0.065    0.759    0.948    0.998    FALSE 1 1.000 15000
# s[2,1,2]       0.948  0.013    0.920    0.949    0.971    FALSE 1 1.000  8594
# s[1,2,2]       0.843  0.097    0.612    0.861    0.979    FALSE 1 1.000 11434
# s[2,2,2]       0.970  0.016    0.932    0.973    0.994    FALSE 1 1.000 15000
# s[1,1,3]       0.878  0.055    0.759    0.884    0.969    FALSE 1 1.000 14762
# s[2,1,3]       0.927  0.022    0.879    0.929    0.964    FALSE 1 1.000 15000
# s[1,2,3]       0.863  0.063    0.728    0.868    0.966    FALSE 1 1.001  3878
# s[2,2,3]       0.915  0.036    0.834    0.920    0.972    FALSE 1 1.000  5561
# s[1,1,4]       0.937  0.027    0.881    0.939    0.983    FALSE 1 1.000  6933
# s[2,1,4]       0.967  0.011    0.943    0.968    0.984    FALSE 1 1.000  8364
# s[1,2,4]       0.927  0.032    0.858    0.930    0.980    FALSE 1 1.000 15000
# s[2,2,4]       0.956  0.019    0.912    0.959    0.985    FALSE 1 1.000  5847
# ann.s[1,1]     0.490  0.071    0.359    0.489    0.634    FALSE 1 1.000 10877
# ann.s[2,1]     0.551  0.045    0.462    0.550    0.637    FALSE 1 1.000 15000
# ann.s[1,2]     0.378  0.068    0.253    0.376    0.519    FALSE 1 1.000  7952
# ann.s[2,2]     0.476  0.066    0.347    0.476    0.603    FALSE 1 1.000 14173
# mean.rho       2.175  0.163    1.867    2.170    2.514    FALSE 1 1.000 15000
# gamma          0.493  0.021    0.451    0.493    0.534    FALSE 1 1.000 15000
# mean.p         0.430  0.038    0.361    0.428    0.510    FALSE 1 1.000 15000
# sigma.rho      0.298  0.065    0.194    0.290    0.447    FALSE 1 1.000  8608
# sigma.p        0.376  0.112    0.187    0.364    0.620    FALSE 1 1.001  2673
# rho[1]         1.963  0.252    1.500    1.950    2.484    FALSE 1 1.000 14867
# [... output truncated ...]
# rho[20]        2.092  0.248    1.634    2.082    2.602    FALSE 1 1.000 10745
# p[1]           0.408  0.059    0.298    0.405    0.531    FALSE 1 1.000  7604
# [... output truncated ...]
# p[20]          0.406  0.057    0.304    0.403    0.526    FALSE 1 1.000 15000
# Nsp[1,1,1]    38.809 23.047    2.000   38.000   84.000    FALSE 1 1.000  3816
# Nsp[2,1,1]    38.726 23.279    2.000   37.000   84.000    FALSE 1 1.000 14971
# Nsp[1,2,1]    19.692 11.514    1.000   19.000   41.000    FALSE 1 1.001  3196
# Nsp[2,2,1]    19.750 11.504    1.000   19.000   41.000    FALSE 1 1.001  2832
# Nsp[1,1,2]    32.980  7.503   20.000   32.000   49.000    FALSE 1 1.000 12507
# [... output truncated ...]
# Nsp[1,1,20]   45.686 10.614   27.000   45.000   69.000    FALSE 1 1.000 15000
# Nsp[2,1,20]   42.317  7.023   30.000   42.000   57.000    FALSE 1 1.000 15000
# Nsp[1,2,20]   41.317  7.747   27.000   41.000   57.000    FALSE 1 1.000 11306
# Nsp[2,2,20]   31.499  5.899   20.000   31.000   43.000    FALSE 1 1.000  7807
# Nls[1,1,1]    63.130 12.022   42.000   62.000   89.000    FALSE 1 1.000 15000
# Nls[2,1,1]    63.566  9.722   46.000   63.000   84.000    FALSE 1 1.000  5830
# Nls[1,2,1]    60.773 12.075   40.000   60.000   87.000    FALSE 1 1.000 15000
# Nls[2,2,1]    32.500  4.567   24.000   32.000   42.000    FALSE 1 1.000 15000
# [... output truncated ...]
# Nls[1,1,20]   74.401 17.089   45.000   73.000  111.000    FALSE 1 1.000 15000
# Nls[2,1,20]   72.206 10.858   53.000   72.000   95.000    FALSE 1 1.000 15000
# Nls[1,2,20]   76.579 17.259   47.000   75.000  114.000    FALSE 1 1.000 15000
# Nls[2,2,20]   59.468  6.765   47.000   59.000   73.000    FALSE 1 1.000 11886



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
    ylab="Population size", xlab=NA, ylim=c(0, 220), col="blue",
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
    ylim=c(0, 220), col=co[1], axes=FALSE, main="Late summer")
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
    ylab="Female population size", xlab=NA, ylim=c(40, 220), col=co, axes=FALSE)
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
