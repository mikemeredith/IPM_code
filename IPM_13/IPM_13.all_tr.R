# Schaub & Kéry (2022) Integrated Population Models
# Chapter 13 : Horseshoe bat
# --------------------------

# Run time for test script 40 mins, full run 26 hrs.

# 13.4 Single data likelihoods
# ============================

# Load greater horseshoe bat data and produce data overview
library(IPMbook); library(jagsUI)
data(bats)
str(bats)
# List of 7
# $ ch     : int [1:574, 1:29] 0 0 0 0 0 0 0 0 0 0 ...
# $ age    : num [1:574] 1 1 1 1 1 1 1 1 1 1 ...
# $ sex    : num [1:574] 2 2 2 1 2 1 1 2 2 1 ...
# $ J.count: num [1:29] 13 11 11 16 16 16 18 20 17 16 ...
# $ Jm     : num 293
# $ Jf     : num 289
# $ A.count: num [1:29] 27 33 NA 27 42 43 41 40 44 33 ...


# 13.4.1 Capture-recapture data
# '''''''''''''''''''''''''''''

# Compute the age matrix x
f <- getFirst(bats$ch)
x <- matrix(NA, ncol=ncol(bats$ch)-1, nrow=nrow(bats$ch))
index <- cbind(1:nrow(bats$ch), f)
x[index] <- bats$age
for (t in 2:ncol(x)){
  x[is.na(x[,t]),t] <- x[is.na(x[,t]),t-1] + 1
}
x[x>3] <- 3

k <- c(rep(12, 14), 1, 2, 3, 4, 5, 12, 12, 6, 7, 8, 12, 11, 9, 10)

# 13.4.2 Juvenile and population count data (no code)

# 13.5 The integrated population models
# =====================================

# Bundle data and produce a data overview
jags.data <- with(bats, list(ch=ch, f=getFirst(ch), nind=nrow(ch), n.occasions=ncol(ch), x=x,
    sex=sex, k=k, J=J.count, Jf=Jf, Jm=Jm, C=A.count, pinit1=dUnif(1, 50), pinit2=dUnif(1, 100),
    beta.priors=as.numeric(getBeta2Par(0.4, 0.15)), rpres=which(k<11), nrpres=which(k>10)))
str(jags.data)
# List of 16
# $ ch         : int [1:574, 1:29] 0 0 0 0 0 0 0 0 0 0 ...
# $ f          : int [1:574] 4 4 4 4 4 4 4 4 4 4 ...
# $ nind       : int 574
# $ n.occasions: int 29
# $ x          : num [1:574, 1:28] NA NA NA NA NA NA NA NA NA NA ...
# $ sex        : num [1:574] 2 2 2 1 2 1 1 2 2 1 ...
# $ k          : num [1:28] 12 12 12 12 12 12 12 12 12 12 ...
# $ J          : num [1:29] 13 11 11 16 16 16 18 20 17 16 ...
# $ Jf         : num 289
# $ Jm         : num 293
# $ C          : num [1:29] 27 33 NA 27 42 43 41 40 44 33 ...
# $ pinit1     : num [1:50] 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 ...
# $ pinit2     : num [1:100] 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 ...
# $ beta.priors: num [1:2] 3.87 5.8
# $ rpres      : int [1:10] 15 16 17 18 19 22 23 24 27 28
# $ nrpres     : int [1:18] 1 2 3 4 5 6 7 8 9 10 ...

# Write JAGS model file
cat(file="model1.txt", "
model {
  # Priors and linear models
  # Survival and recapture
  for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- sur[sex[i],x[i,t],t]
      p[i,t] <- pt[sex[i],x[i,t],k[t]]
    } #t
  } #i

  for (t in 1:(n.occasions-1)){
    logit(sur[1,1,t]) <- eps[1,t]
    logit(sur[1,2,t]) <- eps[2,t]
    logit(sur[1,3,t]) <- eps[3,t]
    logit(sur[2,1,t]) <- eps[4,t]
    logit(sur[2,2,t]) <- eps[5,t]
    logit(sur[2,3,t]) <- eps[6,t]
    eps[1:6,t] ~ dmnorm.vcov(mu.s[], Sigma[,])
  }

  # Recaptures
  for (s in 1:2){
    for (a in 1:3){
      for (t in 1:10){
        logit.pt[s,a,t] ~ dnorm(mu.p[s,a], tau.p[s,a])
        pt[s,a,t] <- ilogit(logit.pt[s,a,t])
      } #t
      pt[s,a,11] ~ dunif(0, 1)                    # rate in 2015, incomplete capture
      pt[s,a,12] <- 0                             # fix to 0 in years without capture
    } #a
  } #s

  mu.s[1] <- logit(mean.s[1,1])
  mu.s[2] <- logit(mean.s[1,2])
  mu.s[3] <- logit(mean.s[1,3])
  mu.s[4] <- logit(mean.s[2,1])
  mu.s[5] <- logit(mean.s[2,2])
  mu.s[6] <- logit(mean.s[2,3])

  for (j in 1:6){
    sd.s[j] ~ dunif(0.001, 10)                    # Prior of the temporal SD
    Sigma[j,j] <- pow(sd.s[j], 2)
  }

  for (m in 1:5){
    for (j in (m+1):6){
      Sigma[m,j] <- zeta[m,j] * pow(Sigma[m,m] * Sigma[j,j], 0.5)
      Sigma[j,m] <- Sigma[m,j]
      zeta[m,j] ~ dunif(-1, 1)                    # Priors for correlation coefficients
    } #j
  } #m

  for (s in 1:2){
    for (a in 1:3){
      mean.s[s,a] ~ dunif(0, 1)
      mean.p[s,a] ~ dunif(0, 1)
      mu.p[s,a] <- logit(mean.p[s,a])
      sigma.p[s,a] ~ dunif(0, 5)
      tau.p[s,a] <- pow(sigma.p[s,a],-2)
    } #a
  } #s

  # Productivity
  for (t in 1:n.occasions){
    rho[1,t] <- mean.rho[1]                         # of 2y females
    logit(rho[2,t]) <- mu.rho + eps.rho[t]          # of >2y females
    eps.rho[t] ~ dnorm(0, tau.rho)
  }

  # Priors for productivity of 2y females
  # To run IPMs 2-4, manually un-hash one line and hash out the others
  mean.rho[1] <- 0
  # mean.rho[1] <- 0.99                           # does not work if = 1 due to initial values
  # mean.rho[1] ~ dunif(0, 1)
  # mean.rho[1] ~ dbeta(beta.priors[1], beta.priors[2])

  # Priors for productivity of >2y females
  mean.rho[2] ~ dunif(0,1)
  mu.rho <- logit(mean.rho[2])
  sigma.rho ~ dunif(0, 10)
  tau.rho <- pow(sigma.rho, -2)

  # Variance of the observation error
  sigma.y ~ dunif(0.5, 20)
  tau.y <- pow(sigma.y, -2)

  # Sex ratio of newborn individuals
  xi ~ dbeta(Jf, Jm)

  # Compute the proportion of bats present in the colony
  for (s in 1:2){
    for (a in 1:3){
      lkappa[s,a,1] ~ dnorm(mu.p[s,a], tau.p[s,a])
      kappa[s,a,1] <- ilogit(lkappa[s,a,1])
      for (t in 1:length(nrpres)){
        lkappa[s,a,nrpres[t]+1] ~ dnorm(mu.p[s,a], tau.p[s,a])
        kappa[s,a,nrpres[t]+1] <- ilogit(lkappa[s,a,nrpres[t]+1])
      } #t
      for (t in 1:length(rpres)){
        kappa[s,a,rpres[t]+1] <- pt[s,a,t]
      } #t
    } #a
  } #s

  # Population count data (state-space model)
  # Model for the initial population size: uniform priors
  N[1,1,1] ~ dcat(pinit1)
  N[1,2,1] ~ dcat(pinit1)
  N[1,3,1] ~ dcat(pinit2)
  N[2,1,1] ~ dcat(pinit1)
  N[2,2,1] ~ dcat(pinit1)
  N[2,3,1] ~ dcat(pinit2)

  # Process model: population sizes of females and males
  for(t in 2:n.occasions){
    N[1,1,t] ~ dbin(sur[1,1,t-1], NB[1,t-1])
    N[1,2,t] ~ dbin(sur[1,2,t-1], N[1,1,t-1])
    N[1,3,t] ~ dbin(sur[1,3,t-1], N[1,2,t-1] + N[1,3,t-1])
    N[2,1,t] ~ dbin(sur[2,1,t-1], NB[2,t-1])
    N[2,2,t] ~ dbin(sur[2,2,t-1], N[2,1,t-1])
    N[2,3,t] ~ dbin(sur[2,3,t-1], N[2,2,t-1] + N[2,3,t-1])
  }

  # Observation model
  for (t in 1:n.occasions){
    Ns[t] <- sum(N[,,t])                          # total population size
    U[t] <- N[1,1,t] * kappa[1,1,t] + N[1,2,t] * kappa[1,2,t] + N[1,3,t] * kappa[1,3,t] + N[2,1,t] *
        kappa[2,1,t] + N[2,2,t] * kappa[2,2,t] + N[2,3,t] * kappa[2,3,t]
    C[t] ~ dnorm(U[t], tau.y)
  }

  # Capture-recapture data (state-space model)
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(phi[i,t-1] * z[i,t-1])
      # Observation process
      ch[i,t] ~ dbern(p[i,t-1] * z[i,t])
    } #t
  } #i

  # Newborn count data (binomial models)
  for (t in 1:n.occasions){
    G[1,t] ~ dbin(rho[1,t], N[1,2,t])             # Newborn by 2y females
    G[2,t] ~ dbin(rho[2,t], N[1,3,t])             # Newborn by >2y females
    J[t] ~ dsum(G[1,t], G[2,t])                   # Total newborn
    NB[1,t] ~ dbin(xi, J[t])                      # Newborn females
    NB[2,t] <- J[t] - NB[1,t]                     # Newborn males
  }
}
")

# Initial values
Ninit <- array(NA, dim=c(2, 3, 29))
Ninit[1,3,] <- 60
G <- matrix(0, nrow=2, ncol=29)
G[2,] <- bats$J.count
inits <- function(){list(z=zInit(bats$ch), N=Ninit, G=G)}

# Parameters monitored
parameters <- c("mean.s", "Sigma", "mean.rho", "sigma.rho", "xi", "zeta", "mean.p", "sigma.p", "sigma.y",
    "sur", "pt", "rho", "kappa", "N", "U", "Ns", "NB", "G")

# MCMC settings
# ni <- 500000; nb <- 200000; nc <- 3; nt <- 100; na <- 5000
ni <- 5000; nb <- 2000; nc <- 3; nt <- 1; na <- 5000  # ~~~ for testing, 8 mins

# Call JAGS (ART 600 min) and check convergence
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out1)


# ~~~~ code to run the other IPMs ~~~~

# IPM2
# Write JAGS model file
cat(file="model2.txt", "
model {
  # Priors and linear models
  # Survival and recapture
  for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- sur[sex[i],x[i,t],t]
      p[i,t] <- pt[sex[i],x[i,t],k[t]]
    } #t
  } #i

  for (t in 1:(n.occasions-1)){
    logit(sur[1,1,t]) <- eps[1,t]
    logit(sur[1,2,t]) <- eps[2,t]
    logit(sur[1,3,t]) <- eps[3,t]
    logit(sur[2,1,t]) <- eps[4,t]
    logit(sur[2,2,t]) <- eps[5,t]
    logit(sur[2,3,t]) <- eps[6,t]
    eps[1:6,t] ~ dmnorm.vcov(mu.s[], Sigma[,])
  }

  # Recaptures
  for (s in 1:2){
    for (a in 1:3){
      for (t in 1:10){
        logit.pt[s,a,t] ~ dnorm(mu.p[s,a], tau.p[s,a])
        pt[s,a,t] <- ilogit(logit.pt[s,a,t])
      } #t
      pt[s,a,11] ~ dunif(0, 1)   # rate in 2015, incomplete capture
      pt[s,a,12] <- 0            # fix to 0 in years without capture
    } #a
  } #s

  mu.s[1] <- logit(mean.s[1,1])
  mu.s[2] <- logit(mean.s[1,2])
  mu.s[3] <- logit(mean.s[1,3])
  mu.s[4] <- logit(mean.s[2,1])
  mu.s[5] <- logit(mean.s[2,2])
  mu.s[6] <- logit(mean.s[2,3])

  for (j in 1:6){
    sd.s[j] ~ dunif(0.001, 10)  # Prior of the temporal SD
    Sigma[j,j] <- pow(sd.s[j], 2)
  }

  for (m in 1:5){
    for (j in (m+1):6){
      Sigma[m,j] <- zeta[m,j] * pow(Sigma[m,m] * Sigma[j,j], 0.5)
      Sigma[j,m] <- Sigma[m,j]
      zeta[m,j] ~ dunif(-1, 1)     # Priors for correlation coefficients
    } #j
  } #m

  for (s in 1:2){
    for (a in 1:3){
      mean.s[s,a] ~ dunif(0, 1)
      mean.p[s,a] ~ dunif(0, 1)
      mu.p[s,a] <- logit(mean.p[s,a])
      sigma.p[s,a] ~ dunif(0, 5)
      tau.p[s,a] <- pow(sigma.p[s,a],-2)
    } #a
  } #s

  # Productivity
  for (t in 1:n.occasions){
    rho[1,t] <- mean.rho[1]                  # of 2y females
    logit(rho[2,t]) <- mu.rho + eps.rho[t]   # of >2y females
    eps.rho[t] ~ dnorm(0, tau.rho)
  }
  # Priors for productivity of 2y females
  #mean.rho[1] <- 0
  mean.rho[1] <- 0.99       # does not work if = 1 due to initial values
  #mean.rho[1] ~ dunif(0, 1)
  #mean.rho[1] ~ dbeta(beta.priors[1], beta.priors[2])

  # Priors for productivity of >2y females
  mean.rho[2] ~ dunif(0,1)
  mu.rho <- logit(mean.rho[2])
  sigma.rho ~ dunif(0, 10)
  tau.rho <- pow(sigma.rho, -2)

  # Variance of the observation error
  sigma.y ~ dunif(0.5, 20)
  tau.y <- pow(sigma.y, -2)

  # Sex ratio of newborn individuals
  xi ~ dbeta(Jf, Jm)

  # Compute the proportion of bats present in the colony
  for (s in 1:2){
    for (a in 1:3){
      lkappa[s,a,1] ~ dnorm(mu.p[s,a], tau.p[s,a])
      kappa[s,a,1] <- ilogit(lkappa[s,a,1])
      for (t in 1:length(nrpres)){
        lkappa[s,a,nrpres[t]+1] ~ dnorm(mu.p[s,a], tau.p[s,a])
        kappa[s,a,nrpres[t]+1] <- ilogit(lkappa[s,a,nrpres[t]+1])
      } #t
      for (t in 1:length(rpres)){
        kappa[s,a,rpres[t]+1] <- pt[s,a,t]
      } #t
    } #a
  } #s

  # Population count data (state-space model)
  # Model for the initial population size: uniform priors
  N[1,1,1] ~ dcat(pinit1)
  N[1,2,1] ~ dcat(pinit1)
  N[1,3,1] ~ dcat(pinit2)
  N[2,1,1] ~ dcat(pinit1)
  N[2,2,1] ~ dcat(pinit1)
  N[2,3,1] ~ dcat(pinit2)

  # Process model: population sizes of females and males
  for(t in 2:n.occasions){
    N[1,1,t] ~ dbin(sur[1,1,t-1], NB[1,t-1])
    N[1,2,t] ~ dbin(sur[1,2,t-1], N[1,1,t-1])
    N[1,3,t] ~ dbin(sur[1,3,t-1], N[1,2,t-1] + N[1,3,t-1])
    N[2,1,t] ~ dbin(sur[2,1,t-1], NB[2,t-1])
    N[2,2,t] ~ dbin(sur[2,2,t-1], N[2,1,t-1])
    N[2,3,t] ~ dbin(sur[2,3,t-1], N[2,2,t-1] + N[2,3,t-1])
  }

  # Observation model
  for (t in 1:n.occasions){
    Ns[t] <- sum(N[,,t])       # total population size
    U[t] <- N[1,1,t] * kappa[1,1,t] + N[1,2,t] * kappa[1,2,t] + N[1,3,t] * kappa[1,3,t]  + N[2,1,t] * kappa[2,1,t]  + N[2,2,t] * kappa[2,2,t] + N[2,3,t] * kappa[2,3,t]
    C[t] ~ dnorm(U[t], tau.y)
  }

  # Capture-recapture data (state-space model)
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(phi[i,t-1] * z[i,t-1])
      # Observation process
      ch[i,t] ~ dbern(p[i,t-1] * z[i,t])
    } #t
  } #i

  # Newborn count data (binomial models)
  for (t in 1:n.occasions){
    G[1,t] ~ dbin(rho[1,t], N[1,2,t])   # Newborn by 2y females
    G[2,t] ~ dbin(rho[2,t], N[1,3,t])   # Newborn by >2y females
    J[t] ~ dsum(G[1,t], G[2,t])         # Total newborn
    NB[1,t] ~ dbin(xi, J[t])            # Newborn females
    NB[2,t] <- J[t] - NB[1,t]           # Newborn males
  }
}
")

# Initial values
Ninit <- array(NA, dim=c(2, 3, 29))
Ninit[1,3,] <- 60
G <- matrix(0, nrow=2, ncol=29)
G[2,] <- bats$J.count
inits <- function(){list(z=zInit(bats$ch), N=Ninit, G=G)}

# Parameters monitored
parameters <- c("mean.s", "Sigma", "mean.rho", "sigma.rho", "xi", "zeta", "mean.p",
    "sigma.p", "sigma.y", "sur", "pt", "rho", "kappa", "N", "U", "Ns", "NB", "G")

# MCMC settings
# ni <- 500000; nb <- 200000; nc <- 3; nt <- 100; na <- 5000
ni <- 5000; nb <- 2000; nc <- 3; nt <- 1; na <- 5000  # ~~~ for testing

# Call JAGS (ART 554 min) and check convergence
out2 <- jags(jags.data, inits, parameters, "model2.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out2)


# IPM3
# Write JAGS model file
cat(file="model3.txt", "
model {
  # Priors and linear models
  # Survival and recapture
  for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- sur[sex[i],x[i,t],t]
      p[i,t] <- pt[sex[i],x[i,t],k[t]]
    } #t
  } #i

  for (t in 1:(n.occasions-1)){
    logit(sur[1,1,t]) <- eps[1,t]
    logit(sur[1,2,t]) <- eps[2,t]
    logit(sur[1,3,t]) <- eps[3,t]
    logit(sur[2,1,t]) <- eps[4,t]
    logit(sur[2,2,t]) <- eps[5,t]
    logit(sur[2,3,t]) <- eps[6,t]
    eps[1:6,t] ~ dmnorm.vcov(mu.s[], Sigma[,])
  }

  # Recaptures
  for (s in 1:2){
    for (a in 1:3){
      for (t in 1:10){
        logit.pt[s,a,t] ~ dnorm(mu.p[s,a], tau.p[s,a])
        pt[s,a,t] <- ilogit(logit.pt[s,a,t])
      } #t
      pt[s,a,11] ~ dunif(0, 1)   # rate in 2015, incomplete capture
      pt[s,a,12] <- 0            # fix to 0 in years without capture
    } #a
  } #s

  mu.s[1] <- logit(mean.s[1,1])
  mu.s[2] <- logit(mean.s[1,2])
  mu.s[3] <- logit(mean.s[1,3])
  mu.s[4] <- logit(mean.s[2,1])
  mu.s[5] <- logit(mean.s[2,2])
  mu.s[6] <- logit(mean.s[2,3])

  for (j in 1:6){
    sd.s[j] ~ dunif(0.001, 10)  # Prior of the temporal SD
    Sigma[j,j] <- pow(sd.s[j], 2)
  }

  for (m in 1:5){
    for (j in (m+1):6){
      Sigma[m,j] <- zeta[m,j] * pow(Sigma[m,m] * Sigma[j,j], 0.5)
      Sigma[j,m] <- Sigma[m,j]
      zeta[m,j] ~ dunif(-1, 1)     # Priors for correlation coefficients
    } #j
  } #m

  for (s in 1:2){
    for (a in 1:3){
      mean.s[s,a] ~ dunif(0, 1)
      mean.p[s,a] ~ dunif(0, 1)
      mu.p[s,a] <- logit(mean.p[s,a])
      sigma.p[s,a] ~ dunif(0, 5)
      tau.p[s,a] <- pow(sigma.p[s,a],-2)
    } #a
  } #s

  # Productivity
  for (t in 1:n.occasions){
    rho[1,t] <- mean.rho[1]                  # of 2y females
    logit(rho[2,t]) <- mu.rho + eps.rho[t]   # of >2y females
    eps.rho[t] ~ dnorm(0, tau.rho)
  }
  # Priors for productivity of 2y females
  #mean.rho[1] <- 0
  #mean.rho[1] <- 0.99       # does not work if = 1 due to initial values
  mean.rho[1] ~ dunif(0, 1)
  #mean.rho[1] ~ dbeta(beta.priors[1], beta.priors[2])

  # Priors for productivity of >2y females
  mean.rho[2] ~ dunif(0,1)
  mu.rho <- logit(mean.rho[2])
  sigma.rho ~ dunif(0, 10)
  tau.rho <- pow(sigma.rho, -2)

  # Variance of the observation error
  sigma.y ~ dunif(0.5, 20)
  tau.y <- pow(sigma.y, -2)

  # Sex ratio of newborn individuals
  xi ~ dbeta(Jf, Jm)

  # Compute the proportion of bats present in the colony
  for (s in 1:2){
    for (a in 1:3){
      lkappa[s,a,1] ~ dnorm(mu.p[s,a], tau.p[s,a])
      kappa[s,a,1] <- ilogit(lkappa[s,a,1])
      for (t in 1:length(nrpres)){
        lkappa[s,a,nrpres[t]+1] ~ dnorm(mu.p[s,a], tau.p[s,a])
        kappa[s,a,nrpres[t]+1] <- ilogit(lkappa[s,a,nrpres[t]+1])
      } #t
      for (t in 1:length(rpres)){
        kappa[s,a,rpres[t]+1] <- pt[s,a,t]
      } #t
    } #a
  } #s

  # Population count data (state-space model)
  # Model for the initial population size: uniform priors
  N[1,1,1] ~ dcat(pinit1)
  N[1,2,1] ~ dcat(pinit1)
  N[1,3,1] ~ dcat(pinit2)
  N[2,1,1] ~ dcat(pinit1)
  N[2,2,1] ~ dcat(pinit1)
  N[2,3,1] ~ dcat(pinit2)

  # Process model: population sizes of females and males
  for(t in 2:n.occasions){
    N[1,1,t] ~ dbin(sur[1,1,t-1], NB[1,t-1])
    N[1,2,t] ~ dbin(sur[1,2,t-1], N[1,1,t-1])
    N[1,3,t] ~ dbin(sur[1,3,t-1], N[1,2,t-1] + N[1,3,t-1])
    N[2,1,t] ~ dbin(sur[2,1,t-1], NB[2,t-1])
    N[2,2,t] ~ dbin(sur[2,2,t-1], N[2,1,t-1])
    N[2,3,t] ~ dbin(sur[2,3,t-1], N[2,2,t-1] + N[2,3,t-1])
  }

  # Observation model
  for (t in 1:n.occasions){
    Ns[t] <- sum(N[,,t])       # total population size
    U[t] <- N[1,1,t] * kappa[1,1,t] + N[1,2,t] * kappa[1,2,t] + N[1,3,t] * kappa[1,3,t]  + N[2,1,t] * kappa[2,1,t]  + N[2,2,t] * kappa[2,2,t] + N[2,3,t] * kappa[2,3,t]
    C[t] ~ dnorm(U[t], tau.y)
  }

  # Capture-recapture data (state-space model)
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(phi[i,t-1] * z[i,t-1])
      # Observation process
      ch[i,t] ~ dbern(p[i,t-1] * z[i,t])
    } #t
  } #i

  # Newborn count data (binomial models)
  for (t in 1:n.occasions){
    G[1,t] ~ dbin(rho[1,t], N[1,2,t])   # Newborn by 2y females
    G[2,t] ~ dbin(rho[2,t], N[1,3,t])   # Newborn by >2y females
    J[t] ~ dsum(G[1,t], G[2,t])         # Total newborn
    NB[1,t] ~ dbin(xi, J[t])            # Newborn females
    NB[2,t] <- J[t] - NB[1,t]           # Newborn males
  }
}
")

# Initial values
Ninit <- array(NA, dim=c(2, 3, 29))
Ninit[1,3,] <- 60
G <- matrix(0, nrow=2, ncol=29)
G[2,] <- bats$J.count
inits <- function(){list(z=zInit(bats$ch), N=Ninit, G=G)}

# Parameters monitored
parameters <- c("mean.s", "Sigma", "mean.rho", "sigma.rho", "xi", "zeta", "mean.p",
    "sigma.p", "sigma.y", "sur", "pt", "rho", "kappa", "N", "U", "Ns", "NB", "G")

# MCMC settings
# ni <- 500000; nb <- 200000; nc <- 3; nt <- 100; na <- 5000
ni <- 5000; nb <- 2000; nc <- 3; nt <- 1; na <- 5000  # ~~~ for testing

# Call JAGS (ART 554 min) and check convergence
out3 <- jags(jags.data, inits, parameters, "model3.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out3)


# IPM4
# Write JAGS model file
cat(file="model4.txt", "
model {
  # Priors and linear models
  # Survival and recapture
  for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- sur[sex[i],x[i,t],t]
      p[i,t] <- pt[sex[i],x[i,t],k[t]]
    } #t
  } #i

  for (t in 1:(n.occasions-1)){
   logit(sur[1,1,t]) <- eps[1,t]
   logit(sur[1,2,t]) <- eps[2,t]
   logit(sur[1,3,t]) <- eps[3,t]
   logit(sur[2,1,t]) <- eps[4,t]
   logit(sur[2,2,t]) <- eps[5,t]
   logit(sur[2,3,t]) <- eps[6,t]
   eps[1:6,t] ~ dmnorm.vcov(mu.s[], Sigma[,])
  }

  # Recaptures
  for (s in 1:2){
    for (a in 1:3){
      for (t in 1:10){
        logit.pt[s,a,t] ~ dnorm(mu.p[s,a], tau.p[s,a])
        pt[s,a,t] <- ilogit(logit.pt[s,a,t])
      } #t
      pt[s,a,11] ~ dunif(0, 1)   # rate in 2015, incomplete capture
      pt[s,a,12] <- 0            # fix to 0 in years without capture
    } #a
  } #s

  mu.s[1] <- logit(mean.s[1,1])
  mu.s[2] <- logit(mean.s[1,2])
  mu.s[3] <- logit(mean.s[1,3])
  mu.s[4] <- logit(mean.s[2,1])
  mu.s[5] <- logit(mean.s[2,2])
  mu.s[6] <- logit(mean.s[2,3])

  for (j in 1:6){
    sd.s[j] ~ dunif(0.001, 10)  # Prior of the temporal SD
    Sigma[j,j] <- pow(sd.s[j], 2)
  }

  for (m in 1:5){
    for (j in (m+1):6){
      Sigma[m,j] <- zeta[m,j] * pow(Sigma[m,m] * Sigma[j,j], 0.5)
      Sigma[j,m] <- Sigma[m,j]
      zeta[m,j] ~ dunif(-1, 1)     # Priors for correlation coefficients
    } #j
  } #m

  for (s in 1:2){
    for (a in 1:3){
      mean.s[s,a] ~ dunif(0, 1)
      mean.p[s,a] ~ dunif(0, 1)
      mu.p[s,a] <- logit(mean.p[s,a])
      sigma.p[s,a] ~ dunif(0, 5)
      tau.p[s,a] <- pow(sigma.p[s,a],-2)
    } #a
  } #s

  # Productivity
  for (t in 1:n.occasions){
    rho[1,t] <- mean.rho[1]                  # of 2y females
    logit(rho[2,t]) <- mu.rho + eps.rho[t]   # of >2y females
    eps.rho[t] ~ dnorm(0, tau.rho)
  }
  # Priors for productivity of 2y females
  #mean.rho[1] <- 0
  #mean.rho[1] <- 0.99       # does not work if = 1 due to initial values
  #mean.rho[1] ~ dunif(0, 1)
  mean.rho[1] ~ dbeta(beta.priors[1], beta.priors[2])

  # Priors for productivity of >2y females
  mean.rho[2] ~ dunif(0,1)
  mu.rho <- logit(mean.rho[2])
  sigma.rho ~ dunif(0, 10)
  tau.rho <- pow(sigma.rho, -2)

  # Variance of the observation error
  sigma.y ~ dunif(0.5, 20)
  tau.y <- pow(sigma.y, -2)

  # Sex ratio of newborn individuals
  xi ~ dbeta(Jf, Jm)

  # Compute the proportion of bats present in the colony
  for (s in 1:2){
    for (a in 1:3){
      lkappa[s,a,1] ~ dnorm(mu.p[s,a], tau.p[s,a])
      kappa[s,a,1] <- ilogit(lkappa[s,a,1])
      for (t in 1:length(nrpres)){
        lkappa[s,a,nrpres[t]+1] ~ dnorm(mu.p[s,a], tau.p[s,a])
        kappa[s,a,nrpres[t]+1] <- ilogit(lkappa[s,a,nrpres[t]+1])
      } #t
      for (t in 1:length(rpres)){
        kappa[s,a,rpres[t]+1] <- pt[s,a,t]
      } #t
    } #a
  } #s

  # Population count data (state-space model)
  # Model for the initial population size: uniform priors
  N[1,1,1] ~ dcat(pinit1)
  N[1,2,1] ~ dcat(pinit1)
  N[1,3,1] ~ dcat(pinit2)
  N[2,1,1] ~ dcat(pinit1)
  N[2,2,1] ~ dcat(pinit1)
  N[2,3,1] ~ dcat(pinit2)

  # Process model: population sizes of females and males
  for(t in 2:n.occasions){
    N[1,1,t] ~ dbin(sur[1,1,t-1], NB[1,t-1])
    N[1,2,t] ~ dbin(sur[1,2,t-1], N[1,1,t-1])
    N[1,3,t] ~ dbin(sur[1,3,t-1], N[1,2,t-1] + N[1,3,t-1])
    N[2,1,t] ~ dbin(sur[2,1,t-1], NB[2,t-1])
    N[2,2,t] ~ dbin(sur[2,2,t-1], N[2,1,t-1])
    N[2,3,t] ~ dbin(sur[2,3,t-1], N[2,2,t-1] + N[2,3,t-1])
  }

  # Observation model
  for (t in 1:n.occasions){
    Ns[t] <- sum(N[,,t])       # total population size
    U[t] <- N[1,1,t] * kappa[1,1,t] + N[1,2,t] * kappa[1,2,t] + N[1,3,t] * kappa[1,3,t]  + N[2,1,t] * kappa[2,1,t]  + N[2,2,t] * kappa[2,2,t] + N[2,3,t] * kappa[2,3,t]
    C[t] ~ dnorm(U[t], tau.y)
  }

  # Capture-recapture data (state-space model)
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(phi[i,t-1] * z[i,t-1])
      # Observation process
      ch[i,t] ~ dbern(p[i,t-1] * z[i,t])
    } #t
  } #i

  # Newborn count data (binomial models)
  for (t in 1:n.occasions){
    G[1,t] ~ dbin(rho[1,t], N[1,2,t])   # Newborn by 2y females
    G[2,t] ~ dbin(rho[2,t], N[1,3,t])   # Newborn by >2y females
    J[t] ~ dsum(G[1,t], G[2,t])         # Total newborn
    NB[1,t] ~ dbin(xi, J[t])            # Newborn females
    NB[2,t] <- J[t] - NB[1,t]           # Newborn males
  }
}
")

# Initial values
Ninit <- array(NA, dim=c(2, 3, 29))
Ninit[1,3,] <- 60
G <- matrix(0, nrow=2, ncol=29)
G[2,] <- bats$J.count
inits <- function(){list(z=zInit(bats$ch), N=Ninit, G=G)}

# Parameters monitored
parameters <- c("mean.s", "Sigma", "mean.rho", "sigma.rho", "xi", "zeta", "mean.p",
    "sigma.p", "sigma.y", "sur", "pt", "rho", "kappa", "N", "U", "Ns", "NB", "G")

# MCMC settings
# ni <- 500000; nb <- 200000; nc <- 3; nt <- 100; na <- 5000
ni <- 5000; nb <- 2000; nc <- 3; nt <- 1; na <- 5000  # ~~~ for testing

# Call JAGS (ART 554 min) and check convergence
out4 <- jags(jags.data, inits, parameters, "model4.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out4)

save(out1, out2, out3, out4, file="BatResults.Rdata")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 13.6 Results
# ============

# ~~~~ Fig. 13.3 ~~~~
library(scales)
cl <- viridis_pal(option='E')(10)[c(9,6,1)]
years <- 1989:2017

# Survival
# Females
nyears <- length(years)-1
op <- par(mfrow=c(3,2), mar=c(1, 4, 2, 1))
d <- 0.2
plot(y=out1$mean$sur[1,1,], x=(1:nyears)+0.5, type="b", pch=16, ylim=c(0, 1), axes=FALSE,
    xlab=NA, ylab=expression(paste("Survival probability (", phi,")")), col=cl[1])
segments((1:nyears)+0.5, out1$q2.5$sur[1,1,], (1:nyears)+0.5, out1$q97.5$sur[1,1,], col=cl[1])
points(y=out1$mean$sur[1,2,], x=(1:nyears)+0.5-d, pch=16, type="b", col=cl[2])
segments((1:nyears)+0.5-d, out1$q2.5$sur[1,2,], (1:nyears)+0.5-d, out1$q97.5$sur[1,2,], col=cl[2])
points(y=out1$mean$sur[1,3,], x=(1:nyears)+0.5+d, pch=16, type="b", col=cl[3])
segments((1:nyears)+0.5+d, out1$q2.5$sur[1,3,], (1:nyears)+0.5+d, out1$q97.5$sur[1,3,], col=cl[3])
nyears <- length(years)
axis(1, at=1:nyears, tcl=-0.25, labels=NA)
axis(1, at=seq(1, nyears, by=2), tcl=-0.5, labels=NA)
axis(2, las=1)
legend("bottomleft", pch=rep(16, 3), col=cl, legend=c("Juvenile", "1y", "Adult (2y & >2y)"), bty="n")
mtext('Females', side=3, line=0.5)

# Males
nyears <- length(years)-1
par(mar=c(1, 4, 2, 1))
d <- 0.2
plot(y=out1$mean$sur[2,1,], x=(1:nyears)+0.5, type="b", pch=16, ylim=c(0, 1), axes=FALSE,
    xlab=NA, ylab=NA, col=cl[1])
segments((1:nyears)+0.5, out1$q2.5$sur[2,1,], (1:nyears)+0.5, out1$q97.5$sur[2,1,], col=cl[1])
points(y=out1$mean$sur[2,2,], x=(1:nyears)+0.5-d, pch=16, type="b", col=cl[2])
segments((1:nyears)+0.5-d, out1$q2.5$sur[2,2,], (1:nyears)+0.5-d, out1$q97.5$sur[2,2,], col=cl[2])
points(y=out1$mean$sur[2,3,], x=(1:nyears)+0.5+d, pch=16, type="b", col=cl[3])
segments((1:nyears)+0.5+d, out1$q2.5$sur[2,3,], (1:nyears)+0.5+d, out1$q97.5$sur[2,3,], col=cl[3])
nyears <- length(years)
axis(1, at=1:nyears, tcl=-0.25, labels=NA)
axis(1, at=seq(1, nyears, by=2), tcl=-0.5, labels=NA)
axis(2, las=1, labels=NA)
mtext('Males', side=3, line=0.5)


# Probability to be present at the colony
# Females
nyears <- length(years)
par(mar=c(1, 4, 2, 1))
d <- 0.2
plot(out1$mean$kappa[1,1,], type="b", pch=16, ylim=c(0, 1), axes=FALSE, xlab=NA,
    ylab=expression(paste("Probability to be present (", kappa, ")")), col=cl[1])
segments(1:nyears, out1$q2.5$kappa[1,1,], 1:nyears, out1$q97.5$kappa[1,1,], col=cl[1])
points(y=out1$mean$kappa[1,2,], x=(1:nyears)-d, pch=16, type="b", col=cl[2])
segments((1:nyears)-d, out1$q2.5$kappa[1,2,], (1:nyears)-d, out1$q97.5$kappa[1,2,], col=cl[2])
points(y=out1$mean$kappa[1,3,], x=(1:nyears)+d, pch=16, type="b", col=cl[3])
segments((1:nyears)+d, out1$q2.5$kappa[1,3,], (1:nyears)+d, out1$q97.5$kappa[1,3,], col=cl[3])
axis(1, at=1:nyears, tcl=-0.25, labels=NA)
axis(1, at=seq(1, nyears, by=2), tcl=-0.5, labels=NA)
axis(2, las=1)
legend("bottomleft", pch=rep(16, 3), col=cl, legend=c("1y", "2y", ">2y"), bty="n")

# Males
par(mar=c(1, 4, 2, 1))
plot(out1$mean$kappa[2,1,], type="b", pch=16, ylim=c(0, 1), axes=FALSE, xlab=NA, ylab=NA, col=cl[1])
segments(1:nyears, out1$q2.5$kappa[2,1,], 1:nyears, out1$q97.5$kappa[2,1,], col=cl[1])
points(y=out1$mean$kappa[2,2,], x=(1:nyears)-d, pch=16, type="b", col=cl[2])
segments((1:nyears)-d, out1$q2.5$kappa[2,2,], (1:nyears)-d, out1$q97.5$kappa[2,2,], col=cl[2])
points(y=out1$mean$kappa[2,3,], x=(1:nyears)+d, pch=16, type="b", col=cl[3])
segments((1:nyears)+d, out1$q2.5$kappa[2,3,], (1:nyears)+d, out1$q97.5$kappa[2,3,], col=cl[3])
axis(1, at=1:nyears, tcl=-0.25, labels=NA)
axis(1, at=seq(1, nyears, by=2), tcl=-0.5, labels=NA)
axis(2, las=1, labels=NA)

# Population sizes
# Females
d <- 0.15
par(mar=c(3, 4, 2, 1))
plot(y=out1$mean$N[1,1,], x=(1:nyears)-0.5*d, type="b", pch=16, ylim=c(0, 50), axes=FALSE,
    xlab=NA, ylab="Numbers", col=cl[1])
segments((1:nyears)-0.5*d, out1$q2.5$N[1,1,], (1:nyears)-0.5*d, out1$q97.5$N[1,1,], col=cl[1])
points(y=out1$mean$N[1,2,], x=(1:nyears)-1.5*d, pch=16, type="b", col=cl[2])
segments((1:nyears)-1.5*d, out1$q2.5$N[1,2,], (1:nyears)-1.5*d, out1$q97.5$N[1,2,], col=cl[2])
points(y=out1$mean$N[1,3,], x=(1:nyears)+0.5*d, pch=16, type="b", col=cl[3])
segments((1:nyears)+0.5*d, out1$q2.5$N[1,3,], (1:nyears)+0.5*d, out1$q97.5$N[1,3,], col=cl[3])
axis(1, at=1:nyears, tcl=-0.25, labels=NA)
axis(1, at=seq(1, nyears, by=2), tcl=-0.5, labels=years[seq(1, nyears, by=2)])
axis(2, las=1)
legend("topleft", pch=rep(16, 3), col=cl, legend=c("1y", "2y", ">2y"), bty="n")

# Males
par(mar=c(3, 4, 2, 1))
plot(y=out1$mean$N[2,1,], x=(1:nyears)-0.5*d, type="b", pch=16, ylim=c(0, 50), axes=FALSE,
    xlab=NA, ylab=NA, col=cl[1])
segments((1:nyears)-0.5*d, out1$q2.5$N[2,1,], (1:nyears)-0.5*d, out1$q97.5$N[2,1,], col=cl[1])
points(y=out1$mean$N[2,2,], x=(1:nyears)-1.5*d, pch=16, type="b", col=cl[2])
segments((1:nyears)-1.5*d, out1$q2.5$N[2,2,], (1:nyears)-1.5*d, out1$q97.5$N[2,2,], col=cl[2])
points(y=out1$mean$N[2,3,], x=(1:nyears)+0.5*d, pch=16, type="b", col=cl[3])
segments((1:nyears)+0.5*d, out1$q2.5$N[2,3,], (1:nyears)+0.5*d, out1$q97.5$N[2,3,], col=cl[3])
axis(1, at=1:nyears, tcl=-0.25, labels=NA)
axis(1, at=seq(1, nyears, by=2), tcl=-0.5, labels=years[seq(1, nyears, by=2)])
axis(2, las=1, labels=NA)
par(op)

# ~~~~ Fig. 13.4 ~~~~
library(scales)
cl <- viridis_pal(option='E')(10)[c(1,6)]

years <- 1989:2017
nyears <- length(years)

# Productivity
op <- par(mfrow=c(2, 1), mar=c(1, 4, 2, 1))
plot(out1$mean$rho[2,], type="b", pch=16, ylim=c(0.7, 1), axes=FALSE, xlab=NA,
    ylab=expression(paste('Productivity  (', rho[2], ')')))
segments(1:nyears, out1$q2.5$rho[2,], 1:nyears, out1$q97.5$rho[2,])
axis(1, at=1:nyears, tcl=-0.25, labels=NA)
axis(1, at=seq(1, nyears, by=2), tcl=-0.5, labels=NA)
axis(2, las=1, at=c(0.75, 0.85, 0.95), labels=NA, tcl=-0.25)
axis(2, las=1, at=c(0.7, 0.8, 0.9, 1), labels=c("0.7", "0.8", "0.9", "1.0"), tcl=-0.5)

# Total population sizes
par(mar=c(3, 4, 2, 1))
d <- 0.25
plot(out1$mean$U, type="b", pch=16, ylim=c(0, 100), axes=FALSE, col=cl[1], xlab=NA, ylab="Number")
segments(1:nyears, out1$q2.5$U, 1:nyears, out1$q97.5$U, col=cl[1])
axis(1, at=1:nyears, tcl=-0.25, labels=NA)
axis(1, at=seq(1, nyears, by=2), tcl=-0.5, labels=years[seq(1, nyears, by=2)])
axis(2, las=1)
points(y=out1$mean$Ns, x=(1:nyears)+d, pch=16, type="b", col=cl[2])
segments((1:nyears)+d, out1$q2.5$Ns, (1:nyears)+d, out1$q97.5$Ns, col=cl[2])
points(bats$A.count, pch=18, cex=1, col=alpha("orange", 0.7))
legend("bottomright", pch=c(16, 16, 18), col=c(rev(cl), alpha("orange", 0.7)),
    legend=c("Superpopulation", "Colony size", "Population count"), bty="n")
par(op)

# 13.7 Prior sensitivity analysis
# ===============================

# ~~~~ Fig 13.5 ~~~~
library(scales)
cl <- viridis_pal(option='E')(20)[c(18,11,5,1)]
lw <- 2

op <- par(mfrow=c(3,2), mar=c(1, 4, 2, 1))
plot(density(out1$sims.list$mean.s[,1,1]), ylim=c(0, 10), xlim=c(0.2, 1), col=cl[1],
    axes=FALSE, ylab="Density", xlab=NA, main=NA, lwd=lw)
mtext(expression(phi[1]), side=1, line=1)
mtext('Females', side=3, line=0.5)
axis(2, las=1)
axis(1, labels=NA)
lines(density(out2$sims.list$mean.s[,1,1]), col=cl[2], lwd=lw)
lines(density(out3$sims.list$mean.s[,1,1]), col=cl[3], lwd=lw)
lines(density(out4$sims.list$mean.s[,1,1]), col=cl[4], lwd=lw)

par(mar=c(1, 1, 2, 1))
plot(density(out1$sims.list$mean.s[,2,1]), ylim=c(0, 10), xlim=c(0.2, 1), col=cl[1],
    axes=FALSE, ylab="Density", xlab=NA, main=NA, lwd=lw)
mtext(expression(phi[1]), side=1, line=1)
mtext('Males', side=3, line=0.5)
axis(2, las=1, labels=NA)
axis(1, labels=NA)
lines(density(out2$sims.list$mean.s[,2,1]), col=cl[2], lwd=lw)
lines(density(out3$sims.list$mean.s[,2,1]), col=cl[3], lwd=lw)
lines(density(out4$sims.list$mean.s[,2,1]), col=cl[4], lwd=lw)

par(mar=c(1, 4, 2, 1))
plot(density(out1$sims.list$mean.s[,1,2]), ylim=c(0, 10), xlim=c(0.2, 1), col=cl[1],
    axes=FALSE, ylab="Density", xlab=NA, main=NA, lwd=lw)
mtext(expression(phi[2]), side=1, line=1)
axis(2, las=1)
axis(1, labels=NA)
lines(density(out2$sims.list$mean.s[,1,2]), col=cl[2], lwd=lw)
lines(density(out3$sims.list$mean.s[,1,2]), col=cl[3], lwd=lw)
lines(density(out4$sims.list$mean.s[,1,2]), col=cl[4], lwd=lw)

par(mar=c(1, 1, 2, 1))
plot(density(out1$sims.list$mean.s[,2,2]), ylim=c(0, 10), xlim=c(0.2, 1), col=cl[1],
    axes=FALSE, ylab="Density", xlab=NA, main=NA, lwd=lw)
mtext(expression(phi[2]), side=1, line=1)
axis(2, las=1, labels=NA)
axis(1, labels=NA)
lines(density(out2$sims.list$mean.s[,2,2]), col=cl[2], lwd=lw)
lines(density(out3$sims.list$mean.s[,2,2]), col=cl[3], lwd=lw)
lines(density(out4$sims.list$mean.s[,2,2]), col=cl[4], lwd=lw)

par(mar=c(4, 4, 2, 1))
plot(density(out1$sims.list$mean.s[,1,3]), ylim=c(0, 20), xlim=c(0.2, 1), col=cl[1],
    axes=FALSE, ylab="Density", xlab=NA, main=NA, lwd=lw)
mtext(expression(phi[3]), side=1, line=2.5)
axis(2, las=1)
axis(1)
lines(density(out2$sims.list$mean.s[,1,3]), col=cl[2], lwd=lw)
lines(density(out3$sims.list$mean.s[,1,3]), col=cl[3], lwd=lw)
lines(density(out4$sims.list$mean.s[,1,3]), col=cl[4], lwd=lw)

par(mar=c(4, 1, 2, 1))
plot(density(out1$sims.list$mean.s[,2,3]), ylim=c(0, 20), xlim=c(0.2, 1), col=cl[1],
    axes=FALSE, ylab="Density", xlab=NA, main=NA, lwd=lw)
mtext(expression(phi[3]), side=1, line=2.5)
axis(2, las=1, labels=NA)
axis(1)
lines(density(out2$sims.list$mean.s[,2,3]), col=cl[2], lwd=lw)
lines(density(out3$sims.list$mean.s[,2,3]), col=cl[3], lwd=lw)
lines(density(out4$sims.list$mean.s[,2,3]), col=cl[4], lwd=lw)
legend("topleft", lwd = rep(lw,4), col = cl, bty = "n",
    legend = c(expression(paste('IPM'[1],': ',rho[1],' = 0')),
        expression(paste('IPM'[2],': ',rho[1],' = 0.99')),
        expression(paste('IPM'[3],': ',rho[1],' ~ U(0,1)')),
        expression(paste('IPM'[4],': ',rho[1],' ~ Beta(3.9,5.8)'))))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Fig. 13.6 ~~~~
library(scales)
cl <- viridis_pal(option='E')(20)[c(18,11,5,1)]
lw <- 2

op <- par(mfrow=c(2,2), mar=c(3.5, 4, 2, 1))
plot(density(out1$sims.list$mean.rho[,2]), ylim=c(0, 15), xlim=c(0.6, 1), col=cl[1],
    axes=FALSE, ylab="Density", xlab=NA, main=NA, lwd=lw)
mtext(expression(rho[2]), side=1, line=2.5)
axis(2, las=1)
axis(1)
lines(density(out2$sims.list$mean.rho[,2]), col=cl[2], lwd=lw)
lines(density(out3$sims.list$mean.rho[,2]), col=cl[3], lwd=lw)
lines(density(out4$sims.list$mean.rho[,2]), col=cl[4], lwd=lw)
legend("topleft", lwd=rep(1, 4), col=cl, bty="n",
    legend=c(expression(paste('IPM'[1],': ',rho[1],' = 0')),
        expression(paste('IPM'[2],': ',rho[1],' = 0.99')),
        expression(paste('IPM'[3],': ',rho[1],' ~ U(0,1)')),
        expression(paste('IPM'[4],': ',rho[1],' ~ Beta(3.9,5.8)'))))

par(mar=c(3.5, 4, 2, 1))
plot(density(out1$sims.list$U[,29]), ylim=c(0, 0.3), xlim=c(35, 65), col=cl[1],
    axes=FALSE, ylab="", xlab=NA, main=NA, lwd=lw)
mtext(expression(italic(U)[29]), side=1, line=2.5)
axis(2, labels=NA, tcl=-0.25)
axis(2, las=1, at=c(0, 0.1, 0.2, 0.3), labels=c(0, 0.1, 0.2, 0.3), tcl=-0.5)
axis(1)
lines(density(out2$sims.list$U[,29]), col=cl[2], lwd=lw)
lines(density(out3$sims.list$U[,29]), col=cl[3], lwd=lw)
lines(density(out4$sims.list$U[,29]), col=cl[4], lwd=lw)

par(mar=c(3.5, 4, 2, 1))
plot(density(out1$sims.list$sigma.y), ylim=c(0, 0.3), xlim=c(0, 12), col=cl[1],
    axes=FALSE, ylab="Density", xlab=NA, main=NA, lwd=lw)
mtext(expression(sigma), side=1, line=2.5)
axis(2, labels=NA, tcl=-0.25)
axis(2, las=1, at=c(0, 0.1, 0.2, 0.3), labels=c(0, 0.1, 0.2, 0.3), tcl=-0.5)
axis(1)
lines(density(out2$sims.list$sigma.y), col=cl[2], lwd=lw)
lines(density(out3$sims.list$sigma.y), col=cl[3], lwd=lw)
lines(density(out4$sims.list$sigma.y), col=cl[4], lwd=lw)

ipm1 <- table(out1$sims.list$N[,1,3,29])
ipm2 <- table(out2$sims.list$N[,1,3,29])
ipm3 <- table(out3$sims.list$N[,1,3,29])
ipm4 <- table(out4$sims.list$N[,1,3,29])
par(mar=c(3.5, 4, 2, 1))
plot(ipm1, type="l", ylim=c(0, 3500), xlim=c(18, 42), col=cl[1], axes=FALSE,
    ylab="", xlab=NA, main=NA, lwd=lw, frame.plot=FALSE)
mtext(expression(italic(N)['3,29']), side=1, line=2.5)
segments(31, 0, 31, 3451, col=cl[1])
axis(2, labels=NA, tcl=-0.25)
axis(2, las=1, at=c(0, 1000, 2000, 3000), labels=c(0, 1000, 2000, 3000), tcl=-0.5)
axis(1)
lines(ipm2, type="l", col=cl[2], lwd=lw)
lines(ipm3, type="l", col=cl[3], lwd=lw)
lines(ipm4, type="l", col=cl[4], lwd=lw)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Fig. 13.7 ~~~~
library(scales)
library(wiqid)
cl <- viridis_pal(option='E')(20)[c(11,1)]
lw <- 2

op <- par(mar=c(4, 4, 2, 1))
plot(densityFolded (out3$sims.list$mean.rho[,1]), ylim=c(0, 3), xlim=c(0, 1), col=cl[1],
    axes=FALSE, ylab="Density", xlab=NA, main=NA, lwd=lw)
mtext(expression(rho[1]), side=1, line=2.5)
axis(2, las=1)
axis(1)
lines(densityFolded (out4$sims.list$mean.rho[,1]), col=cl[2], lwd=lw)
lines(densityFolded (runif(100000, 0, 1)), col=cl[1], lwd=lw, lty=2)
lines(densityFolded(rbeta2(100000, 0.4, 0.15)), col=cl[2], lwd=lw, lty=2)
legend("topright", lwd = rep(2,4), lty = c(1,2,1,2), col =cl[c(1,1,2,2)], bty = "n",
    legend = c(expression('Posterior under IPM'[3]),
        expression(paste('Prior in IPM'[3],' prior')),
        expression('Posterior under IPM'[4]),
        expression(paste('Prior in IPM'[4],' prior'))))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
