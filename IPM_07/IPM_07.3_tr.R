# Schaub & Kéry (2022) Integrated Population Models
# Chapter 7 : Assessment of integrated population models
# ------------------------------------------------------
# Code from final MS.

# Run time for test script 10 mins, full run 19 hrs

library(IPMbook) ; library(jagsUI)

# This script uses the files "model8.txt" and "model9.txt" created
#   in section 7.2.3


# 7.3 Under- and overfitting
# ==========================

# ~~~ Code for the simulations and the figures of section 7.3 ~~~

# Analysing models

# IPMs

# 1. all 3 data sets
# A. Model with constant f(.) and sa(.): model8.txt (same as before)

# B. Model with f(.) and time-dependent sa(t)
cat(file="model20.txt", "
model {
  # Priors and constraints
  mean.sj ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)
  mean.sa <- pow(prod(sa), 1 /(n.occasions-1))     # mean adult survival

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] ~ dunif(0, 1)
    p[t] <- mean.p
  }

  for (t in 1:n.occasions){
    f[t] <- mean.f
  }

  sigma.obs ~ dunif(0.5, upper.sigma.obs)
  tau.obs <- pow(sigma.obs, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dcat(pNinit)
  N[2,1] ~ dcat(pNinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] ~ dpois((N[1,t] + N[2,t]) * f[t] / 2 * sj[t])
    N[2,t+1] ~ dbin(sa[t], N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau.obs)
  }

  # Assessing the fit of the state-space model
  for (t in 1:n.occasions){
    C.exp[t] <- N[1,t] + N[2,t]                    # Expected counts
    Dssm.obs[t] <- abs((C[t] - C.exp[t]) / C[t])   # Discrepancy measure

    C.rep[t] ~ dnorm(N[1,t] + N[2,t], tau.obs)     # Generate replicate data
    Dssm.rep[t] <- abs((C.rep[t] - C.exp[t]) / C.rep[t]) # Discrepancy measure
  }
  Dmape.obs <- sum(Dssm.obs)
  Dmape.rep <- sum(Dssm.rep)

  # Productivity data (Poisson regression model)
  for (t in 1:n.occasions){
    J[t] ~ dpois(B[t] * f[t])

    J.exp[t] <- B[t] * f[t]              # Expected data
    D.obs[t] <- J[t] * log(J[t]/J.exp[t]) - (J[t] - J.exp[t])

    J.rep[t] ~ dpois(B[t] * f[t])
    D.rep[t] <- J.rep[t] * log(J.rep[t]/J.exp[t]) - (J.rep[t] - J.exp[t])
  }
  Dd.obs <- sum(D.obs)
  Dd.rep <- sum(D.rep)

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

  # Assessing the fit of the capture-recapture model (Freeman-Tukey)
  for (t in 1:(n.occasions-1)){
    marr.j.rep[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])   # Generate replicate data
    marr.a.rep[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
    for (j in 1:n.occasions){
      marr.j.exp[t,j] <- pr.j[t,j] * rel.j[t]  # Expected values
      marr.a.exp[t,j] <- pr.a[t,j] * rel.a[t]  # Expected values
      Dcjs.obs[t,j] <- pow(pow(marr.j[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.obs[t+n.occasions-1,j] <- pow(pow(marr.a[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)

      Dcjs.rep[t,j] <- pow(pow(marr.j.rep[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.rep[t+n.occasions-1,j] <- pow(pow(marr.a.rep[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)
    } #j
  } #t
  DFT.obs <- sum(Dcjs.obs)
  DFT.rep <- sum(Dcjs.rep)

  # Derived parameters
  # Mean population growth rate (geometric mean)
  geom.rate <- pow(Ntot[n.occasions] / Ntot[1], 1 / (n.occasions-1))
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
  }
}
")

# C. Model with time-dependent f(t) and constant sa(.)
cat(file="model21.txt", "
model {
  # Priors and constraints
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f <- mean(f)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }

  for (t in 1:n.occasions){
    f[t] ~ dunif(0, 10)
  }

  sigma.obs ~ dunif(0.5, upper.sigma.obs)
  tau.obs <- pow(sigma.obs, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dcat(pNinit)
  N[2,1] ~ dcat(pNinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] ~ dpois((N[1,t] + N[2,t]) * f[t] / 2 * sj[t])
    N[2,t+1] ~ dbin(sa[t], N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau.obs)
  }

  # Assessing the fit of the state-space model
  for (t in 1:n.occasions){
    C.exp[t] <- N[1,t] + N[2,t]                    # Expected counts
    Dssm.obs[t] <- abs((C[t] - C.exp[t]) / C[t])   # Discrepancy measure

    C.rep[t] ~ dnorm(N[1,t] + N[2,t], tau.obs)     # Generate replicate data
    Dssm.rep[t] <- abs((C.rep[t] - C.exp[t]) / C.rep[t]) # Discrepancy measure
  }
  Dmape.obs <- sum(Dssm.obs)
  Dmape.rep <- sum(Dssm.rep)

  # Productivity data (Poisson regression model)
  for (t in 1:n.occasions){
    J[t] ~ dpois(B[t] * f[t])

    J.exp[t] <- B[t] * f[t]              # Expected data
    D.obs[t] <- J[t] * log(J[t]/J.exp[t]) - (J[t] - J.exp[t])

    J.rep[t] ~ dpois(B[t] * f[t])
    D.rep[t] <- J.rep[t] * log(J.rep[t]/J.exp[t]) - (J.rep[t] - J.exp[t])
  }
  Dd.obs <- sum(D.obs)
  Dd.rep <- sum(D.rep)

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

  # Assessing the fit of the capture-recapture model (Freeman-Tukey)
  for (t in 1:(n.occasions-1)){
    marr.j.rep[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])   # Generate replicate data
    marr.a.rep[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
    for (j in 1:n.occasions){
      marr.j.exp[t,j] <- pr.j[t,j] * rel.j[t]  # Expected values
      marr.a.exp[t,j] <- pr.a[t,j] * rel.a[t]  # Expected values
      Dcjs.obs[t,j] <- pow(pow(marr.j[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.obs[t+n.occasions-1,j] <- pow(pow(marr.a[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)

      Dcjs.rep[t,j] <- pow(pow(marr.j.rep[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.rep[t+n.occasions-1,j] <- pow(pow(marr.a.rep[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)
    } #j
  } #t
  DFT.obs <- sum(Dcjs.obs)
  DFT.rep <- sum(Dcjs.rep)

  # Derived parameters
  # Mean population growth rate (geometric mean)
  geom.rate <- pow(Ntot[n.occasions] / Ntot[1], 1 / (n.occasions-1))
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
  }
}
")

# D. Model with time-dependent f(t) and sa(t)
cat(file="model22.txt", "
model {
  # Priors and constraints
  mean.sj ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.sa <- pow(prod(sa), 1 /(n.occasions-1))     # mean adult survival
  mean.f <- mean(f)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] ~ dunif(0, 1)
    p[t] <- mean.p
  }

  for (t in 1:n.occasions){
    f[t] ~ dunif(0, 10)
  }

  sigma.obs ~ dunif(0.5, upper.sigma.obs)
  tau.obs <- pow(sigma.obs, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dcat(pNinit)
  N[2,1] ~ dcat(pNinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] ~ dpois((N[1,t] + N[2,t]) * f[t] / 2 * sj[t])
    N[2,t+1] ~ dbin(sa[t], N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau.obs)
  }

  # Assessing the fit of the state-space model
  for (t in 1:n.occasions){
    C.exp[t] <- N[1,t] + N[2,t]                    # Expected counts
    Dssm.obs[t] <- abs((C[t] - C.exp[t]) / C[t])   # Discrepancy measure

    C.rep[t] ~ dnorm(N[1,t] + N[2,t], tau.obs)     # Generate replicate data
    Dssm.rep[t] <- abs((C.rep[t] - C.exp[t]) / C.rep[t]) # Discrepancy measure
  }
  Dmape.obs <- sum(Dssm.obs)
  Dmape.rep <- sum(Dssm.rep)

  # Productivity data (Poisson regression model)
  for (t in 1:n.occasions){
    J[t] ~ dpois(B[t] * f[t])

    J.exp[t] <- B[t] * f[t]              # Expected data
    D.obs[t] <- J[t] * log(J[t]/J.exp[t]) - (J[t] - J.exp[t])

    J.rep[t] ~ dpois(B[t] * f[t])
    D.rep[t] <- J.rep[t] * log(J.rep[t]/J.exp[t]) - (J.rep[t] - J.exp[t])
  }
  Dd.obs <- sum(D.obs)
  Dd.rep <- sum(D.rep)

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

  # Assessing the fit of the capture-recapture model (Freeman-Tukey)
  for (t in 1:(n.occasions-1)){
    marr.j.rep[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])   # Generate replicate data
    marr.a.rep[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
    for (j in 1:n.occasions){
      marr.j.exp[t,j] <- pr.j[t,j] * rel.j[t]  # Expected values
      marr.a.exp[t,j] <- pr.a[t,j] * rel.a[t]  # Expected values
      Dcjs.obs[t,j] <- pow(pow(marr.j[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.obs[t+n.occasions-1,j] <- pow(pow(marr.a[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)

      Dcjs.rep[t,j] <- pow(pow(marr.j.rep[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.rep[t+n.occasions-1,j] <- pow(pow(marr.a.rep[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)
    } #j
  } #t
  DFT.obs <- sum(Dcjs.obs)
  DFT.rep <- sum(Dcjs.rep)

  # Derived parameters
  # Mean population growth rate (geometric mean)
  geom.rate <- pow(Ntot[n.occasions] / Ntot[1], 1 / (n.occasions-1))
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
  }
}
")

# 2. IPM with counts and CMR only
# A. Model with constant f(.) and sa(.): model9.txt (same is before)

# B. Model with constant f(.) and time-dependent sa(t)
cat(file="model23.txt", "
model {
  # Priors and constraints
  mean.sj ~ dunif(0, 1)
  mean.sa <- pow(prod(sa), 1 /(n.occasions-1))     # mean adult survival
  mean.p ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] ~ dunif(0, 1)
    p[t] <- mean.p
  }

  for (t in 1:n.occasions){
    f[t] <- mean.f
  }

  sigma.obs ~ dunif(0.5, upper.sigma.obs)
  tau.obs <- pow(sigma.obs, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dcat(pNinit)
  N[2,1] ~ dcat(pNinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] ~ dpois((N[1,t] + N[2,t]) * f[t] / 2 * sj[t])
    N[2,t+1] ~ dbin(sa[t], N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau.obs)
  }

  # Assessing the fit of the state-space model
  for (t in 1:n.occasions){
    C.exp[t] <- N[1,t] + N[2,t]                    # Expected counts
    Dssm.obs[t] <- abs((C[t] - C.exp[t]) / C[t])   # Discrepancy measure

    C.rep[t] ~ dnorm(N[1,t] + N[2,t], tau.obs)     # Generate replicate data
    Dssm.rep[t] <- abs((C.rep[t] - C.exp[t]) / C.rep[t]) # Discrepancy measure
  }
  Dmape.obs <- sum(Dssm.obs)
  Dmape.rep <- sum(Dssm.rep)

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

  # Assessing the fit of the capture-recapture model (Freeman-Tukey)
  for (t in 1:(n.occasions-1)){
    marr.j.rep[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])   # Generate replicate data
    marr.a.rep[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
    for (j in 1:n.occasions){
      marr.j.exp[t,j] <- pr.j[t,j] * rel.j[t]  # Expected values
      marr.a.exp[t,j] <- pr.a[t,j] * rel.a[t]  # Expected values
      Dcjs.obs[t,j] <- pow(pow(marr.j[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.obs[t+n.occasions-1,j] <- pow(pow(marr.a[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)

      Dcjs.rep[t,j] <- pow(pow(marr.j.rep[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.rep[t+n.occasions-1,j] <- pow(pow(marr.a.rep[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)
    } #j
  } #t
  DFT.obs <- sum(Dcjs.obs)
  DFT.rep <- sum(Dcjs.rep)

  # Derived parameters
  # Mean population growth rate (geometric mean)
  geom.rate <- pow(Ntot[n.occasions] / Ntot[1], 1 / (n.occasions-1))
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
  }
}
")

# C. Model with time-dependent f(t) and constant sa(.)
cat(file="model24.txt", "
model {
  # Priors and constraints
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f <- mean(f)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }

  for (t in 1:n.occasions){
    f[t] ~ dunif(0, 10)
  }

  sigma.obs ~ dunif(0.5, upper.sigma.obs)
  tau.obs <- pow(sigma.obs, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dcat(pNinit)
  N[2,1] ~ dcat(pNinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] ~ dpois((N[1,t] + N[2,t]) * f[t] / 2 * sj[t])
    N[2,t+1] ~ dbin(sa[t], N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau.obs)
  }

  # Assessing the fit of the state-space model
  for (t in 1:n.occasions){
    C.exp[t] <- N[1,t] + N[2,t]                    # Expected counts
    Dssm.obs[t] <- abs((C[t] - C.exp[t]) / C[t])   # Discrepancy measure

    C.rep[t] ~ dnorm(N[1,t] + N[2,t], tau.obs)     # Generate replicate data
    Dssm.rep[t] <- abs((C.rep[t] - C.exp[t]) / C.rep[t]) # Discrepancy measure
  }
  Dmape.obs <- sum(Dssm.obs)
  Dmape.rep <- sum(Dssm.rep)

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

  # Assessing the fit of the capture-recapture model (Freeman-Tukey)
  for (t in 1:(n.occasions-1)){
    marr.j.rep[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])   # Generate replicate data
    marr.a.rep[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
    for (j in 1:n.occasions){
      marr.j.exp[t,j] <- pr.j[t,j] * rel.j[t]  # Expected values
      marr.a.exp[t,j] <- pr.a[t,j] * rel.a[t]  # Expected values
      Dcjs.obs[t,j] <- pow(pow(marr.j[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.obs[t+n.occasions-1,j] <- pow(pow(marr.a[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)

      Dcjs.rep[t,j] <- pow(pow(marr.j.rep[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.rep[t+n.occasions-1,j] <- pow(pow(marr.a.rep[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)
    } #j
  } #t
  DFT.obs <- sum(Dcjs.obs)
  DFT.rep <- sum(Dcjs.rep)

  # Derived parameters
  # Mean population growth rate (geometric mean)
  geom.rate <- pow(Ntot[n.occasions] / Ntot[1], 1 / (n.occasions-1))
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
  }
}
")

# D. Model with time-dependent f(t) and sa(t)
cat(file="model25.txt", "
model {
  # Priors and constraints
  mean.sj ~ dunif(0, 1)
  mean.sa <- pow(prod(sa), 1 /(n.occasions-1))     # mean adult survival
  mean.p ~ dunif(0, 1)
  mean.f <- mean(f)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] ~ dunif(0, 1)
    p[t] <- mean.p
  }

  for (t in 1:n.occasions){
    f[t] ~ dunif(0, 10)
  }

  sigma.obs ~ dunif(0.5, upper.sigma.obs)
  tau.obs <- pow(sigma.obs, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dcat(pNinit)
  N[2,1] ~ dcat(pNinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] ~ dpois((N[1,t] + N[2,t]) * f[t] / 2 * sj[t])
    N[2,t+1] ~ dbin(sa[t], N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau.obs)
  }

  # Assessing the fit of the state-space model
  for (t in 1:n.occasions){
    C.exp[t] <- N[1,t] + N[2,t]                    # Expected counts
    Dssm.obs[t] <- abs((C[t] - C.exp[t]) / C[t])   # Discrepancy measure

    C.rep[t] ~ dnorm(N[1,t] + N[2,t], tau.obs)     # Generate replicate data
    Dssm.rep[t] <- abs((C.rep[t] - C.exp[t]) / C.rep[t]) # Discrepancy measure
  }
  Dmape.obs <- sum(Dssm.obs)
  Dmape.rep <- sum(Dssm.rep)

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

  # Assessing the fit of the capture-recapture model (Freeman-Tukey)
  for (t in 1:(n.occasions-1)){
    marr.j.rep[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])   # Generate replicate data
    marr.a.rep[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
    for (j in 1:n.occasions){
      marr.j.exp[t,j] <- pr.j[t,j] * rel.j[t]  # Expected values
      marr.a.exp[t,j] <- pr.a[t,j] * rel.a[t]  # Expected values
      Dcjs.obs[t,j] <- pow(pow(marr.j[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.obs[t+n.occasions-1,j] <- pow(pow(marr.a[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)

      Dcjs.rep[t,j] <- pow(pow(marr.j.rep[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.rep[t+n.occasions-1,j] <- pow(pow(marr.a.rep[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)
    } #j
  } #t
  DFT.obs <- sum(Dcjs.obs)
  DFT.rep <- sum(Dcjs.rep)

  # Derived parameters
  # Mean population growth rate (geometric mean)
  geom.rate <- pow(Ntot[n.occasions] / Ntot[1], 1 / (n.occasions-1))
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
  }
}
")


# Put model names into 2 vectors
modfile3 <- c("model8.txt", "model20.txt", "model21.txt", "model22.txt")
modfile2 <- c("model9.txt", "model23.txt", "model24.txt", "model25.txt")


# Simulations

# Number of simulations
# nsim <- 600  # 15.5 hrs
nsim <- 6  # ~~~ for testing

# Number of years
T <- 10

# Observation error for the population survey
sigma <- 10

# Capture and recapture probabilities
cap <- 0.4                        # initial capture probability
recap <- 0.6                      # recapture probability

# Probability to find a brood whose reproductive ouput is recorded
pprod <- 0.5

# Initial population size per age class
Ni <- c(50, 50)

# Demographic rates
# Survival
sj <- 0.3
sa <- 0.55                          # when constant
sat <- seq(0.45, 0.65, by=0.01)   # possible values of adult survival when time dependent

# Productivity
f <- 3.1                            # when constant
ft <- seq(2.1, 4.1, by=0.05)      # possible values of productivity when time dependent


# Define arrays to store results
# res3 contains the summary of each of the 16 model fits for 3 data sets, res2 for 2 data sets
res3 <- array(NA, dim=c(62, 11, nsim, 4, 4)) # vars x stats x sims x data x model
res2 <- array(NA, dim=c(60, 11, nsim, 4, 4))
# bP3 contains the GOF p-values for 3 data set fits, bP2 for 2 data sets
bP3 <- array(NA, c(nsim, 3, 4, 4))  # sims x stats x data x model
bP2 <- array(NA, c(nsim, 2, 4, 4))
# sa.true and f.true true adult survival and fecundity
sa.true <- matrix(NA, nrow=nsim, ncol=T-1)
f.true <- matrix(NA, nrow=nsim, ncol=T)

# Initial values
inits <- function(){
  N <- matrix(NA, nrow=2, ncol=T)
  N[1,] <- round(runif(T, 40, 60))
  N[2,] <- round(runif(T, 40, 60))
  ini <- list(mean.sj=runif(1, 0, 0.5), N=N)
  return(ini)
}

# Parameters monitored
parameters1 <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "sa", "f", "N",
    "sigma.obs", "geom.rate", "Ntot", "Dmape.obs", "Dmape.rep",
    "DFT.obs", "DFT.rep", "Dd.obs", "Dd.rep")
parameters2 <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "sa", "f", "N",
    "sigma.obs", "geom.rate", "Ntot", "Dmape.obs", "Dmape.rep", "DFT.obs", "DFT.rep")

# MCMC settings
ni <- 4000; nb <- 2000; nc <- 3; nt <- 4; na <- 2000


# Start simulations
system.time(
for (s in 1:nsim){
  set.seed(s)

 # Prepare 8 data bundles for JAGS: 4 scenarios x 2 data sets (all 3 vs only 2 (no productivity))
 # -----------------------------------------------------------------------
 # Put them into 2 lists, one for all-3, one for just-2
 jagsdata3 <- jagsdata2 <- vector("list", 4)

  # Scenario 1: all rates constant
  # Create populations
  pop1 <- simPop(Ni=Ni, phi=c(sj, sa), f=f, nYears=T)
  pop2 <- simPop(Ni=Ni, phi=c(sj, sa), f=f, nYears=T)
  pop3 <- simPop(Ni=Ni, phi=c(sj, sa), f=f, nYears=T)

  # Create capture histories & m-arrays from population 1
  ch <- simCapHist(state=pop1$state, cap=cap, recap=recap, maxAge=2)
  marr <- marrayAge(ch$ch, ch$age)

  # Create productivity data from population 2
  pro <- simProd(reprod=pop2$reprod, pInclude=pprod)
  # Aggregate productivity data to make the model run faster
  J <- pro$prod.agg[,1]
  B <- pro$prod.agg[,2]

  # Create the population survey data from population 3
  count <- simCountNorm(N=pop3$totB, sigma=sigma)$count

  # Bundle data
  # All 3 data sets
  jagsdata3[[1]] <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), J=J, B=B,
      C=count, pNinit=dUnif(1, 100), upper.sigma.obs=100)

  # No productivity data
  jagsdata2[[1]] <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]),
      C=count, pNinit=dUnif(1, 100), upper.sigma.obs=100)


  # Scenario 2: adult survival time-dependent
  # Create populations
  sas <- sample(sat, (T-1), replace=TRUE)
  pop1 <- simPop(Ni=Ni, phi=matrix(c(rep(sj,T-1), sas), nrow=2, byrow=TRUE), f=f, nYears=T)
  pop2 <- simPop(Ni=Ni, phi=matrix(c(rep(sj,T-1), sas), nrow=2, byrow=TRUE), f=f, nYears=T)
  pop3 <- simPop(Ni=Ni, phi=matrix(c(rep(sj,T-1), sas), nrow=2, byrow=TRUE), f=f, nYears=T)

  # Create capture histories & m-arrays from population 1
  ch <- simCapHist(state=pop1$state, cap=cap, recap=recap, maxAge=2)
  marr <- marrayAge(ch$ch, ch$age)

  # Create productivity data from population 2
  pro <- simProd(reprod=pop2$reprod, pInclude=pprod)
  # Aggregate productivity data to make the model run faster
  J <- pro$prod.agg[,1]
  B <- pro$prod.agg[,2]

  # Create the population survey data from population 3
  count <- simCountNorm(N=pop3$totB, sigma=sigma)$count

  # Bundle data
  # All 3 data sets
  jagsdata3[[2]] <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), J=J, B=B,
      C=count, pNinit=dUnif(1, 100), upper.sigma.obs=100)

  # No productivity data
  jagsdata2[[2]] <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]),
      C=count, pNinit=dUnif(1, 100), upper.sigma.obs=100)


  # Scenario 3: productivity time-dependent
  # Create populations
  fs <- sample(ft, T, replace=TRUE)
  pop1 <- simPop(Ni=Ni, phi=c(sj, sa), f=matrix(c(fs, fs), nrow=2, byrow=TRUE), nYears=T)
  pop2 <- simPop(Ni=Ni, phi=c(sj, sa), f=matrix(c(fs, fs), nrow=2, byrow=TRUE), nYears=T)
  pop3 <- simPop(Ni=Ni, phi=c(sj, sa), f=matrix(c(fs, fs), nrow=2, byrow=TRUE), nYears=T)

  # Create capture histories & m-arrays from population 1
  ch <- simCapHist(state=pop1$state, cap=cap, recap=recap, maxAge=2)
  marr <- marrayAge(ch$ch, ch$age)

  # Create productivity data from population 2
  pro <- simProd(reprod=pop2$reprod, pInclude=pprod)
  # Aggregate productivity data to make the model run faster
  J <- pro$prod.agg[,1]
  B <- pro$prod.agg[,2]

  # Create the population survey data from population 3
  count <- simCountNorm(N=pop3$totB, sigma=sigma)$count

  # Bundle data
  # All 3 data sets
  jagsdata3[[3]] <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), J=J, B=B,
      C=count, pNinit=dUnif(1, 100), upper.sigma.obs=100)

  # No productivity data
  jagsdata2[[3]] <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]),
      C=count, pNinit=dUnif(1, 100), upper.sigma.obs=100)


  # Scenario 4: productivity and adult survival time-dependent
  # Create populations
  pop1 <- simPop(Ni=Ni, phi=matrix(c(rep(sj,T-1), sas), nrow=2, byrow=TRUE),
      f=matrix(c(fs, fs), nrow=2, byrow=TRUE), nYears=T)
  pop2 <- simPop(Ni=Ni, phi=matrix(c(rep(sj,T-1), sas), nrow=2, byrow=TRUE),
      f=matrix(c(fs, fs), nrow=2, byrow=TRUE), nYears=T)
  pop3 <- simPop(Ni=Ni, phi=matrix(c(rep(sj,T-1), sas), nrow=2, byrow=TRUE),
      f=matrix(c(fs, fs), nrow=2, byrow=TRUE), nYears=T)

  # Create capture histories & m-arrays from population 1
  ch <- simCapHist(state=pop1$state, cap=cap, recap=recap, maxAge=2)
  marr <- marrayAge(ch$ch, ch$age)

  # Create productivity data from population 2
  pro <- simProd(reprod=pop2$reprod, pInclude=pprod)
  # Aggregate productivity data to make the model run faster
  J <- pro$prod.agg[,1]
  B <- pro$prod.agg[,2]

  # Create the population survey data from population 3
  count <- simCountNorm(N=pop3$totB, sigma=sigma)$count

  # Bundle data
  # All 3 data sets
  jagsdata3[[4]] <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), J=J, B=B,
      C=count, pNinit=dUnif(1, 100), upper.sigma.obs=100)

  # No productivity data
  jagsdata2[[4]] <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]),
      C=count, pNinit=dUnif(1, 100), upper.sigma.obs=100)



  # Data analyses - call JAGS from R (jagsUI)
  # All 3 data sets

  # 32 fits = 2 types of data sets x 4 scenarios x 4 estimation models

  for (dat in 1:4){
    for (mod in 1:4){
      ni <- 4000
      out1 <- try(jags(jagsdata3[[dat]], inits, parameters1, modfile3[mod],
          n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, verbose=FALSE, parallel=TRUE))
      if(!inherits(out1, "try-error")){
        res3[,,s,dat,mod] <- out1$summary
        bP3[s,1,dat,mod] <- mean(out1$sims.list$Dmape.rep > out1$sims.list$Dmape.obs)
        bP3[s,2,dat,mod] <- mean(out1$sims.list$DFT.rep > out1$sims.list$DFT.obs)
        bP3[s,3,dat,mod] <- mean(out1$sims.list$Dd.rep > out1$sims.list$Dd.obs)
      } # if

      # just 2 data sets
      ni <- c(4000, 4000, 6000, 6000)[mod]
      out2 <- try(jags(jagsdata2[[dat]], inits, parameters2, modfile2[mod],
          n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, verbose=FALSE, parallel=TRUE))
      if(!inherits(out2, "try-error")){
        res2[,,s,dat,mod] <- out2$summary
        bP2[s,1,dat,mod] <- mean(out2$sims.list$Dmape.rep > out2$sims.list$Dmape.obs)
        bP2[s,2,dat,mod] <- mean(out2$sims.list$DFT.rep > out2$sims.list$DFT.obs)
      } # if
    } # mod
  } # dat

  # Store to time variable true values for ad. survival and productivity
  sa.true[s,] <- sas
  f.true[s,] <- fs
  print(s)
} ) # s


save(res3, res2, sj, sa, sat, f, ft, sa.true, f.true, sigma, recap, bP3, bP2,
    file="ResultsChapter7.3.Rdata")


# Produce graphs

# Things needed for all the plots

# get relative bias
# est is parameters x stats x sims, true is sims x params
rbias <- function(est, true){
  means <- t(est[,1,])
  return((means - true) / true)
}

# get RMSE
rmse <- function(est, true){
  means <- t(est[,1,])
  sds <- t(est[,2,])
  mse <- (means - true)^2 + sds^2
  sqrt(mse)
}

ord <- c(1,10,2,11,3,12,4,13,5,14,6,15,7,16,8,17,9,18)
loc <- c(0.8, 1.2, 1.8, 2.2, 2.8, 3.2, 3.8, 4.2, 4.8, 5.2, 5.8, 6.2, 6.8, 7.2, 7.8, 8.2, 8.8, 9.2)

# Panel colours
library(scales)
co1 <- viridis_pal(option='E')(20)[c(4, 13, 20)] # blue-grey, grey, yellow
co2 <- matrix(c(
    NA,1, 1, 1,
    2, NA,3, 1,
    2, 3, NA,1,
    2, 2, 2, NA), 4, 4, byrow=TRUE)
panelco <- matrix(co1[co2], 4, 4)

# Labels
topLabels <- c(expression(bold('Data: '*italic(f)*'(.), '*italic(s)[a]*'(.)')),
    expression(bold('Data: '*italic(f)*'(.), '*italic(s)[a]*'(t)')),
    expression(bold('Data: '*italic(f)*'(t), '*italic(s)[a]*'(.)')),
    expression(bold('Data: '*italic(f)*'(t), '*italic(s)[a]*'(t)')) )
leftLabels <- c(expression(bold('Model: '*italic(f)*'(.), '*italic(s)[a]*'(.)')),
    expression(bold('Model: '*italic(f)*'(.), '*italic(s)[a]*'(t)')),
    expression(bold('Model: '*italic(f)*'(t), '*italic(s)[a]*'(.)')),
    expression(bold('Model: '*italic(f)*'(t), '*italic(s)[a]*'(t)')) )

# Figure function
plotmystuff <- function(ests1, ests2, truth, STAT, ord, loc, ylim, panelco, al,
    ylab, xlab, topLabels, leftLabels){
  op <- par("mfrow", "las", "mar")
  layout(matrix(1:16, 4, 4, byrow=TRUE), widths=c(2, 1.25, 1.25, 1.25),
      heights=c(1.1, 1, 1, 1.25), TRUE)
  for (row in 1:4){     # rows are models
    for (col in 1:4){  # columns are data scenarios
      # fix margins
      mymars <- c(1,1,1,1)
      if(col == 1) mymars[2] <- 7
      if(row == 1) mymars[3] <- 2
      if(row == 4) mymars[1] <- 3
      par(las=1, mar=mymars)
      # do a blank plot, add coloured rectangle
      plot(NA, ylim=ylim, xlim=c(0, 10), ylab="", xlab=NA, axes=FALSE)
      if(!is.na(panelco[row,col])){
        rect(0, ylim[1], 10, ylim[2], border=NA, col=alpha(panelco[row,col], al))
      }

      boxplot(cbind(STAT(ests1[,,,col,row], truth[,,col]),
          STAT(ests2[,,,col,row], truth[,,col]))[,ord], ylab=NA, outline=FALSE,
          xlab="", col=rep(c("grey", "white"), 9), axes=FALSE, boxwex=0.3,
          ylim=ylim, at=loc, add=TRUE)

       # axes with tick-marks
       axis(1, at=1:9, tcl=-0.25, labels=FALSE)
       axis(1, at=c(1, 3, 5, 7, 9), tcl=-0.5, labels = (row == 4))
       axis(2, labels = (col==1), las=1)
       # axis labels
       if(col == 1) title(ylab=ylab)
       if(row == 4) mtext(xlab, side=1, line=2)
       # outer margin labels
       if(row == 1) mtext(topLabels[col], side=3, line=0.5)
       if(col == 1) mtext(leftLabels[1], side=2, line=5, las=3)
    } # col
  } # row
  par(op)
}


# Load data
load("ResultsChapter7.3.Rdata")

# Select only converged simulations, based on Rhat for params 1:23
r.crit <- 1.1

rhats3 <- res3[1:23,8,,,] < r.crit
good3 <- apply(rhats3, 2, all)
rhats2 <- res2[1:23,8,,,] < r.crit
good2 <- apply(rhats2, 2, all)
incl <- which(good3 & good2)
incl <- incl[1:500]


# Plots for productivity
params <- 14:22      # f[1:9]
toplot3 <- res3[params,,incl,,]
toplot2 <- res2[params,,incl,,]

# Take a look
str(toplot3[,,,1,1]) # 9 x 11 x 500

truef <- abind::abind(matrix(f, length(incl), 9), matrix(f, length(incl), 9),
    f.true[incl,1:9], f.true[incl,1:9], along=3)
# first 2 columns, f constant
# last 2 columns, f varies, different for each sim

# Fig. 7.12
plotmystuff(ests1 = toplot3, ests2 = toplot2, truth = truef, STAT = rbias,
  ord = ord, loc=loc,
  ylim = c(-1.2,1.2), panelco = panelco, al = 0.5,
  ylab = expression('Relative bias ('*italic(f)*')'), xlab = "Time",
  topLabels = topLabels, leftLabels = leftLabels)

# Fig 7.13
plotmystuff(ests1 = toplot3, ests2 = toplot2, truth = truef, STAT = rmse,
  ord = ord, loc=loc,
  ylim = c(0,3.5), panelco = panelco, al = 0.5,
  ylab = expression('RMSE ('*italic(f)*')'), xlab = "Time",
  topLabels = topLabels, leftLabels = leftLabels)


# Plots for adult survival
# New 'toplot' objects
params <- 5:13  # sa[1:9]
toplot3 <- res3[params,,incl,,]
toplot2 <- res2[params,,incl,,]

# New truth array
str(sa.true)
truesa <- abind::abind(matrix(sa, length(incl), 9), sa.true[incl, ], matrix(sa, length(incl), 9),
    sa.true[incl, ], along=3)
# columns 1 & 3 constant, 2 & 4 varying
str(truesa)

# Fig 7.14
plotmystuff(ests1=toplot3, ests2=toplot2, truth=truesa, STAT=rbias, ord=ord,
    loc=loc, ylim=c(-0.4, 0.4), panelco=panelco, al=0.5,
    ylab=expression('Relative bias ('*italic(s)[a]*')'), xlab="Time",
    topLabels=topLabels, leftLabels=leftLabels)

# Fig 7.15
plotmystuff(ests1=toplot3, ests2=toplot2, truth=truesa, STAT=rmse, ord=ord,
    loc=loc, ylim=c(0, 0.25), panelco=panelco, al=0.5,
    ylab=expression('RMSE ('*italic(s)[a]*')'), xlab="Time",
    topLabels=topLabels, leftLabels=leftLabels)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
