# Schaub & Kéry (2022) Integrated Population Models
# Chapter 7 : Assessment of integrated population models
# ------------------------------------------------------
# Code from final MS.

# Run time for test script 4 mins, full run 9 hrs

library(IPMbook) ; library(jagsUI)


# 7.2 Assumptions of integrated population models
# ===============================================

# 7.2.3 The common demography assumption
# --------------------------------------

# ~~~~ Code for the simulations ~~~~

# Analysing models
# 1. IPM with all 3 data sets
cat(file="model8.txt", "
model {
  # Priors and constraints
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
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
    C.exp[t] <- N[1,t] + N[2,t]                            # Expected counts
    Dssm.obs[t] <- abs((C[t] - C.exp[t]) / C[t])           # Discrepancy measure

    C.rep[t] ~ dnorm(N[1,t] + N[2,t], tau.obs)             # Generate replicate data
    Dssm.rep[t] <- abs((C.rep[t] - C.exp[t]) / C.rep[t])   # Discrepancy measure
  }
  Dmape.obs <- sum(Dssm.obs)
  Dmape.rep <- sum(Dssm.rep)

  # Productivity data (Poisson regression model)
  for (t in 1:n.occasions){
    J[t] ~ dpois(B[t] * f[t])

    J.exp[t] <- B[t] * f[t]                                # Expected data
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
    q[t] <- 1-p[t]                                         # Probability of non-recapture
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
cat(file="model9.txt", "
model {
  # Priors and constraints
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
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
    C.exp[t] <- N[1,t] + N[2,t]                            # Expected counts
    Dssm.obs[t] <- abs((C[t] - C.exp[t]) / C[t])           # Discrepancy measure

    C.rep[t] ~ dnorm(N[1,t] + N[2,t], tau.obs)             # Generate replicate data
    Dssm.rep[t] <- abs((C.rep[t] - C.exp[t]) / C.rep[t])   # Discrepancy measure
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


# Simulations
# 1. Small sample size

# Number of simulations
# nsim <- 600  # 2.5 hrs
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
# Populations 1, 2 and 3
  # Age specific survival probabilities (juv, adult)
  phi1 <- c(0.3, 0.55)
  # Fecundity rate (females)
  f1 <- 3.1

# Population 4
  # Age specific survival probabilities (juv, adult)
  phi2 <- c(0.3, 0.64)
  # Fecundity rate (females)
  f2 <- 3.1

# Population 5
  # Age specific survival probabilities (juv, adult)
  phi3 <- c(0.3, 0.55)
  # Fecundity rate (females)
  f3 <- 3.7

# Define matrices to store the result
res1 <- res2 <- res3 <- array(NA, dim=c(43, 11, nsim))
res4 <- res5 <- res6 <- array(NA, dim=c(41, 11, nsim))
bP1 <- bP2 <- bP3 <- matrix(NA, ncol=3, nrow=nsim)
bP4 <- bP5 <- bP6 <- matrix(NA, ncol=2, nrow=nsim)

# Initial values
inits <- function(){
  N <- matrix(NA, nrow=2, ncol=T)
  N[1,] <- round(runif(T, 40, 60))
  N[2,] <- round(runif(T, 40, 60))
  ini <- list(mean.sj=runif(1, 0, 0.5), N=N)
  return(ini)
}

# Parameters monitored
parameters1 <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma.obs",
    "geom.rate", "Ntot", "Dmape.obs", "Dmape.rep", "DFT.obs", "DFT.rep", "Dd.obs", "Dd.rep")
parameters2 <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma.obs",
    "geom.rate", "Ntot", "Dmape.obs", "Dmape.rep", "DFT.obs", "DFT.rep")

# MCMC settings
ni <- 4000; nb <- 2000; nc <- 3; nt <- 4; na <- 2000

# Start simulations
system.time(
for (s in 1:nsim){
  set.seed(s)

  # Create populations
  pop1 <- simPop(Ni=Ni, phi=phi1, f=f1, nYears=T)
  pop2 <- simPop(Ni=Ni, phi=phi1, f=f1, nYears=T)
  pop3 <- simPop(Ni=Ni, phi=phi1, f=f1, nYears=T)
  pop4 <- simPop(Ni=Ni, phi=phi2, f=f2, nYears=T)
  pop5 <- simPop(Ni=Ni, phi=phi3, f=f3, nYears=T)

  # Create capture histories & m-arrays from population 1
  ch <- simCapHist(state=pop1$state, cap=cap, recap=recap, maxAge=2)
  marr <- marrayAge(ch$ch, ch$age)

  # Create productivity data from population 2
  pro <- simProd(reprod=pop2$reprod, pInclude=pprod)
  # Aggregate productivity data to make the model run faster
  J <- pro$prod.agg[,1]
  B <- pro$prod.agg[,2]

  # Create the population survey data from population 3, 4 and 5
  count1 <- simCountNorm(N=pop3$totB, sigma=sigma)$count
  count2 <- simCountNorm(N=pop4$totB, sigma=sigma)$count
  count3 <- simCountNorm(N=pop5$totB, sigma=sigma)$count

  # Bundle data
  # A: all 3 data sets
  # CMR from pop 1, productivity from pop 2, counts from pop 3 (same dynamics)
  jags.data1 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), J=J, B=B,
      C=count1, pNinit=dUnif(1, 100), upper.sigma.obs=100)

  # CMR from pop 1, productivity from pop 2, counts from pop 4
  jags.data2 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), J=J, B=B,
      C=count2, pNinit=dUnif(1, 100), upper.sigma.obs=100)

  # CMR from pop 1, productivity from pop 2, counts from pop 5
  jags.data3 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), J=J, B=B,
      C=count3, pNinit=dUnif(1, 100), upper.sigma.obs=100)

  # B: Counts and capture-recapture data only
  # CMR from pop 1, counts from pop 3 (same dynamics)
  jags.data4 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]),
      C=count1, pNinit=dUnif(1, 100), upper.sigma.obs=100)

  # CMR from pop 1, counts from pop 4
  jags.data5 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]),
      C=count2, pNinit=dUnif(1, 100), upper.sigma.obs=100)

  # CMR from pop 1, counts from pop 5
  jags.data6 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]),
      C=count3, pNinit=dUnif(1, 100), upper.sigma.obs=100)


  # Call JAGS from R (jagsUI)
  # CMR 1, P 2, count 3, same dynamics
  out15 <- try(jags(jags.data1, inits, parameters1, "model8.txt",
      n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=TRUE))
  if(!inherits(out15, "try-error")){
    res1[,,s] <- out15$summary
    bP1[s,1] <- mean(out15$sims.list$Dmape.rep > out15$sims.list$Dmape.obs)
    bP1[s,2] <- mean(out15$sims.list$DFT.rep > out15$sims.list$DFT.obs)
    bP1[s,3] <- mean(out15$sims.list$Dd.rep > out15$sims.list$Dd.obs)
  }

  # CMR 1, P 2, count 4
  out16 <- try(jags(jags.data2, inits, parameters1, "model8.txt",
      n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=TRUE))
  if(!inherits(out16, "try-error")){
    res2[,,s] <- out16$summary
    bP2[s,1] <- mean(out16$sims.list$Dmape.rep > out16$sims.list$Dmape.obs)
    bP2[s,2] <- mean(out16$sims.list$DFT.rep > out16$sims.list$DFT.obs)
    bP2[s,3] <- mean(out16$sims.list$Dd.rep > out16$sims.list$Dd.obs)
  }

  # CMR 1, P 2, count 5
  out17 <- try(jags(jags.data3, inits, parameters1, "model8.txt",
      n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=TRUE))
  if(!inherits(out17, "try-error")){
    res3[,,s] <- out17$summary
    bP3[s,1] <- mean(out17$sims.list$Dmape.rep > out17$sims.list$Dmape.obs)
    bP3[s,2] <- mean(out17$sims.list$DFT.rep > out17$sims.list$DFT.obs)
    bP3[s,3] <- mean(out17$sims.list$Dd.rep > out17$sims.list$Dd.obs)
  }

  # CMR 1, count 3, same dynamics
  out18 <- try(jags(jags.data4, inits, parameters2, "model9.txt",
      n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=TRUE))
  if(!inherits(out18, "try-error")){
    res4[,,s] <- out18$summary
    bP4[s,1] <- mean(out18$sims.list$Dmape.rep > out18$sims.list$Dmape.obs)
    bP4[s,2] <- mean(out18$sims.list$DFT.rep > out18$sims.list$DFT.obs)
  }

  # CMR 1, count 4
  out19 <- try(jags(jags.data5, inits, parameters2, "model9.txt",
      n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=TRUE))
  if(!inherits(out19, "try-error")){
    res5[,,s] <- out19$summary
    bP5[s,1] <- mean(out19$sims.list$Dmape.rep > out19$sims.list$Dmape.obs)
    bP5[s,2] <- mean(out19$sims.list$DFT.rep > out19$sims.list$DFT.obs)
  }

  # CMR 1, count 5
  out20 <- try(jags(jags.data6, inits, parameters2, "model9.txt",
      n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=TRUE))
  if(!inherits(out20, "try-error")){
    res6[,,s] <- out20$summary
    bP6[s,1] <- mean(out20$sims.list$Dmape.rep > out20$sims.list$Dmape.obs)
    bP6[s,2] <- mean(out20$sims.list$DFT.rep > out20$sims.list$DFT.obs)
  }
print(s)
} ) #s

save(res1, res2, res3, res4, res5, res6, phi1, f1, phi2, f2, phi3, f3, sigma,
    recap, bP1, bP2, bP3, bP4, bP5, bP6, file="ResultsChapter7.2.3.Rdata")

##############################

# 2. Large sample size (five times larger)

# Number of simulations
# nsim <- 600  # 4.7 hrs
nsim <- 6  # ~~~ for testing

# Number of years
T <- 10

# Observation error for the population survey
sigma <- 10 * 5

# Capture and recapture probabilities
cap <- 0.4                        # initial capture probability
recap <- 0.6                      # recapture probability

# Probability to find a brood whose reproductive ouput is recorded
pprod <- 0.5

# Initial population size per age class
Ni <- c(50 * 5, 50 * 5)

# Demographic rates
# Populations 1, 2 and 3
  # Age specific survival probabilities (juv, adult)
  phi1 <- c(0.3, 0.55)
  # Fecundity rate (females)
  f1 <- 3.1

# Population 4
  # Age specific survival probabilities (juv, adult)
  phi2 <- c(0.3, 0.64)
  # Fecundity rate (females)
  f2 <- 3.1

# Population 5
  # Age specific survival probabilities (juv, adult)
  phi3 <- c(0.3, 0.55)
  # Fecundity rate (females)
  f3 <- 3.7


# Define matrices to store the result
res1 <- res2 <- res3 <- array(NA, dim=c(43, 11, nsim))
res4 <- res5 <- res6 <- array(NA, dim=c(41, 11, nsim))
bP1 <- bP2 <- bP3 <- matrix(NA, ncol=3, nrow=nsim)
bP4 <- bP5 <- bP6 <- matrix(NA, ncol=2, nrow=nsim)


# Initial values
inits <- function(){
  N <- matrix(NA, nrow=2, ncol=T)
  N[1,] <- round(runif(T, 200, 300))
  N[2,] <- round(runif(T, 200, 300))
  ini <- list(mean.sj=runif(1, 0, 0.5), N=N)
  return(ini)
}

# Parameters monitored
parameters1 <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma.obs",
    "geom.rate", "Ntot", "Dmape.obs", "Dmape.rep", "DFT.obs", "DFT.rep", "Dd.obs", "Dd.rep")
parameters2 <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma.obs",
    "geom.rate", "Ntot", "Dmape.obs", "Dmape.rep", "DFT.obs", "DFT.rep")

# MCMC settings
ni <- 4000; nb <- 2000; nc <- 3; nt <- 4; na <- 2000

# Start simulations
system.time(
for (s in 1:nsim){
  set.seed(s)

  # Create populations
  pop1 <- simPop(Ni=Ni, phi=phi1, f=f1, nYears=T)
  pop2 <- simPop(Ni=Ni, phi=phi1, f=f1, nYears=T)
  pop3 <- simPop(Ni=Ni, phi=phi1, f=f1, nYears=T)
  pop4 <- simPop(Ni=Ni, phi=phi2, f=f2, nYears=T)
  pop5 <- simPop(Ni=Ni, phi=phi3, f=f3, nYears=T)

  # Create the capture histories & m-arrays from population 1
  ch <- simCapHist(state=pop1$state, cap=cap, recap=recap, maxAge=2)
  marr <- marrayAge(ch$ch, ch$age)

  # Create productivity data from population 2
  pro <- simProd(reprod=pop2$reprod, pInclude=pprod)
  # Aggregate productivity data to make the model run faster
  J <- pro$prod.agg[,1]
  B <- pro$prod.agg[,2]

  # Create the population survey data from population 3, 4 and 5
  count1 <- simCountNorm(N=pop3$totB, sigma=sigma)$count
  count2 <- simCountNorm(N=pop4$totB, sigma=sigma)$count
  count3 <- simCountNorm(N=pop5$totB, sigma=sigma)$count


  # Bundle data
  # A: all 3 data sets
  # CMR from pop 1, productivity from pop 2, counts from pop 3 (same dynamics)
  jags.data1 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), J=J, B=B, C=count1,
      pNinit=dUnif(150, 350), upper.sigma.obs=300)

  # CMR from pop 1, productivity from pop 2, counts from pop 4
  jags.data2 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), J=J, B=B, C=count2,
      pNinit=dUnif(150, 350), upper.sigma.obs=300)

  # CMR from pop 1, productivity from pop 2, counts from pop 5
  jags.data3 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), J=J, B=B, C=count3,
      pNinit=dUnif(150, 350), upper.sigma.obs=300)

  # B: Counts and capture-recapture data only
  # CMR from pop 1, counts from pop 3 (same dynamics)
  jags.data4 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), C=count1,
      pNinit=dUnif(150, 350), upper.sigma.obs=300)

  # CMR from pop 1, counts from pop 4
  jags.data5 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), C=count2,
      pNinit=dUnif(150, 350), upper.sigma.obs=300)

  # CMR from pop 1, counts from pop 5
  jags.data6 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), C=count3,
      pNinit=dUnif(150, 350), upper.sigma.obs=300)


  # Call JAGS from R (jagsUI)
  # CMR 1, P 2, count 3, same dynamics
  out15 <- try(jags(jags.data1, inits, parameters1, "model8.txt",
      n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=TRUE))
  if(!inherits(out15, "try-error")){
    res1[,,s] <- out15$summary
    bP1[s,1] <- mean(out15$sims.list$Dmape.rep > out15$sims.list$Dmape.obs)
    bP1[s,2] <- mean(out15$sims.list$DFT.rep > out15$sims.list$DFT.obs)
    bP1[s,3] <- mean(out15$sims.list$Dd.rep > out15$sims.list$Dd.obs)
  }

  # CMR 1, P 2, count 4
  out16 <- try(jags(jags.data2, inits, parameters1, "model8.txt",
      n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=TRUE))
  if(!inherits(out16, "try-error")){
    res2[,,s] <- out16$summary
    bP2[s,1] <- mean(out16$sims.list$Dmape.rep > out16$sims.list$Dmape.obs)
    bP2[s,2] <- mean(out16$sims.list$DFT.rep > out16$sims.list$DFT.obs)
    bP2[s,3] <- mean(out16$sims.list$Dd.rep > out16$sims.list$Dd.obs)
  }

  # CMR 1, P 2, count 5
  out17 <- try(jags(jags.data3, inits, parameters1, "model8.txt",
      n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=TRUE))
  if(!inherits(out17, "try-error")){
    res3[,,s] <- out17$summary
    bP3[s,1] <- mean(out17$sims.list$Dmape.rep > out17$sims.list$Dmape.obs)
    bP3[s,2] <- mean(out17$sims.list$DFT.rep > out17$sims.list$DFT.obs)
    bP3[s,3] <- mean(out17$sims.list$Dd.rep > out17$sims.list$Dd.obs)
  }

  # CMR 1, count 3, same dynamics
  out18 <- try(jags(jags.data4, inits, parameters2, "model9.txt",
      n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=TRUE))
  if(!inherits(out18, "try-error")){
    res4[,,s] <- out18$summary
    bP4[s,1] <- mean(out18$sims.list$Dmape.rep > out18$sims.list$Dmape.obs)
    bP4[s,2] <- mean(out18$sims.list$DFT.rep > out18$sims.list$DFT.obs)
  }

  # CMR 1, count 4
  out19 <- try(jags(jags.data5, inits, parameters2, "model9.txt",
      n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=TRUE))
  if(!inherits(out19, "try-error")){
    res5[,,s] <- out19$summary
    bP5[s,1] <- mean(out19$sims.list$Dmape.rep > out19$sims.list$Dmape.obs)
    bP5[s,2] <- mean(out19$sims.list$DFT.rep > out19$sims.list$DFT.obs)
  }

  # CMR 1, count 5
  out20 <- try(jags(jags.data6, inits, parameters2, "model9.txt",
      n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=TRUE))
  if(!inherits(out20, "try-error")){
    res6[,,s] <- out20$summary
    bP6[s,1] <- mean(out20$sims.list$Dmape.rep > out20$sims.list$Dmape.obs)
    bP6[s,2] <- mean(out20$sims.list$DFT.rep > out20$sims.list$DFT.obs)
  }
  print(s)
} ) #s

save(res1, res2, res3, res4, res5, res6, phi1, f1, phi2, f2, phi3, f3, sigma,
    recap, bP1, bP2, bP3, bP4, bP5, bP6, file="ResultsChapter7.2.3large.Rdata")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code for Figure 7.10 ~~~~

load("ResultsChapter7.2.3.Rdata")

# Select 500 simulations that have converged
r.crit <- 1.05
incl <- which(res1[1,8,]<r.crit & res1[2,8,]<r.crit & res1[4,8,]<r.crit & res1[26,8,]<r.crit & res1[25,8,]<r.crit &
               res2[1,8,]<r.crit & res2[2,8,]<r.crit & res2[4,8,]<r.crit & res2[26,8,]<r.crit & res2[25,8,]<r.crit &
               res3[1,8,]<r.crit & res3[2,8,]<r.crit & res3[4,8,]<r.crit & res3[26,8,]<r.crit & res3[25,8,]<r.crit &
               res4[1,8,]<r.crit & res4[2,8,]<r.crit & res4[4,8,]<r.crit & res4[26,8,]<r.crit & res4[25,8,]<r.crit &
               res5[1,8,]<r.crit & res5[2,8,]<r.crit & res5[4,8,]<r.crit & res5[26,8,]<r.crit & res5[25,8,]<r.crit &
               res6[1,8,]<r.crit & res6[2,8,]<r.crit & res6[4,8,]<r.crit & res6[26,8,]<r.crit & res6[25,8,]<r.crit)
incl <- incl[1:500]

library(scales)
co <- viridis_pal(option='E')(20)[c(5, 11, 16)]

op <- par(mfrow=c(2,2), las =1, mar=c(2.5, 4.5, 2.5, 1))
lab <- expression('Juvenile survival ('*italic(s)[italic(j)]*')')
boxplot(cbind(res1[1,1,incl], res2[1,1,incl], res3[1,1,incl], res4[1,1,incl],
    res5[1,1,incl], res6[1,1,incl]), ylab=lab, outline=FALSE,
    col=c(co, rep("white", 3)), border=c(rep("black", 3), co),
    lwd=c(rep(1,3), rep(2,3)), axes=FALSE, boxwex=0.6)
abline(h=phi1[1], lty=1)
axis(2)
mtext('No hidden parameter', side=3, at=2, line=1)
mtext('Productivity hidden', side=3, at=5, line=1)

lab <- expression('Adult survival ('*italic(s)[italic(a)]*')')
boxplot(cbind(res1[2,1,incl], res2[2,1,incl], res3[2,1,incl], res4[2,1,incl],
    res5[2,1,incl], res6[2,1,incl]), ylab=lab, outline=FALSE,
    col=c(co, rep("white", 3)), border=c(rep("black", 3), co),
    lwd=c(rep(1,3), rep(2,3)), axes=FALSE, ylim=c(0.48, 0.65), boxwex=0.6)
abline(h=phi1[2], lty=1)
segments(1.5, phi2[2], 2.5, phi2[2], lty=2)
segments(4.5, phi2[2], 5.5, phi2[2], lty=2)
axis(2)
mtext('No hidden parameter', side=3, at=2, line=1)
mtext('Productivity hidden', side=3, at=5, line=1)

lab <- expression('Productivity ('*italic(f)*')')
boxplot(cbind(res1[4,1,incl], res2[4,1,incl], res3[4,1,incl], res4[4,1,incl],
    res5[4,1,incl], res6[4,1,incl]), ylab=lab, outline=FALSE,
    col=c(co, rep("white", 3)), border=c(rep("black", 3), co),
    lwd=c(rep(1,3), rep(2,3)), axes=FALSE, boxwex=0.6)
abline(h=f1, lty=1)
segments(2.5, f3, 3.5, f3, lty=2)
segments(5.5, f3, 6.5, f3, lty=2)
axis(2)
legend('topleft', pch=rep(15, 3), col=co,
    legend=c('Fulfilled', expression('Violated: '*italic(s)[italic(a)]*' different'),
    expression('Violated: '*italic(f)*' different')), bty='n', inset=0.025)
mtext('Common demography assumption', side=3, at=0.5, line=-0.5, adj=0)

lab <- expression('Population growth rate ('*lambda*')')
boxplot(cbind(res1[26,1,incl], res2[26,1,incl], res3[26,1,incl], res4[26,1,incl],
    res5[26,1,incl], res6[26,1,incl]), ylab=lab,  outline=FALSE,
    col=c(co, rep("white", 3)), border=c(rep("black", 3), co),
    lwd=c(rep(1,3), rep(2,3)), axes=FALSE, boxwex=0.6)
abline(h=phi1[1] * f1 / 2 + phi1[2], lty=1)
segments(1.5, phi2[1] * f2 / 2 + phi2[2], 2.5, phi2[1] * f2 / 2 + phi2[2], lty=2)
segments(2.5, phi3[1] * f3 / 2 + phi3[2], 3.5, phi3[1] * f3 / 2 + phi3[2], lty=2)
segments(4.5, phi2[1] * f2 / 2 + phi2[2], 5.5, phi2[1] * f2 / 2 + phi2[2], lty=2)
segments(5.5, phi3[1] * f3 / 2 + phi3[2], 6.5, phi3[1] * f3 / 2 + phi3[2], lty=2)
axis(2)
par(op)

# Fig. 7.11

load("ResultsChapter7.2.3.Rdata")

op <- par(mfrow=c(3, 1), las =1, mar=c(3, 4, 3, 1), cex=1.05)

lab <- expression('Bayesian '*italic(p)*'-value (SSM)')
boxplot(cbind(bP1[incl,1], bP2[incl,1], bP3[incl,1], bP4[incl,1], bP5[incl,1],
    bP6[incl,1]), ylab=lab,  outline=FALSE, col=c(co, rep("white", 3)),
    border=c(rep("black", 3), co), lwd=c(rep(1,3), rep(2,3)), axes=FALSE,
    boxwex=0.6, ylim=c(0, 1))
axis(2)
mtext('No hidden parameter', side=3, at=2, line=1)
mtext('Productivity hidden', side=3, at=5, line=1)
abline(h=0.5, lty=2)

par(mar=c(3, 4, 1, 1))
lab <- expression('Bayesian '*italic(p)*'-value (CMR)')
boxplot(cbind(bP1[incl,2], bP2[incl,2], bP3[incl,2], bP4[incl,2], bP5[incl,2],
    bP6[incl,2]), ylab=lab,  outline=FALSE, col=c(co, rep("white", 3)),
    border=c(rep("black", 3), co), lwd=c(rep(1,3), rep(2,3)), axes=FALSE, boxwex=0.6)
axis(2)
abline(h=0.5, lty=2)

lab <- expression('Bayesian '*italic(p)*'-value (Prod)')
boxplot(cbind(bP1[incl,3], bP2[incl,3], bP3[incl,3], NA, NA, NA), ylab=lab,
    outline=FALSE, col=c(co, rep("white", 3)), border=c(rep("black", 3), co),
    lwd=c(rep(1,3), rep(2,3)), axes=FALSE, boxwex=0.6)
axis(2)
legend(x=4, y=0.85, pch=rep(15, 3), col=co,
    legend=c('Fulfilled', expression('Violated: '*italic(s)[italic(a)]*' different'),
    expression('Violated: '*italic(f)*' different')), bty='n', inset=0.025)
text(x=4, y=0.925, 'Common demography\nassumption', adj=0)
abline(h=0.5, lty=2)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
