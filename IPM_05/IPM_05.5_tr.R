# Schaub & Kéry (2022) Integrated Population Models
# Chapter 5 : Introduction to integrated population models
# --------------------------------------------------------

# Run time for test script 3 mins, full run 15 hrs

library(IPMbook) ; library(jagsUI)

# 5.5 Simulation assessment of a simple IPM
# =========================================

# 5.5.1 Simulating data under an integrated population model
# ----------------------------------------------------------

# Pick values for the function arguments
T <- 20                                 # Number of years
phi <- c(0.3, 0.55)                     # Age specific survival probabilities (juv, adult)
f <- c(2.6, 3.6)                        # Age-specific productivity (1y, older)
Ni <- c(50, 50)                         # Initial pop. size for each age class (1y, older)

# Apply the function and produce data overview
set.seed(111167)                        # To initialize the RNGs at the same place
pop <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)
str(pop)

# List of 14
# $ Ni         : num [1:2] 50 50
# $ phi        : num [1:3, 1:19] 0.3 0.55 0.55 0.3 0.55 0.55 0.3 0.55 ...
# $ f          : num [1:2, 1:20] 2.6 3.6 2.6 3.6 2.6 3.6 2.6 3.6 2.6 3.6 ...
# $ pBreed     : num [1:2, 1:20] 1 1 1 1 1 1 1 1 1 1 ...
# $ sex.ratio  : num [1:20] 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ...
# $ Im         : num [1:20] 0 0 0 0 0 0 0 0 0 0 ...
# $ ageOfIm    : num [1:20] 1 1 1 1 1 1 1 1 1 1 ...
# $ state      : num [1:3492, 1:20] 1 1 1 1 1 1 1 1 1 1 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : NULL
# .. ..$ : chr [1:20] "Y1" "Y2" "Y3" "Y4" ...
# $ imYear     : logi [1:3492] NA NA NA NA NA NA ...
# $ reprod     : num [1:3492, 1:20, 1:3] 5 2 2 1 0 1 1 1 2 2 ...
# ..- attr(*, "dimnames")=List of 3
# .. ..$ : NULL
# .. ..$ : chr [1:20] "Y1" "Y2" "Y3" "Y4" ...
# .. ..$ : chr [1:3] "F" "M" "Age"
# $ N          : num [1:6, 1:20] 50 50 100 158 336 0 49 52 101 165 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:6] "1-Year" "2-Year" "totAdults" "BornF" ...
# .. ..$ : chr [1:20] "Y1" "Y2" "Y3" "Y4" ...
# $ breeders   : num [1:3, 1:20] 50 50 100 49 52 101 46 56 102 50 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:3] "1-Year" "2-Year" "totBreeders"
# .. ..$ : chr [1:20] "Y1" "Y2" "Y3" "Y4" ...
# $ totAdults  : Named num [1:20] 100 101 102 108 109 116 115 121 124 ...
# ..- attr(*, "names")= chr [1:20] "Y1" "Y2" "Y3" "Y4" ...
# $ totBreeders: Named num [1:20] 100 101 102 108 109 116 115 121 124 ...
# ..- attr(*, "names")= chr [1:20] "Y1" "Y2" "Y3" "Y4" ...

pop$state[488, 1:10]
# Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y9 Y10
# NA NA  0  1  2  2 -1 NA NA  NA

pop$reprod[488,1:10,]
#      F  M Age
# Y1  NA NA  NA
# Y2  NA NA  NA
# Y3  NA NA  NA
# Y4   3  2   1
# Y5   5  3   2
# Y6   2  0   2
# Y7  NA NA  NA
# Y8  NA NA  NA
# Y9  NA NA  NA
# Y10 NA NA  NA

pop$N[,1:10]
#            Y1  Y2  Y3  Y4  Y5  Y6  Y7  Y8  Y9 Y10
# 1-Year     50  49  46  50  53  57  53  61  61  51
# 2-Year     50  52  56  58  56  59  62  60  63  62
# totAdults 100 101 102 108 109 116 115 121 124 113
# BornF     158 165 167 164 171 167 193 193 192 207
# BornT     336 346 321 337 328 356 380 369 395 389
# Im          0   0   0   0   0   0   0   0   0   0

pop$breeders[,1:10]
#              Y1  Y2  Y3  Y4  Y5  Y6  Y7  Y8  Y9 Y10
# 1-Year       50  49  46  50  53  57  53  61  61  51
# 2-Year       50  52  56  58  56  59  62  60  63  62
# totBreeders 100 101 102 108 109 116 115 121 124 113

pop1 <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)
pop2 <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)
pop3 <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)

# Pick a value of the observation error (SD) for the population survey
sigma <- 10

# Create the population survey data and produce data overview
count <- simCountNorm(N=pop1$totB, sigma=sigma)
str(count)
# List of 2
# $ sigma: num [1:20] 10 10 10 10 10 10 10 10 10 10 ...
# $ count: num [1:20] 111 84 85 94 116 132 115 89 82 61 ...

# Pick values for capture and recapture probabilities
cap <- 0.4                              # Initial capture probability (same for juv. and adults)
recap <- 0.6                            # Recapture probability

# Create the capture histories and produce data overview
ch <- simCapHist(state=pop2$state, cap=cap, recap=recap, maxAge=2)
str(ch)
# List of 5
# $ cap   : num [1:3, 1:20] 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 ...
# $ recap : num [1:2, 1:19] 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 ...
# $ maxAge: num 2
# $ ch    : num [1:2645, 1:20] 1 1 1 1 1 1 0 1 0 1 ...
# $ age   : num [1:2645] 2 2 2 2 2 2 2 2 2 2 ...

# Create m-arrays
marr <- marrayAge(ch$ch, ch$age)

# Pick probability to find a brood
pprod <- 0.3

# Create productivity data and produce data overview
pro <- simProd(reprod=pop3$reprod, pInclude=pprod)
str(pro)
# List of 4
# $ pInclude    : num [1:20] 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 ...
# $ females.only: logi FALSE
# $ prod.ind    : num [1:705, 1:3] 2 2 0 3 5 3 5 4 5 3 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : NULL
# .. ..$ : chr [1:3] "Productivity" "Year" "Age of mother"
# $ prod.agg    : num [1:20, 1:2] 117 98 89 99 103 96 89 131 96 129 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:20] "1" "2" "3" "4" ...
# .. ..$ : chr [1:2] "Juveniles" "Surveyed broods"


# 5.5.2 Simulation results
# ------------------------

# ~~~~ code to run the simulations ~~~
# Define the 3 IPMs
# IPM 1: completely deterministic model

# Write JAGS model file
cat(file="model5.txt", "
model {
  # Priors and linear models
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f[1] ~ dunif(0, 10)
  mean.f[2] ~ dunif(0, 10)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }

  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- (N[1,t] * mean.f[1] / 2 + N[2,t] * mean.f[2] / 2) * mean.sj
    N[2,t+1] <- (N[1,t] + N[2,t]) * mean.sa
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Productivity data (Poisson regression model)
  sJ1 ~ dpois(nJ1 * mean.f[1])
  sJ2 ~ dpois(nJ2 * mean.f[2])

  # Capture-recapture data (multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.occasions-1)){
    # Main diagonal
    q[t] <- 1 - p[t]           # Probability of non-recapture
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
    resN[t] <- Ntot[t] - C[t]
  }
}
")

# IPM 2: include demographic stochasticity
# Write JAGS model file
cat(file="model6.txt", "
model {
  # Priors and linear models
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f[1] ~ dunif(0, 10)
  mean.f[2] ~ dunif(0, 10)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }

  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dcat(pinit)
  N[2,1] ~ dcat(pinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] ~ dpois((N[1,t] * mean.f[1] / 2 + N[2,t] * mean.f[2] / 2) * mean.sj)
    N[2,t+1] ~ dbin(mean.sa, N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Productivity data (Poisson regression model)
  sJ1 ~ dpois(nJ1 * mean.f[1])
  sJ2 ~ dpois(nJ2 * mean.f[2])

  # Capture-recapture data (multinomial likelihood)
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

  # Derived parameters
  # Annual population growth rate (added 0.001 to avoid possible division by 0)
  for (t in 1:(n.occasions-1)){
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t] + 0.001)
  }
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
    resN[t] <- Ntot[t] - C[t]
  }
  }
")

# IPM 3: include environmental and demographic stochasticity
# Write JAGS model file
cat(file="model7.txt", "
model {
  # Priors and linear models
  for (t in 1:(n.occasions-1)){
    logit.sj[t] ~ dnorm(logit.mean.sj, tau.sj)
    sj[t] <- ilogit(logit.sj[t])      # Back-transformation from logit scale
    logit.sa[t] ~ dnorm(logit.mean.sa, tau.sa)
    sa[t] <- ilogit(logit.sa[t])      # Back-transformation from logit scale
    p[t] <- mean.p
  }

  for (t in 1:n.occasions){
    log.f[1,t] ~ dnorm(log.mean.f[1], tau.f[1])
    f[1,t] <- exp(log.f[1,t])         # Back-transformation from log scale
    log.f[2,t] ~ dnorm(log.mean.f[2], tau.f[2])
    f[2,t] <- exp(log.f[2,t])         # Back-transformation from log scale
  }

  mean.sj ~ dunif(0, 1)
  logit.mean.sj <- logit(mean.sj)      # Logit transformation
  mean.sa ~ dunif(0, 1)
  logit.mean.sa <- logit(mean.sa)      # Logit transformation
  sigma.sj ~ dunif(0, 3)
  tau.sj <- pow(sigma.sj, -2)
  sigma.sa ~ dunif(0, 3)
  tau.sa <- pow(sigma.sa, -2)

  for (j in 1:2){
    mean.f[j] ~ dunif(0, 10)
    log.mean.f[j] <- log(mean.f[j])   # Log transformation
    sigma.f[j] ~ dunif(0, 3)
    tau.f[j] <- pow(sigma.f[j], -2)
  }

  mean.p ~ dunif(0, 1)

  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dcat(pinit)
  N[2,1] ~ dcat(pinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] ~ dpois((N[1,t] * f[1,t] / 2 + N[2,t] * f[2,t] / 2) * sj[t])
    N[2,t+1] ~ dbin(sa[t], N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Productivity data (Poisson regression model)
  for (i in 1:length(J)){
    J[i] ~ dpois(f[age[i],year[i]])
  }

  # Capture-recapture data (multinomial likelihood)
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

  # Derived parameters
  # Annual population growth rate (added 0.001 to avoid possible division by 0)
  for (t in 1:(n.occasions-1)){
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t] + 0.001)
  }
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
    resN[t] <- Ntot[t] - C[t]
  }
}
")


# Simulation parameters

# Age specific survival probabilities (juv, adult)
phi <- c(0.3, 0.55)

# Age-specific productivity (1y, older)
f <- c(2.6, 3.6)

# Initial population size per age class
Ni <- c(50, 50)

# Number of years
T <- 20

# Observation error for the population survey
sigma <- 10

# Capture and recapture probabilities
cap <- 0.4             # Initial capture probability
recap <- 0.6           # Recapture probability

# Probability to find a brood whose reproductive output is recorded
pprod <- 0.3

# Number of simulations
# nsim <- 1200
nsim <- 4  # ~~~~ for testing

# Define matrices to store the result
res1 <- res2 <- array(NA, dim=c(106, 11, nsim))
res3 <- array(NA, dim=c(188, 11, nsim))
lam <- matrix(NA, nrow=T-1, ncol=nsim)

# Initial values
inits <- function(){list(mean.sj=runif(1, 0, 0.5))}

# Parameters monitored
parameters1 <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma",
    "ann.growth.rate", "Ntot", "resN")

parameters2 <- c("mean.sj", "mean.sa", "mean.f", "mean.p", "sigma.sj", "sigma.sa",
    "sigma.f", "sj", "sa", "f", "N", "sigma", "ann.growth.rate", "Ntot", "resN")

# MCMC settings
ni <- 10000; nb <- 5000; nc <- 3; nt <- 1; na <- 1000


# Start simulations
system.time(
for (s in 1:nsim){
  set.seed(s)

  # Simulate 3 populations such that independent data sets are analyzed
  pop1 <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)
  pop2 <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)
  pop3 <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)
  # Calculate the annual population growth rates
  lam[,s] <- pop1$totA[-1] / pop1$totA[-T]

  # Simulate the population survey data
  count <- simCountNorm(N=pop1$totB, sigma=sigma)$count

  # Simulate the capture histories and the corresponding m-arrays
  ch <- simCapHist(state=pop2$state, cap=cap, recap=recap, maxAge=2)
  marr <- marrayAge(ch$ch, ch$age)

  # Simulate productivity data
  P <- simProd(reprod=pop3$reprod, pInclude=pprod)

  # Aggregate productivity data to make the models 1 and 2 run faster
  J1 <- P$prod.ind[P$prod.ind[,3]==1,1]
  J2 <- P$prod.ind[P$prod.ind[,3]==2,1]
  sJ1 <- sum(J1)
  nJ1 <- length(J1)
  sJ2 <- sum(J2)
  nJ2 <- length(J2)

  # Bundle data
  jags.data1 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), sJ1=sJ1, nJ1=nJ1,
      sJ2=sJ2, nJ2=nJ2, C=count)

  jags.data2 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), sJ1=sJ1, nJ1=nJ1,
      sJ2=sJ2, nJ2=nJ2, C=count, pinit=dUnif(1, 300))

  jags.data3 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), J=P$prod.ind[,1],
      year=P$prod.ind[,2], age=P$prod.ind[,3], C=count, pinit=dUnif(1, 300))

  # Call JAGS
  out5 <- try(jags(jags.data1, inits, parameters1, "model5.txt",
      n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=TRUE))
  if(!inherits(out5, "try-error"))
    res1[,,s] <- out5$summary

  out6 <- try(jags(jags.data2, inits, parameters1, "model6.txt",
      n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=TRUE))
  if(!inherits(out6, "try-error"))
    res2[,,s] <- out6$summary

  out7 <- try(jags(jags.data3, inits, parameters2, "model7.txt",
      n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=TRUE))
  if(!inherits(out7, "try-error"))
    res3[,,s] <- out7$summary

  print(s)
} ) #s

save(res1, res2, res3, lam, phi, f, sigma, recap,
    file="IPM Simulation chapter 5.5.Rdata")

load("IPM Simulation chapter 5.5.Rdata")
library(RColorBrewer)
co <- brewer.pal(n=8, name='Blues')[c(7,5,3)]

# Function to ensure that only converged estimates are included
r.incl <- function(res, param, r.crit=1.1){
  h <- which(res[param, 8, ] < r.crit, arr.ind=TRUE)
  u <- as.numeric(which(table(h[,2])==length(param)))
  return(u)
}

# Number of simulations of converged runs that are included
sample.size <- 1000

i1 <- r.incl(res1, param=c(1:5, 46:65))[1:sample.size]
i2 <- r.incl(res2, param=c(1:5, 46:65))[1:sample.size]
i3 <- r.incl(res3, param=c(1:5, 128:147))[1:sample.size]

# Calculate the deviations of the realized population growth rates from population 1
#   and the estimated annual growth rates from the IPM
diff.ipm1 <- matrix(NA, ncol=length(i1), nrow=nrow(lam))
diff.ipm2 <- matrix(NA, ncol=length(i2), nrow=nrow(lam))
diff.ipm3 <- matrix(NA, ncol=length(i3), nrow=nrow(lam))

for (s in 1:length(i1)){
  diff.ipm1[,s] <- res1[47:65,1,i1[s]] - lam[,i1[s]]
}
for (s in 1:length(i2)){
  diff.ipm2[,s] <- res2[47:65,1,i2[s]] - lam[,i2[s]]
}
for (s in 1:length(i3)){
  diff.ipm3[,s] <- res3[129:147,1,i3[s]] - lam[,i3[s]]
}


# Code for figure 5.5
mag <- 1
cex.tif <- mag * 1.25
lwd.tif <- 1.25*mag
op <- par(mar=c(3, 4.5, 4, 0.5), las=1, cex=cex.tif, "lwd", "mfrow")
layout(matrix(1:6, 2, 3, byrow=TRUE), widths=c(1, 1, 1), heights=c(1, 1.1), TRUE)

labx <- c(expression(IPM[1]), expression(IPM[2]), expression(IPM[3]))

laby <- expression('Juvenile survival ('*phi[j]*')')
boxplot(cbind(res1[1,1,i1], res2[1,1,i2], res3[1,1,i3]), outline=FALSE,
    ylab=laby, axes=FALSE, col=co, lwd=lwd.tif)
axis(2, las=1, lwd=lwd.tif)
axis(1, at=c(1, 2, 3), labels=NA, lwd=lwd.tif)
abline(h=phi[1], col="red", lwd=lwd.tif)
mtext("A", at=0.5, line=0.5, cex=1.5 * cex.tif)

laby <- expression('Adult survival ('*phi[a]*')')
boxplot(cbind(res1[2,1,i1], res2[2,1,i2], res3[2,1,i3]), outline=FALSE,
    ylab=laby, axes=FALSE, col=co, lwd=lwd.tif)
axis(2, las=1, lwd=lwd.tif)
axis(1, at=c(1, 2, 3), labels=NA, lwd=lwd.tif)
abline(h=phi[2], col="red", lwd=lwd.tif)
mtext("B", at=0.5, line=0.5, cex=1.5 * cex.tif)

laby <- expression('Fecundity (f'[1]*')')
boxplot(cbind(res1[4,1,i1], res2[4,1,i2], res3[3,1,i3]), outline=FALSE,
    ylab=laby, axes=FALSE, col=co, lwd=lwd.tif)
axis(2, las=1, lwd=lwd.tif)
axis(1, at=c(1, 2, 3), labels=NA, lwd=lwd.tif)
abline(h=f[1], col="red", lwd=lwd.tif)
mtext("C", at=0.5, line=0.5, cex=1.5 * cex.tif)

par(mar=c(4, 4.5, 3, 0.5), las=1, cex=cex.tif, lwd=4)
laby <- expression('Fecundity (f'[a]*')')
boxplot(cbind(res1[5,1,i1], res2[5,1,i2], res3[4,1,i3]), outline=FALSE,
    ylab=laby, axes=FALSE, col=co, lwd=lwd.tif)
axis(2, las=1, lwd=lwd.tif)
axis(1, at=c(1, 2, 3), labels=labx, lwd=lwd.tif)
abline(h=f[2], col="red", lwd=lwd.tif)
mtext("D", at=0.5, line=0.5, cex=1.5 * cex.tif)

laby <- expression('Bias in population growth rate ('*lambda*')')
k1 <- c(as.vector(diff.ipm1), as.vector(diff.ipm2), as.vector(diff.ipm3))
k2 <- c(rep(1, length(as.vector(diff.ipm1))), rep(2, length(as.vector(diff.ipm2))),
    rep(3, length(as.vector(diff.ipm3))))
boxplot(k1 ~ k2, outline=FALSE, ylab=laby, xlab=NA, axes=FALSE, col=co, lwd=lwd.tif)
axis(2, las=1, lwd=lwd.tif)
axis(1, at=c(1, 2, 3), labels=labx, lwd=lwd.tif)
abline(h=0, col="red", lwd=lwd.tif)
mtext("E", at=0.5, line=0.5, cex=1.5 * cex.tif)

laby <- expression('Residual error ('*sigma*')')
boxplot(cbind(res1[46,1,i1], res2[46,1,i2], res3[128,1,i3]), outline=FALSE,
    ylab=laby, axes=FALSE, col=co, lwd=lwd.tif)
axis(2, las=1, lwd=lwd.tif)
axis(1, at=c(1, 2, 3), labels=labx, lwd=lwd.tif)
abline(h=sigma, col="red", lwd=lwd.tif)
mtext("F", at=0.5, line=0.5, cex=1.5 * cex.tif)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
