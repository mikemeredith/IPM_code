# Schaub & Kéry (2022) Integrated Population Models
# Chapter 4 : Components of integrated population models
# ------------------------------------------------------

# Run time without 'browser()' approx. 7 mins

library(IPMbook) ; library(jagsUI)

# 4.5 Models for survival surveys
# ===============================

# 4.5.1 Cormack-Jolly-Seber (CJS) model for capture-recapture data
# ----------------------------------------------------------------

# 4.5.1.1 State-space formulation
# '''''''''''''''''''''''''''''''

nmarked <- 10                           # Number of marked individuals at each occasion
nyears <- 11                            # Number of years
phi <- 0.8                              # Constant apparent survival probability
p <- 0.4                                # Constant recapture probability

# Determine occasion when an individual first captured and marked
f <- rep(1:(nyears-1), each=nmarked)
nind <- length(f)                       # Total number of marked individuals

# State or ecological process
z <- array(NA, dim=c(nind, nyears))     # Empty alive/dead matrix

# Initial conditions: all individuals alive at f(i)
for (i in 1:nind){
  z[i,f[i]] <- 1
}

set.seed(1)                             # Initialize the RNGs in R
# Propagate alive/dead process forwards via transition rule:
# Alive individuals survive with probability phi
for (i in 1:nind){
  for (t in (f[i]+1):nyears){
    z[i,t] <- rbinom(1, 1, z[i,t-1] * phi)
  } #t
} #i
head(z); tail(z)                        # Not shown: look at start and end of z

# Observation process: simulate observations
y <- array(0, dim=c(nind, nyears))
for (i in 1:nind){
  y[i,f[i]] <- 1
  for(t in (f[i]+1):nyears){
    y[i,t] <- rbinom(1, 1, z[i,t] * p)
  } #t
} #i

y                                       # Complete simulated capture-recapture data set (not shown)
for (i in 1:10){                        # Look at true and observed states of first 10 individuals
  print(rbind("True state (z)" = z[i,], "Observed state (y)" = y[i,]))
  # browser()  # ~~~ take out for testing
}

# Data bundle
jags.data <- list(y=y, f=f, nind=nind, nyears=ncol(y))
str(jags.data)
# List of 4
# $ y     : num [1:100, 1:11] 1 1 1 1 1 1 1 1 1 1 ...
# $ f     : int [1:100] 1 1 1 1 1 1 1 1 1 1 ...
# $ nind  : num 100
# $ nyears: num 11

# Write JAGS model file
cat(file="model14.txt", "
model {
  # Priors and linear models
  phi.const ~ dunif(0, 1)               # Vague prior for constant phi
  p.const ~ dunif(0, 1)                 # Vague prior for constant p

  for (i in 1:nind){                    # Loop over individuals
    for (t in f[i]:(nyears-1)){         # Loop over time intervals/occasions
      phi[i,t] <- phi.const             # Here we model pattern in phi ...
      p[i,t] <- p.const                 # ... and p
    } #t
  } #i

  # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):nyears){
      # State process
      z[i,t] ~ dbern(z[i,t-1] * phi[i,t-1])
      # Observation process
      y[i,t] ~ dbern(z[i,t] * p[i,t-1])
    } #t
  } #i
}
")

# Initial values
inits <- function(){list(z=zInit(y))}

# Parameters monitored
parameters <- c("phi.const", "p.const") # Could also add "z"

# MCMC settings
ni <- 3000; nb <- 1000; nc <- 3; nt <- 1; na <- 1000

# Call JAGS from R (ART < 1 min) and check convergence
out17 <- jags(jags.data, inits, parameters, "model14.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out17) # Not shown
print(out17, 3)
#              mean     sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# phi.const   0.830  0.028   0.774   0.830   0.884    FALSE 1 1.002   834
# p.const     0.414  0.038   0.341   0.414   0.490    FALSE 1 1.001  1567

library(IPMbook); library(jagsUI)
data(woodchat5)
str(woodchat5)
# List of 4
# $ ch   : num [1:1902, 1:20] 1 1 1 1 1 1 1 1 1 1 ...
# $ age  : num [1:1902] 2 2 2 2 2 2 2 2 2 2 ...
# $ repro: num [1:929, 1:3] 6 2 2 5 3 5 3 2 3 2 ...
# $ count: num [1:20] 91 119 131 88 139 145 148 116 112 106 ...

f <- getFirst(woodchat5$ch)
last <- which(f==ncol(woodchat5$ch))
f <- f[-last]
ch <- woodchat5$ch[-last,]
age <- woodchat5$age[-last]

x <- createAge(f, age, 20, 2)

x[50:51,]                               # Most columns omitted
#      [,1] [,2] [,3] [,4] ...
# [1,]    2    2    2    2 ...
# [2,]    1    2    2    2 ...

# Data bundle
jags.data <- list(y=woodchat5$ch[-last,], f=f, nind=length(f), nyears=ncol(woodchat5$ch), x=x)
str(jags.data)
# List of 5
# $ y     : num [1:1799, 1:20] 1 1 1 1 1 1 1 1 1 1 ...
# $ f     : int [1:1799] 1 1 1 1 1 1 1 1 1 1 ...
# $ nind  : int 1799
# $ nyears: int 20
# $ x     : num [1:1799, 1:19] 2 2 2 2 2 2 2 2 2 2 ...

# Write JAGS model file
cat(file="model15.txt", "
model {
  # Priors and linear models
  for (j in 1:2){
    beta[j] ~ dunif(0, 1)               # Vague priors for age-dependent survival
  }
  p.const ~ dunif(0, 1)                 # Vague prior for constant p

  for (i in 1:nind){
    for (t in f[i]:(nyears-1)){
      phi[i,t] <- beta[x[i,t]]          # Linear model for phi ...
      p[i,t] <- p.const                 # ... and p
    } #t
  } #i

  # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):nyears){
      # State process
      z[i,t] ~ dbern(z[i,t-1] * phi[i,t-1])
      # Observation process
      y[i,t] ~ dbern(z[i,t] * p[i,t-1])
    } #t
  } #i
}
")

# Initial values
inits <- function(){list(z=zInit(jags.data$y))}

# Parameters monitored
parameters <- c("beta", "p.const")

# MCMC settings
ni <- 3000; nb <- 1000; nc <- 3; nt <- 1; na <- 1000

# Call JAGS from R (ART 5 min) and check convergence
out18 <- jags(jags.data, inits, parameters, "model15.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out18) # Not shown
print(out18)
#              mean     sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# beta[1]     0.304  0.017    0.273    0.304    0.339    FALSE 1 0.998  6000
# beta[2]     0.542  0.014    0.513    0.542    0.569    FALSE 1 1.002  6000
# p.const     0.603  0.021    0.563    0.604    0.641    FALSE 1 1.003  3450

# 4.5.1.2 Multinomial formulation
# '''''''''''''''''''''''''''''''

marr <- marray(y)
print(marr)
#         recaptured
# released Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y9 Y10 Y11 never
#      Y1   3  2  1  2  0  0  0  0   0   0     2
#      Y2   0  5  1  3  0  0  0  0   0   0     4
#      Y3   0  0  3  4  3  1  1  0   0   0     5
#      Y4   0  0  0  7  2  1  0  0   0   0     5
#      Y5   0  0  0  0 12  1  1  3   0   0     9
#      Y6   0  0  0  0  0 10  3  3   0   0    11
#      Y7   0  0  0  0  0  0  9  5   1   1     7
#      Y8   0  0  0  0  0  0  0  7   6   2     9
#      Y9   0  0  0  0  0  0  0  0   5   6    17
#      Y10  0  0  0  0  0  0  0  0   0   8    14

rel <- rowSums(marr)
print(rel)
# Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y9 Y10
# 10 13 17 15 26 27 23 24 28  22

# Bundle data
jags.data <- list(marr=marr, rel=rel, nyears=ncol(marr))
str(jags.data)
# List of 3
# $ marr  : num [1:10, 1:11] 3 0 0 0 0 0 0 0 0 0 ...
# $ rel   : num [1:10] 10 13 17 15 26 27 23 24 28 22
# $ nyears: int 11

# Write JAGS model file
cat(file="model16.txt", "
model {
  # Priors and linear models
  phi.const ~ dunif(0, 1)
  p.const ~ dunif(0, 1)

  for (t in 1:(nyears-1)){
    phi[t] <- phi.const
    p[t] <- p.const
  }

  # Likelihood
  # Define the multinomial likelihood
  for (t in 1:(nyears-1)){
    marr[t,1:nyears] ~ dmulti(pi[t,], rel[t])
  }
  # Define the cell probabilities of the m-array
  for (t in 1:(nyears-1)){
    # Main diagonal
    q[t] <- 1 - p[t]                    # Probability of non-recapture
    pi[t,t] <- phi[t] * p[t]
    # Above main diagonal
    for (j in (t+1):(nyears-1)){
      pi[t,j] <- prod(phi[t:j]) * prod(q[t:(j-1)]) * p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pi[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(nyears-1)){
    pi[t,nyears] <- 1-sum(pi[t,1:(nyears-1)])
  }
}
")

# Initial values
inits <- function(){list(phi.const=runif(1, 0, 1))}

# Parameters monitored
parameters <- c("phi.const", "p.const")

# MCMC settings
ni <- 3000; nb <- 1000; nc <- 3; nt <- 1; na <- 1000

# Call JAGS (ART <1 min), check convergence and summarize posteriors
out19 <- jags(jags.data, inits, parameters, "model16.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
  n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out19) # Not shown
print(out19)
#              mean    sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# phi.const   0.831 0.027   0.777   0.832   0.881    FALSE 1 1.001  5290
# p.const     0.415 0.037   0.345   0.414   0.488    FALSE 1 1.001  3714


# Analysis of simulated woodchat data
# ...................................

data(woodchat5)
marr <- marrayAge(ch=woodchat5$ch, age=woodchat5$age, mAge=2)

# Bundle data and produce data overview
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], nyears=dim(marr)[2], rel.j=rowSums(marr[,,1]),
    rel.a=rowSums(marr[,,2]))
str(jags.data)
# List of 5
# $ marr.j: num [1:19, 1:20] 8 0 0 0 0 0 0 0 0 0 ...
# $ marr.a: num [1:19, 1:20] 16 0 0 0 0 0 0 0 0 0 ...
# $ nyears: int 20
# $ rel.j : num [1:19] 51 53 55 65 73 66 61 76 65 75 ...
# $ rel.a : num [1:19] 36 39 44 61 61 50 43 61 51 53 ...

# Write JAGS model file
cat(file="model17.txt", "
model {
  # Priors and linear models
  phij.const ~ dunif(0, 1)
  phia.const ~ dunif(0, 1)
  p.const ~ dunif(0, 1)

  for (t in 1:(nyears-1)){
    phij[t] <- phij.const
    phia[t] <- phia.const
    p[t] <- p.const
  }

  # Capture-recapture data (CJS model with multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(nyears-1)){
    marr.j[t,1:nyears] ~ dmulti(pi.j[t,], rel.j[t])
    marr.a[t,1:nyears] ~ dmulti(pi.a[t,], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(nyears-1)){
    # Main diagonal
    q[t] <- 1 - p[t]                                # Probability of non-recapture
    pi.j[t,t] <- phij[t] * p[t]
    pi.a[t,t] <- phia[t] * p[t]
    # Above main diagonal
    for (j in (t+1):(nyears-1)){
      pi.j[t,j] <- phij[t] * prod(phia[(t+1):j]) * prod(q[t:(j-1)]) * p[j]
      pi.a[t,j] <- prod(phia[t:j]) * prod(q[t:(j-1)]) * p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pi.j[t,j] <- 0
      pi.a[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(nyears-1)){
    pi.j[t,nyears] <- 1-sum(pi.j[t,1:(nyears-1)])
    pi.a[t,nyears] <- 1-sum(pi.a[t,1:(nyears-1)])
  }
}
")

# Initial values
inits <- function(){list(phij.const=runif(1, 0, 1))}

# Parameters monitored
parameters <- c("phij.const", "phia.const", "p.const")

# MCMC settings
ni <- 3000; nb <- 1000; nc <- 3; nt <- 1; na <- 1000

# Call JAGS (ART <1 min), check convergence and summarize posteriors
out20 <- jags(jags.data, inits, parameters, "model17.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out20) # Not shown
print(out20)
#               mean    sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# phij.const   0.304 0.016   0.274   0.304   0.337    FALSE 1 1.001  6000
# phia.const   0.542 0.014   0.515   0.542   0.570    FALSE 1 1.002  1302
# p.const      0.603 0.020   0.564   0.603   0.644    FALSE 1 1.000  6000
